from position.singlepoint import SINGLEPOINTPOSITION
from position.common import *
import numpy as np
class RTKPARAM():
    def __init__(self):
        # 基站坐标，ECEF-XYZ
        self.refPoint = [0.0,0.0,0.0,0.0,0.0,0.0]
        # rtk状态
        self.rtkStatus = 0
        # 解状态
        self.sol = solustion.positionsol()

        self.nx = tool.MAXFREF
        self.na = 9
        self.tt = 0.0
        self.x = []
        self.P = []
        self.xa = []
        self.Pa = []
        self.nfix = 0
        self.neb = 0
        self.ssatList = []
        for i in range(tool.MAXSAT):
            ssat = RTKCOMMON.SSAT()
            self.ssatList.append(ssat)

        self.x = tool.zeroMat(self.nx,1)
        self.xa = tool.zeroMat(self.nx, 1)
        self.P = tool.zeroMat(self.nx, self.nx)
        self.Pa = tool.zeroMat(self.nx, self.nx)



    def initRefPoint(self,X,Y,Z):
        self.refPoint = (X,Y,Z)
    def updateSol(self,sol):
        self.sol = sol
    def getrr(self):
        return [self.x[0][0],self.x[1][0],self.x[2][0]]
class zdresParam():
    def __init__(self):
        self.y = []
        self.azel = []
        self.e = []
        self.xp = []
        self.state = 0
        self.vflg = []
        self.Pp = []
        self.xa = []
        self.v = []
        self.H = []
        self.R = []
        self.nv = 0
        self.bias=[]
    def init(self,nf,n):
        self.y = tool.zeroMat(n,nf*2)
        self.azel = tool.zeroMat(n,2)
        self.e = tool.zeroMat(n,3)
        self.vflg = tool.zeroMat(64*2*2+1,1)
        return self
    def init2(self,ny,nx):
        self.Pp = tool.zeroMat(nx,nx)
        self.xa = tool.zeroMat(nx,1)
        self.v = tool.zeroMat(ny,1)
        self.H = tool.zeroMat(nx,ny)
        self.R = tool.zeroMat(ny,ny)
        self.bias = tool.zeroMat(nx,1)
        return self


class RTKPoistion():
    def NP(self,opt):
        return 3 if opt.dynamics == 0 else 9

    def NI(self,opt):
        return 0 if opt.ionomodel != 4 else tool.MAXSAT

    def NT(self,opt):
        return 0 if opt.tropmodel < 3 else 6

    def NL(self,opt):
        return 0 if opt.glomodear != 2 else 2
    def NR(self,opt):
        return self.NP(opt)+self.NI(opt)+self.NT(opt)+self.NL(opt)

    def IB(self,s,f,opt):
        return self.NR(opt)+tool.MAXSAT*f+s-1

    def __init__(self):
        self.opt = OPTION.POSTIONOPTION()

    def exeRTKPositon(self,rtkParam,roverObs,baseObs,navDataList):
        #先计算单点位置，获取首次有效定位、速度
        rtkParam = SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(roverObs, navDataList, rtkParam)

        #baseNavList = []
        #roverNavList = []
        if len(baseObs.obsary) < 4:
            rtkParam.sol.update(rtkParam.sol.X, len(baseObs.obsary), 0, 0)
            return rtkParam.sol
        if len(roverObs.obsary) < 4:
            rtkParam.update(rtkParam.sol.X, len(roverObs.obsary), 0, 0)
            return rtkParam.sol
        #rtkParam.sol.rr = [-2333326.6741502904, 5383652.2228962239, 2492168.6526665115,0.0,0.0,0.0]

        #根据基站位置计算基于基站的卫星位置信息等
        baseNavList = self.getNavList(baseObs,navDataList)
        roverNavList = self.getNavList(roverObs,navDataList)

        nbaseObs = len(baseObs.obsary)
        nroverObs = len(roverObs.obsary)

        allObs = roverObs
        for i in range(nbaseObs):
            allObs.obsary.append(baseObs.obsary[i])

        allNavList = self.getNavList(allObs,navDataList)

        zParam = zdresParam().init(self.opt.nf,nroverObs+nbaseObs)
        #进行基站的非差残差计算
        zParam = self.zdres(0,baseObs,nbaseObs,nbaseObs,nroverObs,baseNavList,rtkParam.refPoint,1,zParam)
        if  zParam.state != 1:
            return rtkParam.sol
        #获取移动站和基站的公共卫星号、对应在allObs数组的位置。
        satComList = self.selectCommonSat(allObs,nroverObs,nbaseObs)

        rtkParam = self.updateState(rtkParam,satComList[0],satComList[1],satComList[2],len(satComList[0]),allNavList,allObs)
        niter = self.opt.iter

        xp = rtkParam.x
        zParam.xp = xp
        ny = len(satComList[0]) * self.opt.nf * 2 + 2
        zParam = zParam.init2(ny, rtkParam.nx)
        for i in range(niter):
            zParam = self.zdres(1,roverObs,nroverObs,nbaseObs,nroverObs,roverNavList,[rtkParam.x[0][0],rtkParam.x[1][0],rtkParam.x[2][0]],1,zParam)
            if zParam.state != 1:
                break
            zParam = self.ddres(rtkParam,allNavList,0,satComList[0],satComList[1],satComList[2],len(satComList[0]),zParam)
            if zParam.state != 1:
                break



        print()


    def getNavList(self,obs,navDataList):
        navList = []
        for k in range(len(obs.obsary)):
            for i in range(len(navDataList)):
                # 计算卫星位置；
                if obs.obsary[k].prn == navDataList[i].prn:
                    # 卫星观测量对应的L1伪距
                    # 整理星历数据打包程nav对象，计算卫星位置pos，计算卫星种差dts，计算卫星位置及时钟方差vare
                    #print("卫星:",navDataList[i].prn,"位置：")
                    navList.append(SATPOS.SatPos().getSatpos(navDataList[i], obs.t_obs,
                                                                 obs.obsary[k].obsfrefList[0].P))
                else:
                    continue
        return navList

    #type ==0是基站
    def zdres(self,type,obs,nobs,nbaseobs,nrover,navList,rr,index,zParam):
        basei = 0
        if type == 0:
            basei = nrover

        #for i in range(len(obs))
        nfre = self.opt.nf
        y = zParam.y
        dant=[]
        # for i in range(nobs):
        #     y_obs = []
        #     dant.append(0)
        #     for j in range(nfre*2):
        #         y_obs.append(0)
        #     y.append(y_obs)
        #print(y[0],y[0][1])
        if tool.dot(rr) <=0 :
            zParam.state =0
            return zParam
        rr_ = rr
        pos = [0.0,0.0,0.0]
        pos = RTKCOMMON.eceftopos(rr_,pos)
        #L1+L2双频计算
        lam = [tool.CLIGHT / tool.FREQ1, tool.CLIGHT / tool.FREQ2]
        for i in range(nobs):
            dant.append(0)
            r = RTKCOMMON.r_corr(navList[i],rr_)
            #print(r)
            if r <= 0 :
                continue

            e = RTKCOMMON.e_corr(navList[i],rr_)
            azel = RTKCOMMON.satAzel(pos,e)
            zParam.e[basei+i] = e
            zParam.azel[basei+i] =azel
            if azel[1] < 10*np.pi/180 :
                continue

            r += (-tool.CLIGHT)* navList[i].dts
            #print(r)
            zhd = self.tropmodel(obs.t_obs,pos,[0,90*np.pi/180],0)
            tropnmf = self.tropmapf_nmf(obs.gpsweek,obs.t_obs,pos,azel,0)
            r += tropnmf*zhd
            #print(r,tropnmf,zhd)
            for j in range(len(obs.obsary[i].obsfrefList)):
                #for k in range():
                if (obs.obsary[i].obsfrefList[j].L > 0.0) :
                    y[basei+i][j] = obs.obsary[i].obsfrefList[j].L*lam[j] - r - dant[i]
                if (obs.obsary[i].obsfrefList[j].P > 0.0):
                    y[basei+i][j+nfre] = obs.obsary[i].obsfrefList[j].P -r - dant[i]
            zParam.y = y
            zParam.state = 1
        return zParam
    def tropmodel(self,time,pos,azel,humi):
        temp0 = 15
        if pos[2] <-100 or pos[2]>1E4 or azel[1] <= 0 :
            return 0
        hgt = (0 if pos[2]<0 else pos[2])

        pres = 1013.25*np.power(1-2.2557E-5*hgt,5.2568)
        temp = temp0 - 6.5E-3*hgt+273.16
        e = 6.108*humi*np.exp((17.15*temp-4684)/(temp-38.45))

        z = np.pi /2 -azel[1]
        trph = 0.0022768 * pres / (1-0.00266*np.cos(2*pos[0])-0.00028*hgt/1E3)/np.cos(z)
        trpw = 0.002277 *(1255/temp+0.05)*e/np.cos(z)
        return trph+trpw

    def tropmapf_nmf(self,gps_week,gps_sec,pos,azel,map):
        coef = [
            [1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3],
            [2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3],
            [62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3],

            [0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5],
            [0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5],
            [0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5],

            [5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4],
            [1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3],
            [4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2]
        ]
        aht = [2.53E-5, 5.49E-3, 1.14E-3]

        ah = [0,0,0]
        aw = [0,0,0]
        el = azel[1]
        lat = pos[0]*180/np.pi
        hgt = pos[2]

        if el <= 0 :
            return 0

        ep = tool.GPST2EPOCH(gps_week,gps_sec)
        y = (tool.dayofyear(ep)-28 )/365.25+(0.5 if lat< 0 else 0)
        cosy = np.cos(2*np.pi*y)
        lat = np.fabs(lat)

        for i in range(3):
            ah[i] = self.interpc(coef[i],lat) - self.interpc(coef[i+3],lat)*cosy
            aw[i] = self.interpc(coef[i+6],lat)

        dm = (1/np.sin(el)-self.mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3

        if map >0 :
            map = self.mapf(el,aw[0],aw[1],aw[2])

        return self.mapf(el,ah[0],ah[1],ah[2])+dm

    def interpc(self,coef,lat):
        i = int(lat/15)
        if i < 1:
            return coef[0]
        elif i > 4:
            return coef[4]
        return coef[i-1]*(1 - lat/15 +i)+coef[i]*(lat/15 -i)
    def mapf(self,el,a,b,c):
        sinel = np.sin(el)
        return (1+a/(1+b/(1+c)))/(sinel+(a/(sinel+b/(sinel+c))))

    def selectCommonSat(self,allObs,nroverObs,nbaseObs):
        satList=[]
        irover = []
        ibase = []
        i = 0
        j = nroverObs
        while i < nroverObs and j < nroverObs + nbaseObs:

            if allObs.obsary[i].prn < allObs.obsary[j].prn :
                j -=1
            elif allObs.obsary[i].prn > allObs.obsary[j].prn:
                i-=1
            elif allObs.obsary[i].prn == allObs.obsary[j].prn:
                satList.append(allObs.obsary[i].prn)
                irover.append(i)
                ibase.append(j)
            i += 1
            j += 1

        # for i in range(0,nroverObs):
        #     for j in range(nroverObs,nbaseObs):
        #         if allObs.obsary[i].prn == allObs.obsary[j].prn:

        return [satList,irover,ibase]
    #self.updateState(rtkParam,satComList[0],satComList[1],satComList[2],len(satComList[0]),allNavList,allObs)
    def updateState(self,rtkParam,satList,irover,ibase,ns,nav,obs):
        tt = np.fabs(rtkParam.tt);
        rtkParam = self.updos(rtkParam,tt)
        rtkParam = self.upbias(rtkParam,tt,satList,irover,ibase,ns,nav,obs)
        return rtkParam

    def initX(self,rtkParam,xi,var,i):
        rtkParam.x[i][0] = xi
        for j in range(rtkParam.nx):
            a = (var if i==j else 0)
            rtkParam.P[j][i] = a
            rtkParam.P[i][j] = a
    def sdobs(self,obs,ir,ib,f):
        pi = obs.obsary[ir].obsfrefList[f].L if f < tool.NFREQ else obs.obsary[ir].obsfrefList[f-tool.NFREQ].P
        pj = obs.obsary[ib].obsfrefList[f].L if f < tool.NFREQ else obs.obsary[ib].obsfrefList[f-tool.NFREQ].P
        return 0 if pi == 0 or pj == 0 else pi-pj

    def gfobs(self,obs,ir,ib):

        pi = self.sdobs(obs,ir,ib,0)*tool.lam[0]
        pj = self.sdobs(obs,ir,ib,1)*tool.lam[1]
        return 0 if pi ==0 or pj==0 else pi-pj
    def detslp_gf(self,rtkParam,obs,ir,ib,nav):
        sat = obs.obsary[ir].prn
        g1 = self.gfobs(obs,ir,ib)
        if self.opt.nf <= 1 or g1 == 0:
            return 0
        g0 = rtkParam.ssatList[sat-1].gf
        rtkParam.ssatList[sat-1].gf = g1
        if g0 != 0 and np.fabs(g1-g0) > self.opt.thresslip :
            rtkParam.ssatList[sat-1].slip[0] |= 1
            rtkParam.ssatList[sat-1].slip[1] |= 1

        return rtkParam
    #detect cycle slip by LLI
    def detslp_ll(self,rtkParam,obs,i,rcv):
        sat = obs.obsary[i].prn
        for f in range(self.opt.nf):
            if obs.obsary[i].obsfrefList[f].L == 0 :
                continue
            LLI1 = (rtkParam.ssatList[sat-1].slip[f] >> 6) & 3
            LLI2 = (rtkParam.ssatList[sat-1].slip[f] >> 4) & 3
            LLI = LLI1 if rcv == 1 else LLI2

            slip = (rtkParam.ssatList[sat-1].slip[f]|obs.obsary[i].obsfrefList[f].LLI)&3
            if obs.obsary[i].obsfrefList[f].LLI & 1 :
                print("error")

            if (((LLI&2) and (not (obs.obsary[i].obsfrefList[f].LLI & 2))) or ((not(LLI & 2)) and ((obs.obsary[i].obsfrefList[f].LLI)&2))):
                slip |= 1

            if(rcv == 1):
                rtkParam.ssatList[sat-1].slip[f] = (obs.obsary[i].obsfrefList[f].LLI << 6)|(LLI2<<4)|slip
            else:
                rtkParam.ssatList[sat - 1].slip[f] = (obs.obsary[i].obsfrefList[f].LLI << 4) | (LLI1 << 6) | slip

        return rtkParam

    def upbias(self,rtkParam,tt,satList,irover,ibase,ns,nav,obs):
        #nf = 2
        for i in range(ns):
            for j in range(self.opt.nf):
                rtkParam.ssatList[satList[i]-1].slip[j] &= 0xFC
            self.detslp_ll(rtkParam,obs,irover[i],1)
            self.detslp_ll(rtkParam,obs,ibase[i],2)

            self.detslp_gf(rtkParam,obs,irover[i],ibase[i],nav)

        for f in range(self.opt.nf):
        #     for i in range(1,tool.MAXSAT+1):
        #         a = rtkParam.ssatList[i-1].outc[f] +1
        #         reset = 1 if a>self.opt.maxout else 0
        #         if self.opt.modear == tool.ARMODE_INST
            for i in range(ns):
                j = self.IB(satList[i],f,self.opt)

                rtkParam.P[j][j] += self.opt.prn[0]*self.opt.prn[0]*tt
                slip = rtkParam.ssatList[satList[i]-1].slip[f]
                if self.opt.ionomodel == 2 :
                    slip |= rtkParam.ssatList[satList[i]-1].slip[1]
                if self.opt.modear == 2 or (not (slip & 1)) :
                    continue
                rtkParam.x[j][0]= 0
                rtkParam.ssatList[satList[i]-1].lock[f] = -self.opt.minlock

            bias = tool.zeroMat(ns,1)
            for i in range(ns):
                j =0
                offset = 0.0
                if(self.opt.ionomodel != 2):
                    cp = self.sdobs(obs,irover[i],ibase[i],f)
                    pr = self.sdobs(obs,irover[i],ibase[i],f+tool.NFREQ)
                    lami = tool.lam[f]
                    if( cp == 0  or pr == 0 or lami == 0) :
                        continue
                    bias[i][0] = cp- (pr /lami)
                else:
                    print("")

                ib = self.IB(satList[i],f,self.opt)
                if( rtkParam.x[ib][0]  != 0):
                    offset += bias[i][0]-rtkParam.x[ib][0]
                    j+=1
            if j > 0 :
                for i in range(tool.MAXSAT):
                    if(rtkParam.x[self.IB(i,f,self.opt)][0] != 0) :
                        rtkParam.x[self.IB(i,f,self.opt)][0] += (offset/f)

            for i in range(ns):
                ib = self.IB(satList[i],f,self.opt)
                if (bias[i][0] == 0.0) or (rtkParam.x[ib][0] !=0.0) :
                    continue
                self.initX(rtkParam,bias[i][0],self.opt.std[0]*self.opt.std[0],self.IB(satList[i],f,self.opt))
        return rtkParam


    def updos(self,rtkParam,tt):
        Q = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        pos = [0.0,0.0,0.0]

        if(tool.dot(rtkParam.getrr()) <= 0.0):
            for i in range(3):
                self.initX(rtkParam,rtkParam.sol.rr[i],30*30,i)

            for i in range(3,6):
                self.initX(rtkParam,rtkParam.sol.rr[i],10*10,i)

            for i in range(6,9):
                self.initX(rtkParam,1E-6,10*10,i)
        variance = 0
        for i in range(3):
            variance += rtkParam.P[i][i]
        variance /= 3

        if variance > 30*30 :
            print()

        F = tool.eyeMat(rtkParam.nx)
        for i in range(6):
            F[i+3][i] = tt
        #FP = tool.zeroMat(rtkParam.nx,rtkParam.nx)
        #xp = tool.zeroMat(rtkParam.nx,1)

        for i in range(6):
            F[i][i+3] = tt
        #print(rtkParam.x)
        xp = tool.mulmatirix(rtkParam.nx,1,rtkParam.nx,F,rtkParam.x,"TN")
        FP = tool.mulmatirix(rtkParam.nx,rtkParam.nx,rtkParam.nx,F,rtkParam.P,"TN")
        rtkParam.P = tool.mulmatirix(rtkParam.nx,rtkParam.nx,rtkParam.nx,FP,F,"TT")
        rtkParam.x = xp
        Q[0] = Q[4] = self.opt.prn[3]*self.opt.prn[3]
        Q[8] = self.opt.prn[4]*self.opt.prn[4]
        Qv=[]
        pos = RTKCOMMON.eceftopos(rtkParam.getrr(),pos)
        #print(pos)
        #pos = [0.40407093076627232, 1.9797683819037855, 20.058078320696950]
        Qv = self.convEcef(pos,Q)
        for i in range(3):
            for j in range(3):
                rtkParam.P[i+6][j+6] += Qv[i+j*3]
        return rtkParam

    def convEcef(self,pos,Q):
        Qv = []
        E = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        EQ= []
        E = RTKCOMMON.xyz2enu(pos,E)
        #E^T *Q,
        for i in range(3):
            #dl = []
            for j in range(3):
                d = 0
                for k in range(3):
                    d += E[k+i*3]*Q[k+j*3]
                    #print("E[",k+i*3,"]*Q[",k+j*3,"]")
                EQ.append(d)
        #EQ^T *E
        for i in range(3):
            for j in range(3):
                d = 0
                for k in range(3):
                    d += EQ[k+i*3]*E[k+j*3]
                Qv.append(d)
        return Qv
            #EQ.append()

    def baseline(self,xu,xr):
        dr = []
        for i in range(3):
            dr.append(xu[i] - xr[i])
        return tool.norm(dr,3)
    def validobs(self,ir,ib,f,nf,y):

        a = y[ir][f] != 0
        b = y[ib][f] != 0
        c = (True if f < nf else (y[ir][f-nf] != 0 and y[ib][f-nf]!=0))
        return a and b and c

    def ddres(self,rtkParam,nav,dt,satList,iRover,iBase,ns,zParam):
        nf = self.opt.nf
        nv = 0
        posu = [0.0,0.0,0.0]
        posr = [0.0, 0.0, 0.0]
        xu = [zParam.xp[0][0],zParam.xp[1][0],zParam.xp[2][0]]
        xr = rtkParam.refPoint
        bl = self.baseline(xu,xr)
        posu = RTKCOMMON.eceftopos(xu,posu)
        posr = RTKCOMMON.eceftopos(xr,posr)
        Ri = tool.zeroMat(ns*self.opt.nf*2+2,1)
        Rj = tool.zeroMat(ns * self.opt.nf * 2 + 2, 1)
        im = tool.zeroMat(ns,1)
        tropu = tool.zeroMat(ns,1)
        tropr = tool.zeroMat(ns,1)
        dtdxu = tool.zeroMat(ns,3)
        dtdxr = tool.zeroMat(ns,3)
        for i in range(tool.MAXSAT):
            for j in range(3):
                rtkParam.ssatList[i].resp[j] = 0
                rtkParam.ssatList[i].resc[j] = 0
        #计算电离层、对流层延迟；
        #for i in range(ns):
        for m in range(2):
            for f in range(nf*2):
                i = -1
                for j in range(ns):
                    sysi = rtkParam.ssatList[satList[j]-1].sys
                    #0 = gps/sbas/qzss，1=glo
                    if((m == 0 and sysi == 1) or (m == 1 and sysi != 1)):
                        continue
                    if(not self.validobs(iRover[j],iBase[j],f,nf,zParam.y)) :
                        continue
                    if(i<0) :
                        i = j
                    elif(zParam.azel[iRover[j]][1] >= zParam.azel[iRover[i]][1]):
                        i = j

                if i < 0 : continue
                for j in range(ns):
                    sysi = rtkParam.ssatList[satList[i]-1].sys
                    sysj = rtkParam.ssatList[satList[j]-1].sys
                    if((m == 0 and sysj == 1) or (m == 1 and sysj != 1)) : continue
                    if(not(self.validobs(iRover[j],iBase[j],f,nf,zParam.y))) : continue
                    ff = np.mod(f,nf)
                    lami = tool.lam[ff]
                    lamj = tool.lam[ff]
                    if lami <= 0 or lamj <= 0 : continue
                    if zParam.H :
                        Hi = zParam.H[nv*rtkParam.nx]
                    zParam.v[nv][0] = (zParam.y[iRover[i]][f]-zParam.y[iBase[i]][f])\
                                      -(zParam.y[iRover[j]][f] - zParam.y[iBase[j]][f])
                    if zParam.H :
                        for k in range(3):
                            Hi[k] = -zParam.e[iRover[i]][k]+zParam.e[iRover[j][k]]
                    if f < nf:
                        if self.opt.ionomodel != 2:
                            zParam.v[nv] -= lami*rtkParam.x[self.IB(satList[i],f,self.opt)][0]\
                                            -lamj*rtkParam.x[self.IB(satList[j],f,self.opt)]
                            if zParam.H:
                                Hi[self.IB(satList[i],f,self.opt)] = lami
                                Hi[self.IB(satList[j],f,self.opt)] = -lamj
                        else:
                            print()
                    if f < nf :
                        rtkParam.ssatList[satList[j]-1].resc[f] = zParam.v[nv]
                    else:
                        rtkParam.ssatList[satList[j] - 1].resp[f-nf] = zParam.v[nv]







        return zParam













