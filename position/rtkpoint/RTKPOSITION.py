from position.singlepoint import SINGLEPOINTPOSITION
from position.common import *
import numpy as np

#RTK参数类
class RTKPARAM():
    def __init__(self):
        # 基站坐标，ECEF-XYZ
        self.refPoint = [0.0,0.0,0.0,0.0,0.0,0.0]
        # rtk状态
        self.rtkStatus = 0
        # 解状态
        self.sol = solustion.positionsol()
        #lambda状态
        self.LAMBDAStat = 0
        #最大频率数
        self.nx = tool.MAXFREF
        self.na = 9
        self.tt = 0.0
        #float解数据,前9个分别是位置、速度、加速度参数，后面前maxsat是L1的相位差，后面maxsat是L2 的相位差
        self.x = []
        #float解精度，前9个是卫星、速度、加速度的误差模型参数
        self.P = []
        #固定解数组
        self.xa = []
        #固定解精度
        self.Pa = []
        self.bias = []
        self.nfix = 0
        self.neb = 0
        # 每个卫星对应各自频点的状态列表
        self.ssatList = []
        for i in range(tool.MAXSAT):
            ssat = RTKCOMMON.SSAT()
            self.ssatList.append(ssat)
        self.x = tool.zeroMat(self.nx,1)
        self.xa = tool.zeroMat(self.nx, 1)
        self.bias = tool.zeroMat(self.nx,1)
        self.P = tool.zeroMat(self.nx, self.nx)
        self.Pa = tool.zeroMat(self.nx, self.nx)



    def initRefPoint(self,X,Y,Z):
        self.refPoint = (X,Y,Z)
    def updateSol(self,sol):
        self.sol = sol
    #获取坐标
    def getrr(self):
        return [self.x[0][0],self.x[1][0],self.x[2][0]]

#残差参数类
class zdresParam():
    def __init__(self):
        # y数组的结构：
        # 卫星的L1的载波，L2的载波，L1的伪距，L2的伪距的残差值
        self.y = []
        # azel是记录每个卫星的方位及截止角
        self.azel = []
        # 记录接收机到卫星的向量
        self.e = []
        self.xp = []
        self.state = 0
        self.vflg = []
        self.Pp = []
        self.v = []
        # H矩阵的 1-3是保存观测卫星的向量- 参考卫星的向量
        self.H = []
        self.R = []
        self.nv = 0
        self.filterinfo = 0
    def init(self,nf,n):
        self.y = tool.zeroMat(n,nf*2)
        self.azel = tool.zeroMat(n,2)
        self.e = tool.zeroMat(n,3)
        self.vflg = tool.zeroMat(64*2*2+1,1)
        return self
    def init2(self,ny,nx):
        self.v = tool.zeroMat(ny,1)
        self.H = tool.zeroMat(nx,ny)
        self.R = tool.zeroMat(ny,ny)
        return self


class RTKPoistion():
    def NP(self,opt):
        return 3 if opt.dynamics == 0 else 9  #9

    def NI(self,opt):
        return 0 if opt.ionomodel != 4 else tool.MAXSAT  #0

    def NT(self,opt):
        return 0 if opt.tropmodel < 3 else 6 #0

    def NL(self,opt):
        return 0 if opt.glomodear != 2 else 2  #0
    def NR(self,opt):
        return self.NP(opt)+self.NI(opt)+self.NT(opt)+self.NL(opt)  #9

    #IB是获取一个卫星一个频点在矩阵中的位置；此option中，前9位分别是位置、速度、加速度的参数。
    # 所以IB输出是NR+最大卫星数*频率数2+卫星号-1（数据从0开始）
    def IB(self,s,f,opt):
        return self.NR(opt)+tool.MAXSAT*f+s-1   #9+s-1+117

    def __init__(self):
        self.opt = OPTION.POSTIONOPTION()

    def exeRTKPositon(self,rtkParam,roverObs,baseObs,navDataList):
        #上一次有效解算时间
        prev_epocht = rtkParam.sol.gpssec

        #先计算单点位置，获取首次有效定位、速度
        rtkParam = SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(roverObs, navDataList, rtkParam)
        rtkParam.rtkStatus = self.opt.mode
        #baseNavList = []
        #roverNavList = []
        if len(baseObs.obsary) < 4:
            rtkParam.update(rtkParam.sol.X, len(baseObs.obsary), 0, 0)
            return rtkParam.sol
        if len(roverObs.obsary) < 4:
            rtkParam.update(rtkParam.sol.X, len(roverObs.obsary), 0, 0)
            return rtkParam.sol

        #获取当前历元时间与上一历元时间差：
        if prev_epocht != 0 : rtkParam.tt = rtkParam.sol.gpssec - prev_epocht

        #rtkParam.sol.rr = [-2333326.6741502904, 5383652.2228962239, 2492168.6526665115,0.0,0.0,0.0]

        #根据基站位置计算基于基站的卫星位置信息等
        baseNavList = self.getNavList(baseObs,navDataList)
        #获取移动在的卫星位置信息
        roverNavList = self.getNavList(roverObs,navDataList)
        #计算基站、移动站在该时刻的观测量数
        nbaseObs = len(baseObs.obsary)
        nroverObs = len(roverObs.obsary)
        #将基站、移动站的观测量数据进行合并，在移动站观测量数据后面累加基站数据
        allObs = roverObs
        for i in range(nbaseObs):
            allObs.obsary.append(baseObs.obsary[i])
        #获取新的导航数据数组
        allNavList = self.getNavList(allObs,navDataList)
        #初始化残差参数类
        zParam = zdresParam().init(self.opt.nf,nroverObs+nbaseObs)

        #进行基站的非差残差计算
        zParam = self.zdres(0,baseObs,nbaseObs,nbaseObs,nroverObs,baseNavList,rtkParam.refPoint,1,zParam)
        if  zParam.state != 1:
            return rtkParam.sol
        #获取移动站和基站的公共卫星号，对应在allObs数组的位置赋值给对应数组。
        satComList = self.selectCommonSat(allObs,nroverObs,nbaseObs)
        satList = satComList[0]
        irover = satComList[1]
        ibase = satComList[2]
        ns = len(satList)
        #对状态进行更新
        rtkParam = self.updateState(rtkParam,satList,irover,ibase,ns,allNavList,allObs)
        niter = self.opt.iter
        print("single point : " ,rtkParam.sol.rr,"  qr : ",rtkParam.sol.qr)
        #xp =
        zParam.xp = rtkParam.x
        ny = len(satComList[0]) * self.opt.nf * 2 + 2
        zParam = zParam.init2(ny, rtkParam.nx)
        for i in range(niter):
            #计算移动站的非差残差
            zParam = self.zdres(1,roverObs,nroverObs,nbaseObs,nroverObs,roverNavList,[rtkParam.x[0][0],rtkParam.x[1][0],rtkParam.x[2][0]],1,zParam)
            if zParam.state != 1:
                break
            #计算双差
            zParam = self.ddres(rtkParam,allNavList,0,satList,irover,ibase,ns,zParam,1)
            if zParam.nv <=0 :
                break

            zParam.Pp = rtkParam.P
            zParam = self.filter(zParam, rtkParam.nx, zParam.nv)
            if zParam.filterinfo != 1:
                break
            zParam = self.ddres(rtkParam,allNavList,0,satList,irover,ibase,ns,zParam,0)
            if self.valpos(rtkParam,zParam,4) == 1:
                rtkParam.x = zParam.xp
                rtkParam.P = zParam.Pp
                rtkParam.sol.ns = 0
                for j in range(ns):
                    for f in range(self.opt.nf):
                        if rtkParam.ssatList[satList[i]-1].vsat[f] !=0 : continue
                        rtkParam.ssatList[satList[i]-1].lock[f] += 1
                        rtkParam.ssatList[satList[i]-1].outc[f] = 0
                        if f == 0 :
                            rtkParam.sol.ns +=1
                if rtkParam.sol.ns < 4 : rtkParam.rtkStatus = 0
            else:
                rtkParam.rtkStatus = 0




        #if rtkParam.rtkStatus != 0 and self.LAMBDA(rtkParam) > 1:
            #print()

        if rtkParam.rtkStatus == 4:
            for i in range(3):
                rtkParam.sol.rr[i] = rtkParam.xa[i][0]
                rtkParam.sol.qr[i] = rtkParam.Pa[i][i]
            rtkParam.sol.qr[3] = rtkParam.Pa[0][1]
            rtkParam.sol.qr[4] = rtkParam.Pa[2][1]
            rtkParam.sol.qr[5] = rtkParam.Pa[0][2]
        else:
            for i in range(3):
                rtkParam.sol.rr[i] = rtkParam.x[i][0]
                rtkParam.sol.qr[i] = rtkParam.P[i][i]
            rtkParam.sol.qr[3] = rtkParam.P[0][1]
            rtkParam.sol.qr[4] = rtkParam.P[2][1]
            rtkParam.sol.qr[5] = rtkParam.P[0][2]
        print("float point : ", rtkParam.sol.rr, "  qr : ", rtkParam.sol.qr)




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
    #计算载波相位或伪距的非差残差。残差是指通过已知参数计算出来的值，再被实际记录的伪距/载波减去，得到的差值。
    #根据载波方程 ψ
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
            #计算接收机到卫星的几何距离
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
            #几何距离减去上时钟偏差*光速
            r += (-tool.CLIGHT)* navList[i].dts
            #print(r)
            zhd = self.tropmodel(obs.t_obs,pos,[0,90*np.pi/180],0)
            tropnmf = self.tropmapf_nmf(obs.gpsweek,obs.t_obs,pos,azel,0)
            #几何距离 加上 对流层延时。
            r += tropnmf*zhd
            #print(r,tropnmf,zhd)
            for j in range(len(obs.obsary[i].obsfrefList)):
                #计算残差：如果有载波，则载波*频率 - 几何距离,否则用伪距-几何距离
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
        #计算基站与移动站的同一个卫星伪距或载波差，再乘以对应的频率
        pi = self.sdobs(obs,ir,ib,0)*tool.lam[0]
        pj = self.sdobs(obs,ir,ib,1)*tool.lam[1]
        return 0 if pi ==0 or pj==0 else pi-pj
    def detslp_gf(self,rtkParam,obs,ir,ib,nav):
        sat = obs.obsary[ir].prn
        g1 = self.gfobs(obs,ir,ib)
        if self.opt.nf <= 1 or g1 == 0:
            return rtkParam
        g0 = rtkParam.ssatList[sat-1].gf
        rtkParam.ssatList[sat-1].gf = g1
        if g0 != 0 and np.fabs(g1-g0) > self.opt.thresslip :
            rtkParam.ssatList[sat-1].slip[0] |= 1
            rtkParam.ssatList[sat-1].slip[1] |= 1

        return rtkParam
    #detect cycle slip by LLI
    #LLI 说明：
    # 当LLI为0或空表示确定周跳或未知；
    # 其中bit0:为1时，表示当前记录与上一次出现失锁，有整周周跳的可能；
    # 其中Bit1：为1时，表示有半周模糊度或周跳可能
    def detslp_ll(self,rtkParam,obs,i,rcv):
        sat = obs.obsary[i].prn
        for f in range(self.opt.nf):
            if obs.obsary[i].obsfrefList[f].L == 0 :
                continue
            #
            LLI1 = (rtkParam.ssatList[sat-1].slip[f] >> 6) & 3
            LLI2 = (rtkParam.ssatList[sat-1].slip[f] >> 4) & 3
            LLI = LLI1 if rcv == 1 else LLI2
            #3 = 0011 检测该颗卫星在该频点（L1/L2）下的LLI的第0、1两位。
            slip = (rtkParam.ssatList[sat-1].slip[f]|obs.obsary[i].obsfrefList[f].LLI)&3
            if obs.obsary[i].obsfrefList[f].LLI & 1 :
                print("error")
            #对于第一次
            #LLI&2 ==1 && 卫星的LLI&2不等于0  or LLI&2 ==0
            if (((LLI&2) and (not (obs.obsary[i].obsfrefList[f].LLI & 2))) or ((not(LLI & 2)) and ((obs.obsary[i].obsfrefList[f].LLI)&2))):
                slip |= 1

            if(rcv == 1):
                #LLI = 0000 0001,《《6 0100 0000 ,LLI2 = 0000 0001  <<4 0001 0000
                #slip最终表示为8位二进制，其中前8、7位是移动站LLI，6、5是基站的LLI，4-1为是周跳
                #这里是将slip赋值到ssatlist中，也就是说，记录了这个卫星的在基站、移动站的LLI、SLIP情况
                rtkParam.ssatList[sat-1].slip[f] = (obs.obsary[i].obsfrefList[f].LLI << 6)|(LLI2<<4)|slip
            else:
                rtkParam.ssatList[sat - 1].slip[f] = (obs.obsary[i].obsfrefList[f].LLI << 4) | (LLI1 << 6) | slip

        return rtkParam

    def upbias(self,rtkParam,tt,satList,irover,ibase,ns,nav,obs):
        #nf = 2
        for i in range(ns):
            for j in range(self.opt.nf):
                #ssatlist的中的slip赋值为0，基站、移动在的LLI不变（因为本身就是数据中包含的）
                rtkParam.ssatList[satList[i]-1].slip[j] &= 0xFC# 11111100
            #检测基站、移动的LLI、slip，并重新赋值到ssatlist中。
            self.detslp_ll(rtkParam,obs,irover[i],1)
            self.detslp_ll(rtkParam,obs,ibase[i],2)

            self.detslp_gf(rtkParam,obs,irover[i],ibase[i],nav)

        for f in range(self.opt.nf):
        #     for i in range(1,tool.MAXSAT+1):
        #         a = rtkParam.ssatList[i-1].outc[f] +1
        #         reset = 1 if a>self.opt.maxout else 0
        #         if self.opt.modear == tool.ARMODE_INST
            #重置P、x数据
            for i in range(ns):
                #P数组是最大
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
                    #通过同一个卫星的基站、移动的相位差 - （码差/频率） 获取近似估计的相位偏差
                    cp = self.sdobs(obs,irover[i],ibase[i],f)#相位差
                    pr = self.sdobs(obs,irover[i],ibase[i],f+tool.NFREQ)#码差
                    lami = tool.lam[f]
                    if( cp == 0  or pr == 0 or lami == 0) :
                        continue
                    #波长lami = 光速C / 频点频率fi
                    #载波相位（carrier-phase）
                    #载波观测量（相位矩、phase-range） = carrier-phase * lami
                    #载波观测量  = 载波相位*波长  =  伪距观测量 + 波长*相位差（bias）+误差
                    #相位差 = 载波相位 -  伪距观测量 / 波长
                    #zParam.v[nv][0] -= lami*rtkParam.x[self.IB(satList[i],f,self.opt)][0] -lamj*rtkParam.x[self.IB(satList[j],f,self.opt)][0]
                    #zParam.v[nv][0] = y(双差)- (lami* cp_u_i- pr_u_i -lami *cp_r_j - pr_r_j)
                    #这里通过移动站、基站的同一个卫星的载波差、伪距差，
                    #cp1 - pr1/lami - cp2 - pr2 / lami ,得出两者的相位差的差。
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

            #设置初始的相位偏差（phase-bias）状态
            for i in range(ns):
                ib = self.IB(satList[i],f,self.opt)
                if (bias[i][0] == 0.0) or (rtkParam.x[ib][0] !=0.0) :
                    continue
                self.initX(rtkParam,bias[i][0],self.opt.std[0]*self.opt.std[0],self.IB(satList[i],f,self.opt))
        return rtkParam

    #
    def updos(self,rtkParam,tt):
        tt = 0.01
        Q = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        pos = [0.0,0.0,0.0]

        #如果是第一次定位，给rtk参数类赋值x数组、P数组
        #x是状态矩阵，P是误差模型矩阵
        #x =(位置、速度、加速度，各个移动站卫星的L1的相位差，各个移动站卫星的L2的相位差)= （x,y,z,vx,vy,vz,accx,accy,accz,b1i。。。b1n,b2i...b2n,b5i..b5n）
        #预测噪声P = [位置误差模型3*3
        #                    速度误差模型3*3
        #                                  加速度误差模型3*3
        #                                                移动站卫星数的L1的误差模型n*n
        #                                                                         移动站卫星数的L1的误差模型m*m（m<=n）]
        #当第二次定位时，这里的rtkparam中的x就是上一次最优解。
        if(tool.dot(rtkParam.getrr()) <= 0.0):
            for i in range(3):
                #赋值位置
                self.initX(rtkParam,rtkParam.sol.rr[i],30*30,i)

            for i in range(3,6):
                #赋值速度
                self.initX(rtkParam,rtkParam.sol.rr[i],10*10,i)

            for i in range(6,9):
                #赋值加速度x:6-9:1E-6,P: 7:7~ 9:9 = 100
                self.initX(rtkParam,1E-6,10*10,i)
        variance = 0
        for i in range(3):
            variance += rtkParam.P[i][i]
        variance /= 3

        if variance > 30*30 :
            print()

        F = tool.eyeMat(rtkParam.nx)
        #给F数组中的41、52、63、74、85、96,14、25、36、47、58、69赋值tt时间差。tt为0即不管时间差
        #for i in range(6):
            #F[i+3][i] = tt
        #FP = tool.zeroMat(rtkParam.nx,rtkParam.nx)
        #xp = tool.zeroMat(rtkParam.nx,1)
        for i in range(6):
            F[i][i+3] = tt
        #print(rtkParam.x)

        #数学模型：P（下一个时刻，更新转态前的协方差矩阵） =
        #    F（历元从t到t+1的系统噪声转换矩阵） *P（上一个时刻，更新转态后的协方差矩阵） *F_t（历元从t到t+1的系统噪声转换矩阵的倒置）
        #    + Q（历元从t到t+1的系统噪声协方差矩阵）
        #
        #    x（状态矩阵，历元t+1，未更新状态前） = F（历元从t到t+1的系统噪声转换矩阵） * x（历元t，更新状态后）
        #x矩阵乘以F矩阵，更新时间差tt，这里tt是0，所以xp就=x
        xp = tool.mulmatirix(rtkParam.nx,1,rtkParam.nx,F,rtkParam.x,1,0,"NN",[])
        #与上同理，FP = P
        FP = tool.mulmatirix(rtkParam.nx,rtkParam.nx,rtkParam.nx,F,rtkParam.P,1,0,"NN",[])
        #FP*F 同理，P= FP
        rtkParam.P = tool.mulmatirix(rtkParam.nx,rtkParam.nx,rtkParam.nx,FP,F,1,0,"NT",[])
        #重新赋值x
        rtkParam.x = xp
        #prn[3]\[4]分别是平面加速度、垂直加速度的噪声
        #Q就是噪声的协方差矩阵:
        # 0 3 6
        # 1 4 5
        # 2 5 8
        Q[0] = Q[4] = self.opt.prn[3]*self.opt.prn[3]
        Q[8] = self.opt.prn[4]*self.opt.prn[4]
        Qv=[]
        pos = RTKCOMMON.eceftopos(rtkParam.getrr(),pos)
        #print(pos)
        #pos = [0.40407093076627232, 1.9797683819037855, 20.058078320696950]
        #Qv = [ 0 0 0
        #       0 0 0
        #       0 0 0
        #             Q(3*3)]
        Qv = self.convEcef(pos,Q)
        for i in range(3):
            for j in range(3):
                rtkParam.P[i+6][j+6] += Qv[i+j*3]
        return rtkParam

    def convEcef(self,pos,Q):
        # 0 1 2
        # 3 4 5
        # 6 7 8
        Qv = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        E = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        EQ= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        E = RTKCOMMON.xyz2enu(pos,E)
        #E_T *Q,
        for i in range(3):
            #dl = []
            for j in range(3):
                d = 0
                for k in range(3):
                    #E ： 0 1 2 * Q 0 1 2/0 1 2 * 3 4 5
                    d += E[k+i*3]*Q[k+j*3]
                    #print("E[",k+i*3,"]*Q[",k+j*3,"]")
                EQ[i + j*3] = d
        #EQ :     = EQ :
        #0  3  6
        #1  4  7
        #2  5  8
        #EQ *E
        for i in range(3):
            for j in range(3):
                d = 0
                for k in range(3):
                    d += EQ[i+k*3]*E[k+j*3]
                Qv[i+j*3] = d
        return Qv
            #EQ.append()

    def baseline(self,xu,xr):
        dr = []
        for i in range(3):
            dr.append(xu[i] - xr[i])
        return tool.norm(dr,3)

    #如果是多频解算，必须满足L1、L2的伪距、载波都有
    def validobs(self,ir,ib,f,nf,y):
        a = y[ir][f] != 0
        b = y[ib][f] != 0
        c = (True if f < nf else (y[ir][f-nf] != 0 and y[ib][f-nf]!=0))
        return a and b and c

    def ddres(self,rtkParam,nav,dt,satList,iRover,iBase,ns,zParam,HStatus):
        nv = 0
        nf = self.opt.nf
        b = 0
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
        nb = []
        for i in range(tool.NFREQ*2*2+2):
            nb.append(0)
        for i in range(tool.MAXSAT):
            for j in range(3):
                #将伪距残差、载波残差重置为0
                rtkParam.ssatList[i].resp[j] = 0
                rtkParam.ssatList[i].resc[j] = 0
        #计算电离层、对流层延迟，暂时忽略

        #这一步生成H矩阵，其中，h（x） = （h（L1载波），h（L2载波），h（L5载波），h（L1伪距），h（L2伪距），h（L5伪距））
        #h（L1载波） = p+lam1*B1 。。。
        #h（L1伪距） = p =sqrt((xs-xu)^2+(ys-yu)^2+(zs-zu)^2)
        #h(rb,L1载波) = h（r,L1载波）-h（b,L1载波）。。。。
        #H(x)是对h（x）对X0（x0,y0,z0）求偏导数，
        #H（x） = (ex,ey,ez），及用户机到卫星的几何距离的三个方向的向量。
        #H = -D*E = [ ex1(L1C),ex2(L1C)....exn(L1C),ex1(L2C)...exn(L2C),ex1(L1P)...exn(L1P),ex1(L2P)..exn(L2P)
        #      ey1(L1C) ...
        #      ez1(L1C) ....
        #D矩阵是单差矩阵，列都为的，表示该位置的卫星号作为单差参考卫星，这里是截止角最高的卫星
        #D = [ 1 -1  0 ...
        #      1  0 -1 ...
        #      ....
        #      1  0 ....-1 ]
        for m in range(2):
            for f in range(nf*2): #f :0 - 3 ,nf : 0 - 1
                i = -1
                #查找参考卫星钟截止角最高的卫星号i作为用于做双差的参考卫星
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
                    #不需要计算截止角最高的卫星
                    if i == j : continue

                    sysi = rtkParam.ssatList[satList[i]-1].sys
                    sysj = rtkParam.ssatList[satList[j]-1].sys
                    if((m == 0 and sysj == 1) or (m == 1 and sysj != 1)) : continue
                    if(not(self.validobs(iRover[j],iBase[j],f,nf,zParam.y))) : continue
                    ff = np.mod(f,nf)
                    lami = tool.lam[ff]
                    lamj = tool.lam[ff]
                    if lami <= 0 or lamj <= 0 : continue
                    # Hi = []
                    # if zParam.H :
                    #     for h in range(rtkParam.nx):
                    #         Hi.append(zParam.H[h][zParam.nv])
                    #     #Hi = zParam.H[zParam.nv*rtkParam.nx]
                    #计算双差，
                    # 将计算参考卫星i的单差：（移动站截止角最高卫星的残差-基站对应的卫星残差）
                    # 再减去其它观测卫星j的单差 （移动站其它卫星的残差-基站对应卫星的残差）
                    zParam.v[nv][0] = (zParam.y[iRover[i]][f]-zParam.y[iBase[i]][f])\
                                      -(zParam.y[iRover[j]][f] - zParam.y[iBase[j]][f])
                    if HStatus == 1 :
                        for k in range(3):
                            #Hi[k] = -zParam.e[iRover[i]][k]+zParam.e[iRover[j]][k]
                            # 计算参考卫星的与观测卫星的向量差（实际是H = h（x）对X（x，y，z）=X0(x0，y0，z0)=x0求偏导）。
                            zParam.H[k][nv] = -zParam.e[iRover[i]][k] + zParam.e[iRover[j]][k]
                    if f < nf:
                        if self.opt.ionomodel != 2:
                            #计算双差相位差项，x数组在第10位开始是相位偏差，
                            # V是双差 - （最高截止角的相位偏差*频率 - 当前卫星的相位偏差*频率）
                            zParam.v[nv][0] -= lami*rtkParam.x[self.IB(satList[i],f,self.opt)][0]\
                                            -lamj*rtkParam.x[self.IB(satList[j],f,self.opt)][0]
                            #更新H数组中对应的位置的频率，当前卫星的频率为-
                            if HStatus == 1:

                                zParam.H[self.IB(satList[i],f,self.opt)][nv] = lami
                                zParam.H[self.IB(satList[j],f,self.opt)][nv] = -lamj
                        else:
                            print()
                    #V在循环中一直更新，对应L1、L2的伪距和载波，f = 0、1表示载波，2、3表示伪距
                    #此处更新卫星列表中对应卫星的L1、L2的伪距、载波的双差项。
                    if f < nf :
                        rtkParam.ssatList[satList[j]-1].resc[f] = zParam.v[nv][0]
                    else:
                        rtkParam.ssatList[satList[j] - 1].resp[f-nf] = zParam.v[nv][0]

                    if(self.opt.maxinno  > 0 and np.fabs(zParam.v[nv][0])> self.opt.maxinno):
                        if f < nf:
                            rtkParam.ssatList[satList[i]-1].rejc[f] += 1
                            rtkParam.ssatList[satList[j]-1].rejc[f] += 1
                        continue
                    #单差测量误差方差矩阵初始化
                    #Ri是截止角最高的卫星的协方差，Rj表示除截止角最高的卫星外的其它卫星的协方差，
                    #其保存参数与nv相关联，满足f从0-3，也就是L1、L2载波，L1、L2伪距的顺序，结合nb数组可以提取出来。
                    Ri[nv][0] = self.varerr(satList[i],sysi,zParam.azel[iRover[i]][1],bl,dt,f,self.opt)
                    Rj[nv][0] = self.varerr(satList[j],sysj,zParam.azel[iRover[j]][1],bl,dt,f,self.opt)
                    #设置有效卫星位、及flag
                    if self.opt.mode > 1 :
                        if f < nf :
                            rtkParam.ssatList[satList[i]-1].vsat[f] = 1
                            rtkParam.ssatList[satList[j]-1].vsat[f] = 1
                        else:
                            rtkParam.ssatList[satList[i] - 1].vsat[f-nf] = 1
                            rtkParam.ssatList[satList[j] - 1].vsat[f-nf] = 1
                    zParam.vflg[nv][0] = (satList[i]<<16)|(satList[j]<<8)|((0 if f < nf else 1)<<4)|(np.mod(f,nf))
                    #计算出有效的载波、伪距卫星数，b对应0-3，表示L1、L2的载波、L1、L2的伪距有效计数
                    nb[b] += 1
                    #nv是移动站中L1、L2对应基站有效的双差的计数
                    nv+=1
                zParam.nv = nv

                if f < nf:
                    #计算载波残差的平均值，再各自减去平均值
                    s = 0
                    for j in range(tool.MAXSAT):
                        s += rtkParam.ssatList[j].resc[f]
                    # 有效计数+1表示总卫星数
                    s /= (nb[b] +1)
                    for j in range(tool.MAXSAT):
                        if j == (satList[i]-1) or rtkParam.ssatList[j].resc[f] != 0 :
                            rtkParam.ssatList[j].resc[f] -= s
                else:
                    s = 0
                    for j in range(tool.MAXSAT):
                        s += rtkParam.ssatList[j].resp[f-nf]
                    s /= (nb[b] + 1)
                    for j in range(tool.MAXSAT):
                        if j == (satList[i] - 1) or rtkParam.ssatList[j].resp[f-nf] != 0:
                            rtkParam.ssatList[j].resp[f-nf] -= s
                b +=1
            zParam = self.ddcov(nb,b,Ri,Rj,zParam)
        return zParam
    #将参考卫星i的L1、L2的载波和伪距的方差转换为协方差矩阵，
    #R = [ D*Ri_L1C*D_T
    #                  D * Ri_L2C * D_T
    #                                  D*Ri_L1P*D_T
    #                                               D*Ri_L2P*D_T
    #D= [1 -1  0  0  0  0
    #    1  0 -1 ...
    #    ...
    #    1              -1
    # ]
    #Ri = [ Ri  0   0  0 0 0
    #       0   Rj1 0
    #
    #         ......       Rj5]
    # ]
    def ddcov(self,nb,n,Ri,Rj,zParam):
        for i in range(zParam.nv):
            for j in range(zParam.nv):
                zParam.R[i][j] = 0
        k = 0
        for b in range(n):
            for i in range(nb[b]):
                for j in range(nb[b]):
                    #生成标准矩阵，从数据上表达，就是先生成列，
                    #nb[0]*nb[0] .....
                    #0 ....nb[0] nb[1]*nb[1]
                    #0 ....nb[0] ......nb[1] .....
                    #0 .....
                    #0............................nb[n]*nb[n]
                    zParam.R[k+j][k+i] = (Ri[k+i][0]+(Rj[k+i][0] if i == j else 0))
            k+=nb[b]

        return zParam
    #通过设置的标准差参数，计算出方差。
    #a是载波相位的误差标准差的基本项，a^2是方差
    #b是载波相位的误差标准差的仰角相关项，(b/sinel)^2就是方差
    #计算公式；
    #这里是输出2倍σ，2 * eratio*（a^2+(b/sinel)^2）
    def varerr(self,satNo,sys,el,bl,dt,f,opt):
        a = b = c = opt.err[3]*bl/1E4
        d = tool.CLIGHT*opt.sclkstab*dt
        fact = 1
        nf = opt.nf
        sinel = np.sin(el)

        if(f >= nf and opt.ena[0] > 0):
            print()
        elif(f < nf and opt.ena[1] > 0):
            print()
        else:
            #伪距的fact标度因子为100，载波为1
            if(f >= nf) : fact = opt.eratio[f - nf]
            if(fact <= 0) : fact = opt.eratio[0]
            fact *= 1
            a = fact * opt.err[1]#0.3
            b = fact * opt.err[2]#0.3
        return 2.0* (3 if opt.ionomodel == 2 else 1)*(a*a +b*b/sinel/sinel +c*c) +d*d

    def filter(self,zParam,n,m):
        ix = tool.zeroMat(n,1)
        k = 0
        for i in range(n):
            if zParam.xp[i][0] != 0 and zParam.Pp[i][i] > 0 :
                ix[k][0] = i
                k += 1
        x_ = tool.zeroMat(k,1)
        xp_ = tool.zeroMat(k,1)
        P_ = tool.zeroMat(k,k)
        Pp_ = tool.zeroMat(k,k)
        H_ = tool.zeroMat(k,m)
        I = tool.eyeMat(k)
        for i in range(k):
            x_[i][0] = zParam.xp[ix[i][0]][0]
            for j in range(k):
                P_[i][j] = zParam.Pp[ix[i][0]][ix[j][0]]
            for j in range(m):
                #H矩阵是每一列为一个卫星的数据。
                H_[i][j] = zParam.H[ix[i][0]][j]
        Q = zParam.R
        xp_ = x_
        #K = P*H*(H'*P*H + R)^-1   xp = x + K* v   Pp = (I - K* H')*P
        #P * H 就是每个卫星的每个参数都乘对应的误差模型参数
        F = tool.mulmatirix(k,m,k,P_,H_,1,0,"NN",[])#F:k行m列
        #Q = H' * P * H
        Q = tool.mulmatirix(m,m,k,H_,F,1,1,"TN",Q)#Q：k行K列


        #try:
        Q_1 = np.linalg.inv(Q)
        #except np.linalg.LinAlgError:
        # K = F * Q_1_T
        K = tool.mulmatirix(k,m,m,F,Q_1,1,0,"NN",[])#K k col m row
        #xp = K*v +xp_
        xp_ = tool.mulmatirix(k,1,m,K,zParam.v,1,1,"NN",xp_)  #xp: k col,1 rows
        #I = I - K*H_t
        I = tool.mulmatirix(k,k,m,K,H_,-1,1,"NT",I) # I: k col k row
        #Pp = I*P_
        Pp_ = tool.mulmatirix(k,k,k,I,P_,1,0,"NN",[])

        for i in range(k):
            zParam.xp[ix[i][0]][0] = xp_[i][0]
            for j in range(k):
                zParam.Pp[ix[i][0]][ix[j][0]] = Pp_[i][j]
        zParam.filterinfo = 1
        return zParam

    def valpos(self,rtkParam,zParam,thres):
        fact = thres*thres
        stat =1
        for i in range(zParam.nv):
            if zParam.v[i][0] * zParam.v[i][0] <= fact * zParam.R[i][i] : continue

            sat1 = (zParam.vflg[i][0] >> 16) & 0xff
            sat2 = (zParam.vflg[i][0] >> 8) & 0xff
            type = (zParam.vflg[i][0] >> 4) & 0xff

            freq = zParam.vflg[i][0] & 0xff
            stype = "L" if type == 0 else ("L" if type == 1 else "C")
            print("Large residual (sat",sat1,"- ",sat2," ",stype,freq+1," v= ",zParam.v[i][0],"sig = ",np.sqrt(zParam.R[i][i]) ," )")
        return stat

    def ddmat(self,rtkParam,D):
        return [0,D]
    def resamb_LAMBDA(self,rtkParam):
        nx = rtkParam.nx
        na = rtkParam.na
        s = [0,0]
        rtkParam.sol.ratio = 0

        if self.opt.mode <= 1:
            rtkParam.LAMBDAStat = 0
            return rtkParam
        D = tool.zeroMat(nx, nx)
        ddmatParam = self.ddmat(rtkParam,D)
        nb = ddmatParam[0]
        D = ddmatParam[1]
        if(nb <= 0):
            rtkParam.LAMBDAStat = 0
            return rtkParam
        ny = na + nb
        y = tool.zeroMat(ny,1)
        Qy = tool.zeroMat(ny,ny)
        DP = tool.zeroMat(ny,nx)
        b = tool.zeroMat(nb,2)
        db = tool.zeroMat(nb,1)
        Qb = tool.zeroMat(nb,nb)
        Qab = tool.zeroMat(na,nb)
        QQ = tool.zeroMat(na,nb)

        y = tool.mulmatirix(ny,1,nx,D,rtkParam.x,1,0,"TN",[])
        DP = tool.mulmatirix(ny,nx,nx,D,rtkParam.P,1,0,"TN",[])
        Qy = tool.mulmatirix(ny,ny,nx,DP,D,1,0,"NN",[])

        y_nb = []
        y_na = []
        for i in range(na):
            y_na.append(y[i][0])
        for i in range(na, ny):
            y_nb.append(y[i][0])
        for i in range(nb):
            for j in range(nb):
                Qb[i][j] = Qy[na+i][na+j]
        for i in range(na):
            for j in range(nb):
                Qab[i][j] = Qy[i][na+j]
        LAMBDAParam = self.LAMBDA(nb,2,y_nb,Qb,b,s)
        b = LAMBDAParam[1]
        s = LAMBDAParam[2]

        if LAMBDAParam[0] != 0:
            rtkParam.sol.ratio = s[0]/s[1] if s[0] > 0 else 0
            if rtkParam.sol.ratio > 999.9 : rtkParam.sol.ratio = 999.9

            if s[0] <= 0 or s[1]/s[0] >= self.opt.thresar[0] :
                for i in range(na):
                    rtkParam.xa[i][0] = rtkParam.x[i][0]
                    for j in range(na):
                        rtkParam.Pa[i][j] = rtkParam.P[i][j]
                for i in range(nb):
                    rtkParam.bias[i][0] = b[i][0]
                    y[na+i][0] -= b[i][0]
                Qb_1 = np.linalg.inv(Qb)

                db = tool.mulmatirix(nb,1,nb,Qb_1,y_nb,1,0,"NN",[])
                rtkParam.xa = tool.mulmatirix(na,1,nb,Qab,db,-1,1,"NN",rtkParam.xa)

                QQ = tool.mulmatirix(na,nb,nb,Qab,Qb,1,0,"NN",[])
                rtkParam.Pa = tool.mulmatirix(na,na,nb,QQ,Qab,-1,1,"NT",rtkParam.Pa)

                rtkParam = self.restamb(rtkParam)
            else:
                nb = 0
        return rtkParam

    def restamb(self,rtkParam):
        return rtkParam
    def LAMBDA(self,n,m,a,Q,b,s):

        print()



















