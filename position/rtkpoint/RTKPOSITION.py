from position.singlepoint import SINGLEPOINTPOSITION
from position.common import *
import numpy as np
class RTKPARAM():
    #基站坐标，ECEF-XYZ
    refPoint = []
    #rtk状态
    rtkStatus = 0
    #解状态
    sol = solustion.positionsol()

    def initRefPoint(self,X,Y,Z):
        self.refPoint = (X,Y,Z)
    def updateSol(self,sol):
        self.sol = sol


class RTKPoistion():
    def exeRTKPositon(self,rtkParam,roverObs,baseObs,navDataList):
        #先计算单点位置，获取首次有效定位、速度
        rtkParam.sol = SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(roverObs, navDataList, rtkParam.sol)

        #baseNavList = []
        #roverNavList = []
        if len(baseObs.obsary) < 4:
            rtkParam.sol.update(rtkParam.sol.X, len(baseObs.obsary), 0, 0)
            return rtkParam.sol
        if len(roverObs.obsary) < 4:
            rtkParam.update(rtkParam.sol.X, len(roverObs.obsary), 0, 0)
            return rtkParam.sol

        #根据基站位置计算基于基站的卫星位置信息等
        baseNavList = self.getNavList(baseObs,navDataList)
        nbaseObs = len(baseObs.obsary)
        nroverObs = len(roverObs.obsary)

        #roverNavList = self.getNavList(roverObs,navDataList)

        #进行基站的非差残差计算
        if self.zdres(0,baseObs,len(baseObs.obsary),baseNavList,rtkParam.refPoint,1) != 1:
            return rtkParam.sol
        #获取移动站和基站的公共卫星号、对应各自obs数组的位置。
        self.selectCommonSat(roverObs,baseObs,nroverObs,nbaseObs)




    def getNavList(self,obs,navDataList):
        navList = []
        for k in range(len(obs.obsary)):
            for i in range(len(navDataList)):
                # 计算卫星位置；
                if obs.obsary[k].prn == navDataList[i].prn:
                    # 卫星观测量对应的L1伪距
                    # 整理星历数据打包程nav对象，计算卫星位置pos，计算卫星种差dts，计算卫星位置及时钟方差vare
                    print("卫星:",navDataList[i].prn,"位置：")
                    navList.append(SATPOS.SatPos().getSatpos(navDataList[i], obs.t_obs,
                                                                 obs.obsary[k].obsfrefList[0].P))
                else:
                    continue
        return navList
    def zdres(self,type,obs,nobs,navList,rr,index):

        #for i in range(len(obs))
        nfre = 2
        y=[]
        dant=[]
        for i in range(nobs):
            y_obs = []
            dant.append(0)
            for j in range(nfre*2):
                y_obs.append(0)
            y.append(y_obs)
        print(y[0],y[0][1])
        if tool.dot(rr) <=0 :
            return 0
        rr_ = rr
        pos = [0.0,0.0,0.0]
        pos = RTKCOMMON.eceftopos(rr_,pos)
        #L1+L2双频计算
        lam = [tool.CLIGHT / tool.FREQ1, tool.CLIGHT / tool.FREQ2]
        for i in range(nobs):
            r = RTKCOMMON.r_corr(navList[i],rr_)
            #print(r)
            if r <= 0 :
                continue

            e = RTKCOMMON.e_corr(navList[i],rr_)
            azel = RTKCOMMON.satAzel(pos,e)
            if azel[1] < 10*np.pi/180 :
                continue

            r += (-tool.CLIGHT)* navList[i].dts
            print(r)
            zhd = self.tropmodel(obs.t_obs,pos,[0,90*np.pi/180],0)
            tropnmf = self.tropmapf_nmf(obs.gpsweek,obs.t_obs,pos,azel,0)
            r += tropnmf*zhd
            print(r,tropnmf,zhd)
            for j in range(len(obs.obsary[i].obsfrefList)):
                #for k in range():
                if (obs.obsary[i].obsfrefList[j].L > 0.0) :
                    y[i][j] = obs.obsary[i].obsfrefList[j].L*lam[j] - r - dant[i]
                if (obs.obsary[i].obsfrefList[j].P > 0.0):
                    y[i][j+nfre] = obs.obsary[i].obsfrefList[j].P -r - dant[i]

        return 1
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

    def selectCommonSat(self,roverObs,baseObs,nroverObs,nbaseObs):
        satList=[]
        irover = []
        ibase = []
        for i in range(nroverObs):
            for j in range(nbaseObs):
                if roverObs.obsary[i].prn == baseObs.obsary[j].prn:
                    satList.append(roverObs.obsary[i].prn)
                    irover.append(i)
                    ibase.append(j)

        print(satList,irover,ibase)







