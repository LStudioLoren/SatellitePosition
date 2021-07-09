from position.singlepoint import *
from position.common import *
from position.rtkpoint import *
import time

if __name__ == '__main__':
    baseObs = OBSDATA.ObsData().initObsData2("D:\program\python\SPPositionWithGPSandGLO\data\\base.21O",2)
    roverObs = OBSDATA.ObsData().initObsData2("D:\program\python\SPPositionWithGPSandGLO\data\\rover.21O",1)
    navDataList = NAVDATA.NavData().initNavData("D:\program\python\SPPositionWithGPSandGLO\data\\nav.21N")
    nRover = len(roverObs)
    nBase = len(baseObs)
    rtkParam = RTKPOSITION.RTKPARAM()
    rtkParam.initRefPoint(-2333329.0760,5383350.5903,2492875.2940)
    RTKPOSITION.RTKPoistion().exeRTKPositon(rtkParam,roverObs[0],baseObs[0],navDataList)
    # for i in range(len(roverObs)):
    #     print("开始处理时间：", time.time())
    #     rtkParam.sol = SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(roverObs[i], navDataList,rtkParam.sol)
    #     print("GPS week = ",rtkParam.sol.gpsweek," GPS sec = ",rtkParam.sol.gpssec," POS = ",rtkParam.sol.X," NS = ",rtkParam.sol.ns,"AGE = ",rtkParam.sol.age," POS TYPE = ",rtkParam.sol.pos_type)
    #     print("结束处理时间：", time.time())
    # for i in range(len(baseObs)):
    #     print(baseObs[i].t_obs)
    #     for j in range(len(baseObs[i].obsary)):
    #         print("    ",baseObs[i].obsary[j].prn)
    #         for n in range(len(baseObs[i].obsary[j].obsfrefList)):
    #             print("        ",baseObs[i].obsary[j].obsfrefList[n].P)
    # for i in range(len(roverObs)):
    #     print(roverObs[i].t_obs)
    #     for j in range(len(roverObs[i].obsary)):
    #         print("    ",roverObs[i].obsary[j].prn)
    #         for n in range(len(roverObs[i].obsary[j].obsfrefList)):
    #             print("        ",roverObs[i].obsary[j].obsfrefList[n].P)
    # for i in range(len(navDataList)):
    #     print(navDataList[i].prn,navDataList[i].t_toc)