from position.singlepoint import *
from position.common import *


if __name__ == '__main__':
    baseObs = OBSDATA.ObsData().initObsData("D:\program\python\SPPositionWithGPSandGLO\data\\base.21O")
    roverObs = OBSDATA.ObsData().initObsData("D:\program\python\SPPositionWithGPSandGLO\data\\rover.21O")
    navDataList = NAVDATA.NavData().initNavData("D:\program\python\SPPositionWithGPSandGLO\data\\nav.21N")




    SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(baseObs, navDataList)
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