import time
from position.singlepoint import *
from position.common import *
if __name__ == '__main__':
    #读取广播星历数据
    nav_data_list = NAVDATA.NavData().initNavData("D:\program\python\python-project\project1\data\\nav_gps.21N")
    #读取观测量数据
    OBS_DATA = OBSDATA.ObsData().initObsData("D:\program\python\python-project\project1\data\\Obs_gps.21O")
    sol = solustion.positionsol()
    for i in range(len(OBS_DATA)):
        print("开始处理时间：", time.time())
        sol = SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(OBS_DATA[i], nav_data_list,sol)
        print("GPS week = ", sol.gpsweek, " GPS sec = ", sol.gpssec, " POS = ", sol.X,
              " NS = ", sol.ns, "AGE = ", sol.age, " POS TYPE = ", sol.pos_type)
        print("结束处理时间：", time.time())





