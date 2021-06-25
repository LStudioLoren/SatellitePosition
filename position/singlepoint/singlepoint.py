
from position.singlepoint import *
from position.common import *
if __name__ == '__main__':
    #读取广播星历数据
    nav_data_list = NAVDATA.NavData().initNavData("D:\program\python\python-project\project1\data\\nav_gps.21N")
    #读取观测量数据
    OBS_DATA = OBSDATA.ObsData().initObsData("D:\program\python\python-project\project1\data\\Obs_gps.21O")
    SINGLEPOINTPOSITION.SinglePointPosition().exesinglepoint(OBS_DATA, nav_data_list)





