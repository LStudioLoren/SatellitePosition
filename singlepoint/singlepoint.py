import time

import NAVDATA
import OBSDATA
import SATPOS
import SINGLEPOINTPOSITION as spp

if __name__ == '__main__':
    #读取广播星历数据
    nav_data_list = NAVDATA.NavData().initNavData("D:\program\python\python-project\project1\data\\nav_gps.21N")
    #读取观测量数据
    OBS_DATA = OBSDATA.ObsData().initObsData("D:\program\python\python-project\project1\data\\Obs_gps.21O")

    X = [0.0, 0.0, 0.0, 0.0]
    for j in range(len(OBS_DATA)):
        OBS_P = []
        nav_list = []
        for k in range(len(OBS_DATA[j].obsary)):
            for i in range(len(nav_data_list)):
                # 计算卫星位置；
                if OBS_DATA[j].obsary[k].prn == nav_data_list[i].prn :
                    # 卫星观测量对应的L1伪距
                    # 整理星历数据打包程nav对象，计算卫星位置pos，计算卫星种差dts，计算卫星位置及时钟方差vare
                    nav_list.append(SATPOS.SatPos().getSatpos(nav_data_list[i],OBS_DATA[j].t_obs,OBS_DATA[j].obsary[k].p1))
                else:
                    continue

            #     print("prn :", navData.nav.prn, " x = ", navData.nav.x, " y=", navData.nav.y, " z= ", navData.nav.z, " r = ",
            #       np.sqrt(dot([navData.nav.x, navData.nav.y, navData.nav.z])))
            # print("t_obs = ",OBS_ALL[j].t_obs," prn = ",OBS_ALL[j].obsary[k].prn," p1 = ", OBS_ALL[j].obsary[k].p1)
            OBS_P.append(OBS_DATA[j].obsary[k].p1)
        print("开始处理时间：",time.time())
        X = spp.SinglePointPosition().estpos(nav_list, OBS_P,X)

        print("time = ",OBS_DATA[j].t_obs,"  X = ",X)
        print("结束处理时间：", time.time())



