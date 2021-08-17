import position.common.tool as tool
import numpy as np
import math as m
#根据星历、时间、该卫星对应该时刻的L1伪距，计算该颗卫星的位置
class SatPos():
    def nav2pos(self,t_obs,navData):
        #navData.t_toc = navData.toe - tttt
        # print(navData.nav.t_toc)
        # 根据开普勒第三定律
        n0 = np.sqrt(tool.GM) / (navData.Sqrt_A * navData.Sqrt_A * navData.Sqrt_A)

        n = n0 + navData.Delta_n
        # print("参考时刻TOE的平均角速度n0 = ", n0)
        # print("计算卫星运动的平均角速度n = n0 + delta_n：", n,"  delte_n = ",Delta_n)
        #delta_t = navData.a0 + navData.a1 * (t_obs - navData.t_toc) + navData.a2 * (t_obs - navData.t_toc) * (
        #        t_obs - navData.t_toc)
        # print("delta_t = ",delta_t)
        tk = t_obs - navData.toe
        #print("tk = ",tk)
        Mk = navData.M0 + n * tk
        # print("信号发射时卫星的平近点角Mk：", Mk,"  M0  = ",M0)
        ed = Mk
        for i in range(4):
            # Ek = ed
            ed = ed - (ed - navData.e * np.sin(ed) - Mk) / (1 - navData.e * np.cos(ed))

        # print("ed = ",Ek -ed )
        Ek = ed

        Vk = m.atan2(np.sqrt(1 - navData.e * navData.e) * np.sin(Ek), (np.cos(Ek) - navData.e))
        # print("Vk = ",Vk)
        u = navData.omega + Vk
        COS2U = np.cos(2 * u)
        SIN2U = np.sin(2 * u)
        delta_u = navData.Cuc * COS2U + navData.Cus * SIN2U
        delta_r = navData.Crc * COS2U + navData.Crs * SIN2U
        delta_i = navData.Cic * COS2U + navData.Cis * SIN2U
        # print("delta_u = ",delta_u)
        # print("delta_r = ", delta_r)
        # print("delta_i = ", delta_i)
        uk = u + delta_u
        rk = (navData.Sqrt_A * navData.Sqrt_A) * (1 - navData.e * np.cos(Ek)) + delta_r
        ik = navData.i0 + delta_i + navData.IDOT * tk
        #print(uk,rk,ik)
        x = rk * np.cos(uk)
        y = rk * np.sin(uk)
        #print("x = ",x, "y = ", y)

        L = navData.OMEGA + (navData.OMEGA_DOT - tool.EARTH_RAD) * tk - tool.EARTH_RAD * navData.toe

        # print("L = ",L,"  ik =  ",ik)
        COSL = np.cos(L)
        SINL = np.sin(L)
        COSiK = np.cos(ik)
        SINiK = np.sin(ik)

        Xs = x * COSL - y * COSiK * SINL
        Ys = x * SINL + y * COSiK * COSL
        Zs = y * SINiK
        tk_dts = t_obs - navData.t_toc
        # dts = self.eph2clk(t_obs,navData_t_toc,navData)
        dts = navData.a0 + navData.a1 * tk_dts + navData.a2 * tk_dts * tk_dts
        # dts卫星的钟差，钟偏{bias,drift},
        # vare用户测距精度，从卫星星历钟的SAV中获取
        dts = dts - 2 * np.sqrt(tool.GM * navData.Sqrt_A * navData.Sqrt_A) * navData.e * np.sin(Ek) / (
                tool.CLIGHT * tool.CLIGHT)
        vare = navData.SAV * navData.SAV

        #print()
        #navData.updateParam(Xs, Ys, Zs, dts, vare)
        # navData.updataDtsAndVare()
        # navData.updataVare(vare)
        return [Xs, Ys, Zs, dts, vare]

    def getSatpos(self,navData, t_obs, OBS_P):
        # 根据观测量的伪距，除以光速，算出时差
        tttt = OBS_P / tool.CLIGHT
        t_sol = t_obs
        # print(tttt)
        t_obs = t_obs - tttt
        #计算广播星历的钟差
        dt = self.ephclk(t_obs,t_sol,navData)
        # 更新星历的toc时间,减去钟差
        t_obs -=dt

        navData2 = navData

        navData1 = self.nav2pos(t_obs,navData)
        t_obs += 1E-3
        navData2 = self.nav2pos(t_obs, navData)
        # += 0.1
        #navData2 = self.nav2pos(t_obs,navDataTemp)
        print("卫星：", navData.prn, "  卫星坐标： ", navData1[0],navData1[1],navData1[2]," 卫星钟差：", dt)
        navData.updateParam(navData1[0],navData1[1],navData1[2],((navData2[0] - navData1[0])/1E-3),((navData2[1] - navData1[1])/1E-3),((navData2[2] - navData1[2])/1E-3),navData1[3],navData1[4])
        return navData

    def eph2clk(self,ts_obs,navData):
        t = ts = ts_obs - navData.t_toc
        for i in range(3):
            t = ts - (navData.a0 + navData.a1 * t + navData.a2 * t * t)
        #print(ts_obs,navData.t_toc,t)
        return navData.a0 + navData.a1 * t + navData.a2 * t * t

    def seleph(self):
        return 0

    def ephclk(self,t_obs_PC,t_obs,navData):
        #寻找时间最接近的星历数据；
        self.seleph()
        #计算广播星历中的卫星种差（clock bias）
        dt = self.eph2clk(t_obs_PC,navData)
        #print(dt)
        return dt
