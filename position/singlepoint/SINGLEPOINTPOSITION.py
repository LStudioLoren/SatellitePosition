from position.common import *
import numpy as np

#单点定位计算
'''
单点定位的原理：ts卫星时间，tu接收机时间，dts是卫星时间与gps时差，dtu接收机与gps时差,P是伪距
1）P = C*(ts - tu)+C*(dts-dtu)
2）P(s2u) = C*（ts-tu）是卫星与接收机的几何距离
3）真实伪距P = P(s2u)+C（dts-dtu）+I(t)+T(t)
4）P（stu）(t) = sqrt（（Xs- Xu）^2+(Ys-Yu)^2+(Zs -Zu)^2），(Xs,Ys,Zs)是卫星的位置，（Xu，Yu，Zu）是接收机的位置
5）进行泰勒级数展开并取一次项后：
P（stu）(t) = (P（stu）(t))0 + (-1*(Xs-Xu)/(P（stu）(t))0 )*dX+(-1*(Ys-Yu)/(P（stu）(t))0 )*dY+(-1*(Zs-Zu)/(P（stu）(t))0 )*dZ
令(-1*(Xs-Xu)/(P（stu）(t))0 ) = -k
(-1*(Ys-Yu)/(P（stu）(t))0 ) = -l
(-1*(Zs-Zu)/(P（stu）(t))0 ) = -m
C（dts-dtu） = CdT
(P（stu）(t))0  = r   r是卫星与接收机的几何距离
6)P = (P（stu）(t))0  -k(t)*dX-l(t)*dY-m(t)*dZ +C*dT +I(t)+T(t)
...
单系统情况下，有几个卫星就有几个上述方程
7)进行组装矩阵：
P = [P1,P2,P3...Pn]
r = [r1+I(t)+T(t),r2+I(t)+T(t),r3+I(t)+T(t)...rn+I(t)+T(t)]

e = [-k1,-l1,-m1,C
     -k2,-l2,-m2,C
     ...
     -kn,-ln,-mn,C]
未知参数DX = [dX,dY,dZ,dT]
最后就变成P = r + e*DX
8)根据最小二乘法
DX = (eT*e)^-1*eT（P-r）
     
'''
#LSP参数类
class SPPParm():
    #参与解算卫星数据
    vn = 0
    #几何矩阵
    e_matrix = []
    #伪距矩阵
    V = []
    #加权矩阵
    Var = []

    def init(self,vn,e,var,v):
        self.vn = vn
        self.e_matrix = e
        self.Var = var
        self.V = v

#单点定位算法类
class SinglePointPosition():







    #
    '''
    对接收机位置进行估算
    流程：
    1）先将假设接收机位置未X[dx,dy,dz,dt]
    2）将X赋值给rr，rr记录的是上一次计算出来的X
    3）钟差dtr =X[3]
    4)将rr从ecef坐标系转换为大地坐标系pos（B、L、H）
    5）将同一时刻的星历数据、观测量数据、ecef坐标的接收机位置、BLH坐标的接收机位置，误差参数、钟差传入rescode方法，
        进行伪距残差、几何矩阵、伪距矩阵、参与解算卫星数等计算，并加载到sppparm中；
    6）对几何矩阵e,伪距矩阵V，进行加权，权重为Var
    7）进行LSP计算，得出dX[dx,dy,dz,dt],并且X+=dX
    8）将LSP输出的X进行统计dot(dX)，如果小于1E-4，则返回X为最后结果，否则重复2-7
    '''
    def estpos(self,nav_list, OBS_P, sol):
        # V是伪距残差数列
        V = []
        # Var是误差方差（与vare关联）
        Var = []
        for i in range(len(nav_list)):
            V.append(0)
            Var.append(0)
        # 误差参数
        err = [100, 0.003, 0.003]

        rr = [0.0, 0.0, 0.0]
        # X是x,y,z,dtr位置
        X = sol.X
        # dx是x,y,z误差
        dx = []

        count = 0
        while (1):
            pos = [0, 0, 0]
            for i in range(3):
                # rr是上一次计算出来的接收机的x、y、z
                rr[i] = X[i]
            # dtr是
            dtr = X[3]
            # rr是x，Y,z，pos是blh的位置
            pos = RTKCOMMON.eceftopos(rr, pos)
            spparm = self.rescode(nav_list, OBS_P, rr, pos, err, dtr, V, Var)
            # print("e = ",e_matrix)
            n = spparm.vn
            m = 4
            # print(e_matrix)
            # print(V)
            for i in range(n):
                # 对矩阵进行加权，权重值为每个卫星对应的Var。
                sig = np.sqrt(spparm.Var[i])
                spparm.V[i] /= sig
                for j in range(m):
                    spparm.e_matrix[i][j] /= sig
            # print("e= ", e_matrix)
            # print("V=", V)
            # print(Var)
            n = 4
            k = 4
            m = spparm.vn
            # r = [satpos]*dx[dx,dy,dz,dtr]
            dx = tool.LSP(n, k, m, spparm.e_matrix, spparm.V, dx)
            for i in range(4):
                X[i] += dx[i]
            # print("X=", X)
            # print("dx=", dx)
            count += 1
            if np.sqrt(tool.dot(dx) < 1E-4):
                print("LSP处理次数：", count)
                break
        sol.update(X,spparm.vn,0,1)
        return sol
    #计算出伪距矩阵V、接收机相对于卫星位置的矩阵H（乘以-1）、解算卫星数vn、卫星权重矩阵Var
    def rescode(self,nav_list, OBS_P, rr, pos, err, dtr, V, Var):
        #vn是参与解算卫星数
        vn = 0
        e_matrix = []
        for i in range(len(nav_list)):
            # 伪距修正，根据卫星位置，减去估算的用户位置
            r = RTKCOMMON.r_corr(nav_list[i], rr)  # [nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z]
            # 计算接收机相对于卫星位置的向量数组
            e = RTKCOMMON.e_corr(nav_list[i], rr)
            # 计算卫星的方位角和截止角

            azel = RTKCOMMON.satAzel(pos, e)

            # 小于10度截止的的就不参与解算。
            if azel[1] < 0.173:
                continue
            #得到伪距矩阵V
            #V[vn] = OBS_P[i] - (r + dtr - tool.CLIGHT * nav_list[i].dts)
            #dion电离层误差，dtrp对流层误差，未考虑两者的误差
            dion = 0
            dtrp = 0
            V[vn] = OBS_P[i] - (r + dtr - tool.CLIGHT * nav_list[i].dts + dion + dtrp)
            # 将e数组变乘以-1，主要是用于后续最小二乘法中去。并合并到大数组中，得到LSP需要的矩阵H
            for j in range(3):
                e[j] *= -1
            e_matrix.append(e)

            # Var 计算伪距权重。vion是电离层的方差，vtrp是对流层的方差
            vion = 0.09
            vtrp = 0.07438
            Var[vn] = 1 * err[0] * err[0] * (err[1] * err[1] + err[2] * err[2]) / np.sin(azel[1]) + nav_list[
                i].vare + vion + vtrp
            # print("r = ", r, "e = ",              e)  # nav_list[i].nav.prn,nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z,nav_list[i].dts,nav_list[i].vare,V[i],r)
            vn += 1
            s = SPPParm()
            s.init(vn, e_matrix, Var, V)
        return s

    def exesinglepoint(self,OBS_DATA,nav_data_list,sol):
        OBS_P = []
        nav_list = []
        if len(OBS_DATA.obsary) <4 :
            sol.update(sol.X,len(OBS_DATA.obsary),0,0)
            return  sol
        for k in range(len(OBS_DATA.obsary)):
            for i in range(len(nav_data_list)):
                # 计算卫星位置；
                if OBS_DATA.obsary[k].prn == nav_data_list[i].prn:
                    # 卫星观测量对应的L1伪距
                    # 整理星历数据打包程nav对象，计算卫星位置pos，计算卫星种差dts，计算卫星位置及时钟方差vare
                    nav_list.append(SATPOS.SatPos().getSatpos(nav_data_list[i], OBS_DATA.t_obs,
                                                              OBS_DATA.obsary[k].obsfrefList[0].P))
                else:
                    continue

                #     print("prn :", navData.nav.prn, " x = ", navData.nav.x, " y=", navData.nav.y, " z= ", navData.nav.z, " r = ",
                #       np.sqrt(dot([navData.nav.x, navData.nav.y, navData.nav.z])))
                # print("t_obs = ",OBS_ALL[j].t_obs," prn = ",OBS_ALL[j].obsary[k].prn," p1 = ", OBS_ALL[j].obsary[k].p1)
            OBS_P.append(OBS_DATA.obsary[k].obsfrefList[0].P)

        sol.init(OBS_DATA.gpsweek,OBS_DATA.t_obs)
        sol = self.estpos(nav_list, OBS_P,sol)
        #print("week=",sol.gpsweek,"  time = ", sol.gpssec, "  X = ", sol.X," ns = ",sol.ns,"  age = ",sol.age)

        return sol
