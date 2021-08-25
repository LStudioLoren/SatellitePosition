from position.common import *
import numpy as np
# 单点定位计算
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


# LSP参数类
class SPPParm():
    def __init__(self):
        # 参与解算卫星数据
        self.vn = 0
        # 几何矩阵
        self.e_matrix = []
        # 伪距修正矩阵
        self.V = []
        # 加权矩阵
        self.Var = []
        # 卫星方位、截止角
        self.azel = []
        # 可用卫星
        self.vsat = []
        # 伪距残差
        self.resp = []

        self.sol = solustion.positionsol()

    def init(self, vn, e, var, v,azelList ,vsatList,respList):
        self.vn = vn
        self.e_matrix = e
        self.Var = var
        self.V = v
        self.azel = azelList
        self.vsat = vsatList
        self.resp = respList

    def reV(self):
        Vl = []

        for i in range(len(self.V)):
            Vl.append([self.V[i]])
        self.V = Vl
        return self


# 单点定位算法类
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

    def estpos(self, nav_list, OBS_DATA, sol):
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
        print("天线位置：X = ", X[0], " Y = ", X[1], " Z = ", X[2])
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
            #reV是将V数组从一维数组，重整为n行1列的二维数组
            spparm = self.rescode(nav_list, OBS_DATA, rr, pos, err, dtr, V, Var,sol.gpssec).reV()
            # print("e = ",e_matrix)
            n = spparm.vn
            m = 4
            # print(e_matrix)
            # print(V)
            for i in range(n):
                # 对矩阵进行加权，权重值为每个卫星对应的Var。
                sig = np.sqrt(spparm.Var[i])
                spparm.V[i][0] /= sig
                for j in range(m):
                    spparm.e_matrix[i][j] /= sig
            # print("e= ", e_matrix)
            # print("V=", V)
            # print(Var)
            n = 4
            k = 4
            m = spparm.vn
            #print("H = ",spparm.e_matrix)
            #print("V = ",spparm.V)
            #n列，m行矩阵，乘以K列M行
            dx = tool.LSQ(n, k, m, spparm.e_matrix, spparm.V, dx)
            for i in range(4):
                X[i] += dx[i]

            count += 1
            print("迭代次数：",count)
            print("天线位置：X = ", X[0], " Y = ", X[1], " Z = ", X[2]," dt = ",X[3])
            if np.sqrt(tool.dot(dx) < 1E-4):
                #print("处理次数：", count)
                break
        sol.update(X, spparm.vn, 0, 1)
        spparm.sol = sol
        #rtkParam.updateSol(sol)

        return spparm

    # 计算出伪距矩阵V、接收机相对于卫星位置的矩阵H（乘以-1）、解算卫星数vn、卫星权重矩阵Var
    def rescode(self, nav_list, OBS_DATA, rr, pos, err, dtr, V, Var,t_obs):
        # vn是参与解算卫星数
        vn = 0
        e_matrix = []
        vsatList = []
        azelList = []
        respList = []
        for i in range(len(nav_list)):
            vsatList.append(0)
            azelList.append([0.0,0.0])
            respList.append(0)
            # 伪距修正，根据卫星位置，减去估算的用户位置
            r = RTKCOMMON.r_corr(nav_list[i], rr)  # [nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z]
            #print("before r",r,"  rr = ",rr," RS = ",nav_list[i].x,nav_list[i].y,nav_list[i].z)
            # 计算接收机相对于卫星位置的向量数组
            e = RTKCOMMON.e_corr(nav_list[i], rr)
            # 计算卫星的方位角和截止角

            azelList[i] = RTKCOMMON.satAzel(pos, e)

            # 小于10度截止的的就不参与解算。
            if azelList[i][1] < 0.173:
                continue
            # 得到伪距矩阵V
            # V[vn] = OBS_P[i] - (r + dtr - tool.CLIGHT * nav_list[i].dts)
            # dion电离层误差，dtrp对流层误差，未考虑两者的误差
            #        *ion=ionmodel(time,nav->ion_gps,pos,azel);
            #*var=SQR(*ion*ERR_BRDCI);
            #return 1;
            ionCorr = self.ionocorr(t_obs,pos,azelList[i],1)
            tropCorr = self.tropcorr(t_obs,pos,azelList[i],1)
            pCorr = self.prange(OBS_DATA,nav_list,azelList[i],i)
            #dtrp = 0
            print("卫星号 = ",OBS_DATA.obsary[i].prn,"，接收机输出的伪距 =",pCorr[0], "  r=  ",r,"  "," 接收机钟差 = ",dtr," 卫星钟差= ",nav_list[i].dts," 电离层误差= ",ionCorr[0]," 对流层误差= ",tropCorr[0])
            V[vn] = pCorr[0] - (r + dtr - tool.CLIGHT * nav_list[i].dts + ionCorr[0] + tropCorr[0])
            # 将e数组变乘以-1，主要是用于后续最小二乘法中去。并合并到大数组中，得到LSP需要的矩阵H
            for j in range(3):
                e[j] *= -1
            e_matrix.append(e)

            vsatList[i] = 1
            respList[i] = V[vn]

            # Var 计算伪距权重。vion是电离层的方差，vtrp是对流层的方差
            #vion = 0.09
            #vtrp = 0.07438
            Var[vn] = 1 * err[0] * err[0] * (err[1] * err[1] + err[2] * err[2]/ np.sin(azelList[i][1])) + nav_list[
                i].vare +pCorr[1] + ionCorr[1] + tropCorr[1]
            # print("r = ", r, "e = ",              e)  # nav_list[i].nav.prn,nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z,nav_list[i].dts,nav_list[i].vare,V[i],r)
            vn += 1
            s = SPPParm()
            s.init(vn, e_matrix, Var, V,azelList,vsatList,respList)
        #print("V[nv]",V)
        return s
    def prange(self,obs,nav,azel,iter):
        P1 = obs.obsary[iter].obsfrefList[0].P
        #P2 = obs.obsary[iter].obsfrefList[1].P
        #lam = [tool.CLIGHT / tool.FREQ1, tool.CLIGHT / tool.FREQ2]
        #gamma
        var = 0.3*0.3

        b1= nav[iter].TGD *tool.CLIGHT #/ * TGD (m) * /
        return [P1 - b1,var]

    def tropcorr(self,t_obs,pos,azel,troptype):

        trop =0
        var = 0
        if troptype == 1:
            temp0 = 15.0 #/ *temparature at sea level * /
            humi = 0.7


            if (pos[2] < -100.0 or 1E4 < pos[2] or azel[1] <= 0) :
                trop = 0
            else:
                #/ *standard atmosphere * /
                hgt = 0.0 if pos[2] < 0.0 else pos[2]

                pres = 1013.25 * np.power(1.0 - 2.2557E-5 * hgt, 5.2568)
                temp = temp0 - 6.5E-3 * hgt + 273.16;
                e = 6.108 * humi * np.exp((17.15 * temp - 4684.0) / (temp - 38.45))

                #/ *saastamoninen model * /
                z = np.pi / 2.0 - azel[1]
                trph = 0.0022768 * pres / (1.0 - 0.00266 * np.cos(2.0 * pos[0]) - 0.00028 * hgt / 1E3) / np.cos(z);
                trpw = 0.002277 * (1255.0 / temp + 0.05) * e / np.cos(z);
                trop = trph + trpw
        var = (0.3 / (np.sin(azel[1]) +0.1))*(0.3 / (np.sin(azel[1]) +0.1))
        return [trop,var]
    def ionocorr(self,t_obs,pos,azel,iontype):
        ion = 0
        var = 0
        if iontype == 1:
            #2004/1/1
            ion2 = [0,0,0,0,0,0,0,0]
            ion_default= [0.1118E-07, -0.7451E-08, -0.5961E-07, 0.1192E-06,
                            0.1167E+06, -0.2294E+06, -0.1311E+06, 0.1049E+07]
            tt = f = psi = phi = lam = amp = per = x = week = 0

            if (pos[2] < -1E3 or azel[1] <= 0):
                ion = 0
            else:
                if (tool.norm(ion2, 8) <= 0.0):
                    ion2=ion_default

                #earth centered angle(semi - circle) * /
                psi = 0.0137 / (azel[1] / np.pi + 0.11) - 0.022

                #/ *subionospheric latitude / longitude(semi - circle) * /
                phi = pos[0] / np.pi + psi * np.cos(azel[0])
                if (phi > 0.416) :
                    phi = 0.416
                elif (phi < -0.416):
                    phi = -0.416

                lam = pos[1] / np.pi + psi * np.sin(azel[0]) / np.cos(phi * np.pi)

                #/ *geomagnetic latitude(semi - circle) * /
                phi += 0.064 * np.cos((lam - 1.617) * np.pi)
                #/ *local time(s) * /
                tt = 43200.0 * lam + t_obs
                tt -= np.floor(tt / 86400.0) * 86400.0 # / *0 <= tt < 86400 * /

                #/ *slant factor * /
                f = 1.0 + 16.0 * pow(0.53 - azel[1] / np.pi, 3.0)

                #/ *ionospheric delay * /
                amp = ion2[0] + phi * (ion2[1] + phi * (ion2[2] + phi * ion2[3]));
                per = ion2[4] + phi * (ion2[5] + phi * (ion2[6] + phi * ion2[7]));
                amp = 0.0 if amp < 0.0 else  amp
                per = 72000.0 if per < 72000.0 else  per
                x = 2.0 * np.pi * (tt - 50400.0) / per
                ion = tool.CLIGHT * f * ( (5E-9+amp * (1.0+x * x * (-0.5+x * x / 24.0))) if np.fabs(x) < 1.57 else 5E-9)
        var = ion *0.5 *ion*0.5
        #print("[ion,var]",ion,var)
        return [ion,var]
    def resdop(self,obs,nav,sol,vsat,azel,x):
        vn = 0
        H = []
        v = []
        a = [[0.0],[0.0],[0.0]]
        vs = [0,0,0]
        pos = [0.0,0.0,0.0]
        E9 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        pos = RTKCOMMON.eceftopos( sol.rr,pos)
        E9 = RTKCOMMON.xyz2enu(pos,E9)
        E = []
        for i in range(3):
            Ei = []
            for j in range(3):
                Ei.append(E9[j+i*3])
            E.append(Ei)

        for i in range(len(obs.obsary)):
            lam = tool.lam[0]
            if obs.obsary[i].obsfrefList[0].D == 0.0 or lam == 0.0 or vsat[i] == 0 or tool.norm(nav[i].getVel(),3) <= 0.0 :
                continue

            cosel = np.cos(azel[i][1])
            a[0][0] = np.sin(azel[i][0])*cosel
            a[1][0] = np.cos(azel[i][0])*cosel
            a[2][0] = np.sin(azel[i][1])

            e = tool.mulmatirix(3,1,3,E,a,1,0,"NN",[])

            for j in range(3):
                vs[j] = nav[i].getVel()[j] -x[j]
            rs = nav[i].getRS()
            rate = tool.dot_n(vs,[e[0][0],e[1][0],e[2][0]],3)+(tool.OMGE / tool.CLIGHT) * (rs[4]*sol.rr[0] +rs[1]*x[0] -rs[3]*sol.rr[1] - rs[0]*x[1])
            vi = -lam * obs.obsary[i].obsfrefList[0].D - (rate + x[3] -tool.CLIGHT*0)#dts[1],钟飘==0)
            v.append([vi])
            Hi = []
            for j in range(4):
                Hi.append(-e[j][0] if j < 3 else 1.0)
            H.append(Hi)
            vn += 1

        #print(vn)
        #print("vel H:",H)
        #print(("vel v ",v))
        return [vn,H,v]
    def estvel(self,obs,nav,sol,vsat,azel):
        x = [0.0,0.0,0.0,0.0]
        for i in range(10):
            velparam = self.resdop(obs,nav,sol,vsat,azel,x)
            if velparam[0] <4 :
                break

            dx = tool.LSQ(4,4,velparam[0],velparam[1],velparam[2],[])
            for j in range(4):
                x[j] += dx[j]

            if(tool.norm(dx,4) < 1E-6):
                for j in range(3):
                    sol.rr[3+j] = x[j]
                break
        return sol


    def exesinglepoint(self, OBS_DATA, nav_data_list, rtkParam):
        sol = rtkParam.sol
        OBS_P = []
        nav_list = []
        if len(OBS_DATA.obsary) < 4:
            sol.update(sol.X, len(OBS_DATA.obsary), 0, 0)
            return sol
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
            #OBS_P.append(OBS_DATA.obsary[k].obsfrefList[0].P)
        #sppparm = SPPParm()
        sol.init(OBS_DATA.gpsweek, OBS_DATA.t_obs)
        #rtkParam.updateSol(sol)
        #输入卫星坐标、接收机的伪距观测值，计算天线位置
        sppparm = self.estpos(nav_list, OBS_DATA, sol)

        self.estvel(OBS_DATA,nav_list,sol,sppparm.vsat,sppparm.azel)
        #初始化rtkparam中的ssatlist数据，以最大卫星数创建对应每个卫星的参数（移动站的）。
        for i in range(tool.MAXSAT):
            rtkParam.ssatList[i].vs = 0
            rtkParam.ssatList[i].azel = [0.0,0.0]
            rtkParam.ssatList[i].resp = [0.0,0.0,0.0]
            rtkParam.ssatList[i].resc = [0.0,0.0,0.0]
            rtkParam.ssatList[i].snr = [0,0,0]
        for i in range(len(OBS_DATA.obsary)):
            rtkParam.ssatList[OBS_DATA.obsary[i].prn-1].azel = sppparm.azel[i]
            rtkParam.ssatList[OBS_DATA.obsary[i].prn - 1].snr[0] =OBS_DATA.obsary[i].obsfrefList[0].S
            if sppparm.vsat[i] != 1 :
                continue
            rtkParam.ssatList[OBS_DATA.obsary[i].prn - 1].vs = 1
            rtkParam.ssatList[OBS_DATA.obsary[i].prn - 1].resp[0] = sppparm.resp[i]
        # print("week=",sol.gpsweek,"  time = ", sol.gpssec, "  X = ", sol.X," ns = ",sol.ns,"  age = ",sol.age)
        rtkParam.updateSol(sppparm.sol)
        return rtkParam
