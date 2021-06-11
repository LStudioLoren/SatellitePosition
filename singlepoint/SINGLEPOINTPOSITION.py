import tool
import numpy as np
#单点定位计算
class SinglePointPosition():
    # 计算卫星的对应测量点的方位角；这是固定算法，需要将
    def satAzel(self,pos, e):
        # print("satAzel e = ",e)
        az = 0
        el = np.pi / 2
        # enu = []
        if pos[2] > -tool.RE_WGS84:
            enu = self.xyz2enu(pos, e)
            # print("enu = ",enu)
            if enu[0] * enu[0] + enu[1] * enu[1] < 1E-12:
                az = 0
            else:
                az = np.arctan2(enu[0], enu[1])
            if az < 0.0:
                az += 2 * np.pi
            el = np.arcsin(enu[2])
            # print("az = ",az,"el",el,enu[0],enu[1],enu[2])
        return [az, el]
    #将卫星的x、y、z坐标系转换为enu坐标系，公式是固定的；
    def xyz2enu(self,pos, e):
        enu = []
        E9 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        # xyz2enu
        SINP = np.sin(pos[0])
        COSP = np.cos(pos[0])
        SINL = np.sin(pos[1])
        COSL = np.cos(pos[1])
        E9[0] = -SINL
        E9[1] = -SINP * COSL
        E9[2] = COSP * COSL
        E9[3] = COSL
        E9[4] = -SINP * SINL
        E9[5] = COSP * SINL
        E9[6] = 0.0
        E9[7] = COSP
        E9[8] = SINP

        n = 3
        k = 1
        m = 3

        for i in range(n):
            d = 0
            for j in range(k):
                for x in range(m):
                    # print("E9[",x+i*n,"] = ",E9[x+i*n],"e[",x,"] = ",e[x])
                    d += E9[i + x * n] * e[x]  # E1*e[x]/(dot(e_x,e,y,e_z)+E2*e[y]/dot+E3*e[
                # print(d)
                enu.append(d)
        # print(enu)
        enu[0] /= np.sqrt(tool.dot(e))
        enu[1] /= np.sqrt(tool.dot(e))
        enu[2] /= np.sqrt(tool.dot(e))
        # print(enu)
        # print("asdfa enu =  [0.9832349746631751, -0.0017021136149081018, 0.18233509647993557]")
        return enu

    # pos是卫星的位置x，y，z
    def eceftopos(self,rr, pos):
        # print("eceftopos: rr = ",rr)
        e2 = tool.FE_WGS84 * (2 - tool.FE_WGS84)
        v = tool.RE_WGS84
        R2 = tool.dot2(rr)
        zk = 0
        z = rr[2]
        # print(R2)
        while (1):

            if np.fabs(z - zk) >= 1E-4:
                zk = z
                sinp = z / np.sqrt(R2 + z * z)
                v = tool.RE_WGS84 / np.sqrt(1 - e2 * sinp * sinp)
                z = rr[2] + v * e2 * sinp
                # print("zk = ",zk," z = ",z," np.fabs(z-zk) = ",np.fabs(z-zk))
            else:
                break;

        # pos1 = pos2 = pos3 = 0

        if R2 > 1E-12:
            # print(np.sqrt(R2),z)
            pos[0] = np.arctan(z / np.sqrt(R2))
        else:
            if R2 > 0:
                pos[0] = np.pi / 2
            else:
                pos[0] = -np.pi / 2

        if R2 > 1E-12:
            pos[1] = np.arctan2(rr[1], rr[0])
        else:
            pos[1] = 0

        pos[2] = np.sqrt(R2 + z * z) - v

        # print(pos)
        return pos

    # 计算伪距残差
    def r_corr(self,nav, rr):

        r = np.sqrt(tool.dot([nav.x - rr[0], nav.y - rr[1], nav.z - rr[2]]))
        # 再加上
        r += tool.OMGE * (nav.x * rr[1] - nav.y * rr[0]) / tool.CLIGHT
        return r

    #nav是星历数据，rr是当前接收机位置
    def e_corr(self,nav, rr):
        e = [0, 0, 0, 1]
        #将卫星的x、y、z分别减去接收机的位置，得出接收机相对于卫星的x、y、z方向距离，即e数组
        e[0] = nav.x - rr[0]
        e[1] = nav.y - rr[1]
        e[2] = nav.z - rr[2]
        # r是e数组的距离
        r = np.sqrt(tool.dot(e))
        # print("r = ", r, "  e_before= ", e)
        for i in range(3):
            # 计算接收机到卫星的向量vector
            # 这里是重点：将e数组中x、y、z除以距离r，e数组重新赋值为接收机相对于卫星的x、y、z的向量值。
            e[i] = e[i] / r
        # print("e_after= ", e)
        return e

    def estpos(self,nav_list, OBS_P, last_X):
        i = 0
        # V是伪距残差数列
        V = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        # 误差参数
        err = [100, 0.003, 0.003]

        rr = [0.0, 0.0, 0.0]
        # X是x,y,z,dtr位置
        X = last_X
        # dx是x,y,z误差
        dx = []
        # Var是误差方差（与vare关联）
        Var = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        count = 0
        while (1):
            # e矩阵
            e_matrix = []
            pos = [0, 0, 0]
            for i in range(3):
                # rr是上一次计算出来的接收机的x、y、z
                rr[i] = X[i]
            # dtr是
            dtr = X[3]
            # rr是x，Y,z，pos是ecef的位置
            pos = self.eceftopos(rr, pos)


            rescode_data = self.rescode(nav_list, OBS_P, rr, pos, err, dtr, V, Var)
            # vn可用卫星数
            vn = rescode_data[0]
            #
            e_matrix = rescode_data[1]
            Var = rescode_data[2]
            V = rescode_data[3]
            # print("e = ",e_matrix)
            i = j = 0
            n = vn
            m = 4
            # print(e_matrix)
            # print(V)
            for i in range(n):
                # 对矩阵进行加权，权重值为每个卫星对应的Var。
                sig = np.sqrt(Var[i])
                V[i] /= sig
                for j in range(m):
                    e_matrix[i][j] /= sig
            # print("e= ", e_matrix)
            # print("V=", V)
            # print(Var)
            n = 4
            k = 4
            m = vn
            # r = [satpos]*dx[dx,dy,dz,dtr]
            dx = tool.LSP(n, k, m, e_matrix, V, dx)
            for i in range(4):
                X[i] += dx[i]
            # print("X=", X)
            # print("dx=", dx)
            count += 1
            if np.sqrt(tool.dot(dx) < 1E-4):
                print("LSP处理次数：", count)
                break
        return X
    #计算出伪距矩阵V、接收机相对于卫星位置的矩阵H（乘以-1）、解算卫星数vn、卫星权重矩阵Var
    def rescode(self,nav_list, OBS_P, rr, pos, err, dtr, V, Var):
        #vn是参与解算卫星数
        vn = 0
        e_matrix = []
        for i in range(len(nav_list)):
            # 伪距修正，根据卫星位置，减去估算的用户位置
            r = self.r_corr(nav_list[i], rr)  # [nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z]
            # 计算接收机相对于卫星位置的向量数组
            e = self.e_corr(nav_list[i], rr)
            # 计算卫星的方位角和截止角

            azel = self.satAzel(pos, e)

            # 小于10度截止的的就不参与解算。
            if azel[1] < 0.173:
                continue
            #得到伪距矩阵V
            #V[vn] = OBS_P[i] - (r + dtr - tool.CLIGHT * nav_list[i].dts)
            V[vn] = OBS_P[i] - (r + dtr - tool.CLIGHT * nav_list[i].dts)
            # 将e数组变乘以-1，主要是用于后续最小二乘法中去。并合并到大数组中，得到LSP需要的矩阵H
            for j in range(3):
                e[j] *= -1
            e_matrix.append(e)

            # iorr
            # dion
            # dtorr
            # dtrp
            # Var 计算伪距权重。
            Var[vn] = 1 * err[0] * err[0] * (err[1] * err[1] + err[2] * err[2]) / np.sin(azel[1]) + nav_list[
                i].vare + 0.09 + 0.07438
            # print("r = ", r, "e = ",              e)  # nav_list[i].nav.prn,nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z,nav_list[i].dts,nav_list[i].vare,V[i],r)
            vn += 1
        return [vn, e_matrix, Var, V]
