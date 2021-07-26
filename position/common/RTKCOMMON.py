from position.common import tool
import numpy as np

class SSAT():
    def __init__(self):
        self.sys = 0 #gps:0,glo:1,bds:2,GAL:3
        self.vs = 0
        self.azel = [0.0,0.0]
        #伪距残差pseudorange
        self.resp = [0.0,0.0,0.0]
        #载波残差carrier-phase
        self.resc = [0.0,0.0,0.0]
        self.vsat = [0.0,0.0,0.0]
        self.snr = [0.0, 0.0, 0.0]
        self.fix = [0.0, 0.0, 0.0]
        self.slip = [0,0,0]
        self.lock = [0.0, 0.0, 0.0]
        self.outc = [0.0, 0.0, 0.0]
        self.slipc = [0.0,0.0,0.0]
        self.rejc = [0.0,0.0,0.0]
        self.gf = 0
        self.gf2 = 0
        self.phw = 0
        self.pt = [[0,0.0],[0,0.0]]
        self.ph = [[0.0,0.0,0.0],[0.0,0.0,0.0]]


#将接收机位置从ecef转成BLH格式）
#
def eceftopos(rr, blhpos):
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
        blhpos[0] = np.arctan(z / np.sqrt(R2))
    else:
        if R2 > 0:
            blhpos[0] = np.pi / 2
        else:
            blhpos[0] = -np.pi / 2

    if R2 > 1E-12:
        blhpos[1] = np.arctan2(rr[1], rr[0])
    else:
        blhpos[1] = 0

    blhpos[2] = np.sqrt(R2 + z * z) - v

    #print("POS = ",pos)
    return blhpos

def ecef2enu(pos,e):
    enu = []
    E9 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    E9 = xyz2enu(pos,E9)
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
#将卫星的x、y、z坐标系转换为enu坐标系，公式是固定的；
    #pos是接收机的BLH，e是接收机相对于卫星的几何向量
def xyz2enu(pos,E9):

    # xyz2enu
    SINP = np.sin(pos[0])
    COSP = np.cos(pos[0])
    SINL = np.sin(pos[1])
    COSL = np.cos(pos[1])
    #0 1 2是x
    # 0 3 6
    # 1 4 7
    # 2 5 8
    E9[0] = -SINL
    E9[1] = -SINP * COSL
    E9[2] = COSP * COSL
    E9[3] = COSL
    E9[4] = -SINP * SINL
    E9[5] = COSP * SINL
    E9[6] = 0.0
    E9[7] = COSP
    E9[8] = SINP
    return E9



# 计算卫星的对应测量点的方位角；这是固定算法，需要将
def satAzel(blhpos, e):
    # print("satAzel e = ",e)
    az = 0
    el = np.pi / 2
    # enu = []
    #print("blh pos = ",pos)
    if blhpos[2] > -tool.RE_WGS84:
        enu = ecef2enu(blhpos, e)
        #print("enu = ",enu)
        if enu[0] * enu[0] + enu[1] * enu[1] < 1E-12:
            az = 0
        else:
            az = np.arctan2(enu[0], enu[1])
        if az < 0.0:
            az += 2 * np.pi
        el = np.arcsin(enu[2])
        # print("az = ",az,"el",el,enu[0],enu[1],enu[2])
    return [az, el]

    # 计算伪距残差
def r_corr(nav, rr):

    r = np.sqrt(tool.dot([nav.x - rr[0], nav.y - rr[1], nav.z - rr[2]]))
    #print("r_corr : " ,r)
    # 再加上
    r += tool.OMGE * (nav.x * rr[1] - nav.y * rr[0]) / tool.CLIGHT
    return r

def e_corr(nav, rr):
    e = [0, 0, 0, 1]
    # 将卫星的x、y、z分别减去接收机的位置，得出接收机相对于卫星的x、y、z方向距离，即e数组
    e[0] = nav.x - rr[0]
    e[1] = nav.y - rr[1]
    e[2] = nav.z - rr[2]
    # r是e数组的距离
    r = np.sqrt(tool.dot(e))
    # print("r = ", r, "  e_before= ", e)
    for i in range(3):
        # 计算接收机到卫星的向量vector，就是x
        # 这里是重点：将e数组中x、y、z除以距离r，e数组重新赋值为接收机相对于卫星的x、y、z的向量值。
        e[i] = e[i] / r
    # print("e_after= ", e)
    return e
