import numpy as np
import math as m
import time
from scipy import linalg
CLIGHT =299792458.0  #光速
OMGE=7.2921151467E-5  #
FE_WGS84 = 1 / 298.257223563 # 地球扁率 1/长半轴-短半轴
RE_WGS84 = 6378137.0   #地球长半轴长度

GM = 3.986005E+14  #引力常数GM
EARTH_RAD = 7.29211567E-5   #地球自转常数EARTH_RAD
class NavData():
    #星历数据
    nav = 0
    #钟差
    dts = 0
    #方差
    vare = 0
    def initNavData(self,nav,dts,vare):
        self.nav = nav
        self.dts = dts
        self.vare = vare
class Obss():
    prn = 0
    p1 = 0
    p2 = 0
    def init_obss(self,prn,p1,p2):
        self.prn = prn
        self.p1 = p1
        self.p2 = p2
class Obs():
    t_obs = 0
    obsary = []
    def init_obs(self,t_obs,obssary):
        self.t_obs = t_obs
        self.obsary = obssary
class Nav():
    t_toc = 0
    sys = 0
    prn = 0
    a0 = 0.489834230393E-03
    a1 = 0.522959453519E-11
    a2 = 0.000000000000E+00
    IODE = 0.680000000000E+02
    Crs = -0.228125000000E+01
    Delta_n = 0.529879214453E-08
    M0 = 0.960630544678E+00
    Cuc = -0.465661287308E-07
    e = 0.127291339450E-01
    Cus = 0.915490090847E-05
    Sqrt_A = 0.515367845154E+04
    toe = 0.518400000000E+06
    Cic = -0.130385160446E-06
    OMEGA = 0.108459207992E+01
    Cis = 0.931322574615E-07
    i0 = 0.926542338573E+00
    Crc = 0.187343750000E+03
    omega = 0.973659372090E+00
    OMEGA_DOT = -0.844606609827E-08
    IDOT = 0.442161274935E-09
    CODE = 0
    WEEK = 0
    FLAG = 0
    SAV = 0
    SVH = 0

    x = 0
    y = 0
    z = 0
    def satPos(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    #user range accurcy  从ICD中获取，固定的。
    def ura_index(self,sva):
        ura_eph = [2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,3072.0,6144.0]
        for i in range(15):
            if ura_eph[i] >= sva:
                break
        return ura_eph[i]

    def initNav(self,nav_str):
        navarray = nav_str
        #print(navarray)
        self.prn = navarray[0]
        self.t_toc = navarray[1]
        self.a0 = navarray[2]
        self.a1 = navarray[3]
        self.a2 = navarray[4]
        self.IODE = navarray[5]
        self.Crs = navarray[6]
        self.Delta_n = navarray[7]
        self.M0 = navarray[8]
        self.Cuc = navarray[9]
        self.e = navarray[10]
        self.Cus = navarray[11]
        self.Sqrt_A = navarray[12]
        self.toe = navarray[13]
        self.Cic = navarray[14]
        self.OMEGA = navarray[15]
        self.Cis = navarray[16]
        self.i0 = navarray[17]
        self.Crc = navarray[18]
        self.omega = navarray[19]
        self.OMEGA_DOT = navarray[20]
        self.IDOT = navarray[21]
        self.CODE = navarray[22]
        self.WEEK = navarray[23]
        self.FLAG = navarray[24]
        self.SAV = self.ura_index(navarray[25])
        self.SVH = navarray[26]

def dot(rr):
    return rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]
def dot2(rr):
    return rr[0]*rr[0]+rr[1]*rr[1]

#计算卫星的对应测量点的方位角；
#pos是卫星的位置x，y，z
#e是是卫星的位置x，y,z与伪距r的比
def satAzel(pos,e):
    #print("satAzel e = ",e)
    az = 0
    el = np.pi/2
    #enu = []
    if pos[2] > -RE_WGS84:
        enu = xyz2enu(pos,e)
        #print("enu = ",enu)
        if enu[0]*enu[0]+enu[1]*enu[1] < 1E-12 :
            az = 0
        else:
            az = np.arctan2(enu[0],enu[1])
        if az < 0.0 :
            az+=2*np.pi
        el = np.arcsin(enu[2])
        #print("az = ",az,"el",el,enu[0],enu[1],enu[2])
    return [az,el]

def xyz2enu(pos,e):
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
                d += E9[i + x * n] * e[x]#  E1*e[x]/(dot(e_x,e,y,e_z)+E2*e[y]/dot+E3*e[
            # print(d)
            enu.append(d)
    #print(enu)
    enu[0] = enu[0]/np.sqrt(dot(e))
    enu[1] = enu[1]/np.sqrt(dot(e))
    enu[2] = enu[2]/np.sqrt(dot(e))
    #print(enu)
    #print("asdfa enu =  [0.9832349746631751, -0.0017021136149081018, 0.18233509647993557]")
    return enu

def eceftopos(rr,pos):
    #print("eceftopos: rr = ",rr)
    e2 = FE_WGS84 * (2 - FE_WGS84)
    z = zk =v =RE_WGS84
    R2 = dot2(rr)
    sinp=0
    zk = 0
    z = rr[2]
    #print(R2)
    while(1):

        if np.fabs(z-zk) >=1E-4:
            zk = z
            sinp = z/np.sqrt(R2+z*z)
            v = RE_WGS84/np.sqrt(1-e2*sinp*sinp)
            z = rr[2]+v*e2*sinp
            #print("zk = ",zk," z = ",z," np.fabs(z-zk) = ",np.fabs(z-zk))
        else:
            break;

    #pos1 = pos2 = pos3 = 0

    if R2 > 1E-12:
        #print(np.sqrt(R2),z)
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

    #print(pos)
    return pos
'''根据星历数据计算每颗卫星的位置'''

def satpos(navData,t_obs,OBS_P):


    #根据观测量的伪距，除以光速，算出时差
    tttt = OBS_P/CLIGHT
    #print(tttt)
    t_obs = t_obs-tttt
    #更新星历的toc时间
    navData.nav.t_toc = navData.nav.toe-tttt
    #print(navData.nav.t_toc)
    #根据开普勒第三定律
    n0 = np.sqrt(GM) / (navData.nav.Sqrt_A * navData.nav.Sqrt_A * navData.nav.Sqrt_A)

    n = n0 + navData.nav.Delta_n
    # print("参考时刻TOE的平均角速度n0 = ", n0)
    # print("计算卫星运动的平均角速度n = n0 + delta_n：", n,"  delte_n = ",Delta_n)
    delta_t = navData.nav.a0 + navData.nav.a1 * (t_obs - navData.nav.t_toc) + navData.nav.a2 * (t_obs - navData.nav.t_toc) * (t_obs - navData.nav.t_toc)
    # print("delta_t = ",delta_t)
    tk = t_obs - navData.nav.toe - delta_t
    #print("tk = t_obs - delta_t - toe",tk)
    Mk = navData.nav.M0 + n * tk
    # print("信号发射时卫星的平近点角Mk：", Mk,"  M0  = ",M0)
    ed = Mk
    for i in range(4):
        #Ek = ed
        ed = ed - (ed - navData.nav.e * np.sin(ed) - Mk) / (1 - navData.nav.e * np.cos(ed))

    # print("ed = ",Ek -ed )
    Ek = ed

    Vk = m.atan2(np.sqrt(1 - navData.nav.e * navData.nav.e) * np.sin(Ek), (np.cos(Ek) - navData.nav.e))
    # print("Vk = ",Vk)
    u = navData.nav.omega + Vk
    COS2U = np.cos(2 * u)
    SIN2U = np.sin(2 * u)
    delta_u = navData.nav.Cuc * COS2U + navData.nav.Cus * SIN2U
    delta_r = navData.nav.Crc * COS2U + navData.nav.Crs * SIN2U
    delta_i = navData.nav.Cic * COS2U + navData.nav.Cis * SIN2U
    # print("delta_u = ",delta_u)
    # print("delta_r = ", delta_r)
    # print("delta_i = ", delta_i)
    uk = u + delta_u
    rk = (navData.nav.Sqrt_A * navData.nav.Sqrt_A) * (1 - navData.nav.e * np.cos(Ek)) + delta_r
    ik = navData.nav.i0 + delta_i + navData.nav.IDOT * tk

    x = rk * np.cos(uk)
    y = rk * np.sin(uk)
    # print("x = ",x, "y = ", y)

    L = navData.nav.OMEGA + (navData.nav.OMEGA_DOT - EARTH_RAD) * tk - EARTH_RAD * navData.nav.toe

    # print("L = ",L,"  ik =  ",ik)
    COSL = np.cos(L)
    SINL = np.sin(L)
    COSiK = np.cos(ik)
    SINiK = np.sin(ik)

    Xs = x * COSL - y * COSiK * SINL
    Ys = x * SINL + y * COSiK * COSL
    Zs = y * SINiK
    tk_dts = t_obs-navData.nav.t_toc
    dts = navData.nav.a0 +navData.nav.a1*tk_dts+navData.nav.a2*tk_dts*tk_dts
    #dts卫星的钟差，钟漂{bias,drift},没有考虑
    # vare用户测距精度，从卫星星历钟的SAV中获取
    dts = dts - 2*np.sqrt(GM*navData.nav.Sqrt_A*navData.nav.Sqrt_A)*navData.nav.e*np.sin(Ek)/(CLIGHT*CLIGHT)
    vare = navData.nav.SAV*navData.nav.SAV
    #print(navData.nav.prn,"dts = ",dts,"vare = ",vare)
    #print("Xs = ", Xs, " Ys = ", Ys, " Zs = ", Zs)
    navData.nav.satPos(Xs,Ys,Zs)
    navData.dts = dts
    navData.vare = vare
    return navData

#计算伪距残差
def r_corr(nav,rr):

    r = np.sqrt(dot([nav.x-rr[0],nav.y-rr[1],nav.z-rr[2]]))
    #再加上
    r += OMGE*(nav.x*rr[1]-nav.y*rr[0])/CLIGHT
    return r
def e_corr(nav,rr):
    i = 0
    e = [0,0,0,1]
    #for i in range(3):
    e[0] = nav.x - rr[0]
    e[1] = nav.y - rr[1]
    e[2] = nav.z - rr[2]
    r = np.sqrt(dot(e))
    #print("r = ", r, "  e_before= ", e)
    for i in range(3):
        #计算接收机到卫星的向量vector
        e[i]= e[i]/ r
    #print("e_after= ", e)
    return e
def mulmatirix_signle(n,k,m,A,B):
    C = []
    for i in range(n):
        d = 0
        for j in range(k):
            for x in range(m):
                d += A[x][i] * B[x]
            C.append(d)
    #print(" = ", C)
    return C
def mulmatirix_multi(n,k,m,A,B):
    AAT = []
    for i in range(n):
        d = []
        for j in range(k):
            # 0 1 00 01 10 11 20 21
            d_i = 0
            for x in range(m):
                d_i += A[x][i] * B[x][j]
            d.append(d_i)
        AAT.append(d)
    #print("Q = A*A'=", AAT)
    return AAT
#laster sqrt p
def LSP(n,k,m,A,y,C):
    #A*dx=y[p-r]  -》A*dx-y = 0
    #根据矩阵
    #Ay = A*y
    Ay = mulmatirix_signle(n,1,m,A,y)
    #Q = A*A'
    Q = mulmatirix_multi(n,4,m,A,A)
    #Q_1 = Q inv,Q-1
    Q_1 = np.linalg.inv(Q)
    C = mulmatirix_signle(n,1,n,Q_1,Ay)
    return C
def estpos(nav_list, OBS_P,last_X):
    i = 0
    # V是伪距残差数列
    V = [0,0,0,0,0,0,0,0,0]
    err = [100, 0.003, 0.003]
    rr = [0.0, 0.0, 0.0]
    # X是x,y,z,dtr位置
    X = last_X
    # dx是x,y,z误差
    dx = []
    # Var是误差方差（与vare关联）
    Var = [0,0,0,0,0,0,0,0,0]
    count = 0
    while (1):
        # e矩阵
        e_matrix = []
        pos = [0, 0, 0]
        for i in range(3):
            rr[i] = X[i]
        dtr = X[3]
        # rr是x，Y,z，pos是ecef的位置
        pos = eceftopos(rr, pos)

        # vn可用卫星数
        # vn = 0

        rescode_data = []
        rescode_data = rescode(nav_list, OBS_P, rr, pos, err,dtr,V,Var)
        vn = rescode_data[0]
        e_matrix = rescode_data[1]
        Var = rescode_data[2]
        V = rescode_data[3]
        #print("e = ",e_matrix)
        i = j = 0
        n = vn
        m = 4
        #print(e_matrix)
        #print(V)
        for i in range(n):
            #对矩阵进行加权，权重值为每个卫星对应的Var。
            sig = np.sqrt(Var[i])
            V[i] /= sig
            for j in range(m):
                e_matrix[i][j] /= sig
        #print("e= ", e_matrix)
        #print("V=", V)
        #print(Var)
        n = 4
        k = 4
        m = vn
        #r = [satpos]*dx[dx,dy,dz,dtr]
        dx = LSP(n, k, m, e_matrix, V, dx)
        i = 0
        for i in range(4):
            X[i] += dx[i]
        #print("X=", X)
        #print("dx=", dx)
        count += 1
        if np.sqrt(dot(dx) < 1E-4):

            print("LSP处理次数：",count)
            break
    return X

def rescode(nav_list,OBS_P,rr,pos,err,dtr,V,Var):
    i = 0
    vn = 0
    e_matrix = []
    for i in range(len(nav_list)):
        #伪距修正，根据卫星位置，减去估算的用户位置
        r = r_corr(nav_list[i].nav, rr)  # [nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z]

        e = e_corr(nav_list[i].nav, rr)
        #计算卫星的方位角和截止角
        azel = satAzel(pos, e)

        #小于10度截止的的就不参与解算。
        if azel[1] < 0.173:
            continue
        V[vn] = OBS_P[i] - (r + dtr - CLIGHT * nav_list[i].dts)
        for j in range(3):
            e[j] = -e[j]
        e_matrix.append(e)

        # iorr
        # dion
        # dtorr
        # dtrp
        #Var 计算伪距权重。
        Var[vn] = 1 * err[0] * err[0] * (err[1] * err[1] + err[2] * err[2]) / np.sin(azel[1]) + nav_list[
            i].vare + 0.09 + 0.07438
        #print("r = ", r, "e = ",              e)  # nav_list[i].nav.prn,nav_list[i].nav.x,nav_list[i].nav.y,nav_list[i].nav.z,nav_list[i].dts,nav_list[i].vare,V[i],r)
        vn += 1
    return [vn,e_matrix,Var,V]

def UTC2GPST(year,month,day,hour,min,sec,leapsec):

    DayOfYear = 0
    DayOfMonth = 0
    for i in range(1980,year):
        if (i % 4 == 0 and i%100 !=0) or i%400 ==0 :
            DayOfYear += 366
        else:
            DayOfYear += 365

    for i in range(1,month):
        if i==1 or i==3 or i==5 or i==7 or i==8 or i==10 or i==12:
            DayOfMonth+=31
        elif i==4 or i==6 or i==9 or i==11:
            DayOfMonth+=30
        else:
            if (year%4 == 0 and year %100 !=0) or year %400 ==0:
                DayOfMonth+=29
            else:
                DayOfMonth+=28
    allday = DayOfMonth+day+DayOfYear-6
    #print(allday)
    GPSWeek = int(allday / 7)
    GPSSec = allday%7*86400+hour*3600+min*60+sec+leapsec
    #print(GPSWeek,GPSSec)
    return [GPSWeek,GPSSec]

def readNavData(filepath):
    fo = open(filepath, mode='r', newline='\r\n')
    while True:
        #fo.read()
        line = fo.readline()
        #print(line)
        if line.find("END OF HEADER") != -1:
            break;
    lines = []
    NAV_DATA = []
    i = 0
    while True:
        line = fo.readline().strip()
        if line =="" :
            break
        i+=1
        line = line.replace("      ","").replace("     ","").replace("  "," ").replace(".","0.").replace("D","E").split(" ")
        if i == 1 :
            prn = int(line[0].replace("G",""))
            year = int(line[1])
            month = int(line[2])
            day = int(line[3])
            hour = int(line[4])
            min = int(line[5])
            sec = float(line[6])
            gpst = UTC2GPST(year,month,day,hour,min,sec,18)
            lines.append(prn)
            #lines.append(gpst[0])
            lines.append(gpst[1])
            for j in range(7,len(line)):
                lines.append(float(line[j]))
        else:
            for j in range(len(line)):
                lines.append(float(line[j]))
        if i == 8 :
            NAV_DATA.append(lines)
            lines = []
            i = 0
    return NAV_DATA

def readObsData(filepath):
    fo = open(filepath, mode='r', newline='\r\n')
    while True:
        # fo.read()
        line = fo.readline()
        # print(line)
        if line.find("END OF HEADER") != -1:
            break;
    #[[GPSsec,[[prn,p1],[prn,p1],...]],...]
    i = 0
    obss = Obss()
    obss_line = []
    obs_line = []
    gpst = []
    while True:
        line = fo.readline().strip()
        if line == "" :
            break
        if line.find(">") != -1 :
            #读取5个历元数据进行测试。
            if i == 6 :
                break
            if i>0 :
                obs = Obs()
                obs.init_obs(gpst[1], obss_line)
                obs_line.append(obs)
                obss_line = []
            line = line.split(" ")
            i += 1
            gpst = UTC2GPST(int(line[1]),int(line[2]),int(line[3]),int(line[4]),int(line[5]),float(line[6]),18)

        else:
            line_len = len(line)
            #print(line_len)
            if line[:3].find("G") == -1 :
                continue
            if line_len <4:
                continue
            elif line_len < 66:
                obss.init_obss(int(line[:3].replace("G","")),float(line[4:17]),0)
            elif line_len <330:
                obss.init_obss(int(line[:3].replace("G","")), float(line[4:17]), float(line[66:82]))
            # elif line_len <200:
            #     continue
            # elif line_len <260:
            #     continue
            # elif line_len <330:
            #     continue
            obss_line.append(obss)
            obss = Obss()
            #print(obss.prn,obss.p1,obss.p2)
    # for i in range(len(obs_line)):
    #     #print(obs.t_obs)
    #     for j in range(len(obs_line[i].obsary)):
    #         print(obs_line[i].obsary[j].prn,obs_line[i].obsary[j].p1,obs_line[i].obsary[j].p2)
    return obs_line
if __name__ == '__main__':
    #读取广播星历数据
    NAV_DATA = readNavData("D:\program\python\python-project\project1\data\\nav_gps.21N")
    #读取观测量数据
    OBS_DATA = readObsData("D:\program\python\python-project\project1\data\\Obs_gps.21O")

    #卫星观测量对应的L1伪距
    #整理星历数据打包程nav对象，计算卫星位置pos，计算卫星种差dts，计算卫星位置及时钟方差vare，这几个值都是固定，只需要执行一次。
    nav_data_list = []
    for i in range(len(NAV_DATA)):
        navData = NavData()
        nav = Nav()
        nav.initNav(NAV_DATA[i])
        navData.nav = nav
        nav_data_list.append(navData)
    X = [0.0, 0.0, 0.0, 0.0]
    for j in range(len(OBS_DATA)):
        OBS_P = []
        nav_list = []
        for k in range(len(OBS_DATA[j].obsary)):
            for i in range(len(nav_data_list)):
                # 计算卫星位置；
                if OBS_DATA[j].obsary[k].prn == nav_data_list[i].nav.prn :
                    nav_list.append(satpos(nav_data_list[i],OBS_DATA[j].t_obs,OBS_DATA[j].obsary[k].p1))
                else:
                    continue

            #     print("prn :", navData.nav.prn, " x = ", navData.nav.x, " y=", navData.nav.y, " z= ", navData.nav.z, " r = ",
            #       np.sqrt(dot([navData.nav.x, navData.nav.y, navData.nav.z])))
            # print("t_obs = ",OBS_ALL[j].t_obs," prn = ",OBS_ALL[j].obsary[k].prn," p1 = ", OBS_ALL[j].obsary[k].p1)
            OBS_P.append(OBS_DATA[j].obsary[k].p1)
        print("开始处理时间：",time.time())
        X = estpos(nav_list, OBS_P,X)

        print("time = ",OBS_DATA[j].t_obs,"  X = ",X)
        print("结束处理时间：", time.time())



