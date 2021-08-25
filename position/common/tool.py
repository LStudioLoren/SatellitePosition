import numpy as np
import math as m
from scipy import linalg
#工具集
MAXGPS = 32#32-1+1
MAXSBAS = 23#142-120+1
MAXCMP = 35#35-1+1
MAXQZSS = 3#195-193+1
MAXGLO = 24#24-1+1
MAXSAT = 117
MAXFREF =243
CLIGHT =299792458.0  #光速
OMGE=7.2921151467E-5  #
FE_WGS84 = 1 / 298.257223563 # 地球扁率 1/长半轴-短半轴
RE_WGS84 = 6378137.0   #地球长半轴长度

GM = 3.986005E+14  #引力常数GM
EARTH_RAD = 7.2921151467E-5   #地球自转常数EARTH_RAD

#L1的频率
FREQ1 = 1.57542E9
#L2的频率
FREQ2 = 1.22760E9


#NFREQ,最大载波频点
NFREQ = 3
#光速除以L1/L2的频率
lam = [CLIGHT / FREQ1,CLIGHT / FREQ2]

ARMODE_INST = 2


def dot_n(a,b,n):
    r = 0
    for i in range(n):
        r+= a[i]*b[i]
    return r
def norm(rr,n):
    return np.sqrt(dot_n(rr,rr,n))


def dot(rr):
    return rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]
def dot2(rr):
    return rr[0]*rr[0]+rr[1]*rr[1]

#创建一个n行，m列的0矩阵
def zeroMat(n,m):
    A =[]
    for i in range(n):
        B = []
        for j in range(m):
            B.append(0)
        A.append(B)
    return A
#创建一个n行，n列的00,11,22为1的矩阵
def eyeMat(n):
    A = zeroMat(n,n)
    for i in range(n):
        A[i][i] = 1
    return A

#输入n列*m行 A矩阵  *  k列*m行 B矩阵，输出N列*K行的C矩阵。
#矩阵是A矩阵的行，乘以B矩阵的列，累加成行、列
def mulmatirix(n,k,m,A,B,alpha,bate,type,Q):
    C = []
    for i in range(n):
        dlist = []
        for j in range(k):
            d = 0
            if type == "NN":
                for x in range(m):
                    d += A[i][x] * B[x][j]
            elif type == "NT":
                for x in range(m):
                    d += A[i][x] * B[j][x]
            elif type == "TN":
                for x in range(m):
                    d += A[x][i] * B[x][j]
            elif type == "TT":
                for x in range(m):
                    d += A[x][i] * B[j][x]

            if bate != 0:
                d = alpha*d + Q[i][j]
            else:
                d = alpha*d
            dlist.append(d)
        C.append(dlist)
    # # NT,k ==n
    #
    #     for i in range(n):
    #         dlist = []
    #         for j in range(k):
    #             d = 0
    #
    #             if bate != 0:
    #                 d += Q[i][j]
    #             dlist.append(d)
    #         C.append(dlist)
    #
    #     for i in range(n):
    #         dlist = []
    #         for j in range(k):
    #             d = 0
    #
    #             if bate != 0:
    #                 d += Q[i][j]
    #             dlist.append(d)
    #         C.append(dlist)
    #
    #     for i in range(n):
    #         dlist = []
    #         for j in range(k):
    #             d = 0
    #
    #             if bate != 0:
    #                 d += Q[i][j]
    #             dlist.append(d)
    #         C.append(dlist)
    return C
# def mulmatirix(n,k,m,A,B,type):
#     Q = []
#     return mulmatirix2(n,k,m,A,B,0,type,Q)




#这里的得出来的是A的倒转*B矩阵，
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
#这里是A倒转*A
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
#Least squares approximation  最小二乘法估算，二乘就是平方的意思。
'''
单点定位核心数学问题就是，已知每个卫星的位置、以及接收机到卫星的伪距，计算接收机的位置
r = P - e*X
-1*(r-p) = -1*(-e)*X
这里的e就是接收机相对于卫星的向量矩阵，且在组装的时候，乘以-1
P-r就是伪距减去卫星的距离，即得出来是负的接收机的位置
eX = P-r
这个公式也是固定的最小二乘法公式，（e*e的倒转矩阵）的逆矩阵 * （接收机距离*e的倒转矩阵）就得出X，其中包含4个结果，分别是接收机的X、Y、Z以及dtr，钟差值。
X = (e*e_t)^-1(P-r)*e_t

'''
#A:
# x1 y1 z1 c1
# x2 y2 z2 c2
# ....
#xn yn zn cn
def LSQ(n,k,m,A,y,C):
    #dx = (A^T*A)^-1  *  (A^T *y)
    #A*dx=y[p-r]  -》A*dx-y = 0
    #根据矩阵
    #Ay = A^T*y,A^T的倒转*y
    #A矩阵输入进来是 M行，N列矩阵，需要倒转
    Ay = mulmatirix(n,1,m,A,y,1,0,"TN",[])
    #Q = A^t*A
    #Q = mulmatirix_multi(n,k,m,A,A)
    Q = mulmatirix(n,k,m,A,A,1,0,"TN",[])
    #print("Q = ",Q)
    #Q_1 = Q inv,Q-1
    #这个是将A^T*A求逆
    Q_1 = np.linalg.inv(Q)
    #这里是求Q^-1 * Ay
    C = mulmatirix(n,1,n,Q_1,Ay,1,0,"NN",[])
    #C数组是 一个 n列*1行的二维数组，因此此次将其重新编辑为一维数组
    C2 = []
    for i in range(len(C)):
        C2.append(C[i][0])
    #print("C",C2)
    return C2

def dayofyear(ep):
    dayofyear = 0
    for i in range(1,ep[1]):
        if i==1 or i==3 or i==5 or i==7 or i==8 or i==10 or i==12:
            dayofyear+=31
        elif i==4 or i==6 or i==9 or i==11:
            dayofyear+=30
        elif i==2 :
            if (ep[0] % 4 == 0 and ep[0] % 100 != 0) or ep[0] % 400 == 0:
                dayofyear += 29
            else:
                dayofyear += 28

    dayofyear = dayofyear + ep[2] + ep[3]/24 + ep[4]/1440 +ep[5]/86400
    return dayofyear

def GPST2EPOCH(gpsweek,gpssec):
    mday = [31,29,31,30,31,30,31,31,30,31,30,31,
           31,28,31,30,31,30,31,31,30,31,30,31,
           31,28,31,30,31,30,31,31,30,31,30,31,
           31,28,31,30,31,30,31,31,30,31,30,31]
    days = (gpsweek*7+6)+int(gpssec/86400)
    sec = np.mod(gpssec,86400)
    day = np.mod(days, 1461)
    ep = [0,0,0,0,0,0]
    for mon in range(48):

        if day >= mday[mon]:
            day -= mday[mon]
        else:
            break
    ep[0] = 1980+int(days/1461)*4+int(mon/12)
    ep[1] = np.mod(mon,12)+1
    ep[2] = day
    ep[3] = int(sec / 3600)
    ep[4] = int(np.mod(sec,3600)/60)
    ep[5] = np.mod(sec,60)
    #print(ep)
    return ep
#将年月日时分秒，转为GPS周+周秒，且GPS从1980年1月6日开始计算。
#当leapsec为0时，表示将rinex中的gps年月日时分秒转为周秒格式
#当leapsec不为0时，表示将UTC时转换为gps周、周秒
def UTC2GPST(year,month,day,hour,min,sec,leapsec):
    DayOfYear = 0
    DayOfMonth = 0
    #计算从1980年开始至当前时间的前一年，天数的总和
    for i in range(1980,year):
        #判断如果是闰年，则年内天数未366，否则未365
        if (i % 4 == 0 and i%100 !=0) or i%400 ==0 :
            DayOfYear += 366
        else:
            DayOfYear += 365
    #计算当前年内天数
    #先计算当前月份前一个月的天数，并判断是否为闰年。是2月份则为29天，否为28天
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
    #最后计算所有天数的总和，并减去6天，因为GPS时从UTC时1980年1月6日，00:00:00开始计算
    allday = DayOfMonth+day+DayOfYear-6
    #print(allday)
    #总天数除以7的整数部分为GPS周
    GPSWeek = int(allday / 7)
    #总天数除以7的余数部分乘以86400加上时分秒的总秒数，再加上闰秒（闰秒从1980年1月6日开始累计），为GPS周内秒
    GPSSec = allday%7*86400+hour*3600+min*60+sec+leapsec
    #print(GPSWeek,GPSSec)
    return [GPSWeek,GPSSec]