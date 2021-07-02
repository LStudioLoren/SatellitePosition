import numpy as np
import math as m
from scipy import linalg
#工具集

CLIGHT =299792458.0  #光速
OMGE=7.2921151467E-5  #
FE_WGS84 = 1 / 298.257223563 # 地球扁率 1/长半轴-短半轴
RE_WGS84 = 6378137.0   #地球长半轴长度

GM = 3.986005E+14  #引力常数GM
EARTH_RAD = 7.29211567E-5   #地球自转常数EARTH_RAD

def dot(rr):
    return rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]
def dot2(rr):
    return rr[0]*rr[0]+rr[1]*rr[1]

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
def LSP(n,k,m,A,y,C):
    #A*dx=y[p-r]  -》A*dx-y = 0
    #根据矩阵
    #Ay = A*y
    Ay = mulmatirix_signle(n,1,m,A,y)
    #Q = A*A'
    Q = mulmatirix_multi(n,k,m,A,A)
    #Q_1 = Q inv,Q-1
    #这个是将A*A的倒转求逆
    Q_1 = np.linalg.inv(Q)
    C = mulmatirix_signle(n,1,n,Q_1,Ay)
    return C
#计算从UTC时转GPS时，且GPS从1980年1月6日开始计算。
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