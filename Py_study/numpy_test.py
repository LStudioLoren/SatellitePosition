import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import linalg
#gps week,gps sec,lat,lon,h,v-e,v-n,v-h,roll,pitch,heading,num_sat;
fig = plt.figure(0)
ax = fig.add_subplot(111, aspect='equal')
e = Ellipse(xy = (0,0), width =  63.78* 2, height = 63.56 * 2, angle = 0)
ax.add_artist(e)




#e.set_facecolor("white")
plt.xlim(-2000, 2000)
plt.ylim(-2000, 2000)
ax.grid(True)
plt.title("50% Probablity Contour - Homework 4.2")

# #plt.show()
# xd = [0,6378137.0-521854.00,6378137.0,6378137.0+521854.00,6378137.0*2]
# yd = [0,0,6356752.31414,0,0]
# #x = np.arra.arange(1, 10)
# y = np.arange(11,20)
# plt.title("Matplotlib demo")
# plt.xlabel("x axis caption")
# plt.ylabel("y axis caption")
# plt.plot(x,y,'bp')
# #plt.show()


class Earth():
    # 建立椭球类，，采用SGCS2000中国大地坐标系数据
    def __init__(self):  # 给定椭球数据
        # 长半轴
        self.a = 6378137.0
        # 短半轴
        self.b = 6356752.31414
        # 扁率
        # (a - b )/ a = 1 / (a / ( a - b ))
        self.f = 1/298.257222101
        # 第一偏心率的平方
        # sqrt(a*a - b*b)/a = sqrt(1-b*b/a*a)  c/a = sin(x)
        self.e1s = m.pow(0.0818191910428, 2)
        # 第二偏心率的平方
        # sqrt(a*a - b*b)/b   c/b = tan(x)
        self.e2s = self.e1s / (1 - self.e1s)
earth = Earth()

e1 = m.sqrt(earth.a*earth.a-earth.b*earth.b)/earth.a
e2 = m.sqrt(earth.a*earth.a-earth.b*earth.b)/earth.b
x = m.asin(e1)
print("e1 = ",e1," e2 = ",e2)
print("x = ",x,"e1 = ",m.sin(x),"e2 = ",m.tan(x))
print("c = ",e1*earth.a)
print("c2 = ",e2*earth.b)


# p1(1,0,1),p2(1,0,0),p3(0,1,0),p4(0,1,1)
# R1 = 1.04,R2 = sqrt(2)+0.02,R3=sqrt(2),R4 = 1
# f1: R1 = sqrt((1-x)^2 + (0-y)^2 + (1-z)^2) = 1
#1^2  = 1-2x+x^2 +y2 + 1 - 2z +z^2
# f2: R2 = sqrt((1-x)^2 + (0-y)^2 + (0-z)^2) = 1+sqrt(2)
#1 +2sqrt(2)+2 = 1-2x+x^2 + y2 +z2
# f3: R3 = sqrt((0-x)^2 + (1-y)^2 + (0-z)^2) = 1+sqrt(2)
# f4: R4 = sqrt((0-x)^2 + (1-y)^2 + (1-z)^2) = 1
#p是卫星坐标，up是假设一个坐标值，r是卫星到up的伪距
def c_s(up,d):
    p = [[1, 0, 1], [1, 0, 0], [0, 1, 1]]  # ,[0,1,1]]
    for i in range(len(up)):
        up[i]+=d[i]

    #print("-----up=",up)
    r = [1, m.sqrt(2), 1]  # ,1]

    ar2 = []
    dpr = []
    for i in range(len(r)):
        ar1 = []
        dp = 0
        for j in range(len(up)):
            ar1.append((p[i][j] - up[j]) / r[i])
            dp += (m.pow((p[i][j] - up[j]),2))
        # ar1.append(-1)
        dpr.append(r[i]-m.sqrt(dp))
        ar2.append(ar1)

    A = np.array(np.linalg.inv(ar2)) # A代表系数矩阵
    b = np.array(dpr)  # b代表常数列
    x =
    x = linalg.solve(A, b)
    #print("ar2=", A)
    #print("dpr=", b)
    print("-----x= ",x,"----up= ",up,"init_std",d)
    result = []
    a=[]
    for i in range(len(x)):
        a.append(x[i])
    result.append(up)
    result.append(x)
    #print(result)
    return result

if __name__ == '__main__':
    A_arr = []
    '''
    
    init_user = [0.04,0.04,100]
    init_std = [1,1,1]
    x = c_s([0.04,0.04,0.06],init_std)#[2,2,3])#[0.1,0.1,0.1])

    #for i in range(100):
    init_user = [0.04, 0.04, 0.06]
    V = pow(x[1][0],2)+pow(x[1][1],2)+pow(x[1][2],2)
    if 0 < V <= 0.001:
        print("V= " ,V,"up=",x[0])
        #break
    else:
        #print("x[0][0]/x[1][0]=",x[0][0]/x[1][0],"x[0][1]/x[1][1]=",x[0][1]/x[1][1],"x[0][2]/x[1][2]=",x[0][2]/x[1][2])
        up = [x[0][0]*x[1][0]/(x[0][0]+x[1][0]), x[0][1]*x[1][1]/(x[0][1]+x[1][1]), x[0][2]*x[1][2]/(x[0][2]+x[1][2])]
        #up = [x[0][0]*0.95,x[0][1]*0.95,x[0][2]*0.93]
        x = c_s(init_user,up)
    print("V=",V, "up = ", x[0])
    init_user = [0.04, 0.04, 0.06]
    V = pow(x[1][0], 2) + pow(x[1][1], 2) + pow(x[1][2], 2)
    if 0 < V <= 0.001:
        print("V= ", V, "up=", x[0])
        # break
    else:
        # print("x[0][0]/x[1][0]=",x[0][0]/x[1][0],"x[0][1]/x[1][1]=",x[0][1]/x[1][1],"x[0][2]/x[1][2]=",x[0][2]/x[1][2])
        #up = [x[0][0] * x[1][0] / (x[0][0] + x[1][0]), x[0][1] * x[1][1] / (x[0][1] + x[1][1]),
              #x[0][2] * x[1][2] / (x[0][2] + x[1][2])]
        up = [x[0][0]*0.95,x[0][1]*0.95,x[0][2]*0.93]
        x = c_s(init_user, up)
    print("V=", V, "up = ", x[0])
    init_user = [0.04, 0.04, 0.06]
    V = pow(x[1][0], 2) + pow(x[1][1], 2) + pow(x[1][2], 2)
    if 0 < V <= 0.001:
        print("V= ", V, "up=", x[0])
        # break
    else:
        # print("x[0][0]/x[1][0]=",x[0][0]/x[1][0],"x[0][1]/x[1][1]=",x[0][1]/x[1][1],"x[0][2]/x[1][2]=",x[0][2]/x[1][2])
        up = [x[0][0] * x[1][0] / (x[0][0] + x[1][0]), x[0][1] * x[1][1] / (x[0][1] + x[1][1]),
              x[0][2] * x[1][2] / (x[0][2] + x[1][2])]
        # up = [x[0][0]*0.95,x[0][1]*0.95,x[0][2]*0.93]
        x = c_s(init_user, up)
    print("V=", V, "up = ", x[0])
    '''
    
    x = c_s([1.04, 1.04, 1.06], [0.1, 0.1, 0.3])
    print("----------------------------------------------------")
    for i in range(100):
        V = pow(x[1][0],2)+pow(x[1][1],2)+pow(x[1][2],2)

        if 0 < V <= 0.001:
            print("V= " ,V,"up=",x[0])
            break
        elif V < pow(0.1,2)+pow(0.1,2)+pow(0.3,2):
            up = [x[0][0]*0.95,x[0][1]*0.95,x[0][2]*0.93]
            x = c_s([1.04,1.04,1.06],up)
        else:
            print("i=",i,"V=",V, "up = ", x[0])

    x2 = c_s([0.04,0.04,0.06],[0.96,0.96,0.94])
    print(x2)
    '''