import math as m


class Earth():
    # 建立椭球类，，采用SGCS2000中国大地坐标系数据
    #a = 6378137.0,b = 6356752.31414
    #f =0.003352810681238110752

    #WGS84:a = 6378137.0000 ,b = 6356752.3142
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


class Transform():
    # 建立转换类
    def __init__(self):  # 导入椭球数据
        self.earth = Earth()

    def BLH_XYZ(self, B, L, H):  # 大地坐标转换为空间直角坐标
        # 将角度转换为弧度
        #pi = zhouchang / 2r
        #1(deg) = (pi / 180) (rad ) zhouchang / 2r /180/1 =zhouchang / 360 r
        #1rad = (180 /pi) * ( 1deg)
        #x rad = 22.8 * (pi /180) (rad)
        #y deg = 0.399 * (180/pi)(deg)
        self.b = m.radians(B)
        self.l = m.radians(L)
        #print("B= ",B,"rad b= ",self.b,"L= ",L," rad l = ",self.l)

        # 计算辅助函数
        #w = sqrt(1- (1-b*b/a*a) * sin(b)*sin(b))
        self.W = m.sqrt(1 - self.earth.e1s * m.pow(m.sin(self.b), 2))
        #print("m.sin(self.b) = ",m.sin(self.b))
        #print("m.pow(m.sin(self.b), 2) = ",m.pow(m.sin(self.b), 2))
        #print("self.earth.e1s * m.pow(m.sin(self.b), 2) = ",self.earth.e1s * m.pow(m.sin(self.b), 2))
        #print("")
        #print("w = " ,self.W)
        # 转换为空间直角坐标

        self.N = self.earth.a / self.W
        #print("N = ",self.N)
        self.X = (self.N + H) * m.cos(self.b) * m.cos(self.l)
        self.Y = (self.N + H) * m.cos(self.b) * m.sin(self.l)
        self.Z = (self.N * (1 - self.earth.e1s) + H) * m.sin(self.b)

    def XYZ_BLH(self, X, Y, Z):  # 空间直角坐标转换为大地坐标
        # 求出大地经度
        self.l = m.atan(Y / X)

        # 求出大地纬度
        self.r = m.sqrt(X * X + Y * Y)
        self.tb1 = Z / self.r
        while True:
            self.tb2 = 1 / self.r * (
                Z + self.earth.a * self.earth.e1s * self.tb1 /
                m.sqrt(1 + self.tb1 * self.tb1 * (1 - self.earth.e1s)))
            if abs(self.tb2 - self.tb1) <= 5e-10:
                break
            self.tb1 = self.tb2
        self.b = m.atan(self.tb2)

        # 求出大地高
        self.W = m.sqrt(1 - self.earth.e1s * m.pow(m.sin(self.b), 2))
        self.N = self.earth.a / self.W
        self.H = self.r / m.cos(self.b) - self.N

        # 将弧度转换为角度
        self.B = m.degrees(self.b)
        self.L = m.degrees(self.l)


class Point():
    # 建立点类
    def BLH(self, B, L, H):  # 建立直角坐标点类
        self.B = B
        self.L = L
        self.H = H
        self.P = Transform()
        self.P.BLH_XYZ(B, L, H)
        self.X = self.P.X
        self.Y = self.P.Y
        self.Z = self.P.Z

    def XYZ(self, X, Y, Z):  # 建立直角坐标点类
        self.X = X
        self.Y = Y
        self.Z = Z
        self.P = Transform()
        self.P.XYZ_BLH(X, Y, Z)
        self.B = self.P.B
        self.L = self.P.L
        self.H = self.P.H

if __name__ == '__main__':
    earth = Earth()
    point = Point()
    #-3869296.6953060967, 3436590.1719097085, 3717376.5186141585
    point.XYZ(-2273738.038, 5390362.523, 2532110.459)
    #point.BLH(22.5,90,10)
    print("x = ",point.B,"y = ",point.L+180,"h = ",point.H)