

class RTKPARAM():
    #基站坐标，ECEF-XYZ
    refPoint = []
    #rtk状态
    rtkStatus = 0

    def initRefPoint(self,X,Y,Z):
        self.refPoint = (X,Y,Z)

class RTKPoistion():
    def exeRTKPositon(self):
        print("start")

