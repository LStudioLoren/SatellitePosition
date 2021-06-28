from position.common import solustion

class RTKPARAM():
    #基站坐标，ECEF-XYZ
    refPoint = []
    #rtk状态
    rtkStatus = 0
    #解状态
    sol = solustion.positionsol()

    def initRefPoint(self,X,Y,Z):
        self.refPoint = (X,Y,Z)
    def updateSol(self,sol):
        self.sol = sol


class RTKPoistion():
    def exeRTKPositon(self):
        print("start")

