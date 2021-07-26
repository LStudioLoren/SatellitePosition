
class positionsol():

    def __init__(self):
        self.gpsweek = 0
        self.X = [0.0,0.0,0.0,0.0]
        self.rr = [0.0,0.0,0.0,0.0,0.0,0.0]
        self.qr = [0.0,0.0,0.0,0.0,0.0,0.0]
        #gpsweek = 0
        self.gpssec = 0
        #X = [0.0]
        self.ns = 0
        self.age = 0
        self.pos_type = 0
        self.dtr = 0
        self.ratio = 0

    def update(self,X,ns,age,pos_type):
        for i in range(3):
            self.rr[i] = X[i]
        self.X = X
        self.dtr = X[3]
        self.ns = ns
        self.age = age
        self.pos_type = pos_type
        return self
    def init(self,gpsweek,gpssec):
        self.gpsweek = gpsweek
        self.gpssec = gpssec


