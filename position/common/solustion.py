
class positionsol():

    def __init__(self):
        self.gpsweek = 0
        self.X = [0.0,0.0,0.0,0.0]
        #gpsweek = 0
        self.gpssec = 0
        #X = [0.0]
        self.ns = 0
        self.age = 0
        self.pos_type = 0

    def update(self,X,ns,age,pos_type):
        self.X = X
        self.ns = ns
        self.age = age
        self.pos_type = pos_type
        return self
    def init(self,gpsweek,gpssec):
        self.gpsweek = gpsweek
        self.gpssec = gpssec


