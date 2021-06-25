
class positionsol():
    gpsweek = 0
    gpssec = 0
    X = []
    ns = 0
    age = 0

    def update(self,X,ns,age):
        self.X = X
        self.ns = ns
        self.age = age
        return self
    def init(self,gpsweek,gpssec):
        self.gpsweek = gpsweek
        self.gpssec = gpssec


