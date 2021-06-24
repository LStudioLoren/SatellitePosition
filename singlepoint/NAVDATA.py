import tool

class Nav():
    t_toc = 0
    sys = 0
    prn = 0
    a0 = 0
    a1 = 0
    a2 = 0
    IODE = 0
    Crs = 0
    Delta_n = 0
    M0 = 0
    Cuc = 0
    e = 0
    Cus = 0
    Sqrt_A = 0
    toe = 0
    Cic = -0
    OMEGA = 0
    Cis = 0
    i0 = 0
    Crc = 0
    omega = 0
    OMEGA_DOT = 0
    IDOT = 0
    CODE = 0
    WEEK = 0
    FLAG = 0
    SAV = 0
    SVH = 0

    x = 0
    y = 0
    z = 0

    # 星历数据
    nav = 0
    # 钟差
    dts = 0
    # 方差
    vare = 0
    def updateParam(self,x,y,z,newDts,newVare):
        self.x = x
        self.y = y
        self.z = z
        self.dts = newDts
        self.vare = newVare

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

class NavData():

    def initNavData(self,filepath):
        NAV_DATA = self.readNavData(filepath)
        nav_data_list = []

        for i in range(len(NAV_DATA)):
            # navData = NavData()
            nav = Nav()
            nav.initNav(NAV_DATA[i])
            # navData.nav = nav
            nav_data_list.append(nav)
        return nav_data_list

    def readNavData(self,filepath):
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
                gpst = tool.UTC2GPST(year,month,day,hour,min,sec,18)
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