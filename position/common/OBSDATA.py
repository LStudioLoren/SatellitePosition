import position.common.tool as tool


# 单个OBS频点
class Obsfref():
    #1:L1,2:L2,3:L5
    SIGNEL_TYPE = 0
    #C/A码
    P = 0
    #LLI
    LLI = 0
    #载波
    L = 0
    #Doppler
    D = 0
    #信号强度
    S = 0
    # CODE
    CODE = 0
    #pseudorange, phase, LLI, Doppler and SNR.
    def init(self, SINGEL_TYPE,P, L,LLI,D, S,CODE):
        self.SIGNEL_TYPE = SINGEL_TYPE
        self.P = P
        self.LLI = LLI & 3
        self.L = L
        self.D = D
        self.S = int(S*4+0.5)
        self.CODE = CODE


# 单个OBS数据
class Obss():
    # 卫星系统
    sys = 0
    # 卫星号
    prn = 0
    # L1伪距
    #p1 = 0
    # L2伪距
    #p2 = 0
    revflag =0
    nFref = 0
    obsfrefList = []

    def init_obss(self, prn, p1, p2):
        self.prn = prn
        self.p1 = p1
        self.p2 = p2

    def init_obss2(self, prn, obsfrefList,revflag):
        self.prn = prn
        self.obsfrefList = obsfrefList
        self.nFref = len(obsfrefList)
        self.revflag = revflag


# 同一时刻的OBS数据集
class Obs():
    gpsweek = 0
    # OBS时刻
    t_obs = 0
    # OBS集
    obsary = []

    def init_obs(self, gpsweek,t_obs, obssary):
        self.gpsweek = gpsweek
        self.t_obs = t_obs
        self.obsary = obssary


# 初始化OBS数据，最后打包成OBS数据集
class ObsData():
    def initObsData(self, filepath):
        return self.initObsData2(filepath,1)
    def initObsData2(self, filepath,flag):
        fo = open(filepath, mode='r', newline='\n')
        while True:
            # fo.read()
            line = fo.readline()
            # print(line)
            if line.find("END OF HEADER") != -1:
                break;
        # [[GPSsec,[[prn,p1],[prn,p1],...]],...]
        i = 0
        obss = Obss()
        obss_line = []
        obs_line = []
        gpst = []
        while True:
            line = fo.readline().strip()

            if line == "":
                break
            if line.find(">") != -1:
                # 读取10个历元数据进行测试。
                if i == 10:
                    break
                if i > 0:
                    obs = Obs()
                    obs.init_obs(gpst[0],gpst[1], self.sortObsary(obss_line))
                    obs_line.append(obs)
                    obss_line = []
                #line = line.split(" ")
                i += 1
                #print(2,(line[2:6]), int(line[2:6]),(line[7:9]), int(line[7:9]),(line[10:12]),int(line[10:12]), (line[13:15]),int(line[13:15]), (line[16:18]),int(line[16:18]), (line[19:29]),float(line[19:29]))

                gpst = tool.UTC2GPST(int(line[2:6]), int(line[7:9]), int(line[10:12]), int(line[13:15]), int(line[16:18]), float(line[19:29]), 0)
            else:
                type = 0
                line_len = len(line)
                # print(line_len)
                if line[:3].find("G") == -1:
                    continue
                if line_len < 4:
                    continue
                elif line_len < 66:
                    n =1
                    type = 1
                elif line_len < 330:
                    n=2

                of_list = []
                for j in range(n):
                    #for j in range([14,18,14,18]):
                    #print(line[4+j*64:17+j*64],line[19+j*64:33+j*64],line[37+j*64:49+j*64],line[50+j*64:67+j*64])
                    #G17  21634967.664 8 113692520.379 8       636.504          50.000    21634970.820 8  88591581.902 8       495.977          50.000    21634971.516 8  88591585.906 8       495.977          49.000
                    #01234567890123456789012345678901234567890123456789
                    #    012345678901234

                    of = self.readObsFref(j+1, float(line[4+j*64:17+j*64]),#4:17、68:81
                                              float(line[19+j*64:33+j*64]),#18:35
                                              (0 if line[33+j*64:34+j*64] == ' ' else int(line[33+j*64:34+j*64])),
                                              float(line[37+j*64:49+j*64]),#36:49
                                              float(line[50+j*64:67+j*64]),0)#50:67
                    of_list.append(of)
                if type == 1 :
                    of_list.append(self.readObsFref(2,0,0,0,0,0,0))
                obss.init_obss2(int(line[:3].replace("G", "")),of_list,flag)
                # elif line_len <200:
                #     continue
                # elif line_len <260:
                #     continue
                # elif line_len <330:
                #     continue
                obss_line.append(obss)
                obss = Obss()
                # print(obss.prn,obss.p1,obss.p2)
        # for i in range(len(obs_line)):
        #     #print(obs.t_obs)
        #     for j in range(len(obs_line[i].obsary)):
        #         print(obs_line[i].obsary[j].prn,obs_line[i].obsary[j].p1,obs_line[i].obsary[j].p2)
        return obs_line
    def sortObsary(self,obsary):
        #templist = []
        length = len(obsary)
        for i in range(len(obsary)-1):
            minidex = i
            for j in range((i+1),(length)):
                if(obsary[j].prn<obsary[minidex].prn):
                    minidex = j
            temp = obsary[i]
            obsary[i] = obsary[minidex]
            obsary[minidex]=temp


            #print(obsary[minidex].prn)
        return obsary

    def readObsFref(self,fref_type,P,L,LLI,D,S,CODE):
        #print(P,L,D,S)
        of = Obsfref()
        of.init(fref_type,P,L,LLI,D,S,CODE)
        return of
