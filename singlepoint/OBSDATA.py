import tool
#单个OBS数据
class Obss():
    #卫星系统
    sys=0
    # 卫星号
    prn = 0
    #L1伪距
    p1 = 0
    #L2伪距
    p2 = 0
    def init_obss(self,prn,p1,p2):
        self.prn = prn
        self.p1 = p1
        self.p2 = p2
#同一时刻的OBS数据集
class Obs():
    #OBS时刻
    t_obs = 0
    #OBS集
    obsary = []
    def init_obs(self,t_obs,obssary):
        self.t_obs = t_obs
        self.obsary = obssary
#初始化OBS数据，最后打包成OBS数据集
class ObsData():
    def initObsData(self,filepath):
        fo = open(filepath, mode='r', newline='\r\n')
        while True:
            # fo.read()
            line = fo.readline()
            # print(line)
            if line.find("END OF HEADER") != -1:
                break;
        #[[GPSsec,[[prn,p1],[prn,p1],...]],...]
        i = 0
        obss = Obss()
        obss_line = []
        obs_line = []
        gpst = []
        while True:
            line = fo.readline().strip()
            if line == "" :
                break
            if line.find(">") != -1 :
                #读取5个历元数据进行测试。
                if i == 6 :
                    break
                if i>0 :
                    obs = Obs()
                    obs.init_obs(gpst[1], obss_line)
                    obs_line.append(obs)
                    obss_line = []
                line = line.split(" ")
                i += 1
                gpst = tool.UTC2GPST(int(line[1]),int(line[2]),int(line[3]),int(line[4]),int(line[5]),float(line[6]),18)

            else:
                line_len = len(line)
                #print(line_len)
                if line[:3].find("G") == -1 :
                    continue
                if line_len <4:
                    continue
                elif line_len < 66:
                    obss.init_obss(int(line[:3].replace("G","")),float(line[4:17]),0)
                elif line_len <330:
                    obss.init_obss(int(line[:3].replace("G","")), float(line[4:17]), float(line[66:82]))
                # elif line_len <200:
                #     continue
                # elif line_len <260:
                #     continue
                # elif line_len <330:
                #     continue
                obss_line.append(obss)
                obss = Obss()
                #print(obss.prn,obss.p1,obss.p2)
        # for i in range(len(obs_line)):
        #     #print(obs.t_obs)
        #     for j in range(len(obs_line[i].obsary)):
        #         print(obs_line[i].obsary[j].prn,obs_line[i].obsary[j].p1,obs_line[i].obsary[j].p2)
        return obs_line