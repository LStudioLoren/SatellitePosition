if __name__ == '__main__':
    fo = open("../data/no_git_folder/1v2-0622.txt",
              mode='r', newline='\r\n')
    fo2 = open("../data/no_git_folder/TLK_raw_data4.txt",
               mode='w')
    fuhao = 0
    str = ""
    while True:
        line = fo.readline().strip()
        if line == "#VALUE!" :
            continue
        if line == "":
            break
        str+=line
        #print("count",str.count("#")+str.count("$"))
        if str.count("#") >= 2:
            #print("--str:",str)
            strary = str.split("#")
            for i in range(len(strary)-1):
                #print("#"+strary[i])
                a = strary[i].replace("..","")
                if a.startswith("#") :
                    a = a
                else:
                    a = "#"+a
                fo2.write(a + "\r")
            str = "#"+strary[len(strary)-1]
            #print("---re str:",str)
        elif str.count("#")>0 and str.count("$") >0:
            #print("------2----")
            #print("str2 ",str)
            strary1 = str.split("$")
            #print(str.find("#"),str.find("$"))
            if str.find("#")<str.find("$"):
                strary1 = str.split("$")
                for i in range(len(strary1)-1):
                    #print("$" +strary1[i])
                    a = strary1[i].replace("..","")
                    if a.startswith("$"):
                        a = a
                    else:
                        a = "$" + a
                    fo2.write(a + "\r")

                str = "$"+strary1[len(strary1)-1]

            elif str.find("#")>str.find("$"):
                strary1 = str.split("#")
                for i in range(len(strary1) - 1):
                    #print("#"+strary1[i])
                    a = strary1[i].replace("..","")

                    if a.startswith("#"):
                        a = a
                    else:
                        a = "#" + a
                    fo2.write(a + "\r")
                str = "#" + strary1[len(strary1) - 1]
        elif str.count("$") >=2:
            strary = str.split("$")
            for i in range(len(strary) - 1):
                #print("$" +strary[i])
                #fo2.write( + "\r")
                a = strary[i].replace("..","")
                if a.startswith("$"):
                    a = a
                else:
                    a = "$" + a
                fo2.write(a + "\r")
                #fo2.write(strary[i] + "\r\n")
            str = "$" + strary[len(strary) - 1]
            #print("---re str:", str)


            # if strary1[0].count("#") <= 0 :
            #     print(strary1[0])
            #     strary2 = strary1[0].split("$")
            #     for i in range(len(strary2)-1):
            #         print(strary2[i])
            #     str = strary1[len(strary2)-1]
            # elif strary1[0].count("#") > 0:
            #     print(strary1[0])
            #     if strary1[1].count("$") > 0:
            #         strary2 = strary1[1].split("$")
            #         for i in range(len(strary2)-1):
            #             print(strary2[i])
            #         str = strary2[len(strary2)-1]

