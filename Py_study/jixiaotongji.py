import numpy as np
import struct
import binascii
if __name__ == '__main__':
    fo2 = open("D:\program\python\python-project\project1\data\\2021技术日志-王静龙-fei.txt",mode='r',newline='\r\n')
    clientList = []
    saleList = []
    hangyeList = []
    productList = []
    dataList = []
    while True:
        line = fo2.readline().strip()
        if line == "" :
            break
        lineSplit = line.split("\t")
        print(lineSplit)
        lineArray = [lineSplit[0],lineSplit[1],lineSplit[2],lineSplit[3],lineSplit[4],lineSplit[5],lineSplit[6],int(lineSplit[7])]
        dataList.append(lineArray)
        clientflag = True
        for i in range(len(clientList)):
            if clientList[i] == lineSplit[3]:
                clientflag = False
                break
        if clientflag:
            clientList.append(lineSplit[3])

        saleflag = True
        for i in range(len(saleList)):
            if saleList[i] == lineSplit[4]:
                saleflag = False
                break
        if saleflag:
            saleList.append(lineSplit[4])

        hangyeflag = True
        for i in range(len(hangyeList)):
            if hangyeList[i] == lineSplit[5]:
                hangyeflag = False
                break
        if hangyeflag:
            hangyeList.append(lineSplit[5])

        productflag = True
        for i in range(len(productList)):
            if productList[i] == lineSplit[6]:
                productflag = False
                break
        if productflag:
            productList.append(lineSplit[6])
        #print(clientList)
    for i in range(len(clientList)):
        clienttime = 0
        for j in range(len(dataList)):
            if clientList[i] == dataList[j][3]:
                clienttime += dataList[j][7]
        print( clientList[i], clienttime)

    for i in range(len(saleList)):
        saletime = 0
        for j in range(len(dataList)):
            if saleList[i] == dataList[j][4]:
                saletime += dataList[j][7]
        print(saleList[i],  saletime)

    for i in range(len(hangyeList)):
        hangyetime = 0
        for j in range(len(dataList)):
            if hangyeList[i] == dataList[j][5]:
                hangyetime += dataList[j][7]
        print(hangyeList[i], hangyetime)

    for i in range(len(productList)):
        producttime = 0
        for j in range(len(dataList)):
            if productList[i] == dataList[j][6]:
                producttime += dataList[j][7]
        print(productList[i], producttime)

    saleclentlist = []
    salehangyelist = []
    saleproductlist = []
    for i in range(len(saleList)):
        for j in range(len(clientList)):
            saleclenttime = 0
            #salehangyetime = 0
            #saleproducttime = 0
            for n in range(len(dataList)):
                if saleList[i] == dataList[n][4] and clientList[j] == dataList[n][3] :
                    saleclenttime+=dataList[n][7]
                #if saleList[i] == dataList[n][4] and
            if saleclenttime > 0 :
                saleclentlist.append([saleList[i], clientList[j], saleclenttime])

    for j in range(len(hangyeList)):
        salehangyetime = 0
        # saleproducttime = 0
        for n in range(len(dataList)):
            if hangyeList[j] == dataList[n][5]:
                salehangyetime += dataList[n][7]
            # if saleList[i] == dataList[n][4] and
        if salehangyetime > 0:
            salehangyelist.append([hangyeList[j], salehangyetime])

    for j in range(len(productList)):
        saleproducttime = 0
        for n in range(len(dataList)):
            if productList[j] == dataList[n][6]:
                saleproducttime += dataList[n][7]
        if saleproducttime > 0:
            saleproductlist.append([productList[j], saleproducttime])

    for i in range( len(saleclentlist)):
        print(saleclentlist[i][0],saleclentlist[i][1],saleclentlist[i][2])
    for i in range(len(salehangyelist)):
        print(salehangyelist[i][0],salehangyelist[i][1])
    for i in range(len(saleproductlist)):
        print(saleproductlist[i][0],saleproductlist[i][1])
    print()




    #print(clientList)
    fo2.close()