import numpy as np
import math as m
import matplotlib.pyplot as plt
def myplot(xAry,yAry,zAry,n):
    t = " group "+str(n)+" result"
    plt.title(t)
    plt.xlabel("time(s)")
    plt.ylabel("error(m)")
    # X = [1,2,3,4,5,6,7,8,9,10,11,12,13
    #     ,14,15,16,17,18,19,20,21,22,23,24
    #     ,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39
    #     ,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56
    #     ,57,58,59,60]
    X = np.arange(0, 60)
    plt.xlim(0, 60)
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
    plt.ylim(0, 0.3)
    plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3])
    # print(xAry)
    plt.plot(X, xAry, label="x error")
    plt.plot(X, yAry, label="y error")
    plt.plot(X, zAry, label="z error")
    plt.legend(loc='upper left')
    plt.show()

    print("plot")
if __name__ == '__main__':
    fo = open("D:\IE_project\daoyuan\\20210525\data\\0423\\100C.txt",mode='r', newline='\r\n')
    xAry = []
    yAry = []
    zAry = []
    i=0
    n=0
    while True:
        line = fo.readline().strip()
        if line == "":
            break
        lineAry = line.split("\t")

        xAry.append(float(lineAry[0]))
        yAry.append(float(lineAry[1]))
        zAry.append(float(lineAry[2]))

        i += 1
        # if i== 60 :
        #     n+=1
        #     i=0
        #     print(xAry)
        #     print(yAry)
        #     print(zAry)
        #
        #     myplot(xAry,yAry,zAry,n)
        #     xAry = []
        #     yAry = []
        #     zAry = []

    t = " group 2 result"
    plt.title(t)
    plt.xlabel("time(s)")
    plt.ylabel("error(m)")
    X = np.arange(0, i)
    plt.xlim(0, i)
    # plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
    plt.ylim(0, 0.3)
    plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3])
    # print(xAry)
    plt.plot(X, xAry, label="x error")
    plt.plot(X, yAry, label="y error")
    plt.plot(X, zAry, label="z error")
    plt.legend(loc='upper left')
    plt.show()
        #novatel_header = convert_novatel_header(fline2[0], 'a')
        #   inspvaxlist.append(convert_inspvax(novatel_header, fline2[1], 'a'))
