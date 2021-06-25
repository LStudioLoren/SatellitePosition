import binascii

import serial

import serial.tools.list_ports
import time
class serialConfig:

    def initPort(self,portNum,baudrate):
        self.ser = serial.Serial()
        self.ser.port = "COM"+portNum
        self.ser.baudrate = baudrate
        self.ser.timeout = 1
        self.ser.open()
        if self.ser.isOpen():
            return self.ser
        else:
            print("open serial error!check port num and baudrate!")

def updateFW(filepath,ser):
    #fo = open("D:\program\python\SPPositionWithGPSandGLO\data\OM7MR0800RN0000.shex", mode='r', newline='\r\n')
    fo = open(filepath, mode='r', newline='\r\n')

    maxlenth = 0
    # while True:
    pack_data = ""
    count = 0

    while True:
    #for i in range(10):
        fline = fo.readline().strip()

        if fline == "":
            break


        if fline.startswith("S0") or fline.startswith("S5") or fline.startswith("S7"):
            pack_data = "SOFTLOADSREC \""+fline+"\"\r\n"
            ser.write(pack_data.encode("utf8"))
            print(pack_data)
            #while True:
            time.sleep(0.2)
            if ser.inWaiting() > 0:
                # ser.write("LOG version\r\n".encode("gbk"))

                data_str = ser.read(ser.inWaiting()).decode('ascii')
                print(data_str)
            #print(ser.readlines())

            #if check.find("OK") >0:
            #break
            #print()
            pack_data = ""
        elif fline.startswith("S3"):
            line = fline.replace("S3", "", 1)
            line_len = int(line[0:2], 16)
            line_add = line[2:10]
            line_data = line[10:line_len * 2]
            # print(len(bytes(line_data, encoding='UTF8')))
            pack_data += line_data
            if len(bytes(pack_data, encoding='UTF8')) + 56 > 4096 or line_len < 33:
                pack_data2 = "SOFTLOADDATA " + pack_data +"\r\n"
                print(count,pack_data2)
                ser.write(bytes(pack_data2))
                #while True:
                    #check = ser.read(15)
                    #if check.find("OK") > 0:
                        #break
                time.sleep(0.2)
                if ser.inWaiting() > 0:
                    #ser.write("LOG version\r\n".encode("gbk"))

                    data_str = ser.read(ser.inWaiting()).decode('ascii')
                    print(data_str)
                pack_data = ""
                count += 1
    ser.write(bytes("SOFTLOADCOMMIT", encoding='ascii'))
    # while True:
    # check = ser.read(15)
    # if check.find("OK") > 0:
    # break
    time.sleep(0.2)
    while True:
        if ser.inWaiting() > 0:
        # ser.write("LOG version\r\n".encode("gbk"))
            data_str = ser.read(ser.inWaiting()).decode('ascii')
            print(data_str)
                # print("SOFTLOADDATA",fline,line_len,line_add,bytes(line_data,encoding='UTF8'))

if __name__ == '__main__':
    ser = serialConfig().initPort("13",230400)
    ser.write("SOFTLOADRESET\r\n".encode("ASCII"))
    updateFW("/data/no_git_folder/OM7MR0800RN0000.shex", ser)
    #for i in range(20):
        #ser.write("log version\r\n".encode("gbk"))
    # while True:
    #     if ser.inWaiting() > 0:
    #         ser.write("LOG version\r\n".encode("gbk"))
    #         time.sleep(0.2)
    #         data_str = ser.read(ser.inWaiting()).decode('ascii')
    #         print(data_str,end='')
    #         #else:
                #break

    #while True:
        #print(ser.read(15))
#print(count)
