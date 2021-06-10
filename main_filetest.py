# 这是一个示例 Python 脚本。

# 按 Shift+F10 执行或将其替换为您的代码。
# 按 双击 Shift 在所有地方搜索类、文件、工具窗口、操作和设置。
import numpy as np
import struct
import binascii
import CoordinateTransformation
import MSG.header as Header


from matplotlib import pyplot as plt
# novatel CRC32↓↓↓↓↓↓↓↓
CRC32_POLYNOMIAL = 0xEDB88320


def CRC32Value(ci):
    ulCRC = ci
    for j in range(8):
        if ulCRC & 1:
            ulCRC = (ulCRC >> 1) ^ CRC32_POLYNOMIAL
        else:
            ulCRC >>= 1
    return ulCRC


def CalculateBlockCRC32(ulCount, ucBuff):
    ulCRC = 0
    for ui in range(ulCount-4):
        ulTemp1 = (ulCRC >> 8) & 0x00FFFFFF
        ulTemp2 = CRC32Value((ulCRC ^ ucBuff[ui]) & 0xFF)
        ulCRC = ulTemp1 ^ ulTemp2
        ui += 1
    return ulCRC

def CalculateBlockCRC32_ascii(ulCount, ucBuff):
    ulCRC = 0
    for ui in range(ulCount):
        ulTemp1 = (ulCRC >> 8) & 0x00FFFFFF
        ulTemp2 = CRC32Value((ulCRC ^ ord(ucBuff[ui])) & 0xFF)
        ulCRC = ulTemp1 ^ ulTemp2
        ui += 1
    return ulCRC

# novatel CRC32↑↑↑↑↑↑↑

#XYZ/BLH/xyH


#XYZ/BLH/xyH

def checkCRC(line,format):
    rawCrc = 0
    crc = 0
    line_len = len(line)
    if format == 'a':
        lines = line.split("#")[1].split("*")
        crc = binascii.b2a_hex(CalculateBlockCRC32_ascii(len(lines[0]), lines[0]).to_bytes(4, byteorder='big'))
        rawCrc = bytes(lines[1], encoding='utf8')
    elif format == 'b':
        crc = CalculateBlockCRC32(line_len, line).to_bytes(4, byteorder='little')
        rawCrc = line[line_len - 4:line_len]
    '''
    to_bytes(4,byteorder='littel')方法可以将数字型的crc顺序调整。
    binascii.b2a_hex方法是将byte数组转成ascii显示的二进制；
    '''

    if  rawCrc == crc:
        # print("raw CRC = ", binascii.b2a_hex(line[line_len - 4:line_len]), "cal CRC= ", binascii.b2a_hex(crc))
        return True
    else:
        print("CRC error:RAW DATA:", line, " raw crc:", binascii.b2a_hex(line[line_len - 4:line_len]), "calculate CRC:",
              crc)
        return False




class NovatelInspvax:
    header = Header.NovatelHeader()
    ins_status = 0
    pos_type = 0
    lat = 0
    lon = 0
    height = 0
    undulation = 0
    v_north = 0
    v_east = 0
    v_up = 0
    roll = 0
    pitch = 0
    azimuth = 0
    lat_std = 0
    lon_std = 0
    height_std = 0
    v_north_std = 0
    v_east_std = 0
    v_up_std = 0
    roll_std = 0
    pitch_std = 0
    azimuth_std = 0
    ext_sol_status = 0
    time_since_update = 0

def print_hi(name):
    # 在下面的代码行中使用断点来调试脚本。
    print(f"Hi, {name}")  # 按 Ctrl+F8 切换断点。


def convert_novatel_header(line,format):
    novatel_header = Header.NovatelHeader()
    if format == 'b':
        novatel_header.header_length = int(line[3])
        novatel_header.msg_length = struct.unpack('<H', line[8:10])[0]
        novatel_header.msg_id = struct.unpack('<H', line[4:6])[0]
        novatel_header.msg_type = int(line[6])
        novatel_header.gps_sec = struct.unpack('<l', line[16:20])[0]
        novatel_header.gps_week = struct.unpack('<H', line[14:16])[0]
    elif format == 'a':
        #INSPVAXA,FILE,0,57.0,FINESTEERING,2076,464696.800,00000800,4e14,14102;
        asciiheader = line.split(",")
        novatel_header.msg_name = asciiheader[0]
        novatel_header.port = asciiheader[1]
        novatel_header.sequence = asciiheader[2]
        novatel_header.idle_time = asciiheader[3]
        novatel_header.time_status = asciiheader[4]
        novatel_header.gps_week = int(asciiheader[5])
        novatel_header.gps_sec = float(asciiheader[6])
        novatel_header.receiver_status = asciiheader[7]
        novatel_header.SW_version = asciiheader[9]

    # print("header:")
    # print("AA:", hex(line[0]))
    # print("44:", hex(line[1]))
    # print("12:", hex(line[2]))
    # headerLength = int(line[3])
    # print("header length:", novatel_header.header_length)
    # print("MSG ID:", novatel_header.msg_id)
    # print("MSG Type:", novatel_header.msg_type)
    # print("Port Add:", int(line[7]))
    # meglength = struct.unpack('<H', line[8:10])[0]
    # print("MSG length:", novatel_header.msg_length)
    # print("Seq:", struct.unpack('<H', line[10:12])[0])
    # print("idle time:", int(line[12]))
    # print("time status:", int(line[13]))
    # print("GPS week:", novatel_header.gps_week)
    # print("GPS sec:", novatel_header.gps_sec)
    # print("Receiver Status:", struct.unpack('<l', line[20:24])[0])
    # print("Reserved:", struct.unpack('<H', line[24:26])[0])
    # print("Receiver S/W Ver:", struct.unpack('<H', line[26:28])[0])
    return novatel_header

def convert_inspvax(header,line,format):
    inspvax = NovatelInspvax()
    inspvax.header = header
    if format == 'a':
        asciiline = line.split(",")
        inspvax.ins_status = asciiline[0]
        inspvax.pos_type = asciiline[1]
        inspvax.lat = float(asciiline[2])
        inspvax.lon = float(asciiline[3])
        inspvax.height = float(asciiline[4])
        inspvax.undulation = asciiline[5]
        inspvax.v_north = asciiline[6]
        inspvax.v_east = asciiline[7]
        inspvax.v_up = asciiline[8]
        inspvax.roll = asciiline[9]
        inspvax.pitch = asciiline[10]
        inspvax.azimuth = asciiline[11]
        inspvax.lat_std = asciiline[12]
        inspvax.lon_std = asciiline[13]
        inspvax.height_std = asciiline[14]
        inspvax.roll_std = asciiline[15]
        inspvax.pitch_std = asciiline[16]
        inspvax.azimuth_std = asciiline[17]
        inspvax.ext_sol_status = asciiline[18]
        inspvax.time_since_update = asciiline[19]
    elif format == 'b':
        print("I don't know")

    return inspvax


def convert_bestpos(header, line):
    # print("bestpos binary data:", len(line))
    bestpos = Header.NovatelBestpos()
    bestpos.header = header
    bestpos.sol_status, bestpos.pos_type, bestpos.lat \
        , bestpos.lon, bestpos.height, bestpos.undulation, \
    bestpos.datum_id, bestpos.lat_std, bestpos.lon_std, \
    bestpos.height_std, bestpos.base_id, bestpos.diff_age, \
    bestpos.sol_age, bestpos.sat_tracked_num, bestpos.sat_used_num, \
    bestpos.sat_used_L1, bestpos.sat_used_L1L2, bestpos.reserved, \
    bestpos.ext_sol_status, bestpos.gal_beidou_mask, bestpos.gps_glo_mask, \
        = struct.unpack('<lldddflffflffBBBBBBBB', line[28:100])
    # print("bestpos.sol_stauts,bestpos.pos_type,bestpos.lat\
    #     ,bestpos.lon,bestpos.height,bestpos.undulation,\
    #     bestpos.datum_id,bestpos.lat_std,bestpos.lon_std,\
    #     bestpos.height_std,bestpos.base_id,bestpos.diff_age,\
    #     bestpos.sol_age,bestpos.sat_tracked_num,bestpos.sat_used_num,\
    #     bestpos.sat_used_L1,bestpos.sat_used_L1L2,bestpos.reserved,\
    #     bestpos.ext_sol_status,bestpos.gal_beidou_mask,bestpos.gps_glo_mask,", \
    #       bestpos.sol_stauts, bestpos.pos_type, bestpos.lat \
    #       , bestpos.lon, bestpos.height, bestpos.undulation, \
    #       bestpos.datum_id, bestpos.lat_std, bestpos.lon_std, \
    #       bestpos.height_std, bestpos.base_id, bestpos.diff_age, \
    #       bestpos.sol_age, bestpos.sat_tracked_num, bestpos.sat_used_num, \
    #       bestpos.sat_used_L1, bestpos.sat_used_L1L2, bestpos.reserved, \
    #       bestpos.ext_sol_status, bestpos.gal_beidou_mask, hex(bestpos.gps_glo_mask), \
    #       )
    return bestpos

    # = struct.unpack('<d', line[36:44])[0]

    # print("lat:", bestpos.lat,"lat_std:",bestpos.lat_std,"sat_num:",bestpos.sat_tracked_num)


def novatelplot(bestposlist,inspvaxlist):
    y_lon = []
    y_lat = []
    x_gpssec = []
    index_i = 0
    for j in range(len(bestposlist)):
        for i in range(index_i,len(inspvaxlist)):
            b_sec = bestposlist[j].header.gps_sec/1000
            if b_sec == inspvaxlist[i].header.gps_sec :
                p1 = CoordinateTransformation.Point()
                p1.BLH(bestposlist[j].lat,bestposlist[j].lon,bestposlist[j].height)
                p2 = CoordinateTransformation.Point()
                p2.BLH(inspvaxlist[i].lat,inspvaxlist[i].lon,inspvaxlist[i].height)
                y_lat.append(p1.X - p2.X)
                y_lon.append(p1.Y - p2.Y)
                x_gpssec.append(b_sec)
                index_i = i
        # print("week:",bestposlist[j].header.gps_week,"time:",,"lat:",bestposlist[j].lat)

    ylon = np.array(y_lon)
    ylat = np.array(y_lat)
    x = np.array(x_gpssec)
    plt.title("Matplotlib demo")
    plt.xlabel("x axis caption")
    plt.ylabel("y axis caption")
    #plt.xticks(y_gpssec)
    #plt.yticks(np.linspace(-1,10,20,endpoint=False))
    plt.plot(x,ylon,'rp')
    plt.plot(x,ylat,'bp')
    plt.show()


if __name__ == '__main__':
    convert_success_count = 0
    convert_unsuccess_count = 0
    # novatel_header = NovatelHeader()
    #fo = open("D:\program\python\python-project\project1\data\BMAW17020017K_2_20191025_083848.BIN", mode='rb')
    fo = open("D:\program\python\python-project\project1\data\cutecom(1).log", mode='rb')

    line_i = 0
    bestposlist = []
    while True:
        i = 0
        fline = fo.read(i + 2048)
        fline_len = len(fline)
        # print("--------------line len = ", fline_len, "-----line_i= ", line_i)
        if fline_len == 0:
            break
        while i < fline_len:
            print(hex(fline[i+0]))
            if fline[i + 0] == 0xaa and fline[i + 1] == 0x44 and fline[i + 2] == 0x12:
                novatel_header = convert_novatel_header(fline[i:i + 28],'b')
                if novatel_header.msg_id == 42:
                    msg_length = novatel_header.msg_length + novatel_header.header_length + 4
                    if msg_length != 104:
                        i += 3
                        convert_unsuccess_count += 3
                        break
                    # print("line len = ",msg_length)
                    if i + msg_length > fline_len:
                        break
                    else:
                        if checkCRC(fline[i:i + msg_length],'b'):
                            convert_success_count += 1
                            bestposlist.append(convert_bestpos(novatel_header, fline[i:i + msg_length]))
                            i += msg_length
                        else:
                            i += 3
                            convert_unsuccess_count += 3

            else:
                i += 1
                convert_unsuccess_count += 1
        line_i += i
        fo.seek(line_i)
    print("convert_success_count = ", convert_success_count, "convert_unsuccess_count = ", convert_unsuccess_count)
    fo.close()
    fo2 = open("D:\program\python\python-project\project1\data\BMAW17020017K_2_20191025_083848-INSPVAX_ASC.ASC",mode='r',newline='\r\n')
    inspvaxlist = []
    while True:
        fo2line = fo2.readline().strip()
        if fo2line == "":
            break
        if checkCRC(fo2line,'a'):
            fline2 = fo2line.split(";")
            novatel_header = convert_novatel_header(fline2[0],'a')
            inspvaxlist.append(convert_inspvax(novatel_header,fline2[1],'a'))

    #for i in range(len(inspvaxlist)):
        #print(inspvaxlist[i].header.gps_sec,inspvaxlist[i].lat,inspvaxlist[i].lon)
    novatelplot(bestposlist,inspvaxlist)
    fo2.close()


# 访问 https://www.jetbrains.com/help/pycharm/ 获取 PyCharm 帮助
