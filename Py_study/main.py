# 这是一个示例 Python 脚本。

# 按 Shift+F10 执行或将其替换为您的代码。
# 按 双击 Shift 在所有地方搜索类、文件、工具窗口、操作和设置。
import struct
import binascii
# import numpy as np
# from matplotlib import pyplot as plt
#
# x = np.arange(1, 11)
# y = 2 * x + 5
# plt.title("Matplotlib demo")
# plt.xlabel("x axis caption")
# plt.ylabel("y axis caption")
# plt.plot(x, y)
# plt.show()
def print_hi(name):
    # 在下面的代码行中使用断点来调试脚本。
    print(f"Hi, {name}")  # 按 Ctrl+F8 切换断点。
'''
8*2+2=18*4=72
00000000 10000000 1B0450B3 F28E4940 16FA6BBE 7C825CC0 0060769F 449F9040
A62A82C1 3D000000 125ACB3F CD9E983F DB664040 00303030 00000000 00000000
0B0B0000 00060003

1B0450B3 F28E4940
0001 1011 0000 0100 0101 0000 1011 0011 
'''
b=bytearray()
# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    #hex1[8] = (0x1B,0x04,0x50,0xB3,0xF2,0x8E,0x49,0x40)
    #b = binascii.b2a_hex(int(0x1B0450B3F28E4940))
    #print(b)
    b = b'aa4412'
    b2 = b'B0450B3'
    print(b[0])
    float1 = 51.11635910984
    b2 = binascii.b2a_hex(struct.pack('d',float1))
    print("b2 len = ", hex(b2[15]))
    print(struct.pack('d',float1))
    print(float1)
    print(b2)
    print((b2[0]))
    print(struct.unpack('d',binascii.a2b_hex(b2))[0])
    print("----2")
    a = 0
    print("a>>8 = ",a>>8& 0x00FFFFFF)
# 访问 https://www.jetbrains.com/help/pycharm/ 获取 PyCharm 帮助
