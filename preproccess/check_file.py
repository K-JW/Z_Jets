#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
FileName: check_file.py
Version: 1.0
License: MIT License

Author: KANG Jin-Wen
E-Mail: kangjinwen@vip.qq.com
Date: 2019-10-17 14:28:11

LastEditors: KANG Jin-Wen
LastEditTime: 2019-11-03 15:12:45
Description: This python script is used to detect if the hepmc file is complete.
'''
import sys, os
import linecache as lc

# 检测 HepMC 文件是否头尾声明完整
def isComplete(file_name):
    head1 = (lc.getline(file_name, 2).strip().split()[0] == "HepMC::Version")
    head2 = (lc.getline(file_name, 3).strip() == "HepMC::IO_GenEvent-START_EVENT_LISTING")
    file = open(file_name, 'r')
    end1 = ( file.readlines()[-1].strip() == "HepMC::IO_GenEvent-END_EVENT_LISTING" ) 
    file.close()
    print(head1)
    print(head2)
    print(end1)
    return ( head1 and head2 and end1 )

# 倒序读行，行数为负数，-1 为倒数第一行
def readline_reverse(file_name, line_num):
    fsize = os.path.getsize(file_name)
    file = open(file_name, 'rb')
    file.seek(fsize - 1)
    assert type(line_num) == int
    if line_num > 0:
        raise IndexError("line number must less than 0")
    cur_pos = file.tell()
    buf = b''
    need_enter_num = abs(line_num) + 1
    get_enter_num = 0
    while get_enter_num < need_enter_num:
        b = file.read(1)
        if b == '\n':
            get_enter_num += 1
        buf = b + buf
        cur_pos -= 1
        if cur_pos < 0: break
        file.seek(cur_pos)
    file.close()
    return buf.split('\n')[1]

# 读取 HepMC 文件内的事件数，默认事件数从0计数，且顺序增一
def getEvent(file_name):
    line_num = -1
    while True:
        line_num -= 1
        line = readline_reverse(file_name, line_num).split(' ')
        if len(line) > 1 and line[0] == "E":
            return int(line[1]) + 1
        if line[0] == "HepMC::IO_GenEvent-START_EVENT_LISTING":
            raise UserWarning("File error, can't find event number")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("    argv error, plese input file name")
        exit
    else:
        event_num = 0
        for i in range(1, len(sys.argv)):
            if not isComplete(sys.argv[i]):
                print("--\033[1;31m ERROR:\033[0m file " + str(sys.argv[i]) + " is not complete!")
            else:
                event_num += getEvent(sys.argv[i])
                print("\033[0;32m Info:\033[0m file " + str(sys.argv[i]) + " is complete!")
        
        print('\n')
        print("\033[0;32m INFO:\033[0m \033[0;33mValid event number is\033[0m " + "\033[1;32m" + str(event_num) + "\033[0m")