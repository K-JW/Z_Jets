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
LastEditTime: 2019-10-17 15:19:26
Description: This python script is used to detect if the hepmc file is complete.
'''
import sys
import linecache as lc

def isComplete(file_name):
    head1 = (lc.getline(file_name, 2).strip().split()[0] == "HepMC::Version")
    head2 = (lc.getline(file_name, 3).strip() == "HepMC::IO_GenEvent-START_EVENT_LISTING")
    file = open(file_name, 'r')
    end1 = ( file.readlines()[-1].strip() == "HepMC::IO_GenEvent-END_EVENT_LISTING" ) 
    file.close()
    return ( head1 and head2 and end1 )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("    argv error, plese input file name")
        exit
    else:
        for i in range(1, len(sys.argv)):
            if not isComplete(sys.argv[i]):
                print("--ERROR: file ", sys.argv[i], " is not complete!")
            else:
                print("file ", sys.argv[i], " is complete!")
