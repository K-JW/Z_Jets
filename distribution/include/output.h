/*
 * FileName: output.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-26 21:24:24
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-26 21:53:52
 * Description: Define method use for output data.
 */

#ifndef IHEPTOOLS_OUTPUT_H
#define IHEPTOOLS_OUTPUT_H

#include <string>
#include <fstream>
#include <iomanip>

#include "histo.h"

using namespace std;

namespace iHepTools {
    
    void WriteDataToText(const string &file_name, const vector<distInfo> &mHistInfo, 
        const string &comments = "") {
        //
        ofstream file;
        file.open(file_name);
        file << comments;
        file << "# xlow\txhigh\txmiddle\tval" << '\n';
        for (size_t i = 0; i < mHistInfo.size(); i++) {
            file << setw(12) << scientific << setprecision(6)
                << mHistInfo[i].region.leftValue << '\t' 
                << mHistInfo[i].region.rightValue << '\t'
                << mHistInfo[i].region.middleValue << '\t'
                << mHistInfo[i].distValue << '\n';
        }
        file.close();
    }
}

#endif // IHEPTOOLS_OUTPUT_H