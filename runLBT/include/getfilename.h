/*
 * FileName: getfilename.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-18 09:39:55
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-18 09:41:37
 * Description: Get file name.
 */

#include <string>

namespace iHepTools {
    std::string get_file_name(std::string filepath) {
        if (!filepath.empty()) {
            int location_point = filepath.find_last_of('.');
            int location_filename = filepath.find_last_of('/');
            std::string file_name = filepath.substr(location_filename + 1, 
                location_point - location_filename - 1);
            return file_name;
        } else {
            return "NULL";
        }
    }
}