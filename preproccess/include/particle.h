/*
 * FileName: particle.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 16:51:18
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-15 16:57:23
 * Description: Define a struct that include some information of particle
 */

#ifndef IHEPTOOLS_PARTILCE_H
#define IHEPTOOLS_PARTILCE_H

namespace iHepTools {
    struct Particle {
        int pdg_code;
        double x0;
        double x1;
        double x2; 
        double x3;
        double mass;
    };
}

#endif // IHEPTOOLS_PARTILCE_H