/*
 * FileName: zboson.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-16 14:42:56
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-16 15:41:53
 * Description: file content
 */

#ifndef IHEPTOOLS_ZBOSON_H
#define IHEPTOOLS_ZBOSON_H

#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

namespace iHepTools {
    //
    // returns true if the GenParticle does not decay
    inline bool isFinal( const HepMC::GenParticle* p ) {
        return !p->end_vertex() && p->status()==1;
    }
}

#endif // IHEPTOOLS_ZBOSON_H