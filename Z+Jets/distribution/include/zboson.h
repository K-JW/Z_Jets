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
 * LastEditTime: 2019-11-03 11:02:21
 * Description: file content
 */

#ifndef IHEPTOOLS_ZBOSON_H
#define IHEPTOOLS_ZBOSON_H


#include <vector>

#include <fastjet/PseudoJet.hh>
// #include <HepMC/IO_GenEvent.h>
// #include <HepMC/GenEvent.h>

using namespace std;
using namespace fastjet;

namespace iHepTools {
    //
    // return bool type and rebuild Z boson
    bool isGoodZBoson(const vector<PseudoJet> &mZ_daughters, PseudoJet &ZBoson, const double &ZBosonPtMin = 60.0) {
        if (mZ_daughters.size() != 2) {
            ZBoson.reset(0., 0., 0., 0.);
            return false;
        } else if (abs(mZ_daughters[0].user_index()) == 11 && 
                abs(mZ_daughters[1].user_index()) == 11 ) {
                //
                if (    (fabs(mZ_daughters[0].eta()) < 1.44 || (
                        fabs(mZ_daughters[0].eta()) > 1.57 &&
                        fabs(mZ_daughters[0].eta()) < 2.5)  // daughter0
                        ) && (
                            fabs(mZ_daughters[1].eta()) < 1.44 || (
                                fabs(mZ_daughters[1].eta()) > 1.57 &&
                                fabs(mZ_daughters[1].eta()) < 2.5
                            )   // daughter1
                        )
                    ) {
                    ZBoson = mZ_daughters[0] + mZ_daughters[1];
                    if (mZ_daughters[0].pt() > 20.0 && mZ_daughters[1].pt() > 20.0 && 
                        ZBoson.m() > 70.0 && ZBoson.m() < 110 && ZBoson.pt() > ZBosonPtMin && 
                        fabs(ZBoson.rap()) < 2.5) {
                        return true;
                    } else {
                        ZBoson.reset(0., 0., 0., 0.);
                        return false;
                    }
                } else {
                    ZBoson.reset(0., 0., 0., 0);
                    return false;
                }
        } else if (abs(mZ_daughters[0].user_index()) == 13 && 
                    abs(mZ_daughters[1].user_index()) == 13 ) {
            //
            if (fabs(mZ_daughters[0].eta()) < 2.4 && fabs(mZ_daughters[1].eta()) < 2.4 && 
                mZ_daughters[0].pt() > 10.0 && mZ_daughters[1].pt() > 10.0) {
                //
                ZBoson = mZ_daughters[0] + mZ_daughters[1];
                if (ZBoson.m() > 70 && ZBoson.m() < 110 && ZBoson.pt() > ZBosonPtMin && 
                    fabs(ZBoson.rap()) < 2.5) {
                    return true;
                } else {
                    ZBoson.reset(0., 0., 0., 0.);
                    return false;
                }
            } else {
                ZBoson.reset(0., 0., 0., 0.);
                return false;
            }
        } else {
            ZBoson.reset(0., 0., 0., 0.);
            return false;
        }
    }

    bool isGoodZBoson(const vector<PseudoJet> &mZ_daughters, PseudoJet &ZBoson, const double &ZBosonPtMin, const double &ZBosonPtMax) {
        PseudoJet tmpZBoson;
        bool tmpBool = isGoodZBoson(mZ_daughters, tmpZBoson, ZBosonPtMin);
        if (tmpBool && tmpZBoson.pt() < ZBosonPtMax) {
            ZBoson = tmpZBoson;
            return true;
        } else {
            ZBoson.reset(0., 0., 0., 0.);
            return false;
        }
    }
}

#endif // IHEPTOOLS_ZBOSON_H