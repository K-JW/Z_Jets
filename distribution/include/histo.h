/*
 * FileName: histo.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 21:50:05
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-16 11:53:03
 * Description: Class histogram
 */

#ifndef IHEPTOOLS_HISTO_H
#define IHEPTOOLS_HISTO_H

#include <vector>

using namespace std;

namespace iHepTools {

    struct Region {
        double leftValue;
        double rightValue;
        double middleValue;
    };

    struct distInfo {
        Region region;
        double distValue;
    };

    bool isAtRegion(const Region &mRegin, const double &value);
    
    class Histo {
        public:
            Histo(double xmin, double xmax, double binWidth);
            Histo(const vector<double> &pointList);
            Histo(const vector<distInfo> &distInfoVec) : mDistInfoVec(distInfoVec) {
                for (size_t i = 0; i < mDistInfoVec.size(); i++) {
                    mDistInfoVec[i].region.middleValue = (
                        mDistInfoVec[i].region.leftValue + 
                        mDistInfoVec[i].region.rightValue
                    ) / 2.0;
                    mDistInfoVec[i].distValue = 0.;
                }
                mNormalisationValue = 0.;
            }
            
            void clear();
            void addEventNum(const double &value, const double &weight);
            inline void addEventNorm(double weight) {
                mNormalisationValue += weight;
            }

            // // 重载加法运算符
            // Histo operator+(const Histo &mHisto);

            // 返回值列表
            vector<distInfo> getHisto();
        
        private:
            vector<distInfo> mDistInfoVec;
            double mNormalisationValue;

            int findAtRegion(double value);
    };
}


#endif // IHEPTOOLS_HISTO_H