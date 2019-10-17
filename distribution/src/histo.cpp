/*
 * FileName: histo.cpp
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-16 09:46:45
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-16 15:53:54
 * Description: Define histo class
 */

#include <iostream>
#include <cmath>
#include "histo.h"

namespace iHepTools {

bool isAtRegion(const Region &mRegin, const double &value) {
    if (value >= mRegin.leftValue && value < mRegin.rightValue) {
        return true;
    } else {
        return false;
    }
}

Histo::Histo(double xmin, double xmax, double binWidth) {
    const static int vec_size = std::floor((xmax - xmin) / binWidth) + 1;
    mDistInfoVec.resize(vec_size);
    mDistInfoVec[0].region.leftValue = xmin;
    mDistInfoVec[0].region.rightValue = xmin + binWidth;
    mDistInfoVec[0].distValue = 0.;
    for (size_t i = 1; i < vec_size; i++) {
        mDistInfoVec[i].region.leftValue = mDistInfoVec[i - 1].region.rightValue;
        mDistInfoVec[i].region.rightValue = mDistInfoVec[i].region.leftValue + binWidth;
        mDistInfoVec[i].region.middleValue = (
            mDistInfoVec[i].region.leftValue + mDistInfoVec[i].region.rightValue
        ) / 2.0;
        mDistInfoVec[i].distValue = 0.;
    }
    mNormalisationValue = 0.;
    if (mDistInfoVec.back().region.rightValue > xmax) 
        mDistInfoVec.back().region.rightValue = xmax;
}

Histo::Histo(const vector<double> &pointList) {
    const static int vec_size = pointList.size() - 1;
    mDistInfoVec.resize(vec_size);
    for (size_t i = 0; i < vec_size; i++) {
        mDistInfoVec[i].region.leftValue = pointList[i];
        mDistInfoVec[i].region.rightValue = pointList[i + 1];
        mDistInfoVec[i].region.middleValue = (
            mDistInfoVec[i].region.leftValue + mDistInfoVec[i].region.rightValue
        ) / 2.0;
        mDistInfoVec[i].distValue = 0.;
    }
    mNormalisationValue = 0.;
}

void Histo::clear() {
    for (auto mDistInfo : mDistInfoVec) {
        mDistInfo.distValue = 0.0;
    }
    mNormalisationValue = 0.0;
}

int Histo::findAtRegion(double value) {
    for (size_t i = 0; i < mDistInfoVec.size(); i++) {
        if (isAtRegion(mDistInfoVec[i].region, value))
            return i;
    }
    if (fabs(value - mDistInfoVec.back().region.rightValue) <= 1.0E-06) {
        return mDistInfoVec.size() - 1;
    } else {
        std::cout << value << " don't belong to any region." << std::endl;
    }
}

void Histo::addEventNum(const double &value, const double &weight) {
    int region_index = findAtRegion(value);
    mDistInfoVec[region_index].distValue += weight;
}

// iHepTools::Histo iHepTools::Histo::operator+(const iHepTools::Histo &mHisto) {
//     for (auto mDistVecComponent : mHisto.mDistInfoVec) {
        
//     }
// }

vector<distInfo> Histo::getHisto() {
    vector<iHepTools::distInfo> temp = mDistInfoVec;
    for (size_t i = 0; i < mDistInfoVec.size(); i++) {
        temp[i].distValue = mDistInfoVec[i].distValue / (
            mDistInfoVec[i].region.rightValue - mDistInfoVec[i].region.leftValue
        ) / mNormalisationValue;
    }
    return temp;
}

}