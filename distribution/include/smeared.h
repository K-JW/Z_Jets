/*
 * FileName: smeared.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-23 17:49:05
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-24 14:20:17
 * Description: smearing
 */

#include <iostream>
#include <random>
#include <cmath>

namespace iHepTools {

    struct CSN {
        double C;
        double S;
        double N;
    };

    inline double Variance(double primitive_value, CSN mCSN) {
        return sqrt(
            pow(mCSN.C, 2) + (
                pow(mCSN.S, 2) / primitive_value
            ) + (
                pow(mCSN.N, 2) / pow(primitive_value, 2)
            )
        );
    }
    
    // 传入一个随机数生成引擎作为第一个参数
    template<class RNG>
    double GaussSmeared(RNG &gen, const double &primitive_value, const double &variance) {
        static_assert(
            std::is_same<RNG, std::default_random_engine>::value || 
            std::is_same<RNG, std::minstd_rand0>::value || 
            std::is_same<RNG, std::minstd_rand>::value ||
            std::is_same<RNG, std::mt19937>::value ||
            std::is_same<RNG, std::mt19937_64>::value ||
            std::is_same<RNG, std::ranlux24_base>::value ||
            std::is_same<RNG, std::ranlux48_base>::value ||
            std::is_same<RNG, std::ranlux24>::value ||
            std::is_same<RNG, std::ranlux48>::value ||
            std::is_same<RNG, std::knuth_b>::value, 
            "RNG is not a supported random engine!"
        );
        double gauss_smeared_value;
        std::normal_distribution<double> norm(1.0, variance);
        gauss_smeared_value = primitive_value * norm(gen);
        return gauss_smeared_value;
    }
}