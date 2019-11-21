/*
 * FileName: detector.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-18 12:16:43
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-18 13:42:27
 * Description: Discern the type of particle, use for pdg code.
 */

#ifndef IHEPTOOLS_DETECTOR_H
#define IHEPTOOLS_DETECTOR_H

namespace iHepTools {

inline bool isUp (const int &pdg_id ){
    return abs(pdg_id) == 2;
}

inline bool isDown (const int &pdg_id ){
    return abs(pdg_id) == 1;
}

inline bool isCharm (const int &pdg_id ){
    return abs(pdg_id) == 4;
}

inline bool isStrange (const int &pdg_id ){
    return abs(pdg_id) == 3;
}

inline bool isTop (const int &pdg_id ){
    return abs(pdg_id) == 6;
}

inline bool isBottom (const int &pdg_id ){
    return abs(pdg_id) == 5;
}

inline bool isElectron (const int &pdg_id ){
    return abs(pdg_id) == 11;
}

inline bool isMuon (const int &pdg_id ){
    return abs(pdg_id) == 13;
}

inline bool isTau (const int &pdg_id ){
    return abs(pdg_id) == 15;
}

inline bool isLepton( const int &pdg_id ) {
    return abs(pdg_id) == 11 || abs(pdg_id) == 13 || abs(pdg_id) == 15;
}

inline bool isNeutrino(const int &pdg_id){
    return abs(pdg_id) == 12 || abs(pdg_id) == 14 || abs(pdg_id) == 16;
}

inline bool isGluon( const int &pdg_id ) {
    return pdg_id == 21;
}

inline bool isPhoton( const int &pdg_id ) {
    return pdg_id == 22;
}

inline bool isZ( const int &pdg_id ) {
    return pdg_id == 23;
}

inline bool isW( const int &pdg_id ) {
    return abs(pdg_id) == 24;
}

inline bool isHiggs( const int &pdg_id ) {
    return pdg_id == 25;
}

//@TODO: Zprime,Wprime,Zprimeprime,H0,A0,H+

inline bool isChargedPion( const int &pdg_id ) {
    return abs(pdg_id) == 211;
}

inline bool isChargino1( const int &pdg_id ) {
    return abs(pdg_id) == 1000024;
}

inline bool isNeutralino2( const int &pdg_id ) {
    return abs(pdg_id) == 1000023;
}

inline bool isNeutralino1( const int &pdg_id ) {
    return abs(pdg_id) == 1000022;
}

inline bool iso1Meson(const int &pdg_id){
    return abs(pdg_id) == 111 ||
           abs(pdg_id) == 211 ||
           abs(pdg_id) == 113 ||
           abs(pdg_id) == 213 ||
           abs(pdg_id) == 115 ||
           abs(pdg_id) == 215 ||
           abs(pdg_id) == 117 ||
           abs(pdg_id) == 217 ||
           abs(pdg_id) == 119 ||
           abs(pdg_id) == 219 ||
           abs(pdg_id) == 9000111 ||
           abs(pdg_id) == 9000211 ||
           abs(pdg_id) == 100111 ||
           abs(pdg_id) == 100211 ||
           abs(pdg_id) == 10111 ||
           abs(pdg_id) == 10211 ||
           abs(pdg_id) == 9010111 ||
           abs(pdg_id) == 9010211 ||
           abs(pdg_id) == 10113 ||
           abs(pdg_id) == 10213 ||
           abs(pdg_id) == 20113 ||
           abs(pdg_id) == 20213 ||
           abs(pdg_id) == 9000113 ||
           abs(pdg_id) == 9000213 ||
           abs(pdg_id) == 9010113 ||
           abs(pdg_id) == 9010213 ||
           abs(pdg_id) == 100113 ||
           abs(pdg_id) == 100213 ||
           abs(pdg_id) == 9010113 ||
           abs(pdg_id) == 9010213 ||
           abs(pdg_id) == 9020113 ||
           abs(pdg_id) == 9020213 ||
           abs(pdg_id) == 30113 ||
           abs(pdg_id) == 30213 ||
           abs(pdg_id) == 9030113 ||
           abs(pdg_id) == 9030213 ||
           abs(pdg_id) == 9040113 ||
           abs(pdg_id) == 9040213 ||
           abs(pdg_id) == 10115 ||
           abs(pdg_id) == 10215 ||
           abs(pdg_id) == 9000115 ||
           abs(pdg_id) == 9000215 ||
           abs(pdg_id) == 9010115 ||
           abs(pdg_id) == 9010215 ||
           abs(pdg_id) == 9000117 ||
           abs(pdg_id) == 9000217 ||
           abs(pdg_id) == 9010117 ||
           abs(pdg_id) == 9010217;
}

inline bool iso0Meson(const int &pdg_id){
    return abs(pdg_id) == 221 ||
           abs(pdg_id) == 331 ||
           abs(pdg_id) == 223 ||
           abs(pdg_id) == 333 ||
           abs(pdg_id) == 225 ||
           abs(pdg_id) == 335 ||
           abs(pdg_id) == 227 ||
           abs(pdg_id) == 337 ||
           abs(pdg_id) == 229 ||
           abs(pdg_id) == 339 ||
           abs(pdg_id) == 339 ||
           abs(pdg_id) == 100221 ||
           abs(pdg_id) == 10221 ||
           abs(pdg_id) == 100331 ||
           abs(pdg_id) == 10331 ||
           abs(pdg_id) == 10223 ||
           abs(pdg_id) == 20223 ||
           abs(pdg_id) == 10333 ||
           abs(pdg_id) == 20333 ||
           abs(pdg_id) == 100223 ||
           abs(pdg_id) == 30223 ||
           abs(pdg_id) == 100333 ||
           abs(pdg_id) == 10225 ||
           abs(pdg_id) == 10335 ||
           abs(pdg_id) == 9000221 ||
           abs(pdg_id) == 9010221 ||
           abs(pdg_id) == 9020221 ||
           abs(pdg_id) == 9030221 ||
           abs(pdg_id) == 9040221 ||
           abs(pdg_id) == 9050221 ||
           abs(pdg_id) == 9060221 ||
           abs(pdg_id) == 9070221 ||
           abs(pdg_id) == 9080221 ||
           abs(pdg_id) == 9000223 ||
           abs(pdg_id) == 9010223 ||
           abs(pdg_id) == 9000225 ||
           abs(pdg_id) == 9010225 ||
           abs(pdg_id) == 9020225 ||
           abs(pdg_id) == 9030225 ||
           abs(pdg_id) == 9040225 ||
           abs(pdg_id) == 9050225 ||
           abs(pdg_id) == 9060225 ||
           abs(pdg_id) == 9070225 ||
           abs(pdg_id) == 9080225 ||
           abs(pdg_id) == 9090225 ||
           abs(pdg_id) == 9000229 ||
           abs(pdg_id) == 9010229;
}

inline bool strangeMeson(const int &pdg_id){
    return  abs(pdg_id) == 130 ||
            abs(pdg_id) == 310 ||
            abs(pdg_id) == 311 ||
            abs(pdg_id) == 321 ||
            abs(pdg_id) == 313 ||
            abs(pdg_id) == 323 ||
            abs(pdg_id) == 315 ||
            abs(pdg_id) == 325 ||
            abs(pdg_id) == 317 ||
            abs(pdg_id) == 327 ||
            abs(pdg_id) == 319 ||
            abs(pdg_id) == 329 ||
            abs(pdg_id) == 9000311 ||
            abs(pdg_id) == 9000321 ||
            abs(pdg_id) == 10311 ||
            abs(pdg_id) == 10321 ||
            abs(pdg_id) == 100311 ||
            abs(pdg_id) == 100321 ||
            abs(pdg_id) == 9010311 ||
            abs(pdg_id) == 9010321 ||
            abs(pdg_id) == 9020311 ||
            abs(pdg_id) == 9020321 ||
            abs(pdg_id) == 10313 ||
            abs(pdg_id) == 10323 ||
            abs(pdg_id) == 20313 ||
            abs(pdg_id) == 20323 ||
            abs(pdg_id) == 100313 ||
            abs(pdg_id) == 100323 ||
            abs(pdg_id) == 9000313 ||
            abs(pdg_id) == 9000323 ||
            abs(pdg_id) == 30313 ||
            abs(pdg_id) == 30323 ||
            abs(pdg_id) == 9000315 ||
            abs(pdg_id) == 9000325 ||
            abs(pdg_id) == 10315 ||
            abs(pdg_id) == 10325 ||
            abs(pdg_id) == 20315 ||
            abs(pdg_id) == 20325 ||
            abs(pdg_id) == 9010315 ||
            abs(pdg_id) == 9010325 ||
            abs(pdg_id) == 9020315 ||
            abs(pdg_id) == 9020325 ||
            abs(pdg_id) == 9010317 ||
            abs(pdg_id) == 9010327 ||
            abs(pdg_id) == 9010319 ||
            abs(pdg_id) == 9010329;
}

inline bool charmedMeson(const int &pdg_id){
    return  abs(pdg_id) == 411 ||
            abs(pdg_id) == 421 ||
            abs(pdg_id) == 413 ||
            abs(pdg_id) == 423 ||
            abs(pdg_id) == 415 ||
            abs(pdg_id) == 425 ||
            abs(pdg_id) == 431 ||
            abs(pdg_id) == 433 ||
            abs(pdg_id) == 435 ||
            abs(pdg_id) == 10411 ||
            abs(pdg_id) == 10421 ||
            abs(pdg_id) == 10413 ||
            abs(pdg_id) == 10423 ||
            abs(pdg_id) == 20413 ||
            abs(pdg_id) == 20423 ||
            abs(pdg_id) == 10431 ||
            abs(pdg_id) == 10433 ||
            abs(pdg_id) == 20433;
}

inline bool bottomMeson(const int &pdg_id){
    return abs(pdg_id) == 511 ||
           abs(pdg_id) == 521 ||
           abs(pdg_id) == 513 ||
           abs(pdg_id) == 523 ||
           abs(pdg_id) == 515 ||
           abs(pdg_id) == 525 ||
           abs(pdg_id) == 531 ||
           abs(pdg_id) == 533 ||
           abs(pdg_id) == 535 ||
           abs(pdg_id) == 541 ||
           abs(pdg_id) == 543 ||
           abs(pdg_id) == 545 ||
           abs(pdg_id) == 10511 ||
           abs(pdg_id) == 10521 ||
           abs(pdg_id) == 10513 ||
           abs(pdg_id) == 10523 ||
           abs(pdg_id) == 20513 ||
           abs(pdg_id) == 20523 ||
           abs(pdg_id) == 10531 ||
           abs(pdg_id) == 10533 ||
           abs(pdg_id) == 20533 ||
           abs(pdg_id) == 10541 ||
           abs(pdg_id) == 10543 ||
           abs(pdg_id) == 20543;
}

inline bool ccMeson(const int &pdg_id){
    return abs(pdg_id) == 441 ||
           abs(pdg_id) == 443 ||
           abs(pdg_id) == 445 ||
           abs(pdg_id) == 10441 ||
           abs(pdg_id) == 100441 ||
           abs(pdg_id) == 10443 ||
           abs(pdg_id) == 100443 ||
           abs(pdg_id) == 20443 ||
           abs(pdg_id) == 30443 ||
           abs(pdg_id) == 9000443 ||
           abs(pdg_id) == 9010443 ||
           abs(pdg_id) == 9020443 ||
           abs(pdg_id) == 100445;
}

inline bool bbMeson(const int &pdg_id){
    return abs(pdg_id) == 551 ||
           abs(pdg_id) == 553 ||
           abs(pdg_id) == 555 ||
           abs(pdg_id) == 557 ||
           abs(pdg_id) == 10551 ||
           abs(pdg_id) == 100551 ||
           abs(pdg_id) == 110551 ||
           abs(pdg_id) == 200551 ||
           abs(pdg_id) == 210551 ||
           abs(pdg_id) == 10553 ||
           abs(pdg_id) == 20553 ||
           abs(pdg_id) == 30553 ||
           abs(pdg_id) == 100553 ||
           abs(pdg_id) == 110553 ||
           abs(pdg_id) == 120553 ||
           abs(pdg_id) == 130553 ||
           abs(pdg_id) == 200553 ||
           abs(pdg_id) == 210553 ||
           abs(pdg_id) == 220553 ||
           abs(pdg_id) == 300553 ||
           abs(pdg_id) == 9000553 ||
           abs(pdg_id) == 9010553 ||
           abs(pdg_id) == 10555 ||
           abs(pdg_id) == 20555 ||
           abs(pdg_id) == 100555 ||
           abs(pdg_id) == 110555 ||
           abs(pdg_id) == 120555 ||
           abs(pdg_id) == 200555 ||
           abs(pdg_id) == 100557;
}

inline bool lightBaryon(const int &pdg_id){
    return abs(pdg_id) == 2212 ||
           abs(pdg_id) == 2112 ||
           abs(pdg_id) == 2224 ||
           abs(pdg_id) == 2214 ||
           abs(pdg_id) == 2114 ||
           abs(pdg_id) == 1114;
}

inline bool strangeBaryon(const int &pdg_id){
    return abs(pdg_id) >= 3000 && abs(pdg_id) < 4000;
}

inline bool charmedBaryon(const int &pdg_id){
    return abs(pdg_id) >= 4000 && abs(pdg_id) < 5000;
}

inline bool bottomBaryon(const int &pdg_id){
    return abs(pdg_id) >= 5000 && abs(pdg_id) < 6000;
}

}

#endif // iHEPTOOLS_DETECTOR_H