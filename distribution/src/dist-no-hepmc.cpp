/*
 * FileName: distribution.cpp
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 21:00:36
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-29 12:56:51
 * Description: Calculate distribution.
 */


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include "cmdline/cmdline.h"

#include "jetSelector.h"
#include "histo.h"
#include "zboson.h"
#include "smeared.h"
#include "output.h"

using namespace std;
using namespace iHepTools;

// random number generator
random_device mRD;
default_random_engine mRNG(mRD());

constexpr double PI = acos(-1);
constexpr double R_jet = 0.3;

// Jet definition
JetDefinition jet_def(antikt_algorithm, R_jet);
Selector select_akt = SelectorAbsEtaMax(1.6) && SelectorPtMin(20.); // 因为要做 smeared, 所以这里要取小一点

int main(int argc, char *argv[]) {
    
    // 参数解析
    cmdline::parser args;
    args.add<string>("delta-phi", 'd', "assigned delta_phi output file name", 
        false, "delta_phi.dat");
    args.add<string>("x-jz", 'x', "assigned x_jZ output file name", 
        false, "x_jZ.dat");
    args.add<string>("mean-x-jz", 'm', "assigned <x_jZ> output file name", 
        false, "mean-x_jZ.dat");
    args.add<string>("R_jZ", 'r', "assigned R_jZ output file name", 
        false, "R_jZ.dat");
    args.add("enable-smeared", 's', "execute smearing for data (default: false)");
    args.add<double>("C-CSN", 'C', "given C's value of CSN", false, 0.061);
    args.add<double>("S-CSN", 'S', "given S's value of CSN", false, 0.95);
    args.add<double>("N-CSN", 'N', "given N's value of CSN", false, 0.001);
    args.add("help", 'h', "print help message");
    args.footer("filename ...");
    args.set_program_name("dist");

    bool isGetArgs = args.parse(argc, argv);
    
    if (argc == 1 || args.exist("help")) {
        cerr << args.usage();
        exit(EXIT_FAILURE);
    }

    if (!isGetArgs) {
        cerr << args.error() << endl << args.usage();
        exit(EXIT_FAILURE);
    }
    
    vector<string> file_name_vec;
    for (size_t i = 0; i < args.rest().size(); i++) {
        file_name_vec.emplace_back(args.rest()[i]);
    }
    bool isSmeared = args.exist("enable-smeared");
    CSN mCSN = {
        args.get<double>("C-CSN"), args.get<double>("S-CSN"), args.get<double>("N-CSN")
    };
    // 

    // 定义 phi histo
    vector<double> mPhiHistoPointList = 
        {0, 0.62, 1.26, 1.74, 2.26, 2.54, 2.86, 2.94, PI};
    Histo mPhiHisto(mPhiHistoPointList);

    // 定义 x_jZ histo
    vector<double> mXjZHistoPointList = 
        {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    Histo mXjZHisto(mXjZHistoPointList);

    // 定义 <x_jZ> histo
    vector<double> mMeanXjZHistoPointList = {
        40.0, 50.0, 60.0, 80.0, 120.0
    };
    Histo mMeanXjzHisto(mMeanXjZHistoPointList);
    // 新建一个 Histogram ，保存不同 p_T 下的 N_Z ，用以求 <x_jZ>
    Histo mMeanXjzDenoHisto(mMeanXjZHistoPointList);

    // 定义 R_jZ histo
    vector<double> mRjZHistoPointList = {
        40.0, 50.0, 60.0, 80.0, 120.0
    };
    Histo mRjZHisto(mRjZHistoPointList);
    // 新建一个 Histogram ，保存不同 p_T 下的 N_Z ，用以求 R_jZ
    Histo mRjZDenoHisto(mRjZHistoPointList);
    
    for (size_t i = 0; i < file_name_vec.size(); i++) {

        ifstream in_file(file_name_vec[i]);
        if (!in_file.is_open()) {
            cout << "  \033[1;31mERROR:\033[0m Can't open file!\n";
            exit(EXIT_FAILURE);
        }

        while ( !in_file.eof() ) {
            
            vector<fastjet::PseudoJet> pseduo_jets;
            vector<fastjet::PseudoJet> Z_daughters;
            fastjet::PseudoJet ZBoson;

            int event_id, pdgcode, partons_num, ZDaugCode;
            double ZD1X, ZD1Y, ZD1Z, ZD1E, ZD2X, ZD2Y, ZD2Z, ZD2E;
            double px, py, pz, energy, temp;
            double weight0, weight1, weight2, weight3, weight4;
            
            in_file >> event_id >> partons_num >> weight0 >> weight1 
                >> weight2 >> weight3; // >> weight4;
            in_file >> ZDaugCode >> ZD1X >> ZD1Y >> ZD1Z >> ZD1E >> temp;
            fastjet::PseudoJet tmpZD1(ZD1X, ZD1Y, ZD1Z, ZD1E);
            tmpZD1.set_user_index(ZDaugCode);
            Z_daughters.emplace_back(tmpZD1);
            in_file >> ZDaugCode >> ZD2X >> ZD2Y >> ZD2Z >> ZD2E >> temp;
            fastjet::PseudoJet tmpZD2(ZD2X, ZD2Y, ZD2Z, ZD2E);
            tmpZD2.set_user_index(ZDaugCode);
            Z_daughters.emplace_back(tmpZD2);
            //cout << "event_id: " << event_id << "\n";

            for (size_t i = 0; i < partons_num - 2; i++) {
                in_file >> pdgcode >> px >> py >> pz >> energy >> temp;
                fastjet::PseudoJet tmpParton(px, py, pz, energy);
                pseduo_jets.emplace_back(tmpParton);
                // cout << "Parton's code: " << pdgcode << '\n';
            } // for (size_t i = 0; i < partons_num - 2; i++)

            if (isGoodZBoson(Z_daughters, ZBoson, 40.0)) {
                // cout << "Here\n";
                if (ZBoson.pt() > 60.0) {
                    mPhiHisto.addEventNorm( weight0 / weight2 );
                    mXjZHisto.addEventNorm( weight0 / weight2 );
                }
                mRjZDenoHisto.addEventNum(ZBoson.pt(), weight0 / weight2 );
                vector<PseudoJet> jets = SelectJet(pseduo_jets, jet_def, select_akt);
                if (jets.size() > 0) {
                    for (const auto &jet : jets) {
                        double delta_phi = fabs(jet.phi() - ZBoson.phi());
                        delta_phi = delta_phi > PI ? 2 * PI - delta_phi : delta_phi;
                        // calc deltaR
                        double deltaPhi0J = fabs(Z_daughters[0].phi() - jet.phi());
                        deltaPhi0J = deltaPhi0J > PI ? 2 * PI - deltaPhi0J : deltaPhi0J;
                        double deltaPhi1J = fabs(Z_daughters[1].phi() - jet.phi());
                        deltaPhi1J = deltaPhi1J > PI ? 2 * PI - deltaPhi1J : deltaPhi1J;
                        double deltaRap0J = fabs(Z_daughters[0].rap() - jet.rap());
                        double deltaRap1J = fabs(Z_daughters[1].rap() - jet.rap());
                        double deltaR0J = sqrt(deltaPhi0J * deltaPhi0J + deltaRap0J * deltaRap0J);
                        double deltaR1J = sqrt(deltaPhi1J * deltaPhi1J + deltaRap1J * deltaRap1J);
                        if (deltaR0J >=0.4 && deltaR1J >= 0.4) {
                            double jet_pt = isSmeared ? GaussSmeared(
                                    mRNG, jet.pt(), StandardDeviation(jet.pt(), mCSN)
                                ) : jet.pt();
                            if (jet_pt > 30.0 && ZBoson.pt() > 60.0) {
                                mPhiHisto.addEventNum(delta_phi, 
                                    weight0 / weight2);
                            }
                            if ( delta_phi > (7 * PI / 8.0) ) {
                                if (jet_pt > 30.0) {
                                    if (ZBoson.pt() > 60.0) {
                                        mXjZHisto.addEventNum(jet_pt / ZBoson.pt(), 
                                            weight0 / weight2);
                                    }
                                    double x_jZ_Weighted_Sum = (
                                        jet_pt / ZBoson.pt()
                                    ) * (weight0 / weight2);
                                    mMeanXjzDenoHisto.addEventNum(ZBoson.pt(), weight0 / weight2 );
                                    mMeanXjzHisto.addEventNum(ZBoson.pt(), x_jZ_Weighted_Sum);
                                    mRjZHisto.addEventNum(ZBoson.pt(), 
                                        weight0 / weight2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    string comments = "";
    if (isSmeared) {
        comments = "# Smeared: True, and\n# C: " + to_string(mCSN.C) 
            + ", S: " + to_string(mCSN.S) + ", N: " 
            + to_string(mCSN.N) + "\n";
    }

    // 输出 delta phi_jZ
    vector<distInfo> mPhiHistoInfo = mPhiHisto.getDHisto();
    WriteDataToText(args.get<string>("delta-phi"), mPhiHistoInfo, comments);
    // 输出 x_jZ
    vector<distInfo> mXjZHistoInfo = mXjZHisto.getDHisto();
    WriteDataToText(args.get<string>("x-jz"), mXjZHistoInfo, comments);
    // 输出 <x_jZ> 
    vector<distInfo> mMeanXjzHistoInfo = mMeanXjzHisto.getHistoNoNorm();
    vector<double> meanXjZDenoVec = mMeanXjzDenoHisto.getBinValues();
    for (size_t i = 0; i < mMeanXjzHistoInfo.size(); i++) {
        mMeanXjzHistoInfo[i].distValue /= meanXjZDenoVec[i];
    }
    WriteDataToText(args.get<string>("mean-x-jz"), mMeanXjzHistoInfo, comments);
    // 输出 R_jZ
    vector<distInfo> mRjZHistoInfo = mRjZHisto.getHistoNoNorm();
    vector<double> RjZDenoVec = mRjZDenoHisto.getBinValues();
    for (size_t i = 0; i < mRjZHistoInfo.size(); i++) {
        mRjZHistoInfo[i].distValue /= RjZDenoVec[i];
    }
    WriteDataToText(args.get<string>("R_jZ"), mRjZHistoInfo, comments);
    
    
}