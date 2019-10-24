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
 * LastEditTime: 2019-10-24 15:32:24
 * Description: Calculate distribution.
 */


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include "cmdline/cmdline.h"

#include "jetSelector.h"
#include "histo.h"
#include "zboson.h"
#include "smeared.h"

using namespace std;
using namespace iHepTools;

// returns true if the GenParticle does not decay
inline bool isFinal( const HepMC::GenParticle* p ) {
    return !p->end_vertex() && p->status()==1;
}

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
    args.add("enable-smeared", 's', "execute smearing for data");
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
    
    for (size_t i = 0; i < file_name_vec.size(); i++) {

        HepMC::IO_GenEvent ascii_in(file_name_vec[i], std::ios::in);
        HepMC::GenEvent* evt = ascii_in.read_next_event();

        while ( evt ) {
            
            vector<fastjet::PseudoJet> pseduo_jets;
            vector<fastjet::PseudoJet> Z_daughters;
            fastjet::PseudoJet ZBoson;

            for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); 
                p != evt->particles_end(); ++p) {
                //
                if (isFinal(*p)) {
                    if ((abs((*p)->pdg_id()) == 11 || abs((*p)->pdg_id()) == 13 )) {
                        fastjet::PseudoJet tmp((*p)->momentum().px(),(*p)->momentum().py(), 
                            (*p)->momentum().pz(), (*p)->momentum().e());
                        tmp.set_user_index((*p)->pdg_id());
                        Z_daughters.emplace_back(tmp);
                    } else {
                        fastjet::PseudoJet tmp((*p)->momentum().px(),(*p)->momentum().py(), 
                            (*p)->momentum().pz(), (*p)->momentum().e());
                        pseduo_jets.emplace_back(tmp);
                    }
                }
            }

            if (isZBoson(Z_daughters, ZBoson)) {
                mPhiHisto.addEventNorm( (evt->weights())[0] / (evt->weights())[2] );
                mXjZHisto.addEventNorm( (evt->weights())[0] / (evt->weights())[2] );
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
                            if (jet.pt() >= 30.0) {
                                mPhiHisto.addEventNum(delta_phi, 
                                    (evt->weights())[0] / (evt->weights())[2]);
                            }
                            if ( delta_phi > (7 * PI / 8.0) ) {
                                double jet_pt = isSmeared ? GaussSmeared(
                                    mRNG, jet.pt(), Variance(jet.pt(), mCSN)
                                ) : jet.pt();
                                if (jet_pt > 30.0) {
                                    mXjZHisto.addEventNum(jet_pt / ZBoson.pt(), 
                                        (evt->weights())[0] / (evt->weights())[2]);
                                }
                            }
                        }
                    }
                }
            }

            delete evt;
            ascii_in >> evt;
            
        }
    }

    vector<distInfo> mPhiHistoInfo = mPhiHisto.getHisto();
    ofstream mDeltaPhiFile;
    mDeltaPhiFile.open(args.get<string>("delta-phi"));
    if (isSmeared) {
        mDeltaPhiFile << "# Smeared: True, and\n"
            << "# C: " << mCSN.C << ", S: " << mCSN.S << ", N: "
            << mCSN.N << '\n';
    }
    mDeltaPhiFile << "# xlow\txhigh\txmiddle\tval" << '\n';
    for (size_t i = 0; i < mPhiHistoInfo.size(); i++) {
        mDeltaPhiFile << setw(12) << scientific << setprecision(6)
            << mPhiHistoInfo[i].region.leftValue << '\t' 
            << mPhiHistoInfo[i].region.rightValue << '\t'
            << mPhiHistoInfo[i].region.middleValue << '\t'
            << mPhiHistoInfo[i].distValue << '\n';
    }

    vector<distInfo> mXjZHistoInfo = mXjZHisto.getHisto();
    ofstream mXjZFile;
    mXjZFile.open(args.get<string>("x-jz"));
    if (isSmeared) {
        mXjZFile << "# # Smeared: True, and\n"
            << "# C: " << mCSN.C << ", S: " << mCSN.S << ", N: "
            << mCSN.N << '\n';
    }
    mXjZFile << "# xlow\txhigh\txmiddle\tval" << '\n';
    for (size_t i = 0; i < mXjZHistoInfo.size(); i++) {
        mXjZFile << setw(12) << scientific << setprecision(6)
            << mXjZHistoInfo[i].region.leftValue << '\t' 
            << mXjZHistoInfo[i].region.rightValue << '\t'
            << mXjZHistoInfo[i].region.middleValue << '\t'
            << mXjZHistoInfo[i].distValue << '\n';
    }
    
}