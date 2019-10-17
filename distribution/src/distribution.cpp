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
 * LastEditTime: 2019-10-17 13:54:19
 * Description: Calculate distribution.
 */


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include "cmdline/cmdline.h"

#include "jetSelector.h"
#include "histo.h"
#include "zboson.h"

using namespace std;
using namespace iHepTools;

// returns true if the GenParticle does not decay
inline bool isFinal( const HepMC::GenParticle* p ) {
    return !p->end_vertex() && p->status()==1;
}

constexpr double PI = acos(-1);
constexpr double R_jet = 0.3;

// Jet definition
JetDefinition jet_def(antikt_algorithm, R_jet);
Selector select_akt = SelectorAbsEtaMax(1.6) && SelectorPtMin(30.);

int main(int argc, char *argv[]) {
    
    // 参数解析
    cmdline::parser args;
    args.add<string>("name", 'n', "assigned output file name", true, "out1.dat");
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
    // 

    // 定义 phi histo
    vector<double> mPhiHistoPointList = 
        {0, 0.31, 0.94, 1.25, 1.72, 2.20, 2.51, 2.82, 2.98, PI};
    Histo mPhiHisto(mPhiHistoPointList);
    
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
                vector<PseudoJet> jets = SelectJet(pseduo_jets, jet_def, select_akt);
                if (jets.size() > 0) {
                    for (const auto jet : jets) {
                        double delta_phi = fabs(jet.phi() - ZBoson.phi());
                        delta_phi = delta_phi > PI ? 2 * PI - delta_phi : delta_phi;
                        mPhiHisto.addEventNum(delta_phi, (evt->weights())[0] / (evt->weights())[2]);
                    }
                }
            }

            delete evt;
            ascii_in >> evt;
            
        }
    }

    vector<distInfo> mPhiHistoInfo = mPhiHisto.getHisto();
    ofstream mDeltaPhiFile;
    mDeltaPhiFile.open(args.get<string>("name"));
    mDeltaPhiFile << "# xlow\txhigh\txmiddle\tval" << '\n';
    for (size_t i = 0; i < mPhiHistoInfo.size(); i++) {
        mDeltaPhiFile << setw(12) << scientific << setprecision(6)
            << mPhiHistoInfo[i].region.leftValue << '\t' 
            << mPhiHistoInfo[i].region.rightValue << '\t'
            << mPhiHistoInfo[i].region.middleValue << '\t'
            << mPhiHistoInfo[i].distValue << '\n';
    }
}