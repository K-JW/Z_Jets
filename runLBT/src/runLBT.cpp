/*
 * FileName: runLBT.cpp
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-18 09:09:56
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-20 14:51:01
 * Description: run LBT model
 */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <mpi.h>

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "cmdline/cmdline.h"
#include "writeHepMC2.h"
#include "getfilename.h"
#include "LBT.h"
#include "detector.h"

using namespace std;
namespace hep = iHepTools;

// returns true if the GenParticle does not decay
inline bool isFinal( const HepMC::GenParticle* p ) {
    return !p->end_vertex() && p->status()==1;
}

int main(int argc, char *argv[]) {

    // 参数解析
    cmdline::parser args;
    args.add<string>("path", 'p', "assigned the path of LBT model data", false, "./");
    args.add<string>("out-dir", 'd', "assigned the path of output data", false, "./Output");
    args.add("help", 'h', "print help message");
    args.footer("filenames...");
    args.set_program_name("runLBT");

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
    const int file_num = file_name_vec.size();
    //

    // 初始化 MPI
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (file_num < size) {
        cerr << "--\033[1;31m ERROR:\033[0m Assigned process number greater than files number, please check it!"
            << endl;
        exit(EXIT_FAILURE);
    }

    // 获取一般 MPI 进程可以分配到的数据量
    const int data_size_general = floor(file_num / size);
    // 获取最后一个 MPI 进程分配到的数据量
    const int data_size_last = file_num % size;
    
    // 申请动态数组使用的内存块
    int* send_count = new int[size]();
    int* global_file_index = new int[file_num]();
    int* locale_file_index = new int[file_num]();
    int* offset_array = new int[size]();
    //
    for (size_t i = 0; i < file_num; i++)
        global_file_index[i] = i;
    for (size_t i = 0; i <size; i++)
        send_count[i] = data_size_general;
    for (size_t i = 0; i < data_size_last; i++)
        send_count[i % size] += 1;
    
    offset_array[0] = 0;
    for (size_t i = 1; i < size; i++)
        offset_array[i] = offset_array[ i - 1 ] + send_count[ i - 1 ];

    // 分发数据
    MPI_Scatterv(global_file_index, send_count, offset_array, MPI_INT, 
        locale_file_index, send_count[rank], MPI_INT, 0, MPI_COMM_WORLD);

    // 初始化 LBT model
    LBT *mLBT = new LBT();
    mLBT->Initialize(args.get<string>("path") + "/LBT-Config.yaml");
    double mJetNumElement, mXNucleonElement, mYNucleonElement, mJetNumTotal = 0.;
    ifstream Geomtery(args.get<string>("path") + "/LBT-tables/HydroProfile/geometry.dat");
    vector<double> mJetNum, XNucleon, YNucleon;
    if (Geomtery.is_open()) {
        while (!Geomtery.eof()) {
            Geomtery >> mJetNumElement >> mXNucleonElement >> mYNucleonElement;
            mJetNumTotal += mJetNumElement;
            mJetNum.emplace_back(mJetNumElement);
            XNucleon.emplace_back(mXNucleonElement);
            YNucleon.emplace_back(mYNucleonElement);
        }
    } else {
        cerr << "--\033[1;31m ERROR:\033[0m Can't open geometry.dat file." << endl;
        exit(EXIT_FAILURE);
    }
    

    for (size_t i = 0; i < send_count[rank]; i++) {
        
        // 待处理文件名
        string file_name = file_name_vec[locale_file_index[i]];
        const char *hepmc_in = file_name.c_str();
        cout << "    Rank: " << rank << ", file: " << file_name << '\n';
        
        // 输出文件名
        string hepmc_out = args.get<string>("out-dir") + "/" + 
            hep::get_file_name(file_name) + "Rank-" + to_string(rank) + 
            "-LBT.hepmc";
        
        // Read HepMC
        HepMC::IO_GenEvent ascii_in(hepmc_in, std::ios::in);
        HepMC::GenEvent* evt = ascii_in.read_next_event();
        // Write HepMC
        HepMC::IO_GenEvent out(hepmc_out, std::ios::out);
        
        int mEvent_ID = -1;

        while (evt) {

            vector<Particle> mPartons;
            vector<Particle> mRetainPartons;
            double mRandomXY = mLBT->MyRandom();
            double R1 = 0., XXX, YYY;
            for (int i = 0; i < XNucleon.size(); i++) {
                R1 += mJetNum[i] / mJetNumTotal;
                if (mRandomXY < R1) {
                    XXX = XNucleon[i];
                    YYY = YNucleon[i];
                    break;
                }
            }

            for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); 
                    p != evt->particles_end(); ++p) {
                //
                double mEta, mP_T;
                Particle mJet;
                FourVector mMomentum, mPosition(0., XXX, YYY, 0.);
                
                if ( isFinal(*p) ) {
                    int pdg_code = (*p)->pdg_id();
                    if ( hep::isDown(pdg_code) || hep::isUp(pdg_code) || 
                        hep::isStrange(pdg_code) || hep::isCharm(pdg_code) || 
                        hep::isBottom(pdg_code) || hep::isGluon(pdg_code) ) {
                        //
                        mMomentum.SetVector((*p)->momentum().e(), (*p)->momentum().px(), 
                            (*p)->momentum().py(), (*p)->momentum().pz());
                        mEta = 1.0 / 2.0 * log( (mMomentum.x0() + mMomentum.x3()) / (
                            mMomentum.x0() - mMomentum.x3()
                            ) 
                        );
                        mP_T = sqrt(pow(mMomentum.x1(), 2) + pow(mMomentum.x2(), 2));

                        if (abs(mEta) <= 2.4 && mP_T >= 0.5) {
                            mJet.SetParticleInfo((*p)->pdg_id(), mMomentum, mPosition);
                            mJet.SetFormationTime(2.0 * mMomentum.x0() / (
                                pow(mMomentum.x1(), 2) + pow(mMomentum.x2(), 2)
                                )
                            );
                            mJet.SetRadiationTime(mJet.GetFormationTime());
                            mPartons.emplace_back(mJet);
                        }
                    } else if ( hep::isElectron(pdg_code) || hep::isMuon(pdg_code) ) {
                        mMomentum.SetVector((*p)->momentum().e(), (*p)->momentum().px(), 
                            (*p)->momentum().py(), (*p)->momentum().pz());
                        mJet.SetParticleInfo((*p)->pdg_id(), mMomentum, mPosition);
                        mRetainPartons.emplace_back(mJet);
                    }
                }
            }

                            
            mLBT->SetEvent(mPartons);

            TimeInfo mTimeInfo = mLBT->GetTimeInfo();
            double mTime = mTimeInfo.Start;
            for (unsigned i = 0; i < mTimeInfo.StepNum; i++) {
                mTime += mTimeInfo.Step;
                mLBT->LinearBoltzmannTransport(mTime);
                if (mLBT->isReachTauEnd()) {
                    i = mTimeInfo.StepNum - 1;
                }
            }

            vector<Particle> mFinalPartons; 
            mFinalPartons = mLBT->GetFinalPartons();
            // output
            if (!mFinalPartons.empty()) {
                mEvent_ID++;
                hep::WritePartonsToHepMC2(out, mEvent_ID, evt->alphaQCD(), 
                    evt->alphaQED(), evt->weights()[0], evt->weights()[1], 
                    evt->weights()[2], evt->weights()[3], evt->weights()[4], 
                    mRetainPartons, mFinalPartons);
            }

            mLBT->clear();
            
            delete evt;
            ascii_in >> evt;
        }
        
    }
    
    // 输出运行成功的信息
    cout << "\033[0;32m Info:\033[0m Tasks at rank " << rank << " is finished!" << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // 回收内存
    delete []send_count;
    delete []global_file_index;
    delete []locale_file_index;
    delete []offset_array;

    // 输出运行成功的信息
    if (rank == 0)
        cout << "\n\033[1;32mAll Tasks are finished!\033[0m" << endl;

    MPI_Finalize();

    return 0;
}