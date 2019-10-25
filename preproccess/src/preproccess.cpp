/*
 * FileName: preproccess.cpp
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 17:06:05
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-25 11:41:50
 * Description: Preproccess HepMC2 file.
 */

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "wHepMC2.h"
#include "particle.h"
#include "cmdline/cmdline.h"

using namespace std;

namespace iHepTools {
    string get_file_name(string filepath);
    
    // returns true if the GenParticle does not decay
    inline bool isFinal( const HepMC::GenParticle* p) {
        return !p->end_vertex() && ( p->status() == 1 );
    }
}

int main(int argc, char *argv[]) {
    
    // 参数解析
    cmdline::parser args;
    args.add<string>("out-dir", 'o', "assigned the path of output data", false, "./Output");
    args.add<double>("Z-mass-Min", 'l', "assigned select Z's max mass", false, 60.0);
    args.add("help", 'h', "print help message");
    args.footer("filenames...");
    args.set_program_name("prep");
    bool isGetArgs = args.parse(argc, argv);
    if (argc == 1 || args.exist("help")) {
        cerr << args.usage();
        cout << "\033[1;31mERROR:\033[0m Please supply HepMC files as argument...\n";
        exit(EXIT_FAILURE);
    }
    if (!isGetArgs) {
        cerr << args.error() << endl << args.usage();
        cout << "\033[1;31mERROR:\033[0m Please supply HepMC files as argument...\n";
        exit(EXIT_FAILURE);
    }
    vector<string> file_name_vec;
    for (size_t i = 0; i < args.rest().size(); i++) {
        file_name_vec.emplace_back(args.rest()[i]);
    }
    const int file_num = file_name_vec.size();
    const double ZMassMin = args.get<double>("Z-mass-Min");
    // if (argc == 1) {
    //     cout << "    " << "ERROR: Please supply HepMC file as arguments..." << endl;
    //     exit(EXIT_FAILURE);
    // }

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
    const int data_size_general = file_num / size;

    // 获取最后一个 NPI 进程分配到的数据量
    const int data_size_last = file_num % size;

    // 申请动态数组使用的内存块 
    int* send_count = new int[size]();
    int* global_file_index = new int[file_num]();
    int* locale_file_index = new int[file_num]();
    int* offset_array = new int[size]();

    for (int i = 0; i < file_num; i++) {
        global_file_index[i] = i;
    }

    for (int i = 0; i < size; i++) {
        send_count[i] = data_size_general;
    }

    for (int i = 0; i < data_size_last; i++) {
        send_count[i % size] += 1;
    }

    offset_array[0] = 0;
    for (int i = 1; i < size; i++) {
        offset_array[i] = offset_array[i - 1] + send_count[i - 1];
    }

    // 分发数据
    MPI_Scatterv(global_file_index, send_count, offset_array, MPI_INT, 
        locale_file_index, send_count[rank], MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < send_count[rank]; i++) {
        
        // 待处理文件名
        string hepmc_in = file_name_vec[locale_file_index[i]];

        std::cout << "Rank: " << rank << ", File name: " 
            << hepmc_in << std::endl;


        // 输出文件名
        std::string output = args.get<string>("out-dir") + "/" +
            iHepTools::get_file_name(static_cast<std::string>(hepmc_in)) + 
            static_cast<std::string>("Rank-") + 
            std::to_string(rank) +
            static_cast<std::string>("-Prep.hepmc");

        HepMC::IO_GenEvent out(output, std::ios::out);

        HepMC::IO_GenEvent ascii_in(hepmc_in, std::ios::in);
        HepMC::GenEvent* evt = ascii_in.read_next_event();

        int fact_event_number = -1;
        while (evt) {

            std::vector<iHepTools::Particle> mParticles;
            std::vector<iHepTools::Particle> mZ_daughter;
            bool isGoodEvent = false;


            for (HepMC::GenEvent::particle_iterator p = evt->particles_begin();
                    p != evt->particles_end(); ++p) {
                
                if ( iHepTools::isFinal(*p) && ( abs((*p)->pdg_id()) == 11 || 
                    abs((*p)->pdg_id()) == 13 ) ) {
                    
                    iHepTools::Particle tmp = {
                        (*p)->pdg_id(), 
                        (*p)->momentum().e(), 
                        (*p)->momentum().px(), 
                        (*p)->momentum().py(), 
                        (*p)->momentum().pz(), 
                        (*p)->generated_mass()
                    };
                    mZ_daughter.emplace_back(tmp);
                } else if ( iHepTools::isFinal(*p)) {
                    iHepTools::Particle tmp = {
                        (*p)->pdg_id(), 
                        (*p)->momentum().e(), 
                        (*p)->momentum().px(), 
                        (*p)->momentum().py(), 
                        (*p)->momentum().pz(), 
                        (*p)->generated_mass()
                    };
                    mParticles.emplace_back(tmp);
                } else {
                    continue;
                }
            }

            if (mZ_daughter.size() == 2) {
                if (sqrt(pow(mZ_daughter[0].x1 + mZ_daughter[1].x1, 2) + 
                    pow(mZ_daughter[0].x2 + mZ_daughter[1].x2, 2)) > ZMassMin) {
                    //
                    fact_event_number++;
                    isGoodEvent = true;
                } else {
                    delete evt;
                    ascii_in >> evt;
                    continue;
                }
            }

            if (isGoodEvent) {
                iHepTools::WriteToHepMC2(out, fact_event_number, 
                    evt->alphaQCD(), evt->alphaQED(), evt->weights()[0], 
                    evt->weights()[1], evt->weights()[2], evt->weights()[3], 
                    evt->weights()[4], mZ_daughter, mParticles);
            }

            delete evt;
            ascii_in >> evt;
        }
    }

    std::cout << "Tasks at rank " << rank << " is finished!" << std::endl;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    delete []send_count;
    delete []global_file_index;
    delete []locale_file_index;
    delete []offset_array;

    MPI_Finalize();

    std::cout << "All Tasks are finished!" << std::endl;

    return 0;

}

string iHepTools::get_file_name(string filepath) {
    if (!filepath.empty()) {
        int location_point = filepath.find_last_of('.');
        int location_filename = filepath.find_last_of('/');
        std::string file_name = filepath.substr(location_filename + 1, 
            location_point - location_filename - 1);
        return file_name;
    } else {
        return "NULL";
    }
}