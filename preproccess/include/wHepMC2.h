/*
 * FileName: wHepMC2.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 16:55:11
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-15 17:47:59
 * Description: Write events' information to the file.
 */


#ifndef IHEPTOOLS_WHEPMC2_H
#define IHEPTOOLS_WHEPMC2_H

#include <iostream>
#include <vector>

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "particle.h"

using namespace std;

namespace iHepTools {
    bool WriteToHepMC2(HepMC::IO_GenEvent &mIO_GenEvent, const int &mEvent_Id, 
        const double &mAlpahQCD, const double &mAlphaQED, const double &weight0, 
        const double &weight1, const double &weight2, const double &weight3, 
        const double &weight4, const vector<Particle> &mRetainPartons, 
        const vector<Particle> &mPartons);
}

// 将事件输出到 HepMC2 文件
bool iHepTools::WriteToHepMC2(HepMC::IO_GenEvent &mIO_GenEvent, const int &mEvent_Id, 
        const double &mAlpahQCD, const double &mAlphaQED, const double &weight0, 
        const double &weight1, const double &weight2, const double &weight3, 
        const double &weight4, const vector<Particle> &mRetainPartons, 
        const vector<Particle> &mPartons) {
    //

    // 1st create the event container
    HepMC::GenEvent* evt = new HepMC::GenEvent();

    // 2nd set event info
    evt->set_event_number(mEvent_Id);
    evt->set_event_scale(0);
    evt->set_alphaQCD(mAlpahQCD);
    evt->set_alphaQED(mAlphaQED);
    evt->weights()["0"] = weight0;
    evt->weights()["1"] = weight1;
    evt->weights()["2"] = weight2;
    evt->weights()["3"] = weight3;
    evt->weights()["4"] = weight4;
    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

    // 3rd create 1 vertex and add to event
    HepMC::GenVertex* vertex = new HepMC::GenVertex();
    evt->add_vertex(vertex);

    // 4th create retain partons and write to vertex
    // here record some particles and that not from outgoing but we need they
    for (unsigned i = 0; i < mRetainPartons.size(); i++) {
        HepMC::GenParticle* retain_particle = new HepMC::GenParticle(
            HepMC::FourVector(
                mRetainPartons[i].x1, 
                mRetainPartons[i].x2, 
                mRetainPartons[i].x3, 
                mRetainPartons[i].x0
            ), mRetainPartons[i].pdg_code, 1
        );
        vertex->add_particle_out(retain_particle);
        // delete retain_particle;
    }

    // 5th create outgoing partons and write to vertex
    for (unsigned i = 0; i < mPartons.size(); i++) {
        HepMC::GenParticle* outgoing_particle = new HepMC::GenParticle(
            HepMC::FourVector(
                mPartons[i].x1, 
                mPartons[i].x2, 
                mPartons[i].x3, 
                mPartons[i].x0 
            ), mPartons[i].pdg_code, 1
        );
        vertex->add_particle_out(outgoing_particle);
    }

    // 6th output event to file
    mIO_GenEvent << evt;

    // end clean event object
    delete evt;
}

#endif // IHEPTOOLS_WHEPMC2_H