/*
 * FileName: jetSelector.h
 * Version: 1.0
 * License: MIT License
 * 
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-16 11:22:47
 * 
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-10-16 15:39:48
 * Description: Define function to search jet
 */

#ifndef IHEPTOOLS_JETSELECTOR_H
#define IHEPTOOLS_JETSELECTOR_H

#include <vector>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>

namespace iHepTools {

    std::vector<fastjet::PseudoJet> SelectJet(const std::vector<fastjet::PseudoJet> 
        &hadrons, fastjet::JetDefinition jet_def, fastjet::Selector select);
}

std::vector<fastjet::PseudoJet> iHepTools::SelectJet(const std::vector<fastjet::PseudoJet> 
        &hadrons, fastjet::JetDefinition jet_def, fastjet::Selector select) {
    //
    fastjet::ClusterSequence cluster(hadrons, jet_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(select(
        cluster.inclusive_jets()));
    return jets;
}

#endif // IHEPTOOLS_JETSELECTOR_H