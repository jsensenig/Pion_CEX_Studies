//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef PI0_MC_STUDY_H
#define PI0_MC_STUDY_H

#include "../utilities/Histograms.hpp"

void run_pi0_mc_study( std::string in_file, Histograms &hists );
double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 );


#endif //PI0_MC_STUDY_H

