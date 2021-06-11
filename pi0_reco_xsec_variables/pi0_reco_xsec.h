//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef PI0_MC_STUDY_H
#define PI0_MC_STUDY_H

#include "../utilities/Histograms.hpp"

void run_pi0_reco_xsec( std::string in_file, Histograms &hists );

std::vector<int> get_reco_gamma_index( std::vector<int> &true_pdg, std::vector<int> &true_parent_pdg );

std::vector<int> get_true_reco_gamma_index( std::vector<int> &true_pdg, std::vector<int> &true_id,
                                            std::vector<int> &true_parent_pdg );

int get_reco_pi0_gamma_count( std::vector<int> &true_pdg, std::vector<int> &true_parent_pdg );

double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 );

double E_gamma_gammma( double E_g1, double E_g2, double theta_gg );

double theta_gamma_gammma( double E_g1, double E_g2, double theta_gg, double theta_g1, double theta_g2 );


#endif //PI0_MC_STUDY_H

