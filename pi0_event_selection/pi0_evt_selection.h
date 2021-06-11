//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef PI0_EVT_SELECTION_H
#define PI0_EVT_SELECTION_H

#include "../utilities/Histograms.hpp"
#include "selection.h"

void run_pi0_evt_selection( std::string in_file, Histograms &hists, Selection &sel );
std::vector<int> get_reco_gamma_index( std::vector<int> &true_pdg, std::vector<int> &true_parent_pdg );
std::vector<int> get_true_reco_gamma_index( std::vector<int> &true_pdg, std::vector<int> &true_id,
                                            std::vector<int> &true_parent_pdg );
double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 );
double E_gamma_gammma( double E_g1, double E_g2, double theta_gg );
double theta_gamma_gammma( double E_g1, double E_g2, double theta_gg, double theta_g1, double theta_g2 );
bool set_selection_parameters( std::string param_file );

/// Selection Counts ///
size_t true_cex_cnt = 0;
size_t reco_cex_cnt = 0;
size_t true_reco_cex = 0;
std::map<std::string, size_t> select_count;
std::map<std::string, size_t> reject_count;

/// Data variables ///
int true_daughter_nPi0;
int true_daughter_nNeutron;
int true_daughter_nProton;
int true_daughter_nPiMinus;
int true_daughter_nPiPlus;
int true_beam_PDG;
double true_beam_endP;
double true_beam_endX, true_beam_endY, true_beam_endZ;
double true_beam_endPx, true_beam_endPy, true_beam_endPz;

bool reco_beam_passes_beam_cuts;
int reco_beam_true_byHits_PDG;
double beam_inst_P;

std::vector<double> *beam_inst_TOF = new std::vector<double>;
std::vector<double> *reco_beam_resRange = new std::vector<double>;
std::vector<double> *reco_beam_calibrated_dEdX = new std::vector<double>;

std::vector<double> *reco_daughter_allTrack_ID = new std::vector<double>;
std::vector<std::vector<double>> *reco_daughter_allTrack_calibrated_dEdX_SCE = new std::vector<std::vector<double>>;
std::vector<std::vector<double>> *reco_daughter_allTrack_resRange_SCE = new std::vector<std::vector<double>>;
std::vector<double> *reco_daughter_allTrack_dR = new std::vector<double>;
std::vector<int> *reco_daughter_PFP_nHits = new std::vector<int>;
std::vector<double> *reco_daughter_allTrack_startX = new std::vector<double>;
std::vector<double> *reco_daughter_allTrack_startY = new std::vector<double>;
std::vector<double> *reco_daughter_allTrack_startZ = new std::vector<double>;
std::vector<double> *reco_daughter_PFP_michelScore_collection = new std::vector<double>;

std::vector<int> *true_beam_daughter_PDG = new std::vector<int>;
std::vector<int> *true_beam_daughter_ID = new std::vector<int>;
std::vector<double> *true_beam_daughter_startP = new std::vector<double>;
std::vector<double> *true_beam_daughter_startPx = new std::vector<double>;
std::vector<double> *true_beam_daughter_startPy = new std::vector<double>;
std::vector<double> *true_beam_daughter_startPz = new std::vector<double>;
std::vector<int> *true_beam_Pi0_decay_PDG = new std::vector<int>;
std::vector<int> *true_beam_Pi0_decay_parID = new std::vector<int>;
std::vector<int> *true_beam_Pi0_decay_ID = new std::vector<int>;
std::vector<double> *true_beam_Pi0_decay_startP = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_startPx = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_startPy = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_startPz = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_startX = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_startY = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_startZ = new std::vector<double>;
std::vector<double> *true_beam_Pi0_decay_len = new std::vector<double>;

double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
double reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ;

std::vector<double> *reco_daughter_PFP_trackScore_collection = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_dirX = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_dirY = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_dirZ = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_energy = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_startX = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_startY = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_startZ = new std::vector<double>;
std::vector<int> *reco_daughter_allShower_ID = new std::vector<int>;
std::vector<std::vector<int>> *true_beam_daughter_reco_byHits_allShower_ID = new std::vector<std::vector<int>>;
std::vector<int> *reco_daughter_PFP_true_byHits_parPDG = new std::vector<int>;
std::vector<int> *reco_daughter_PFP_true_byHits_PDG = new std::vector<int>;
std::vector<int> *reco_daughter_PFP_true_byHits_parID = new std::vector<int>;
std::vector<int> *reco_daughter_PFP_true_byHits_ID = new std::vector<int>;
std::vector<double> *reco_daughter_PFP_true_byHits_endX = new std::vector<double>;
std::vector<double> *reco_daughter_PFP_true_byHits_endY = new std::vector<double>;
std::vector<double> *reco_daughter_PFP_true_byHits_endZ = new std::vector<double>;
std::vector<double> *reco_daughter_PFP_true_byHits_startPx = new std::vector<double>;
std::vector<double> *reco_daughter_PFP_true_byHits_startPy = new std::vector<double>;
std::vector<double> *reco_daughter_PFP_true_byHits_startPz = new std::vector<double>;
std::vector<std::string> *reco_daughter_PFP_true_byHits_endProcess = new std::vector<std::string>;

std::vector<std::vector<int>> *true_beam_Pi0_decay_reco_byHits_allShower_ID = new std::vector<std::vector<int>>;
std::string *true_beam_endProcess = new std::string;


#endif //PI0_EVT_SELECTION_H

