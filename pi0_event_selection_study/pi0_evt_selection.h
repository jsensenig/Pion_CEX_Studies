//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef PI0_EVT_SELECTION_H
#define PI0_EVT_SELECTION_H

#include "../utilities/Histograms.hpp"
#include "../utilities/Chi2PIDAlg.hpp"
#include "selection.h"

/// Functions ///
void run_pi0_evt_selection( const std::string& in_file, Histograms &hists, Selection &sel );

bool BeamTofCut( Histograms& hists, Selection& sel );

bool BeamToTpcCut( Histograms& hists, Selection& sel );

bool BeamEndZCut( Histograms& hists, Selection& sel );

bool DaughterTrackScore( Histograms& hists, Selection& sel, size_t i );

bool DaughterMichelScore( Histograms& hists, Selection& sel, size_t i );

bool DaughterPionCut( Histograms& hists, Selection& sel, pid::Chi2PIDAlg& chi2, size_t i );

bool DaughterDeltaRCut( Histograms& hists, Selection& sel, size_t i );

bool DaughterNhitCut( Histograms& hists, Selection& sel, size_t i );

void CleanMemory();

/// Selection Counts ///
bool true_cex;
size_t true_cex_cnt = 0;
size_t reco_cex_cnt = 0;
size_t true_reco_cex = 0;
// Pair contents: pair< true_cnt, reco_cnt >
std::map<std::string, std::pair<size_t, size_t>> select_count;
std::map<std::string, std::pair<size_t, size_t>> reject_count;

/// Data variables ///
int true_daughter_nPi0;
int true_daughter_nNeutron;
int true_daughter_nProton;
int true_daughter_nPiMinus;
int true_daughter_nPiPlus;

// True beam
int reco_beam_true_byHits_PDG;
auto true_beam_endProcess = new std::string;

// Reco beam
bool reco_beam_passes_beam_cuts;
double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
auto beam_inst_TOF = new std::vector<double>;
auto reco_beam_resRange = new std::vector<double>;
auto reco_beam_calibrated_dEdX = new std::vector<double>;

// True Daughter
auto true_beam_daughter_PDG = new std::vector<int>;
auto reco_daughter_PFP_true_byHits_parPDG = new std::vector<int>;
auto reco_daughter_PFP_true_byHits_PDG = new std::vector<int>;

// Reco Daughter
auto reco_daughter_allTrack_calibrated_dEdX_SCE = new std::vector<std::vector<double>>;
auto reco_daughter_allTrack_resRange_SCE = new std::vector<std::vector<double>>;
auto reco_daughter_allTrack_dR = new std::vector<double>;
auto reco_daughter_PFP_nHits = new std::vector<int>;
auto reco_daughter_PFP_michelScore_collection = new std::vector<double>;
auto reco_daughter_allTrack_Chi2_proton = new std::vector<double>;
auto reco_daughter_allTrack_Chi2_ndof = new std::vector<int>;
auto reco_daughter_PFP_trackScore_collection = new std::vector<double>;
auto reco_daughter_allShower_startX = new std::vector<double>;
auto reco_daughter_allShower_startY = new std::vector<double>;
auto reco_daughter_allShower_startZ = new std::vector<double>;


#endif //PI0_EVT_SELECTION_H

