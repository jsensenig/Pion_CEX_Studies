//
// Created by Jon Sensenig on 6/5/21.
//

#ifndef RECO_VERTEX_EVENTSELECTION_H
#define RECO_VERTEX_EVENTSELECTION_H

#include "TTree.h"
#include "../utilities/Histograms.hpp"
#include "../utilities/Chi2PIDAlg.hpp"
#include "selection.h"

class EventSelection {

public:

  EventSelection();
  virtual ~EventSelection();

  void InitTtree( TTree* tree );

  bool SelectEvent( TTree* tree, size_t evt, Histograms& hists, bool true_cex, int nPi0 );

private:

  pid::Chi2PIDAlg chi2;
  Selection sel;

  /// Selection Functions ///

  bool BeamTofCut();

  bool BeamToTpcCut();

  bool BeamEndZCut();

  bool DaughterTrackScore( size_t i );

  bool DaughterMichelScore( size_t i );

  bool DaughterPionCut( size_t i );

  bool DaughterDeltaRCut( size_t i );

  bool DaughterNhitCut( size_t i );

  void CleanMemory();

  /// Selection Counts ///
  bool true_cex;
  size_t true_cex_cnt = 0;
  size_t reco_cex_cnt = 0;
  size_t true_reco_cex = 0;
  // Pair contents: pair< true_cnt, reco_cnt >
  std::map<std::string, std::pair<size_t, size_t>> select_count;
  std::map<std::string, std::pair<size_t, size_t>> reject_count;

  // True beam
  int true_daughter_nPi0;
  int reco_beam_true_byHits_PDG;
  std::string *true_beam_endProcess = new std::string;

  // Reco beam
  bool reco_beam_passes_beam_cuts;
  double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
  std::vector<double> *beam_inst_TOF = new std::vector<double>;

  // True Daughter
  std::vector<int> *true_beam_daughter_PDG = new std::vector<int>;
  std::vector<int> *reco_daughter_PFP_true_byHits_parPDG = new std::vector<int>;
  std::vector<int> *reco_daughter_PFP_true_byHits_PDG = new std::vector<int>;

  // Reco Daughter
  std::vector<std::vector<double>> *reco_daughter_allTrack_calibrated_dEdX_SCE = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *reco_daughter_allTrack_resRange_SCE = new std::vector<std::vector<double>>;
  std::vector<double> *reco_daughter_allTrack_dR = new std::vector<double>;
  std::vector<int> *reco_daughter_PFP_nHits = new std::vector<int>;
  std::vector<double> *reco_daughter_PFP_michelScore_collection = new std::vector<double>;
  std::vector<double> *reco_daughter_allTrack_Chi2_proton = new std::vector<double>;
  std::vector<int> *reco_daughter_allTrack_Chi2_ndof = new std::vector<int>;
  std::vector<double> *reco_daughter_PFP_trackScore_collection = new std::vector<double>;
  std::vector<double> *reco_daughter_allShower_startX = new std::vector<double>;
  std::vector<double> *reco_daughter_allShower_startY = new std::vector<double>;
  std::vector<double> *reco_daughter_allShower_startZ = new std::vector<double>;

};


#endif //RECO_VERTEX_EVENTSELECTION_H
