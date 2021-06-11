//
// Created by Jon Sensenig on 6/5/21.
//

#include "EventSelection.h"
#include "datatypes/vector_vector.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>

EventSelection::EventSelection()
{ ; }

EventSelection::~EventSelection()
{
  std::cout << "True CEX Event Count: " << true_cex_cnt << " Reco CEX Event Count: " << reco_cex_cnt
            << " True Selected CEX Event Count: " << true_reco_cex << std::endl;

  for( const auto& reject : reject_count ) std::cout << "[Reject] Cut " << reject.first << " Count True/Reco: "
                                                     << reject.second.first << "/" << reject.second.second << std::endl;
  std::cout << std::endl;
  for( const auto& select : select_count ) std::cout << "[Select] Cut " << select.first << " Count True/Reco: "
                                                     << select.second.first << "/" << select.second.second << std::endl;

  std::cout << "Clean memory" << std::endl;
  CleanMemory();
}


void EventSelection::InitTtree( TTree *tree ) {

  // Reco
  tree->SetBranchAddress( "reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection );
  tree->SetBranchAddress( "reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange_SCE );
  tree->SetBranchAddress( "reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE );
  tree->SetBranchAddress( "reco_daughter_allTrack_dR", &reco_daughter_allTrack_dR );
  tree->SetBranchAddress( "reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection );
  tree->SetBranchAddress( "reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton );
  tree->SetBranchAddress( "reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof );
  tree->SetBranchAddress( "reco_daughter_allShower_startX", &reco_daughter_allShower_startX );
  tree->SetBranchAddress( "reco_daughter_allShower_startY", &reco_daughter_allShower_startY );
  tree->SetBranchAddress( "reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ );
  tree->SetBranchAddress( "reco_daughter_PFP_true_byHits_parPDG", &reco_daughter_PFP_true_byHits_parPDG );
  tree->SetBranchAddress( "reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG );

  // Beam-line
  tree->SetBranchAddress( "reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG );
  tree->SetBranchAddress( "reco_beam_passes_beam_cuts", &reco_beam_passes_beam_cuts );
  tree->SetBranchAddress( "beam_inst_TOF", &beam_inst_TOF );
  // Beam in TPC
  tree->SetBranchAddress( "reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits );

  // Beam Reco
  tree->SetBranchAddress( "reco_beam_calo_endX", &reco_beam_calo_endX );
  tree->SetBranchAddress( "reco_beam_calo_endY", &reco_beam_calo_endY );
  tree->SetBranchAddress( "reco_beam_calo_endZ", &reco_beam_calo_endZ );


  // Get number of all events
  size_t nevts = tree->GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  // Get number of true CEX events
  std::string cex_def(
      "true_beam_PDG==211 && true_daughter_nPi0==1 && true_daughter_nPiMinus==0 && true_daughter_nPiPlus==0 && (true_daughter_nProton>0 || true_daughter_nNeutron>0)" );
  size_t ncex = tree->GetEntries( cex_def.c_str());
  std::cout << "True CEX Events: " << ncex << std::endl;

  // Get number of true ABS events
  std::string abs_def( "true_beam_PDG==211 && true_daughter_nPi0==0 && true_daughter_nPiMinus==0 && true_daughter_nPiPlus==0" );
  size_t nabs = tree->GetEntries( abs_def.c_str());
  std::cout << "True ABS Events: " << nabs << std::endl;
}

  // pi0 event selection Jake: https://absuploads.aps.org/presentation.cfm?pid=18257
  // https://indico.fnal.gov/event/22908/contributions/70262/attachments/44197/53241/Pion_Absorption_and_Charge_Exchange_Cross-section_analysis1.pdf

bool EventSelection::SelectEvent( TTree* tree, size_t evt, Histograms& hists, bool cex_event, int nPi0 ) {

    tree->GetEntry( evt );
    true_cex = cex_event;
  true_daughter_nPi0 = nPi0;

    if( evt % 10000 == 0 ) std::cout << "Event: " << evt << std::endl;
    if( true_cex ) true_cex_cnt++;

    //.................................................................................
    //.................................................................................
    // Event Selection Cuts

    // global selection variable
    bool global_cex_select = true;

    // ***************************
    // *  Beam Primary Cuts
    // ***************************

    // FIXME temp
    int s_cnt = 0;
    std::vector<int> parPDG;
    for( size_t t = 0; t < reco_daughter_PFP_trackScore_collection->size(); t++  ) {
      if( reco_daughter_PFP_trackScore_collection->at(t) < 0.3 ) {
        s_cnt++;
        parPDG.push_back(reco_daughter_PFP_true_byHits_parPDG->at(t));
      }
    }
    //for(auto p : parPDG ) hists.th2_hists["hShowerCountPdg"]->Fill(utils::pdg::pdg2string(p).c_str(), s_cnt, 1);
    if( true_daughter_nPi0 == 1 ) hists.th2_hists["hTrueCexShower1Pi0"] -> Fill( s_cnt, true_cex );
    if( true_daughter_nPi0 == 2 ) hists.th2_hists["hTrueCexShower2Pi0"] -> Fill( s_cnt, true_cex );
    if( s_cnt != 2 && s_cnt != 3 ) {
      global_cex_select = false;
      reject_count["shower_count_55"].first += true_cex;
      reject_count["shower_count_55"].second++;
    } else {
      select_count["shower_count_55"].first += true_cex;
      select_count["shower_count_55"].second++;
    }
    if( !global_cex_select ) return false;

    // 1. "beam_tof" ========================== (beam_inst_TOF)
    // Beam-line cuts, valid and TOF
    global_cex_select = BeamTofCut();
    if( !global_cex_select ) return false;

    // 2. "beam_tpc_match" ========================== (reco_beam_passes_beam_cuts)
    // Beam-to-TPC matching cuts (XYZ, theta cuts)
    global_cex_select = BeamToTpcCut();
    if( !global_cex_select ) return false;

    // 3. "beam_endz" ========================== (reco_beam_calo_endZ, reco_beam_calibrated_dEdX, reco_beam_resRange)
    // Beam TPC track cuts (EndZ, track score, etc)
    global_cex_select = BeamEndZCut();
    if( !global_cex_select ) return false;

    // TODO The beam PIDA doesn't have much discrimination power for beam particles
    // this isn't surprising given the beam particles are higher momenta which is where all particles dE/dx:R curves converge
    // Use the average PIDA fit for beam PID (can use Chi2 fit to template, if it does better)
    // "beam_pida"

    // ***************************
    // *     Daughter Cuts
    // ***************************

    bool pion_daughter = false;
    bool dR_candidate = false;
    bool nhit_candidate = false;
    bool muon_parent = false;
    int shower_count = 0;

    for( size_t i = 0; i < reco_daughter_PFP_trackScore_collection->size(); i++ ) {

      // 4. "daughter_cnn_track_shower" ========================== (reco_daughter_PFP_trackScore_collection)
      // Beam daughter track/shower - If passed = track failed = shower
      bool track = DaughterTrackScore( i );

      // !track = shower.. count number of daughter showers in event
      if( !track ) shower_count++;

      if( track ) {
        // 5. "daughter_pida" ========================== (reco_daughter_allTrack_calibrated_dEdX_SCE, reco_daughter_allTrack_resRange_SCE)
        // Beam daughter PID: select pions so we can reject the event, i.e. we don't want events with pions in the final state
        if( DaughterPionCut( i ) ) pion_daughter = true;

        // 6 "daughter_cnn_michel" ========================== (reco_daughter_PFP_michelScore_collection)
        // Beam daughter Michel score, if passsed = "not michel" otherwise "michel"
        if( DaughterMichelScore( i ) ) muon_parent = true;

        // The rest is shower so skip if it's a track
        continue;
      }

      if( track ) std::cout << "BROKEN SELECTION!!!!!!!!!!!!!!!" << std::endl;

      // 7. "daughter_dr" ========================== (reco_beam_calo_end{X,Y,Z}, reco_daughter_allTrack_start{X,Y,Z})
      // Beam daughter Distance to Vertex: reco_beam_calo_end* SCE-corrected i.e. (daughter start) - (primary track end)
      if( DaughterDeltaRCut( i ) ) dR_candidate = true;

      // 8. "daughter_nhit" ========================== (reco_daughter_PFP_nHits)
      // Beam daughter nHits to help discriminate against non-pi0 showers
      if( DaughterNhitCut( i ) ) nhit_candidate = true;

    }

    // Cut 5 pion daughter
    if( pion_daughter ) {
      global_cex_select = false;
      reject_count["pion_daughter_04"].first += true_cex;
      reject_count["pion_daughter_04"].second++;
    } else {
      select_count["pion_daughter_04"].first += true_cex;
      select_count["pion_daughter_04"].second++;
    }
    if( !global_cex_select ) return false;

    // Cut 6 Daughter Michel
    if( muon_parent ) {
      global_cex_select = false;
      reject_count["muon_parent_05"].first += true_cex;
      reject_count["muon_parent_05"].second++;
    } else {
      select_count["muon_parent_05"].first += true_cex;
      select_count["muon_parent_05"].second++;
    }
    if( !global_cex_select ) return false;

  // Cut 6.5 Daughter showers count
  hists.th2_hists["hTrueRecoAllShowerCount"] -> Fill( shower_count, true_daughter_nPi0 );
//  if( shower_count != 2 ) {
//    global_cex_select = false;
//    reject_count["shower_count_55"].first += true_cex;
//    reject_count["shower_count_55"].second++;
//  } else {
//    select_count["shower_count_55"].first += true_cex;
//    select_count["shower_count_55"].second++;
//  }
//  if( !global_cex_select ) return false;

    // Cut 7 dR lower limit
    if( !dR_candidate ) {
      global_cex_select = false;
      reject_count["daughter_dR_06"].first += true_cex;
      reject_count["daughter_dR_06"].second++;
    } else {
      select_count["daughter_dR_06"].first += true_cex;
      select_count["daughter_dR_06"].second++;
    }
    if( !global_cex_select ) return false;

    // Cut 8 nHit lower limit
    if( !nhit_candidate ) {
      global_cex_select = false;
      reject_count["daughter_nhit_07"].first += true_cex;
      reject_count["daughter_nhit_07"].second++;
    } else {
      select_count["daughter_nhit_07"].first += true_cex;
      select_count["daughter_nhit_07"].second++;
    }
    if( !global_cex_select ) return false;

    for(auto p : parPDG ) hists.th2_hists["hShowerCountPdg"]->Fill(utils::pdg::pdg2string(p).c_str(), s_cnt, 1);

    if( global_cex_select ) reco_cex_cnt++;
    if( global_cex_select && true_cex ) true_reco_cex++;


    //.................................................................................
    //.................................................................................
    // End Event Selection Cuts

    return global_cex_select;

}


//..................................................................
bool EventSelection::BeamTofCut() {

  if( beam_inst_TOF->size() < 1 ) return true;

  if( !sel.beam_tof_cut( beam_inst_TOF->at(0) ) ) {
    reject_count["beam_tof_01"].first += true_cex;
    reject_count["beam_tof_01"].second++;
    return false;
  } else {
    select_count["beam_tof_01"].first += true_cex;
    select_count["beam_tof_01"].second++;
  }

  return true;

}

//..................................................................
bool EventSelection::BeamToTpcCut() {

  if( !reco_beam_passes_beam_cuts ) {
    reject_count["beam_tpc_match_02"].first += true_cex;
    reject_count["beam_tpc_match_02"].second++;
    return false;
  } else {
    select_count["beam_tpc_match_02"].first += true_cex;
    select_count["beam_tpc_match_02"].second++;
  }
  return true;
}

//..................................................................
bool EventSelection::BeamEndZCut() {
//  std::cout << "EndZ " << reco_beam_calo_endZ << std::endl;
  if ( !sel.beam_endz_cut( reco_beam_calo_endZ )) {
    reject_count["beam_endz_03"].first += true_cex;
    reject_count["beam_endz_03"].second++;
    return false;
  } else {
    select_count["beam_endz_03"].first += true_cex;
    select_count["beam_endz_03"].second++;
  }
  return true;
}

//..................................................................
bool EventSelection::DaughterTrackScore( size_t i ) {

  return sel.cnn_track_shower_score( reco_daughter_PFP_trackScore_collection->at(i) );

}

////..................................................................
//bool EventSelection::ShowerCountCut( size_t i ) {
//
//  return sel.cnn_track_shower_score( reco_daughter_PFP_trackScore_collection->at( i ));
//
//}

//..................................................................
bool EventSelection::DaughterMichelScore( size_t i ) {

  return !sel.daughter_michel_score(reco_daughter_PFP_michelScore_collection->at(i));

}

//..................................................................
bool EventSelection::DaughterPionCut ( size_t i ) {

  int daughter_pdg = reco_daughter_PFP_true_byHits_PDG->at(i);
  double daughter_chi2 = reco_daughter_allTrack_Chi2_proton->at(i) / reco_daughter_allTrack_Chi2_ndof->at(i);

  double avg_daughter_pida = chi2.PidaFit( reco_daughter_allTrack_calibrated_dEdX_SCE->at( i ),
                                           reco_daughter_allTrack_resRange_SCE->at( i ));
  double pi_fit = chi2.Chi2PID( reco_daughter_allTrack_calibrated_dEdX_SCE->at( i ),
                                reco_daughter_allTrack_resRange_SCE->at( i ), utils::pdg::kPdgPiP );
  double p_fit = chi2.Chi2PID( reco_daughter_allTrack_calibrated_dEdX_SCE->at( i ),
                               reco_daughter_allTrack_resRange_SCE->at( i ), utils::pdg::kPdgProton );
  double mu_fit = chi2.Chi2PID( reco_daughter_allTrack_calibrated_dEdX_SCE->at( i ),
                                reco_daughter_allTrack_resRange_SCE->at( i ), utils::pdg::kPdgMuon );
  double k_fit = chi2.Chi2PID( reco_daughter_allTrack_calibrated_dEdX_SCE->at( i ),
                               reco_daughter_allTrack_resRange_SCE->at( i ), utils::pdg::kPdgKP );
  double min_pid = std::min({pi_fit, p_fit, mu_fit, k_fit});

  bool chi2_is_muon = false;
  if( (min_pid < 9999) && (abs(min_pid - mu_fit) < 1.e3) ) chi2_is_muon = true;

  //if( avg_daughter_pida < 4. && daughter_chi2 > 100. ) {
  if( daughter_chi2 > 80. ) return true;
  else return false;


}

//..................................................................
bool EventSelection::DaughterDeltaRCut( size_t i ) {

  double dR = utils::dR( reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ,
                         reco_daughter_allShower_startX->at(i), reco_daughter_allShower_startY->at(i),
                         reco_daughter_allShower_startZ->at(i) );
  return sel.deltaR_to_vertex( dR );

}

//..................................................................
bool EventSelection::DaughterNhitCut( size_t i ) {

  return sel.daughter_nhits( reco_daughter_PFP_nHits->at(i) );

}

//..................................................................
void EventSelection::CleanMemory() {

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_PDG;
  delete beam_inst_TOF;
  delete reco_daughter_PFP_michelScore_collection;
  delete reco_daughter_allShower_startX;
  delete reco_daughter_allShower_startY;
  delete reco_daughter_allShower_startZ;
  delete reco_daughter_allTrack_dR;
  delete reco_daughter_allTrack_Chi2_proton;
  delete reco_daughter_allTrack_Chi2_ndof;
  delete reco_daughter_allTrack_calibrated_dEdX_SCE;
  delete reco_daughter_allTrack_resRange_SCE;
  delete reco_daughter_PFP_trackScore_collection;
  delete reco_daughter_PFP_true_byHits_parPDG;
  delete reco_daughter_PFP_true_byHits_PDG;

}

