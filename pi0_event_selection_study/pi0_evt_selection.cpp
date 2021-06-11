//
// Created by Jon Sensenig on 4/23/21.
//

#include "pi0_evt_selection.h"
#include "datatypes/vector_vector.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>


void run_pi0_evt_selection( const std::string& in_file, Histograms &hists, Selection &sel ) {

  pid::Chi2PIDAlg chi2;
  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  auto* tree = (TTree*)proc_file -> Get("pionana/beamana");

  // True Daughter Counts
  tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
  tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
  tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);
  tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);

  // Reco
  tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
  tree->SetBranchAddress("reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange_SCE);
  tree->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);
  tree->SetBranchAddress("reco_daughter_allTrack_dR", &reco_daughter_allTrack_dR);
  tree->SetBranchAddress("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);
  tree->SetBranchAddress("reco_daughter_allShower_startX", &reco_daughter_allShower_startX);
  tree->SetBranchAddress("reco_daughter_allShower_startY", &reco_daughter_allShower_startY);
  tree->SetBranchAddress("reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_parPDG", &reco_daughter_PFP_true_byHits_parPDG);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG);

  // Beam-line
  tree->SetBranchAddress("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG);
  tree->SetBranchAddress("reco_beam_passes_beam_cuts", &reco_beam_passes_beam_cuts);
  tree->SetBranchAddress("beam_inst_TOF", &beam_inst_TOF);
  // Beam in TPC
  tree->SetBranchAddress("reco_beam_resRange", &reco_beam_resRange);
  tree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);

  // Beam Reco
  tree->SetBranchAddress("reco_beam_calo_endX", &reco_beam_calo_endX);
  tree->SetBranchAddress("reco_beam_calo_endY", &reco_beam_calo_endY);
  tree->SetBranchAddress("reco_beam_calo_endZ", &reco_beam_calo_endZ);


  // Get number of all events
  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  // Get number of true CEX events
  std::string cex_def("true_beam_endProcess==\"pi+Inelastic\" && true_daughter_nPi0>0 && true_daughter_nPiMinus==0 && true_daughter_nPiPlus==0 && (true_daughter_nProton>0 || true_daughter_nNeutron>0)");
  size_t ncex = tree -> GetEntries( cex_def.c_str() );
  std::cout << "True CEX Events: " << ncex << std::endl;

  // pi0 event selection Jake: https://absuploads.aps.org/presentation.cfm?pid=18257

  for ( size_t evt = 0; evt < nevts; evt++ ) {
    tree->GetEntry( evt );

    if( evt % 10000 == 0 ) std::cout << "Event: " << evt << std::endl;

    // Define true CEX
    true_cex = *true_beam_endProcess == "pi+Inelastic" &&
                   true_daughter_nPi0 > 0 && true_daughter_nPiMinus == 0 &&
                   true_daughter_nPiPlus == 0 && ( true_daughter_nProton > 0 || true_daughter_nNeutron > 0 );

    if( true_cex ) true_cex_cnt++;

    //.................................................................................
    //.................................................................................
    // Event Selection Cuts

    // global selection variable
    bool global_cex_select = true;

    // ***************************
    // *  Beam Primary Cuts
    // ***************************

    // 1. "beam_tof" ========================== (beam_inst_TOF)
    // Beam-line cuts, valid and TOF
    global_cex_select = BeamTofCut( hists, sel );
    if( !global_cex_select ) continue;

    // 2. "beam_tpc_match" ========================== (reco_beam_passes_beam_cuts)
    // Beam-to-TPC matching cuts (XYZ, theta cuts)
    global_cex_select = BeamToTpcCut( hists, sel );
    if( !global_cex_select ) continue;

    // 3. "beam_endz" ========================== (reco_beam_calo_endZ, reco_beam_calibrated_dEdX, reco_beam_resRange)
    // Beam TPC track cuts (EndZ, track score, etc)
    global_cex_select = BeamEndZCut( hists, sel );
    if( !global_cex_select ) continue;

    // 3.5 TODO Do beam particle PID here!
    double avg_beam_pida = chi2.PidaFit( *reco_beam_calibrated_dEdX, *reco_beam_resRange );
    //if( avg_beam_pida < 5. ) hists.th2_hists["fit_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Pion",1);
//    hists.th2_hists["fit_beam_pid"]->Fill( avg_beam_pida, utils::pdg::pdg2string(reco_beam_true_byHits_PDG).c_str(), 1 );

    // TODO The beam PIDA doesn't have much discrimination power for beam particles
    // this isn't surprising given the beam particles are higher momenta which is where all particles dE/dx:R curves converge
    // Use the average PIDA fit for beam PID (can use Chi2 fit to template, if it does better)
    // "beam_pida"
    if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgProton ) hists.th1_hists["beam_avg_pida_p"]->Fill( avg_beam_pida );
    else if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgPiP ) hists.th1_hists["beam_avg_pida_pi"]->Fill( avg_beam_pida );
    else if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgElectron ) hists.th1_hists["beam_avg_pida_e"]->Fill( avg_beam_pida );
    else if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgMuon ) hists.th1_hists["beam_avg_pida_mu"]->Fill( avg_beam_pida );
    else hists.th1_hists["beam_avg_pida_other"]->Fill( avg_beam_pida );

    // ***************************
    // *     Daughter Cuts
    // ***************************

    bool pion_daughter = false;
    bool dR_candidate = false;
    bool nhit_candidate = false;
    bool muon_parent = false;

    for( size_t i = 0; i < reco_daughter_PFP_trackScore_collection->size(); i++ ) {

      // 4. "daughter_cnn_track_shower" ========================== (reco_daughter_PFP_trackScore_collection)
      // Beam daughter track/shower - If passed = track failed = shower
      bool track = DaughterTrackScore( hists, sel, i );

      if( track ) {
        // 5. "daughter_pida" ========================== (reco_daughter_allTrack_calibrated_dEdX_SCE, reco_daughter_allTrack_resRange_SCE)
        // Beam daughter PID: select pions so we can reject the event, i.e. we don't want events with pions in the final state
         if( DaughterPionCut( hists, sel, chi2, i ) ) pion_daughter = true;

        // 6 "daughter_cnn_michel" ========================== (reco_daughter_PFP_michelScore_collection)
        // Beam daughter Michel score, if passsed = "not michel" otherwise "michel"
         if( DaughterMichelScore( hists, sel, i ) ) muon_parent = true;

        // The rest is shower so skip if it's a track
        continue;
      }

      if( track ) std::cout << "BROKEN SELECTION!!!!!!!!!!!!!!!" << std::endl;

      // 7. "daughter_dr" ========================== (reco_beam_calo_end{X,Y,Z}, reco_daughter_allTrack_start{X,Y,Z})
      // Beam daughter Distance to Vertex: reco_beam_calo_end* SCE-corrected i.e. (daughter start) - (primary track end)
      if( DaughterDeltaRCut( hists, sel, i ) ) dR_candidate = true;

      // 8. "daughter_nhit" ========================== (reco_daughter_PFP_nHits)
      // Beam daughter nHits to help discriminate against non-pi0 showers
      if( DaughterNhitCut( hists, sel, i ) ) nhit_candidate = true;

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
    if( !global_cex_select ) continue;

    // Cut 6 Daughter Michel
    if( muon_parent ) {
      global_cex_select = false;
      reject_count["muon_parent_05"].first += true_cex;
      reject_count["muon_parent_05"].second++;
    } else {
      select_count["muon_parent_05"].first += true_cex;
      select_count["muon_parent_05"].second++;
    }
    if( !global_cex_select ) continue;

    // Cut 7 dR lower limit
    if( !dR_candidate ) {
      global_cex_select = false;
      reject_count["daughter_dR_06"].first += true_cex;
      reject_count["daughter_dR_06"].second++;
    } else {
      select_count["daughter_dR_06"].first += true_cex;
      select_count["daughter_dR_06"].second++;
    }
    if( !global_cex_select ) continue;

    // Cut 8 nHit lower limit
    if( !nhit_candidate ) {
      global_cex_select = false;
      reject_count["daughter_nhit_07"].first += true_cex;
      reject_count["daughter_nhit_07"].second++;
    } else {
      select_count["daughter_nhit_07"].first += true_cex;
      select_count["daughter_nhit_07"].second++;
    }
    if( !global_cex_select ) continue;



    if( global_cex_select ) reco_cex_cnt++;
    if( global_cex_select && true_cex ) true_reco_cex++;


    //.................................................................................
    //.................................................................................
    // End Event Selection Cuts

  } // end evt loop

  std::cout << "True CEX Event Count: " << true_cex_cnt << " Reco CEX Event Count: " << reco_cex_cnt
            << " True Selected CEX Event Count: " << true_reco_cex << std::endl;

  for( const auto& reject : reject_count ) std::cout << "[Reject] Cut " << reject.first << " Count True/Reco: "
                                              << reject.second.first << "/" << reject.second.second << std::endl;
  std::cout << std::endl;
  for( const auto& select : select_count ) std::cout << "[Select] Cut " << select.first << " Count True/Reco: "
                                              << select.second.first << "/" << select.second.second << std::endl;

  CleanMemory();
  proc_file -> Close();

}

//..................................................................
bool BeamTofCut( Histograms& hists, Selection &sel ) {

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
bool BeamToTpcCut( Histograms& hists, Selection& sel ) {

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
bool BeamEndZCut( Histograms& hists, Selection& sel ) {

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
bool DaughterTrackScore( Histograms& hists, Selection& sel, size_t i ) {

  return sel.cnn_track_shower_score( reco_daughter_PFP_trackScore_collection->at(i) );

}

//..................................................................
bool DaughterMichelScore( Histograms& hists, Selection& sel, size_t i ) {

  return !sel.daughter_michel_score(reco_daughter_PFP_michelScore_collection->at(i));

}

//..................................................................
bool DaughterPionCut ( Histograms& hists, Selection& sel, pid::Chi2PIDAlg& chi2, size_t i ) {

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

  hists.th2_hists["true_beam_pdg_michel_score"]->Fill(
      reco_daughter_PFP_michelScore_collection->at(i), utils::pdg::pdg2string(reco_beam_true_byHits_PDG).c_str(), 1);

  if( (min_pid < 9999) && (abs(min_pid - pi_fit) < 1.e3) ) hists.th2_hists["chi2_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Pion",1);
  else hists.th2_hists["chi2_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Other",1);

  if( abs(daughter_pdg) == utils::pdg::kPdgProton ) { hists.th1_hists["daughter_avg_pida_p"]->Fill( avg_daughter_pida ); hists.th1_hists["daughter_chi2_p"]->Fill(daughter_chi2); }
  else if( abs(daughter_pdg) == utils::pdg::kPdgPiP ) { hists.th1_hists["daughter_avg_pida_pi"]->Fill( avg_daughter_pida ); hists.th1_hists["daughter_chi2_pip"]->Fill(daughter_chi2); }
  else if( abs(daughter_pdg) == utils::pdg::kPdgElectron ) hists.th1_hists["daughter_avg_pida_e"]->Fill( avg_daughter_pida );
  else if( abs(daughter_pdg) == utils::pdg::kPdgMuon ) hists.th1_hists["daughter_avg_pida_mu"]->Fill( avg_daughter_pida );
  else hists.th1_hists["daughter_avg_pida_other"]->Fill( avg_daughter_pida );

  bool chi2_is_muon = false;
  if( (min_pid < 9999) && (abs(min_pid - mu_fit) < 1.e3) ) chi2_is_muon = true;

  //if( avg_daughter_pida < 4. && daughter_chi2 > 100. ) {
  if( daughter_chi2 > 80. ) {
      hists.th2_hists["fit_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Pion",1);
      return true;
  } else {
    hists.th2_hists["fit_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Other",1);
    return false;
  }

  // True = pion, False = not pion
  //return sel.pida_cut( avg_daughter_pida, utils::pdg::kPdgPiP );

}

//..................................................................
bool DaughterDeltaRCut( Histograms& hists, Selection& sel, size_t i ) {

  double dR = utils::dR( reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ,
                         reco_daughter_allShower_startX->at(i), reco_daughter_allShower_startY->at(i), reco_daughter_allShower_startZ->at(i) );
  return sel.deltaR_to_vertex( dR );

}

//..................................................................
bool DaughterNhitCut( Histograms& hists, Selection& sel, size_t i ) {

  return sel.daughter_nhits( reco_daughter_PFP_nHits->at(i) );

}

//..................................................................
void CleanMemory() {

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_PDG;
  delete beam_inst_TOF;
  delete reco_beam_calibrated_dEdX;
  delete reco_beam_resRange;
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

int main() {

  std::string input_file = "../../../pionana_Prod4_mc_1GeV_1_14_21.root";
  TString output_file = "out.root";
  std::string hists_config = "../hists.json";

  // Configure histograms
  Histograms hists;
  hists.ConfigureHistos( hists_config );

  // Initialize selection parameters
  Selection sel;
  sel.set_selection_parameters( "../selection_params.json" );

  std::cout << "Starting pi0 xsec study!" << std::endl;
  run_pi0_evt_selection( input_file, hists, sel );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}
