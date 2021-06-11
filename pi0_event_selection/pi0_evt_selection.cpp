//
// Created by Jon Sensenig on 4/23/21.
//

#include "pi0_evt_selection.h"
#include "../utilities/Utilities.hpp"
#include "../utilities/Chi2PIDAlg.hpp"
#include "datatypes/vector_vector.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>


void run_pi0_evt_selection( std::string in_file, Histograms &hists, Selection &sel ) {

  pid::Chi2PIDAlg chi2;
  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  TTree* tree = (TTree*)proc_file -> Get("pionana/beamana");

  // True
  tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
  tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
  tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);
  tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID);
  tree->SetBranchAddress("true_beam_daughter_startP", &true_beam_daughter_startP);
  tree->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess);
  tree->SetBranchAddress("true_beam_endX", &true_beam_endX);
  tree->SetBranchAddress("true_beam_endY", &true_beam_endY);
  tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ);
  tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
  tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
  tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);
  tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

  tree->SetBranchAddress("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG);
  tree->SetBranchAddress("true_beam_Pi0_decay_startP", &true_beam_Pi0_decay_startP);
  tree->SetBranchAddress("true_beam_Pi0_decay_parID", &true_beam_Pi0_decay_parID);
  tree->SetBranchAddress("true_beam_Pi0_decay_ID", &true_beam_Pi0_decay_ID);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
  tree->SetBranchAddress("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
  tree->SetBranchAddress("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
  tree->SetBranchAddress("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);
  tree->SetBranchAddress("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);

  // Reco
  tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
  tree->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
  tree->SetBranchAddress("reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange_SCE);
  tree->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);
  tree->SetBranchAddress("reco_daughter_allTrack_dR", &reco_daughter_allTrack_dR);
  tree->SetBranchAddress("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX);
  tree->SetBranchAddress("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY);
  tree->SetBranchAddress("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ);
  tree->SetBranchAddress("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);

  tree->SetBranchAddress("reco_daughter_allShower_dirX", &reco_daughter_allShower_dirX);
  tree->SetBranchAddress("reco_daughter_allShower_dirY", &reco_daughter_allShower_dirY);
  tree->SetBranchAddress("reco_daughter_allShower_dirZ", &reco_daughter_allShower_dirZ);
  tree->SetBranchAddress("reco_daughter_allShower_energy", &reco_daughter_allShower_energy);
  tree->SetBranchAddress("reco_daughter_allShower_startX", &reco_daughter_allShower_startX);
  tree->SetBranchAddress("reco_daughter_allShower_startY", &reco_daughter_allShower_startY);
  tree->SetBranchAddress("reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_startPx", &reco_daughter_PFP_true_byHits_startPx);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_startPy", &reco_daughter_PFP_true_byHits_startPy);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_startPz", &reco_daughter_PFP_true_byHits_startPz);
  tree->SetBranchAddress("reco_daughter_allShower_ID", &reco_daughter_allShower_ID);
  tree->SetBranchAddress("true_beam_daughter_reco_byHits_allShower_ID", &true_beam_daughter_reco_byHits_allShower_ID);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_parPDG", &reco_daughter_PFP_true_byHits_parPDG);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_parID", &reco_daughter_PFP_true_byHits_parID);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_endX", &reco_daughter_PFP_true_byHits_endX);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_endY", &reco_daughter_PFP_true_byHits_endY);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_endZ", &reco_daughter_PFP_true_byHits_endZ);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_endProcess", &reco_daughter_PFP_true_byHits_endProcess);

  // Beam-line
  tree->SetBranchAddress("beam_inst_P", &beam_inst_P);
  tree->SetBranchAddress("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG);
  tree->SetBranchAddress("reco_beam_passes_beam_cuts", &reco_beam_passes_beam_cuts);
  tree->SetBranchAddress("beam_inst_TOF", &beam_inst_TOF);
  // Beam in TPC
  tree->SetBranchAddress("reco_beam_resRange", &reco_beam_resRange);
  tree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);

  tree->SetBranchAddress("true_beam_Pi0_decay_reco_byHits_allShower_ID", &true_beam_Pi0_decay_reco_byHits_allShower_ID);

  // Beam Reco
  tree->SetBranchAddress("reco_beam_calo_endX", &reco_beam_calo_endX);
  tree->SetBranchAddress("reco_beam_calo_endY", &reco_beam_calo_endY);
  tree->SetBranchAddress("reco_beam_calo_endZ", &reco_beam_calo_endZ);
  tree->SetBranchAddress("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
  tree->SetBranchAddress("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
  tree->SetBranchAddress("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);

  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  // pi0 event selection Jake: https://absuploads.aps.org/presentation.cfm?pid=18257

  for ( size_t evt = 0; evt < nevts; evt++ ) {
    tree->GetEntry( evt );

    double pi0_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );
    int pi0_idx = utils::FindIndex<int>( *true_beam_daughter_PDG, utils::pdg::kPdgPi0 );

    // Define true CEX
    bool true_cex = true_beam_PDG == 211 &&
                   true_daughter_nPi0 > 0 && true_daughter_nPiMinus == 0 &&
                   true_daughter_nPiPlus == 0;

    if( true_cex ) true_cex_cnt++;

    //.................................................................................
    //.................................................................................
    // Event Selection Cuts

    // global selection variable
    bool global_cex_select = true;

    // ***************************
    // *     Beam Primary
    // ***************************

    // 1. "beam_tof" ========================== (beam_inst_TOF)
    // Beam-line cuts, valid and TOF
    if( !sel.beam_tof_cut( beam_inst_TOF->at(0) ) ) { global_cex_select = false; reject_count["beam_tof"]++; }
    else { select_count["beam_tof"]++; }

    // 2. "beam_tpc_match" ========================== (reco_beam_passes_beam_cuts)
    // Beam-to-TPC matching cuts (XYZ, theta cuts)
    if( !reco_beam_passes_beam_cuts ) { global_cex_select = false; reject_count["beam_tpc_match"]++; }
    else { select_count["beam_tpc_match"]++; }

    // 3. "beam_endz" ========================== (reco_beam_calo_endZ, reco_beam_calibrated_dEdX, reco_beam_resRange)
    // Beam TPC track cuts (EndZ, track score, etc)
    if( !sel.beam_endz_cut( reco_beam_calo_endZ ) ) { global_cex_select = false; reject_count["beam_endz"]++; }
    else { select_count["beam_endz"]++; }

    // TODO The beam PIDA doesn't have much discrimination power for beam particles
    // this isn't surprising given the beam particles are higher momenta which is where all particles dE/dx:R curves converge
    // Use the average PIDA fit for beam PID (can use Chi2 fit to template, if it does better)
    // "beam_pida"
    double avg_beam_pida = chi2.PidaFit( *reco_beam_calibrated_dEdX, *reco_beam_resRange );
//    if( !sel.pida_cut( avg_beam_pida, utils::pdg::kPdgPiP ) ) { global_cex_select = false; reject_count["beam_pida"]++; }
//    else { select_count["beam_pida"]++; }

    if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgProton ) hists.th1_hists["beam_avg_pida_p"]->Fill( avg_beam_pida );
    else if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgPiP ) hists.th1_hists["beam_avg_pida_pi"]->Fill( avg_beam_pida );
    else if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgElectron ) hists.th1_hists["beam_avg_pida_e"]->Fill( avg_beam_pida );
    else if( abs(reco_beam_true_byHits_PDG) == utils::pdg::kPdgMuon ) hists.th1_hists["beam_avg_pida_mu"]->Fill( avg_beam_pida );
    else hists.th1_hists["beam_avg_pida_other"]->Fill( avg_beam_pida );

    // ***************************
    // *      DAUGHTERS
    // ***************************

    bool pion_daughter = false;
    bool dR_candidate = false;
    bool nhit_candidate = false;
    bool muon_parent = false;

    for( size_t i = 0; i < reco_daughter_PFP_trackScore_collection->size(); i++ ) {

      int daughter_pdg = reco_daughter_PFP_true_byHits_PDG->at(i);

      // 4. "daughter_cnn_track_shower" ========================== (reco_daughter_PFP_trackScore_collection)
      // Beam daughter track/shower - If passed = track failed = shower
      bool track = sel.cnn_track_shower_score( reco_daughter_PFP_trackScore_collection->at(i) );

      if( track ) {
        // 5. "daughter_pida" ========================== (reco_daughter_allTrack_calibrated_dEdX_SCE, reco_daughter_allTrack_resRange_SCE)
        // Beam daughter PID: select pions so we can reject the event, i.e. we don't want events with pions in the final state
        double avg_daughter_pida = chi2.PidaFit( reco_daughter_allTrack_calibrated_dEdX_SCE->at( i ),
                                                 reco_daughter_allTrack_resRange_SCE->at( i ));
        if ( !sel.pida_cut( avg_daughter_pida, utils::pdg::kPdgPiP ) ) pion_daughter = true;

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

        if( !sel.daughter_michel_score(reco_daughter_PFP_michelScore_collection->at(i))  ) muon_parent = true;

        if( (min_pid < 9999) && (abs(min_pid - pi_fit) < 1.e3) ) hists.th2_hists["chi2_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Pion",1);
        else hists.th2_hists["chi2_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Other",1);

        if( abs(daughter_pdg) == utils::pdg::kPdgProton ) hists.th1_hists["daughter_avg_pida_p"]->Fill( avg_daughter_pida );
        else if( abs(daughter_pdg) == utils::pdg::kPdgPiP ) hists.th1_hists["daughter_avg_pida_pi"]->Fill( avg_daughter_pida );
        else if( abs(daughter_pdg) == utils::pdg::kPdgElectron ) hists.th1_hists["daughter_avg_pida_e"]->Fill( avg_daughter_pida );
        else if( abs(daughter_pdg) == utils::pdg::kPdgMuon ) hists.th1_hists["daughter_avg_pida_mu"]->Fill( avg_daughter_pida );
        else hists.th1_hists["daughter_avg_pida_other"]->Fill( avg_daughter_pida );

        bool chi2_is_muon = false;
        if( (min_pid < 9999) && (abs(min_pid - mu_fit) < 1.e3) ) chi2_is_muon = true;
        if( avg_daughter_pida < 5. ) hists.th2_hists["fit_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Pion",1);
        else hists.th2_hists["fit_daughter_pid"]->Fill(utils::pdg::pdg2string(daughter_pdg).c_str(),"Other",1);

        // The rest is shower so skip if it's a track
        continue;
      }

      if( track ) std::cout << "BROKEN SELECTION!!!!!!!!!!!!!!!" << std::endl;

      // 6. "daughter_dr" ========================== (reco_beam_calo_end{X,Y,Z}, reco_daughter_allTrack_start{X,Y,Z})
      // Beam daughter Distance to Vertex: reco_beam_calo_end* SCE-corrected i.e. (daughter start) - (primary track end)
      double dR = utils::dR( reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ,
          reco_daughter_allShower_startX->at(i), reco_daughter_allShower_startY->at(i), reco_daughter_allShower_startZ->at(i) );
      if( sel.deltaR_to_vertex( dR ) ) dR_candidate = true;

      // 7. "daughter_nhit" ========================== (reco_daughter_PFP_nHits)
      // Beam daughter nHits to help discriminate against non-pi0 showers
      if( sel.daughter_nhits( reco_daughter_PFP_nHits->at(i) ) ) nhit_candidate = true;

    }

    if( pion_daughter ) { global_cex_select = false; reject_count["pion_daughter"]++; }
    else { select_count["pion_daughter"]++; }

    if( !dR_candidate ) { global_cex_select = false; reject_count["daughter_dR"]++; }
    else { select_count["daughter_dR"]++; }

    if( !nhit_candidate ) { global_cex_select = false; reject_count["daughter_nhit"]++; }
    else { select_count["daughter_nhit"]++; }

    if( muon_parent ) { global_cex_select = false; reject_count["muon_parent"]++; }
    else { select_count["muon_parent"]++; }

    if( global_cex_select ) reco_cex_cnt++;
    if( global_cex_select && true_cex ) true_reco_cex++;


    //.................................................................................
    //.................................................................................
    // End Event Selection Cuts

//    if( !true_cex ) continue;
//
//    // Loop over pi0 decay gammas
//    for( size_t i = 0; i < true_beam_Pi0_decay_parID->size()/2; i++ ) {
//      size_t idx = i * 2; // Gamma from the same pi0 are adjacent in index, i.e., i and i+1
//      if( true_beam_Pi0_decay_parID -> at( idx ) != true_beam_Pi0_decay_parID -> at( idx+1 ) ) {
//        std::cout << "Something's wrong, decay gamma IDs do not match!" << std::endl;
//        std::cout << "PDG 1/2 " << true_beam_Pi0_decay_PDG->at(idx) << "/" << true_beam_Pi0_decay_PDG->at(idx+1) << std::endl;
//      }
//      // Decay gammas momentum and opening angle
//      double gamma1 = true_beam_Pi0_decay_startP -> at( idx ) * 1.e3;
//      double gamma2 = true_beam_Pi0_decay_startP -> at( idx+1 ) * 1.e3;
//      double open_angle_gg = open_angle( true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
//                                 true_beam_Pi0_decay_startPz->at(idx), true_beam_Pi0_decay_startPx->at(idx+1),
//                                 true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );
//
//      double beam_end_pos = utils::Distance(true_beam_endX, true_beam_endY, true_beam_endZ);
//      double gamma_pos1 = utils::Distance(true_beam_Pi0_decay_startX->at(idx), true_beam_Pi0_decay_startY->at(idx),
//                                          true_beam_Pi0_decay_startZ->at(idx));
//      double gamma_pos2 = utils::Distance(true_beam_Pi0_decay_startX->at(idx+1), true_beam_Pi0_decay_startY->at(idx+1),
//                                          true_beam_Pi0_decay_startZ->at(idx+1));
//
//      hists.th1_hists["hGammaOpenAngle"] -> Fill( open_angle_gg );
//      hists.th1_hists["hLeadGammaP"] -> Fill( std::max( gamma1, gamma2 ) );
//      hists.th1_hists["hSubLeadGammaP"] -> Fill( std::min( gamma1, gamma2 ) );
//      hists.th2_hists["hGammaP"] -> Fill( std::max( gamma1, gamma2 ), std::min( gamma1, gamma2 ) );
//      hists.th2_hists["hGammaR"] -> Fill( gamma_pos1-beam_end_pos, gamma_pos2-beam_end_pos );
//      hists.th2_hists["hGammaLen"] -> Fill( true_beam_Pi0_decay_len->at(idx), true_beam_Pi0_decay_len->at(idx+1) );
//
//      if( true_daughter_nPi0 > 1 ) continue;
//
//      // Angular momentum calculation
//      double angle_deg = open_angle_gg * TMath::RadToDeg();
//      // Momentum from polynomial fit Ref https://arxiv.org/pdf/1511.00941.pdf
//      double p_poly = 2202.3 - 94.9*angle_deg + 2.1*pow(angle_deg, 2) - 0.025*pow(angle_deg, 3) + 0.00017*pow(angle_deg, 4)
//               - 6.0e-7*pow(angle_deg, 5) + 8.5e-10*pow(angle_deg, 6);
//      double p_root = sqrt( pow(gamma1,2) + pow(gamma2,2) + 2*gamma1*gamma2*cos(open_angle_gg) );
//
//      // Energy true and calculated
//      double pi0_true_energy = utils::CalculateE( true_beam_daughter_startP->at( pi0_idx ) * 1.e3, pi0_mass );
//      double pi0_calc_energy = E_gamma_gammma( gamma1, gamma2, open_angle_gg );
//
//      hists.th1_hists["hPolyPi0PError"] -> Fill( p_poly / (true_beam_daughter_startP->at( pi0_idx )*1.e3) );
//      hists.th1_hists["hRootPi0PError"] -> Fill( p_root / (true_beam_daughter_startP->at( pi0_idx )*1.e3) );
//      hists.th2_hists["hPi0TrueCalcEnergy"] -> Fill( pi0_calc_energy, pi0_true_energy );
//
//      // Gamma angle wrt pi0
//      double gamma1_open_angle = open_angle( true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
//                                        true_beam_daughter_startPz->at(pi0_idx), true_beam_Pi0_decay_startPx->at(idx),
//                                        true_beam_Pi0_decay_startPy->at(idx), true_beam_Pi0_decay_startPz->at(idx) );
//      double gamma2_open_angle = open_angle( true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
//                                        true_beam_daughter_startPz->at(pi0_idx), true_beam_Pi0_decay_startPx->at(idx+1),
//                                        true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );
//
//      // Gamma angle wrt incoming pi+
//      double gamma1_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
//                                             true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
//                                             true_beam_Pi0_decay_startPz->at(idx) );
//      double gamma2_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
//                                             true_beam_Pi0_decay_startPx->at(idx+1), true_beam_Pi0_decay_startPy->at(idx+1),
//                                             true_beam_Pi0_decay_startPz->at(idx+1) );
//
//      // pi0 angle wrt incoming pi+
//      double pi0_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
//                                     true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
//                                     true_beam_daughter_startPz->at(pi0_idx) );
//
//      hists.th1_hists["hGammaDiff"] -> Fill( gamma1_open_angle - gamma2_open_angle );
//      hists.th2_hists["hGammaPi0Angle"] -> Fill( gamma1_open_angle*TMath::RadToDeg(), gamma2_open_angle*TMath::RadToDeg() );
//
//      // angles wrt incoming pi+
//      double cos_pi0_angle_calc = theta_gamma_gammma( gamma1, gamma2, open_angle_gg, gamma1_angle, gamma2_angle );
//      double cos_pi0_angle_mc = cos( pi0_angle );
//      hists.th2_hists["hPi0TrueCalc"] -> Fill( cos_pi0_angle_calc, cos_pi0_angle_mc );
//      hists.th1_hists["hPi0CalcAngleDiff"] -> Fill( cos_pi0_angle_calc - cos_pi0_angle_mc );
//    }
//
//    // 1. Loop over the daughter reco pi0 decay gammas
//    // 2. Get reco daughter gamma true (by hits) ID
//    // 3. Find the true gamma (its index) by matching ID
//
//    std::vector<int> gamma_idx = get_reco_gamma_index( *reco_daughter_PFP_true_byHits_PDG, *reco_daughter_PFP_true_byHits_parPDG );
//    hists.th2_hists["hTrueRecoGammaCount"] -> Fill( gamma_idx.size(), true_daughter_nPi0 );
//
//    for( auto &reco_idx : gamma_idx ) {
//      TVector3 gamma( reco_daughter_PFP_true_byHits_endX->at(reco_idx), reco_daughter_PFP_true_byHits_endY->at(reco_idx),
//                      reco_daughter_PFP_true_byHits_endZ->at(reco_idx));
//      TVector3 gammaP( reco_daughter_PFP_true_byHits_startPx->at(reco_idx), reco_daughter_PFP_true_byHits_startPy->at(reco_idx),
//                       reco_daughter_PFP_true_byHits_startPz->at(reco_idx));
//      for( size_t i = 0; i < reco_daughter_allShower_startX->size(); i++ ){
//        TVector3 shower( reco_daughter_allShower_startX->at(i),reco_daughter_allShower_startY->at(i),
//                         reco_daughter_allShower_startZ->at(i));
//        TVector3 showerDir( reco_daughter_allShower_dirX->at(i),reco_daughter_allShower_dirY->at(i),
//                            reco_daughter_allShower_dirZ->at(i));
//      }
//    }
//
//    // Get the beam direction, to be used below
//    TVector3 true_beam_dir( true_beam_endPx, true_beam_endPy, true_beam_endPz );
//    TVector3 reco_beam_dir( reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ );
//
//    // Get the indices of the pi0 decay gammas
//    std::vector<int> true_gamma_idx = get_true_reco_gamma_index( *reco_daughter_PFP_true_byHits_PDG, *reco_daughter_PFP_true_byHits_ID,
//                                                                 *reco_daughter_PFP_true_byHits_parPDG );
//
//    // Create map with a the 2 gammas from pi0 -> gg
//    std::map<int, std::vector<int>> gamma_pair_map; // map = <index, {leading idx, sub-leading idx}>
//    for( int idx : true_gamma_idx ) {
//      int id = reco_daughter_PFP_true_byHits_ID -> at(idx);
//      gamma_pair_map[ id ].push_back( idx );
//    }
//
//    for( auto & gpair : gamma_pair_map ) {
//      int idx1 = gpair.second.at(0); int idx2 = gpair.second.at(1);
//      // If either shower energy is invalid skip
//      if( reco_daughter_allShower_energy->at(idx1) == -999 || reco_daughter_allShower_energy->at(idx2) == -999 ) continue;
//      TVector3 shower_reco_dir1( reco_daughter_allShower_dirX->at(idx1), reco_daughter_allShower_dirY->at(idx1),
//                                reco_daughter_allShower_dirZ->at(idx1) );
//      TVector3 shower_reco_dir2( reco_daughter_allShower_dirX->at(idx2), reco_daughter_allShower_dirY->at(idx2),
//                                 reco_daughter_allShower_dirZ->at(idx2) );
//      TVector3 true_gamma_P1( reco_daughter_PFP_true_byHits_startPx->at(idx1), reco_daughter_PFP_true_byHits_startPy->at(idx1),
//                             reco_daughter_PFP_true_byHits_startPz->at(idx1));
//      TVector3 true_gamma_P2( reco_daughter_PFP_true_byHits_startPx->at(idx2), reco_daughter_PFP_true_byHits_startPy->at(idx2),
//                             reco_daughter_PFP_true_byHits_startPz->at(idx2));
//      // Reco quantities
//      double shower_energy1 = reco_daughter_allShower_energy->at(idx1);
//      double shower_energy2 = reco_daughter_allShower_energy->at(idx2);
//      double reco_theta_gg = shower_reco_dir1.Angle( shower_reco_dir2 ); // theta [rad]
//      double reco_theta1 = reco_beam_dir.Angle( shower_reco_dir1 );
//      double reco_theta2 = reco_beam_dir.Angle( shower_reco_dir2 );
//      // True quantities
//      double true_gamma_energy1 = true_gamma_P1.Mag() * 1.e3;
//      double true_gamma_energy2 = true_gamma_P2.Mag() * 1.e3;
//      double true_theta_gg = true_gamma_P1.Angle( true_gamma_P2 );
//      double true_theta1 = true_beam_dir.Angle( true_gamma_P1 );
//      double true_theta2 = true_beam_dir.Angle( true_gamma_P2 );
//
//      double reco_pi0_energy = E_gamma_gammma( shower_energy1, shower_energy2, reco_theta_gg );
//      double reco_pi0_theta = theta_gamma_gammma( shower_energy1, shower_energy2, reco_theta_gg, reco_theta1, reco_theta2 );
//      double true_pi0_energy = E_gamma_gammma( true_gamma_energy1, true_gamma_energy2, true_theta_gg );
//      double true_pi0_theta = theta_gamma_gammma( true_gamma_energy1, true_gamma_energy2, true_theta_gg, true_theta1, true_theta2 );
//      hists.th2_hists["hpi0TrueRecoEnergy"] -> Fill( reco_pi0_energy, true_pi0_energy );
//      hists.th2_hists["hpi0TrueRecoTheta"] -> Fill( reco_pi0_theta * TMath::RadToDeg(), true_pi0_theta * TMath::RadToDeg() );
//
//    }
//
//    // Loop over showers from pi0s which have both showers reconstructed
//    for( int k : true_gamma_idx ) {
//      TVector3 gamma_true_mom( reco_daughter_PFP_true_byHits_startPx->at(k)*1.e3, reco_daughter_PFP_true_byHits_startPy->at(k)*1.e3,
//                      reco_daughter_PFP_true_byHits_startPz->at(k)*1.e3 );
//      TVector3 shower_reco_dir( reco_daughter_allShower_dirX->at(k), reco_daughter_allShower_dirY->at(k),
//                                reco_daughter_allShower_dirZ->at(k) );
//      hists.th2_hists["hShowerGammaEnergy"] -> Fill( reco_daughter_allShower_energy->at(k), gamma_true_mom.Mag() );
//      // true - reco
//      double true_angle = true_beam_dir.Angle( gamma_true_mom ) * TMath::RadToDeg();
//      double energy_diff = gamma_true_mom.Mag() - reco_daughter_allShower_energy->at(k);
//      double reco_angle_tbeam = true_beam_dir.Angle( shower_reco_dir ) * TMath::RadToDeg();
//      double reco_angle_rbeam = reco_beam_dir.Angle( shower_reco_dir ) * TMath::RadToDeg();
//      double angle_diff = true_angle - reco_angle_tbeam;
//      hists.th1_hists["hShowerEnergyDiff"] -> Fill( energy_diff );
//      hists.th1_hists["hShowerAngleTrueBeamDiff"] -> Fill( true_angle - reco_angle_tbeam );
//      hists.th1_hists["hShowerAngleRecoBeamDiff"] -> Fill( true_angle - reco_angle_rbeam );
//      // 2D energy difference
//      hists.th2_hists["hShowerEnergyDifffEnergy"] -> Fill( gamma_true_mom.Mag(), energy_diff );
//      hists.th2_hists["hShowerEnergyDifffAngle"] -> Fill( true_angle, energy_diff );
//      // 2D angle difference
//      hists.th2_hists["hShowerAngleDifffEnergy"] -> Fill( gamma_true_mom.Mag(), angle_diff );
//      hists.th2_hists["hShowerAngleDifffAngle"] -> Fill( true_angle, angle_diff );
//    }

    if( evt % 10000 == 0 ) std::cout << "Event: " << evt << std::endl;

  } // end evt loop

  std::cout << "True CEX Event Count: " << true_cex_cnt << " Reco CEX Event Count: " << reco_cex_cnt
            << " True Selected CEX Event Count: " << true_reco_cex << std::endl;

  for( auto reject : reject_count ) std::cout << "[Reject] Cut " << reject.first << " Count: " << reject.second << std::endl;
  for( auto select : select_count ) std::cout << "[Select] Cut " << select.first << " Count: " << select.second << std::endl;

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_startP;
  delete true_beam_daughter_PDG;
  delete true_beam_daughter_ID;
  delete true_beam_daughter_startPx;
  delete true_beam_daughter_startPy;
  delete true_beam_daughter_startPz;
  delete beam_inst_TOF;
  delete reco_beam_calibrated_dEdX;
  delete reco_beam_resRange;
  delete reco_daughter_allTrack_startX;
  delete reco_daughter_allTrack_startY;
  delete reco_daughter_allTrack_startZ;
  delete reco_daughter_allTrack_ID;
  delete reco_daughter_PFP_michelScore_collection;
  delete reco_daughter_allShower_dirX;
  delete reco_daughter_allShower_dirY;
  delete reco_daughter_allShower_dirZ;
  delete reco_daughter_allShower_startX;
  delete reco_daughter_allShower_startY;
  delete reco_daughter_allShower_startZ;
  delete reco_daughter_allShower_energy;
  delete reco_daughter_allShower_ID;
  delete reco_daughter_allTrack_dR;
  delete reco_daughter_allTrack_calibrated_dEdX_SCE;
  delete reco_daughter_allTrack_resRange_SCE;
  delete true_beam_daughter_reco_byHits_allShower_ID;
  delete true_beam_Pi0_decay_reco_byHits_allShower_ID;
  delete reco_daughter_PFP_trackScore_collection;
  delete reco_daughter_PFP_true_byHits_parPDG;
  delete reco_daughter_PFP_true_byHits_PDG;
  delete reco_daughter_PFP_true_byHits_parID;
  delete reco_daughter_PFP_true_byHits_ID;
  delete reco_daughter_PFP_true_byHits_endX;
  delete reco_daughter_PFP_true_byHits_endY;
  delete reco_daughter_PFP_true_byHits_endZ;
  delete reco_daughter_PFP_true_byHits_startPx;
  delete reco_daughter_PFP_true_byHits_startPy;
  delete reco_daughter_PFP_true_byHits_startPz;
  delete reco_daughter_PFP_true_byHits_endProcess;
  delete true_beam_Pi0_decay_PDG;
  delete true_beam_Pi0_decay_startP;
  delete true_beam_Pi0_decay_parID;
  delete true_beam_Pi0_decay_ID;
  delete true_beam_Pi0_decay_startPx;
  delete true_beam_Pi0_decay_startPy;
  delete true_beam_Pi0_decay_startPz;
  delete true_beam_Pi0_decay_startX;
  delete true_beam_Pi0_decay_startY;
  delete true_beam_Pi0_decay_startZ;
  delete true_beam_Pi0_decay_len;
  proc_file -> Close();

}

double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 ) {

  TVector3 in( px1, py1, pz1 );
  TVector3 out( px2, py2, pz2 );
  // Angle from dot product definition
  return in.Angle( out );
}

std::vector<int> get_reco_gamma_index( std::vector<int> &true_pdg, std::vector<int> &true_parent_pdg ) {
  std::vector<int> daughter_idx;
  for( size_t j = 0; j < true_pdg.size(); j++ ) {
    if( true_pdg.at(j) == utils::pdg::kPdgGamma || true_parent_pdg.at(j) == utils::pdg::kPdgPi0 ) {
      daughter_idx.push_back( j );
    }
  }
  return daughter_idx;
}

// Return vector of indices corresponding to the reco showers for both reco and true_byHits vectors
std::vector<int> get_true_reco_gamma_index( std::vector<int> &true_pdg, std::vector<int> &true_id,
                                            std::vector<int> &true_parent_pdg ) {
  // True-byHits so some showers adn therefore gamma are missed
  // 1. Must be a gamma with parent pi0
  // 2. Must have 2 same IDs which indicates the gammas are from the same parent (pi0)
  // 3. If meets above save index
  std::vector<int> daughter_idx;
  for ( size_t j = 0; j < true_pdg.size(); j++ ) {
    if ( true_pdg.at( j ) == utils::pdg::kPdgGamma && true_parent_pdg.at( j ) == utils::pdg::kPdgPi0 ) {
      int id_cnt = utils::Count<int>( true_id, true_id.at(j) );
      if( id_cnt == 2 ) daughter_idx.push_back( j );
    }
  }
  return daughter_idx;
}

// Calculate the energy of the parent pi0 given the decay gamma's energy and opening angle
// Momentum calculation `p_gg` from Ref https://arxiv.org/pdf/1511.00941.pdf Eq.(3)
double E_gamma_gammma( double E_g1, double E_g2, double theta_gg ) {
  double p_gg = sqrt( pow(E_g1, 2) + pow(E_g2, 2) + 2 * E_g1 * E_g1 * cos(theta_gg) );
  return sqrt( pow( p_gg, 2 ) + pow( utils::pdg::pdg2mass( utils::pdg::kPdgPi0), 2) );
}

// Calculate the angle of the parent pi0 given the decay gamma's energy, angles wrt beam and opening angle
// Angle calculation `theta_gg` from Ref https://arxiv.org/pdf/1511.00941.pdf Eq.(4)
double theta_gamma_gammma( double E_g1, double E_g2, double theta_gg, double theta_g1, double theta_g2 ) {
  double p_gg = sqrt( pow(E_g1, 2) + pow(E_g2, 2) + 2 * E_g1 * E_g1 * cos(theta_gg) );
  return ( E_g1 * cos(theta_g1) + E_g2 * cos(theta_g2) ) / p_gg;
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
