//
// Created by Jon Sensenig on 4/23/21.
//

#include "pi0_reco_xsec.h"
//#include "../utilities/CrossSection.h"
#include "EventSelection.h"
#include "TTree.h"
#include "TVector3.h"
#include "TGraph.h"
#include <iostream>


void run_pi0_mc_xsec( const std::string& in_file, const std::string& in_tree, Histograms &hists, bool truth_xsec ) {

  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  TTree* tree = (TTree*)proc_file -> Get( in_tree.c_str() );

  /// Truth
  tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
  tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
  tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);
  tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
  tree->SetBranchAddress("true_beam_startP", &true_beam_startP);
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
  tree->SetBranchAddress("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
  tree->SetBranchAddress("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
  tree->SetBranchAddress("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
  tree->SetBranchAddress("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);
  tree->SetBranchAddress("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);

  /// Reco
  tree->SetBranchAddress("beam_inst_P", &beam_inst_P);
  tree->SetBranchAddress("reco_beam_calibrated_dEdX_SCE", &reco_beam_calibrated_dEdX_SCE);
  tree->SetBranchAddress("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
  tree->SetBranchAddress("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
  tree->SetBranchAddress("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);

  tree->SetBranchAddress("reco_daughter_allShower_energy", &reco_daughter_allShower_energy);
  tree->SetBranchAddress("reco_daughter_allShower_len", &reco_daughter_allShower_len);
  tree->SetBranchAddress("reco_daughter_allShower_dirX", &reco_daughter_allShower_dirX);
  tree->SetBranchAddress("reco_daughter_allShower_dirY", &reco_daughter_allShower_dirY);
  tree->SetBranchAddress("reco_daughter_allShower_dirZ", &reco_daughter_allShower_dirZ);


  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  pi0_proxy pi0;
  pi0_proxy true_pi0;
  EventSelection evt_sel;
  evt_sel.InitTtree( tree );
  true_cex_count = 0;
  true_inc_piplus_count = 0;

  for ( size_t evt = 0; evt < nevts; evt++ ) {
    tree->GetEntry( evt );

    pi0.reset();

    // Define true CEX
    bool true_cex = true_beam_PDG == utils::pdg::kPdgPiP && *true_beam_endProcess == "pi+Inelastic" &&
                   true_daughter_nPi0 == 1 && true_daughter_nPiMinus == 0 &&
                   true_daughter_nPiPlus == 0 && ( true_daughter_nProton > 0 || true_daughter_nNeutron > 0 );

    // We only want the beam pions contributing to the in-elastic processes
    true_inc_piplus_count += true_beam_PDG == utils::pdg::kPdgPiP && *true_beam_endProcess == "pi+Inelastic";

    std::map<int, double> pi0_energy_map = daughter_pi0_energy( hists );
    for( size_t i = 0; i < true_beam_Pi0_decay_parID->size()/2; i++ ) {
      size_t idx = i * 2;
      true_pi0.reset();

      // Decay gammas from the same pi0 are adjacent in index, i.e., index and index+1
      // They should also share the same parent ID (obviously)
      if( true_beam_Pi0_decay_parID -> at( idx ) != true_beam_Pi0_decay_parID -> at( idx+1 ) ) {
        std::cout << "Something's wrong, decay gamma IDs do not match!" << std::endl;
        std::cout << "PDG 1/2 " << true_beam_Pi0_decay_PDG->at(idx) << "/" << true_beam_Pi0_decay_PDG->at(idx+1) << std::endl;
      }

      /// These are the 5 variables needed for the cross section calculation ///
      // .........................................................................

      /// Truth values
      // Interaction KE of the beam pi+
      true_pi0.pip_energy = utils::CalculateKE( true_beam_endP * 1.e3, utils::pdg::pdg2mass( utils::pdg::kPdgPiP ));
      // Decay gammas momentum and opening angle
      true_pi0.energy.first = true_beam_Pi0_decay_startP->at( idx ) * 1.e3;
      true_pi0.energy.second = true_beam_Pi0_decay_startP->at( idx + 1 ) * 1.e3;
      true_pi0.open_angle = open_angle( true_beam_Pi0_decay_startPx->at( idx ),
                                        true_beam_Pi0_decay_startPy->at( idx ),
                                        true_beam_Pi0_decay_startPz->at( idx ),
                                        true_beam_Pi0_decay_startPx->at( idx + 1 ),
                                        true_beam_Pi0_decay_startPy->at( idx + 1 ),
                                        true_beam_Pi0_decay_startPz->at( idx + 1 ));

      // Gamma angle wrt incoming pi+
      true_pi0.angle.first = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                         true_beam_Pi0_decay_startPx->at( idx ),
                                         true_beam_Pi0_decay_startPy->at( idx ),
                                         true_beam_Pi0_decay_startPz->at( idx ));
      true_pi0.angle.second = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                          true_beam_Pi0_decay_startPx->at( idx + 1 ),
                                          true_beam_Pi0_decay_startPy->at( idx + 1 ),
                                          true_beam_Pi0_decay_startPz->at( idx + 1 ));

      // .........................................................................

      // Plots for events with multiple pi0
      plot_all_pi0(pi0, pi0_energy_map.at(true_beam_Pi0_decay_parID->at( idx )), hists, idx );

      if( true_daughter_nPi0 > 1 ) continue; // Only look at single pi0 events for now

      // Plots for events with only a single pi0
      plot_single_pi0( pi0, hists, idx );

    } //true pi0 loop

    // If event not selected, ie CEX candidate, then skip
    if( !evt_sel.SelectEvent( tree, evt, hists, true_cex, true_daughter_nPi0 ) && !truth_xsec ) continue;

    if( reco_daughter_allShower_energy->size() > 1 ) {
      /// Reco values
      // Interaction KE of the beam pi+
      pi0.pip_energy = pip_interaction_ke();
      // Decay gammas momentum and opening angle
      pi0.energy.first = reco_daughter_allShower_energy->at(0);
      pi0.energy.second = reco_daughter_allShower_energy->at(1);
      pi0.open_angle = open_angle( reco_daughter_allShower_dirX->at(0), reco_daughter_allShower_dirY->at(0),
                                   reco_daughter_allShower_dirZ->at(0), reco_daughter_allShower_dirX->at(1),
                                   reco_daughter_allShower_dirY->at(1), reco_daughter_allShower_dirZ->at(1) );

      // Gamma angle wrt incoming pi+
      pi0.angle.first = open_angle( reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ,
                                    reco_daughter_allShower_dirX->at(0), reco_daughter_allShower_dirY->at(0),
                                    reco_daughter_allShower_dirZ->at(0));
      pi0.angle.second = open_angle( reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ,
                                     reco_daughter_allShower_dirX->at(1), reco_daughter_allShower_dirY->at(1),
                                     reco_daughter_allShower_dirZ->at(1) );
    } else if( reco_daughter_allShower_energy->size() > 2 ) {
      std::cout << "AllShower size big " << reco_daughter_allShower_energy->size() << std::endl;
    }

    if( truth_xsec && true_cex ) cs.FillTrueXsecHisto( true_pi0.pip_energy, pi0_kinetic_energy( true_pi0 ), pi0_cos_angle( true_pi0 ) );
    else if( !truth_xsec ) cs.FillRecoXsecHisto( pi0.pip_energy, pi0_kinetic_energy( pi0 ), pi0_cos_angle( pi0 ) );

    true_cex_count += true_cex;

  }//evt loop

  std::cout << "Incident pi+: " << true_inc_piplus_count << " True CEX: " << true_cex_count << std::endl;

  TString exsec_file = "energy_xsec.root";
  TString axsec_file = "angle_xsec.root";
  //cs.ExtractXsec( nevts, exsec_file );
  cs.ExtractXsecEnergy( true_inc_piplus_count, exsec_file, truth_xsec ); // sigma as function of pi0 KE
  cs.ExtractXsecAngle( true_inc_piplus_count, axsec_file, truth_xsec );  // sigma as function of pi0 theta

  proc_file -> Close();
  clean_pointers();

}

// ............................................................................
double pip_interaction_ke() {

  double pip_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPiP );

  // Get beam particle energy loss by integrating dE/dx
  double pip_eloss = 0;
  for( auto& eloss : *reco_beam_calibrated_dEdX_SCE ) {
    if ( eloss > 1000. ) continue;
    pip_eloss += eloss;
  }

  // FIXME Prod2 has no beam instrumentation info so use truth for now!
  //double interaction_e = ( utils::CalculateE( beam_inst_P, pip_mass ) - pip_eloss ) - pip_mass;
  double interaction_ke = ( utils::CalculateE( true_beam_startP * 1e3, pip_mass ) - pip_eloss ) - pip_mass;

  if( interaction_ke < 0 ) {
    std::cout << "Negative (" << interaction_ke << ") Pion Energy! Returning KE = 0" << std::endl;
    return 0.;
  }

  return interaction_ke;
}

// ............................................................................
double pi0_cos_angle( pi0_proxy& pi0 ) {

  return ( pi0.energy.first * cos(pi0.angle.first) + pi0.energy.second * cos(pi0.angle.second) ) / pi0_mom( pi0 );

}

// ............................................................................
double pi0_energy( pi0_proxy& pi0 ) {

  double pi0_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );
  double pi0_p = pi0_mom( pi0 );
  return sqrt( pi0_mass*pi0_mass + pi0_p*pi0_p );

}

// ............................................................................
double pi0_kinetic_energy( pi0_proxy& pi0 ) {

  return pi0_energy(pi0) - utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );

}

// ............................................................................
double pi0_mom( pi0_proxy& pi0 ) {

  return sqrt( pow(pi0.energy.first,2) + pow(pi0.energy.second,2) +
                  2 * pi0.energy.first * pi0.energy.second * cos(pi0.open_angle) );

}

// ............................................................................
void plot_single_pi0( pi0_proxy& pi0, Histograms& hists, size_t idx ) {

  double pi0_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );
  int pi0_idx = utils::FindIndex<int>( *true_beam_daughter_PDG, utils::pdg::kPdgPi0 );
  double pi0_true_energy = utils::CalculateE( true_beam_daughter_startP->at( pi0_idx ) * 1.e3, pi0_mass );

  // Momentum from polynomial fit Ref https://arxiv.org/pdf/1511.00941.pdf
  // Convert angle to degrees first
  double angle_deg = pi0.open_angle * TMath::RadToDeg();
  double p_poly = 2202.3 - 94.9*angle_deg + 2.1*pow(angle_deg, 2) - 0.025*pow(angle_deg, 3) + 0.00017*pow(angle_deg, 4)
                  - 6.0e-7*pow(angle_deg, 5) + 8.5e-10*pow(angle_deg, 6);

  hists.th1_hists["hPolyPi0PError"] -> Fill( p_poly / (true_beam_daughter_startP->at( pi0_idx )*1.e3) );
  hists.th1_hists["hRootPi0PError"] -> Fill( pi0_mom( pi0 ) / (true_beam_daughter_startP->at( pi0_idx )*1.e3) );
  hists.th2_hists["hPi0TrueCalcEnergy"] -> Fill( pi0_energy( pi0 ), pi0_true_energy );

  // Gamma angle wrt pi0
  double gamma1_open_angle = open_angle( true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
                                         true_beam_daughter_startPz->at(pi0_idx), true_beam_Pi0_decay_startPx->at(idx),
                                         true_beam_Pi0_decay_startPy->at(idx), true_beam_Pi0_decay_startPz->at(idx) );
  double gamma2_open_angle = open_angle( true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
                                         true_beam_daughter_startPz->at(pi0_idx), true_beam_Pi0_decay_startPx->at(idx+1),
                                         true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );
  hists.th1_hists["hGammaDiff"] -> Fill( gamma1_open_angle - gamma2_open_angle );
  hists.th2_hists["hGammaPi0Angle"] -> Fill( gamma1_open_angle*TMath::RadToDeg(), gamma2_open_angle*TMath::RadToDeg() );

  // pi0 angle wrt incoming pi+
  double true_pi0_angle = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                 true_beam_daughter_startPx->at(pi0_idx), true_beam_daughter_startPy->at(pi0_idx),
                                 true_beam_daughter_startPz->at(pi0_idx) );
  hists.th2_hists["hPi0TrueCalc"] -> Fill( pi0_cos_angle( pi0 ), cos( true_pi0_angle ) );
  hists.th1_hists["hPi0CalcAngleDiff"] -> Fill( pi0_cos_angle( pi0 ) - cos( true_pi0_angle ) );
}

// ............................................................................
void plot_all_pi0( pi0_proxy& pi0, double pi0_energy, Histograms& hists, size_t idx ) {

  double beam_end_pos = utils::Distance(true_beam_endX, true_beam_endY, true_beam_endZ);
  double gamma_pos1 = utils::Distance(true_beam_Pi0_decay_startX->at(idx), true_beam_Pi0_decay_startY->at(idx),
                                      true_beam_Pi0_decay_startZ->at(idx));
  double gamma_pos2 = utils::Distance(true_beam_Pi0_decay_startX->at(idx+1), true_beam_Pi0_decay_startY->at(idx+1),
                                      true_beam_Pi0_decay_startZ->at(idx+1));

  hists.th1_hists["hGammaOpenAngle"] -> Fill( pi0.open_angle );
  hists.th1_hists["hLeadGammaP"] -> Fill( std::max( pi0.energy.first, pi0.energy.second ) );
  hists.th1_hists["hSubLeadGammaP"] -> Fill( std::min( pi0.energy.first, pi0.energy.second ) );
  hists.th2_hists["hGammaR"] -> Fill( gamma_pos1-beam_end_pos, gamma_pos2-beam_end_pos );
  hists.th2_hists["hGammaLen"] -> Fill( true_beam_Pi0_decay_len->at(idx), true_beam_Pi0_decay_len->at(idx+1) );
  hists.th2_hists["hPi0EGammaOpenAngle"] -> Fill( TMath::RadToDeg()*pi0.open_angle, pi0_energy );
  hists.th2_hists["hGammaP"] -> Fill( std::max( pi0.energy.first, pi0.energy.second ),
                                      std::min( pi0.energy.first, pi0.energy.second ) );

}

// ............................................................................
std::map<int, double> daughter_pi0_energy( Histograms& hists ) {

  std::map<int, double> daughter_pi0_energy_map;
  double pi0_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );

  for( size_t i = 0; i < true_beam_daughter_PDG->size(); i++ ) {
    // Only look at pi0
    if( true_beam_daughter_PDG->at(i) != utils::pdg::kPdgPi0 ) continue;

    // Plot the daughter pi0
    double daughter_pi0_energy = utils::CalculateE( true_beam_daughter_startP->at( i ) * 1.e3, pi0_mass );
    daughter_pi0_energy_map[true_beam_daughter_ID->at(i)] = daughter_pi0_energy;

    hists.th1_hists["hPi0E"]->Fill( daughter_pi0_energy );
    hists.th2_hists["hPiPPi0P"] -> Fill( true_beam_endP*1.e3, true_beam_daughter_startP->at( i )*1.e3 );
    if( true_daughter_nPi0 == 1 ) {
      hists.th2_hists["hPiP1Pi0P"] -> Fill( true_beam_endP*1.e3, true_beam_daughter_startP->at( i )*1.e3 );
    }
  }

  return daughter_pi0_energy_map;

}

// ............................................................................
double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 ) {

  TVector3 in( px1, py1, pz1 );
  TVector3 out( px2, py2, pz2 );
  // Angle from dot product definition
  return in.Angle( out );
}

// ............................................................................
void clean_pointers() {

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_startP;
  delete true_beam_daughter_PDG;
  delete true_beam_daughter_ID;
  delete true_beam_daughter_startPx;
  delete true_beam_daughter_startPy;
  delete true_beam_daughter_startPz;
  delete true_beam_Pi0_decay_PDG;
  delete true_beam_Pi0_decay_startP;
  delete true_beam_Pi0_decay_parID;
  delete true_beam_Pi0_decay_startPx;
  delete true_beam_Pi0_decay_startPy;
  delete true_beam_Pi0_decay_startPz;
  delete true_beam_Pi0_decay_startX;
  delete true_beam_Pi0_decay_startY;
  delete true_beam_Pi0_decay_startZ;
  delete true_beam_Pi0_decay_len;
  delete reco_beam_calibrated_dEdX_SCE;
  delete reco_daughter_allShower_energy;
  delete reco_daughter_allShower_len;
  delete reco_daughter_allShower_dirX;
  delete reco_daughter_allShower_dirY;
  delete reco_daughter_allShower_dirZ;

}

int main(int argc, char * argv[]){

  std::cout << "starting" << std::endl;

  if( !parseArgs( argc, argv ) )
    return 0;

  std::string input_file;
  std::string input_tree;

  /// 1 GeV 228k events
  if( erange_arg == 1 ) {
    input_file = "../../../pionana_Prod4_mc_1GeV_1_14_21.root";
    input_tree = "pionana/beamana";
  } else if( erange_arg == 2 ) {
    /// 2GeV 2.6k events
    input_file = "../../../pduneana_2gev_n2590.root";
    input_tree = "pduneana/beamana;2";
  }

  TString output_file = "out.root";
  std::string hists_config = "../hists.json";

  // Configure histograms
  hists.ConfigureHistos( hists_config );

  std::cout << "Starting pi0 xsec study! For " << erange_arg << "GeV Truth=" << truth_xsec_arg << std::endl;

  run_pi0_mc_xsec( input_file, input_tree, hists, truth_xsec_arg );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}

bool parseArgs(int argc, char ** argv) {

  bool found_erange = false, found_sim_select = false;

  for( int i = 1; i < argc; ++i ) {

    if( ( strcmp( argv[i], "--help" )  == 0 ) || ( strcmp( argv[i], "-h" ) == 0 ) ){
           std::cout << "Usage: ./pi0_reco_xsec -e E_pi+ -s MC=1/Reco=0 [options]" << std::endl;
           std::cout << std::endl;
           std::cout << "Options: " << std::endl;
           std::cout << "\t-o <output_file_override>.root" << std::endl;

           return false;
    }
    else if( strcmp( argv[i], "-e" ) == 0 ) {
      erange_arg = atoi( argv[i+1] );
      found_erange = true;
    }
    else if( strcmp( argv[i], "-s" ) == 0 ) {
      truth_xsec_arg = atoi( argv[i + 1] );
      found_sim_select = true;
    }
    else if( strcmp( argv[i], "-o" ) == 0 ){
       output_file_override = argv[i+1];
    }

  }

  if(!( found_erange && found_sim_select )) std::cout << "Missing -e or -s arguement!" << std::endl;
  return found_erange && found_sim_select;
}

