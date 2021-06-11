//
// Created by Jon Sensenig on 4/23/21.
//

#include "pi0_mc_xsec.h"
//#include "../utilities/CrossSection.h"
#include "TTree.h"
#include "TVector3.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>


void run_pi0_mc_xsec( std::string in_file, Histograms &hists ) {

  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  //TTree* tree = (TTree*)proc_file -> Get("pionana/beamana");  // 1GeV
  TTree* tree = (TTree*)proc_file -> Get("pduneana/beamana;2"); // 2GeV

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
  tree->SetBranchAddress("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
  tree->SetBranchAddress("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
  tree->SetBranchAddress("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
  tree->SetBranchAddress("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
  tree->SetBranchAddress("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);
  tree->SetBranchAddress("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);

  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  pi0_proxy pi0;
//  CrossSection cs(20, 10);

  for ( size_t evt = 0; evt < nevts; evt++ ) {
    tree->GetEntry( evt );

    //if( evt > 10000 ) break;

    // Define true CEX
    bool true_cex = *true_beam_endProcess == "pi+Inelastic" &&
                   true_daughter_nPi0 > 0 && true_daughter_nPiMinus == 0 &&
                   true_daughter_nPiPlus == 0 && ( true_daughter_nProton > 0 || true_daughter_nNeutron > 0 );

    if( !true_cex ) continue;

    std::map<int, double> pi0_energy_map = daughter_pi0_energy( hists );

    for( size_t i = 0; i < true_beam_Pi0_decay_parID->size()/2; i++ ) {
      size_t idx = i * 2;
      pi0.reset();

      // Decay gammas from the same pi0 are adjacent in index, i.e., index and index+1
      // They should also share the same parent ID (obviously)
      if( true_beam_Pi0_decay_parID -> at( idx ) != true_beam_Pi0_decay_parID -> at( idx+1 ) ) {
        std::cout << "Something's wrong, decay gamma IDs do not match!" << std::endl;
        std::cout << "PDG 1/2 " << true_beam_Pi0_decay_PDG->at(idx) << "/" << true_beam_Pi0_decay_PDG->at(idx+1) << std::endl;
      }

      /// These are the 5 variables needed for the cross section calculation ///
      // .........................................................................

      // Decay gammas momentum and opening angle
      pi0.energy.first  = true_beam_Pi0_decay_startP -> at( idx ) * 1.e3;
      pi0.energy.second = true_beam_Pi0_decay_startP -> at( idx+1 ) * 1.e3;
      pi0.open_angle = open_angle( true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
                                 true_beam_Pi0_decay_startPz->at(idx), true_beam_Pi0_decay_startPx->at(idx+1),
                                 true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );

      // Gamma angle wrt incoming pi+
      pi0.angle.first = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                    true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
                                    true_beam_Pi0_decay_startPz->at(idx) );
      pi0.angle.second = open_angle( true_beam_endPx, true_beam_endPy, true_beam_endPz,
                                     true_beam_Pi0_decay_startPx->at(idx+1), true_beam_Pi0_decay_startPy->at(idx+1),
                                     true_beam_Pi0_decay_startPz->at(idx+1) );

      // .........................................................................

      // Plots for events with multiple pi0
      plot_all_pi0(pi0, pi0_energy_map.at(true_beam_Pi0_decay_parID->at( idx )), hists, idx );

      if( true_daughter_nPi0 > 1 ) continue; // Only look at single pi0 events for now

      // Plots for events with only a single pi0
      plot_single_pi0( pi0, hists, idx );

      cs.FillXsecHisto( pi0_kinetic_energy(pi0), pi0_cos_angle( pi0 ) );

    }
  }

  TString xsec_file = "xsec2.root";
  cs.ExtractXsec( nevts, xsec_file );

  proc_file -> Close();
  clean_pointers();
  //delete proc_file;

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

}

int main() {

  /// 1 GeV 228k events
//  std::string input_file = "../../../pionana_Prod4_mc_1GeV_1_14_21.root";
  /// 2GeV 2.6k events
  std::string input_file = "../../../pduneana_2gev_n2590.root";
  TString output_file = "out.root";
  std::string hists_config = "../hists.json";

  // Configure histograms
  hists.ConfigureHistos( hists_config );

  std::cout << "Starting pi0 xsec study!" << std::endl;
  run_pi0_mc_xsec( input_file, hists );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}
