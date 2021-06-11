//
// Created by Jon Sensenig on 4/23/21.
//

#include "pi0_mc_study.h"
#include "../utilities/Utilities.hpp"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>


void run_pi0_mc_study( std::string in_file, Histograms &hists ) {

  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  TTree* tree = (TTree*)proc_file -> Get("pionana/beamana");

  int true_daughter_nPi0;
  int true_daughter_nNeutron;
  int true_daughter_nProton;
  int true_daughter_nPiMinus;
  int true_daughter_nPiPlus;
  int true_beam_PDG;
  double true_beam_endP;
  double true_beam_endX, true_beam_endY, true_beam_endZ;

  std::vector<int> *true_beam_daughter_PDG = new std::vector<int>;
  std::vector<int> *true_beam_daughter_ID = new std::vector<int>;
  std::vector<double> *true_beam_daughter_startP = new std::vector<double>;
  std::vector<int> *true_beam_Pi0_decay_PDG = new std::vector<int>;
  std::vector<int> *true_beam_Pi0_decay_parID = new std::vector<int>;
  std::vector<double> *true_beam_Pi0_decay_startP = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startPx = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startPy = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startPz = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startX = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startY = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_startZ = new std::vector<double>;
  std::vector<double> *true_beam_Pi0_decay_len = new std::vector<double>;

  std::string *true_beam_endProcess = new std::string;

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

  for ( size_t i = 0; i < nevts; i++ ) {
    tree->GetEntry( i );

    double pi0_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPi0 );
    double pi_mass = utils::pdg::pdg2mass( utils::pdg::kPdgPiP );

    double pi_ke = utils::CalculateKE( true_beam_endP * 1.e3, pi_mass );
    // All in-elastic pi events
    if( true_beam_endProcess->compare("pi+Inelastic") == 0 ) {
      hists.th1_hists["hPiInElKe"]->Fill( pi_ke );
    }

    bool true_abs = true_beam_endProcess->compare("pi+Inelastic") == 0 &&
                    true_daughter_nPi0 == 0 && true_daughter_nPiMinus == 0 &&
                    true_daughter_nPiPlus == 0;
    if( true_abs ) {
      int ngamma = 0;
      for( size_t i = 0; i < true_beam_daughter_PDG->size(); i++ ) {
        if( true_beam_daughter_PDG->at(i) == utils::pdg::kPdgGamma ) {
          hists.th1_hists["hAbsGammaP"] -> Fill( true_beam_daughter_startP->at(i)*1.e3 );
          ngamma++;
        }
      }
      hists.th1_hists["hAbsNgamma"] -> Fill( ngamma );
    }

    // Define true CEX
    bool true_cex = true_beam_endProcess->compare("pi+Inelastic") == 0 &&
                   true_daughter_nPi0 > 0 && true_daughter_nPiMinus == 0 &&
                   true_daughter_nPiPlus == 0;

    if( !true_cex ) continue;

    hists.th1_hists["hNpi0Cex"] -> Fill( true_daughter_nPi0 );
    hists.th1_hists["hFracCexKe"] -> Fill( pi_ke );
    hists.th2_hists["hNucleonDaughters"] -> Fill( true_daughter_nProton, true_daughter_nNeutron );

    std::map<int,double> pi0_energy_map;
    for( size_t i = 0; i < true_beam_daughter_PDG->size(); i++ ) {

      // Only look at pi0
      if( true_beam_daughter_PDG->at(i) != utils::pdg::kPdgPi0 ) continue;
      // Plot the daughter pi0
      double daughter_pi0_ke = utils::CalculateKE( true_beam_daughter_startP->at( i ) * 1.e3, pi0_mass );
      double daughter_pi0_energy = utils::CalculateE( true_beam_daughter_startP->at( i ) * 1.e3, pi0_mass );
      hists.th2_hists["hPiPPi0P"] -> Fill( true_beam_endP*1.e3, true_beam_daughter_startP->at( i )*1.e3 );
      if( true_daughter_nPi0 == 1 ) {
        hists.th2_hists["hPiP1Pi0P"] -> Fill( true_beam_endP*1.e3, true_beam_daughter_startP->at( i )*1.e3 );
      }
      hists.th1_hists["hPi0Ke"]->Fill( daughter_pi0_ke );
      hists.th1_hists["hPi0E"]->Fill( daughter_pi0_energy );
      pi0_energy_map[true_beam_daughter_ID->at(i)] = daughter_pi0_energy;
    }

    for( size_t i = 0; i < true_beam_Pi0_decay_parID->size()/2; i++ ) {
      size_t idx = i * 2; // Gamma from the same pi0 are adjacent in index, i.e., i and i+1
      if( true_beam_Pi0_decay_parID -> at( idx ) != true_beam_Pi0_decay_parID -> at( idx+1 ) ) {
        std::cout << "Something's wrong, decay gamma IDs do not match!" << std::endl;
        std::cout << "PDG 1/2 " << true_beam_Pi0_decay_PDG->at(idx) << "/" << true_beam_Pi0_decay_PDG->at(idx+1) << std::endl;
      }
      // Decay gammas momentum and opening angle
      double gamma1 = true_beam_Pi0_decay_startP -> at( idx ) * 1.e3;
      double gamma2 = true_beam_Pi0_decay_startP -> at( idx+1 ) * 1.e3;
      double angle = open_angle( true_beam_Pi0_decay_startPx->at(idx), true_beam_Pi0_decay_startPy->at(idx),
                                 true_beam_Pi0_decay_startPz->at(idx), true_beam_Pi0_decay_startPx->at(idx+1),
                                 true_beam_Pi0_decay_startPy->at(idx+1), true_beam_Pi0_decay_startPz->at(idx+1) );

      double beam_end_pos = utils::Distance(true_beam_endX, true_beam_endY, true_beam_endZ);
      double gamma_pos1 = utils::Distance(true_beam_Pi0_decay_startX->at(idx), true_beam_Pi0_decay_startY->at(idx), true_beam_Pi0_decay_startZ->at(idx));
      double gamma_pos2 = utils::Distance(true_beam_Pi0_decay_startX->at(idx+1), true_beam_Pi0_decay_startY->at(idx+1), true_beam_Pi0_decay_startZ->at(idx+1));

      hists.th1_hists["hGammaOpenAngle"] -> Fill( angle );
      hists.th2_hists["hPi0EGammaOpenAngle"] -> Fill( pi0_energy_map.at(true_beam_Pi0_decay_parID->at( idx )), angle );
      hists.th1_hists["hLeadGammaP"] -> Fill( std::max( gamma1, gamma2 ) );
      hists.th1_hists["hSubLeadGammaP"] -> Fill( std::min( gamma1, gamma2 ) );
      hists.th2_hists["hGammaP"] -> Fill( std::max( gamma1, gamma2 ), std::min( gamma1, gamma2 ) );
      hists.th2_hists["hGammaR"] -> Fill( gamma_pos1-beam_end_pos, gamma_pos2-beam_end_pos );
      hists.th2_hists["hGammaLen"] -> Fill( true_beam_Pi0_decay_len->at(idx), true_beam_Pi0_decay_len->at(idx+1) );

    }

  }

  // Fraction of all in-elastic scatters are CEX
  hists.th1_hists["hFracCexKe"] -> Divide( hists.th1_hists["hPiInElKe"].get() );

  // Clean up
  delete true_beam_endProcess;
  delete true_beam_daughter_startP;
  delete true_beam_daughter_PDG;
  delete true_beam_daughter_ID;
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
  proc_file -> Close();

}

double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 ) {

  TVector3 in( px1, py1, pz1 );
  TVector3 out( px2, py2, pz2 );
  // Angle from dot product definition
  return in.Angle( out );
}

int main() {

  std::string input_file = "../../../pionana_Prod4_mc_1GeV_1_14_21.root";
  TString output_file = "out.root";
  std::string hists_config = "../hists.json";

  // Configure histograms
  Histograms hists;
  hists.ConfigureHistos( hists_config );

  std::cout << "Starting pi0 study!" << std::endl;
  run_pi0_mc_study( input_file, hists );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}
