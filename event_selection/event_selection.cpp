//
// Created by Jon Sensenig on 4/23/21.
//

#include "event_selection.h"
#include "../utilities/Utilities.hpp"
#include "TTree.h"
#include <iostream>


void run_event_selection( std::string in_file, Histograms &hists ) {

  TFile *proc_file = TFile::Open( in_file.c_str() );

  if( !proc_file -> IsOpen() ) {
    std::cout << "File " << in_file << " not open!" << std::endl;
    return;
  }

  TTree* tree = (TTree*)proc_file -> Get("pionana/beamana");

  double beam_inst_TOF;
  int true_beam_PDG;
  double true_beam_endX;
  double true_beam_endY;
  double true_beam_endZ;
  double reco_beam_endX;
  double reco_beam_endY;
  double reco_beam_endZ;
  double reco_beam_startX;
  double reco_beam_startY;
  double reco_beam_startZ;
  double reco_beam_len;
  int reco_beam_type;
  std::string *true_beam_endProcess = new std::string;

  tree->SetBranchAddress("beam_inst_TOF", &beam_inst_TOF);
  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);
  tree->SetBranchAddress("true_beam_endX", &true_beam_endX);
  tree->SetBranchAddress("true_beam_endY", &true_beam_endY);
  tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ);
  tree->SetBranchAddress("reco_beam_endX", &reco_beam_endX);
  tree->SetBranchAddress("reco_beam_endY", &reco_beam_endY);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("reco_beam_len", &reco_beam_len);
  tree->SetBranchAddress("reco_beam_type", &reco_beam_type);
  tree->SetBranchAddress("reco_beam_startX", &reco_beam_startX);
  tree->SetBranchAddress("reco_beam_startY", &reco_beam_startY);
  tree->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ);

  tree->SetBranchAddress("true_beam_endProcess", &true_beam_endProcess);

  size_t nevts = tree -> GetEntries();
  std::cout << "Processing " << nevts << " events." << std::endl;

  for ( size_t i = 0; i < nevts; i++ ) {
    tree->GetEntry( i );

    // 1. Beam type (track = 13 or shower = 11)
    // 2. ToF cut to select pi+ primary
    // 3. IsBeamLike (BI to TPC position match + angle)
    // 4. 2 daughter tracks (pi+ and >0 proton) no showers using track/shower CNN score
    // 5. 1 pi+ and >0 proton daughters using dE/dx vs Residual range fitting PID

    if( reco_beam_type != 11 ) {
      continue;
    }
    if( beam_inst_TOF < 96.0 || beam_inst_TOF > 97.0 ) {
      continue;
    }


  }

  // Clean up
  delete true_beam_endProcess;
  proc_file -> Close();

}

int main() {

  std::string input_file = "../../../pionana_Prod4_mc_1GeV_1_14_21.root";
  TString output_file = "out.root";
  std::string hists_config = "../hists.json";

  // Configure histograms
  Histograms hists;
  hists.ConfigureHistos( hists_config );

  std::cout << "Starting event selection!" << std::endl;
  run_event_selection( input_file, hists );

  std::cout << "Writing histograms to " << output_file << std::endl;
  // Write histograms ot file
  hists.WriteHistos( output_file );

  return 0;

}
