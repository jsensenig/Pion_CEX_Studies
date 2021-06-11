//
// Created by Jon Sensenig on 5/24/21.
//

#ifndef RECO_VERTEX_SELECTION_H
#define RECO_VERTEX_SELECTION_H

#include "../utilities/Utilities.hpp"

class Selection {

public:

  Selection() { ; }
  virtual ~Selection() { ; };


  //................................................................
  // Selection functions

  // Upper & lower cuts (not-inclusive) on beam-line ToF cuts to do PID
  inline bool beam_tof_cut( double tof ) const {
    return ( tof > lower_beam_tof_cut && tof < upper_beam_tof_cut );
  }

  // Upper & lower cuts (not-inclusive) on beam particles to help remove muons
  inline bool beam_endz_cut( double endz ) const {
    return endz > lower_beam_endz && endz < upper_beam_endz;
  }

  // Upper & lower cuts (not-inclusive) average particle PIDA
  // pdg parameter is the particle type hypothesis
  inline bool pida_cut( double pida, int pdg ) const {

    switch ( pdg ) {
      case utils::pdg::kPdgPiP:
        return pida > lower_pion_pida && pida < upper_pion_pida;
      case utils::pdg::kPdgProton:
        return pida > lower_proton_pida && pida < upper_proton_pida;
      case utils::pdg::kPdgMuon:
        return pida > lower_muon_pida && pida < upper_muon_pida;
      case utils::pdg::kPdgKP:
        return pida > lower_kaon_pida && pida < upper_kaon_pida;
      default:
        std::cout << "Unknown PDG " << pdg << std::endl;
    }

  }

  // Upper and lower cuts for the CNN track/shower score
  inline bool cnn_track_shower_score( double score ) const {
    return score > lower_cnn_trackshower_score && score < upper_cnn_trackshower_score;
  }

  inline bool deltaR_to_vertex( double dR ) const {
    return dR > lower_deltaR;
  }

  inline bool daughter_nhits( double nhits ) const {
    return nhits > lower_daughter_nhits;
  }

  inline bool daughter_michel_score( double score ) const {
    return score < upper_michel_score;
  }


  //................................................................
  // Load the parameters for event selection from the json configuration file
  bool set_selection_parameters( const std::string& param_file ) {

    std::vector<std::string> paramset_name;
    json paramset = utils::LoadConfig( param_file );

    if ( paramset == nullptr ) return false; // Failed to load config

    try {

      upper_beam_tof_cut = paramset.at("beam_tof").at("upper").get<double>();
      lower_beam_tof_cut = paramset.at("beam_tof").at("lower").get<double>();
      paramset_name.emplace_back("beam_tof");

      upper_beam_endz = paramset.at("beam_endz").at("upper").get<double>();
      lower_beam_endz = paramset.at("beam_endz").at("lower").get<double>();
      paramset_name.emplace_back("beam_endz");

      upper_pion_pida = paramset.at("pida").at("pion_upper").get<double>();
      lower_pion_pida = paramset.at("pida").at("pion_lower").get<double>();

      upper_proton_pida = paramset.at("pida").at("proton_upper").get<double>();
      lower_proton_pida = paramset.at("pida").at("proton_lower").get<double>();

      upper_muon_pida = paramset.at("pida").at("muon_upper").get<double>();
      lower_muon_pida = paramset.at("pida").at("muon_lower").get<double>();

      upper_kaon_pida = paramset.at("pida").at("kaon_upper").get<double>();
      lower_kaon_pida = paramset.at("pida").at("kaon_lower").get<double>();
      paramset_name.emplace_back("pida");

      upper_cnn_trackshower_score = paramset.at("cnn_trackshower").at("upper").get<double>();
      lower_cnn_trackshower_score = paramset.at("cnn_trackshower").at("lower").get<double>();
      paramset_name.emplace_back("cnn_trackshower");

      lower_deltaR = paramset.at("deltaR").at("lower").get<double>();
      paramset_name.emplace_back("deltaR");

      lower_daughter_nhits = paramset.at("daughter_nhits").at("lower").get<int>();
      paramset_name.emplace_back("daughter_nhits");

      upper_michel_score = paramset.at("michel_score").at("upper").get<int>();
      paramset_name.emplace_back("michel_score");

    } catch ( json::exception & ex ) {
      std::cout << " Parse error in file " << param_file << " " << ex.what() << std::endl;
      return false;
    }

    std::cout << "************** Loading Selection Parameters *****************" << std::endl;
    for( auto &str : paramset_name ) std::cout << " Loaded: " << paramset.at(str).dump() << std::endl;
    std::cout << "*************************************************************" << std::endl;

    return true;
  }

private:

  //................................................................
  // Selection parameters

  // Upper & lower cuts on beam-line ToF cuts to do PID
  double upper_beam_tof_cut;
  double lower_beam_tof_cut;

  // Upper & lower beam EndZ cut to help remove muons
  double upper_beam_endz;
  double lower_beam_endz;

  // Upper & lower average PIDA cut for PID pions
  double upper_pion_pida;
  double lower_pion_pida;

  // Upper & lower average PIDA cut for PID protons
  double upper_proton_pida;
  double lower_proton_pida;

  // Upper & lower average PIDA cut for PID muons
  double upper_muon_pida;
  double lower_muon_pida;

  // Upper & lower average PIDA cut for PID kaons
  double upper_kaon_pida;
  double lower_kaon_pida;

  // CNN track/shower score cut
  double upper_cnn_trackshower_score;
  double lower_cnn_trackshower_score;

  // Distance between daughter and primary vertex
  double lower_deltaR;

  // Daughter nhits
  int lower_daughter_nhits;

  // Michel score from Ajib's CNN classifier to help remove beam muons
  double upper_michel_score;

};
#endif //RECO_VERTEX_SELECTION_H
