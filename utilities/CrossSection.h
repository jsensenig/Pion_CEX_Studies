//
// Created by Jon Sensenig on 6/5/21.
//

#ifndef RECO_VERTEX_CROSSSECTION_H
#define RECO_VERTEX_CROSSSECTION_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TGraph.h"
#include <iostream>


class CrossSection {

public:

  CrossSection();
  virtual ~CrossSection();

  void FillRecoXsecHisto( double pip_ke, double energy, double angle );

  void FillTrueXsecHisto( double pip_ke, double energy, double angle );

  void ExtractXsec( int total_events, TString& out_file );

  void ExtractXsecEnergy( int nevts, TString& out_file, bool truth_xsec );

  void ExtractXsecAngle( int nevts, TString& out_file, bool truth_xsec );

private:

  TGraph* SinglePlotXsecEnergy( std::vector<double> &xsec, double angle, std::vector<double> &energy, double beam_energy,
                                              std::vector<double> &xerr, std::vector<double> &yerr, bool true_xsec );

  TGraph* SinglePlotXsecAngle( std::vector<double> &xsec, double energy, std::vector<double> &angle, double beam_energy,
                                             std::vector<double> &xerr, std::vector<double> &yerr, bool true_xsec );

  void WriteXsec( TString & out_file );

  double GetDDXsec( double energy, double angle, double beam );

  /// Histograms for X-section
  int energy_bins;
  int angle_bins;
  std::unique_ptr<TH1D> beam_energy_hist;
  std::unique_ptr<TH1D> energy_hist;
  std::unique_ptr<TH1D> angle_hist;
  std::unique_ptr<TH2D> energy_angle_hist;
  std::unique_ptr<TH3D> beam_energy_angle_hist;
  std::unique_ptr<TH3D> truth_beam_energy_angle_hist;

  /// Get the GEANT xsec
  std::unique_ptr<TFile> geant_xsec_file;
  std::map<double, TH2D*> geant_xsec_map;

};


#endif //RECO_VERTEX_CROSSSECTION_H
