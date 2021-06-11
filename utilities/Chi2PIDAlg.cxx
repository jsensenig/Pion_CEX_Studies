//
// Created by Jon Sensenig on 2/15/21.
// from Larsoft PID algorithms 
//
// LArSoft: Toolkit for Simulation, Reconstruction and Analysis of Liquid Argon TPC Neutrino Detectors
// E.L. Snider (Fermilab), G. Petrillo (Fermilab)
// DOI: 10.1088/1742-6596/898/4/042057
// Published in: J.Phys.Conf.Ser. 898 (2017) 4, 042057
//

#include "Chi2PIDAlg.hpp"
#include "TFile.h"


pid::Chi2PIDAlg::Chi2PIDAlg() {

  Config();
  std::cout << "Attempting to open template file: " << fTemplateFile << std::endl;
  TFile *file = TFile::Open(fTemplateFile.c_str());
  dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

}

pid::Chi2PIDAlg::~Chi2PIDAlg()
{ }

//std::pair<double, int> pid::Chi2PIDAlg::Chi2PID( const std::vector<double> & track_dedx, const std::vector<double> & range, int pdg ) {
double pid::Chi2PIDAlg::Chi2PID( const std::vector<double> & track_dedx, const std::vector<double> & range, int pdg ) {

  double pid_chi2 = 0.;
  int npt = 0;

  // Get appropriate profile
  TProfile * profile = SelectProfile( pdg );

  if( profile == 0 ) {
    std::cout << "Cannot find template for particle " << pdg << std::endl;
    //return std::make_pair(9999., -1);
    return 9999.;
  }

  if( track_dedx.size() < 1 || range.size() < 1 )
    //return std::make_pair(9999., -1);
    return 9999.;

  //Ignore first and last point
  for( size_t i = 1; i < track_dedx.size()-1; ++i ){

    //Skip large pulse heights
    if( track_dedx[i] > 1000. )
      continue;

    int bin = profile->FindBin( range[i] );
    if( bin >= 1 && bin <= profile->GetNbinsX() ){

      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;
      }

      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;
      }

      double dedx_res = 0.04231 + 0.0001783 * track_dedx[i] * track_dedx[i];
      dedx_res *= track_dedx[i];

      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )
      pid_chi2 += ( pow( (track_dedx[i] - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }

  // if( npt == 0 ) return std::make_pair(9999., -1);
  if( npt == 0 ) return 9999.;

  // return std::make_pair(pid_chi2, npt);
  return pid_chi2 / npt;
}

TProfile * pid::Chi2PIDAlg::SelectProfile( int pdg ) {

  switch( pdg ) {

    case utils::pdg::kPdgPiP :
      return dedx_range_pi;

    case utils::pdg::kPdgProton :
      return dedx_range_pro;

    case utils::pdg::kPdgMuon :
      return dedx_range_mu;

    case utils::pdg::kPdgKP :
      return dedx_range_ka;

    default :
      return nullptr;

  }

}

double pid::Chi2PIDAlg::PidaFit( const std::vector<double> & track_dedx, const std::vector<double> & res_range ) {

  if( track_dedx.size() == 0 || res_range.size() == 0 ) return -1;

  // Data has functional form f = A * R^(-0.42)  => A = f * R^(0.42)
  // f = dE/dx
  // R = residual range (exponent extracted in Argonuet)
  // A = fit parameter = PIDA
  std::vector<double> pida_values( res_range.size() );
  double fMinResRange = 0;
  double fMaxResRange = 30;
  double fMaxPIDAValue = 50;
  double fExponentConstant = 0.42;

  for( size_t i_r = 0; i_r < res_range.size(); i_r++ ){
    if( res_range[i_r] > fMaxResRange || res_range[i_r] < fMinResRange ) continue;
    float val = track_dedx[i_r] * std::pow( res_range[i_r], fExponentConstant );
    if( val < fMaxPIDAValue ) pida_values.push_back(val);
  }

  if( pida_values.size() == 0 ) return -1;

  double pida_sum = 0;
  for( auto &pida : pida_values ) pida_sum += pida;

  return pida_sum / pida_values.size(); // return avg PIDA value

}

void pid::Chi2PIDAlg::Config() {

  json conf = utils::LoadConfig( _config_file );
  if ( conf == 0x0 ) return;

  if ( conf.contains("PIDAAvg") ) {
    MinResRange = conf.at("PIDAAvg").at( "MinResRange" ).get<double>();
    MaxResRange = conf.at("PIDAAvg").at( "MaxResRange" ).get<double>();
    MaxPIDAValue = conf.at("PIDAAvg").at( "MaxPIDAValue" ).get<double>();
    ExponentConstant = conf.at("PIDAAvg").at( "ExponentConstant" ).get<double>();
  }

}
