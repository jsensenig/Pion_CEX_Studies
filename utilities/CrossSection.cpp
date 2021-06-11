//
// Created by Jon Sensenig on 6/5/21.
//

#include "CrossSection.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <sstream>


CrossSection::CrossSection()
{

  // Define bin edges 2GeV/c
//  std::vector<double> beam_bins = {1000, 1500.,1750., 2000.};
//  std::vector<double> energy_bins = {0., 200., 400., 600., 800., 1000.};
//  std::vector<double> angle_bins = {-1., -0.5, 0., 0.25, 0.5, 0.75, 1.};

  // Define bin edges 1GeV/c
  std::vector<double> beam_bins = {0, 250., 500.,750., 1000.};
  std::vector<double> energy_bins = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};
  std::vector<double> angle_bins = {-1., -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.};

  size_t bbins = beam_bins.size() - 1;
  size_t ebins = energy_bins.size() - 1;
  size_t abins = angle_bins.size() - 1;

  // Each variable plot
  beam_energy_hist = std::make_unique<TH1D>( "beam_energy", "Beam Pi+ Kinetic Energy;T_{#pi^{+}} [MeV/c];Count", bbins, beam_bins.data() );
  energy_hist = std::make_unique<TH1D>( "pi0_energy", "Pi0 Kinetic Energy;T_{#pi^{0}} [MeV/c];Count", ebins, energy_bins.data() );
  angle_hist = std::make_unique<TH1D>( "pi0_angle", "Pi0 Angle;cos(#theta_{pi^{0}});Count", abins, angle_bins.data() );

  // Differential cross section variables
  std::string energy_angle_title("Pi0 Angle vs Kinetic Energy;T_{#pi^{0}} [MeV/c];cos(#theta_{#pi^{0}})");
  energy_angle_hist = std::make_unique<TH2D>( "pi0_energy_angle", energy_angle_title.c_str(),
                                              ebins, energy_bins.data(), abins,angle_bins.data() );

  // Total 3D histogram. Dimensions are: {beam pi+ interaction KE, pi0 KE, pi0 scattering angle}
  std::string title_3d("#pi0 Angle vs #pi0 KE vs Beam #pi+ KE;T_{#pi^{0}} [MeV/c];cos(#theta_{#pi^{0}});T_{#pi^{+}} [MeV/c]");
  beam_energy_angle_hist = std::make_unique<TH3D>( "beam_ke_pi0_energy_angle", title_3d.c_str(),
                                                   ebins, energy_bins.data(), abins, angle_bins.data(), bbins, beam_bins.data() );

  // Truth total 3D histogram. Dimensions are: {beam pi+ interaction KE, pi0 KE, pi0 scattering angle}
  std::string truth_title("Truth #pi0 Angle vs #pi0 KE vs Beam #pi+ KE;T_{#pi^{0}} [MeV/c];cos(#theta_{#pi^{0}});T_{#pi^{+}} [MeV/c]");
  truth_beam_energy_angle_hist = std::make_unique<TH3D>( "truth_beam_ke_pi0_energy_angle", title_3d.c_str(),
                                                   ebins, energy_bins.data(), abins, angle_bins.data(), bbins, beam_bins.data() );


  geant_xsec_file = std::make_unique<TFile>("/Users/jsen/tmp_fit/cross_section_cex_n50k_T870_newcosonly.root");
  if( !geant_xsec_file->IsOpen() ) std::cout << "Can't find Xsec file!" << std::endl;

  // Get the double differential cross section file
  geant_dd_xsec = (TH2D*)geant_xsec_file->Get("inel_cex_870_MeV");

  std::cout << "Rebinned Geant Xsec histogram! xbins = " << geant_dd_xsec->GetXaxis()->GetNbins() << " ybins = "
            << geant_dd_xsec->GetYaxis()->GetNbins() << std::endl;

}

CrossSection::~CrossSection()
{
  geant_xsec_file->Close();
}

// ............................................................................
void CrossSection::FillRecoXsecHisto( double pip_ke, double energy, double angle ) {

  // Interacting beam pion energy
  beam_energy_hist -> Fill( pip_ke );
  // Bin the xsec variables
  energy_hist -> Fill( energy );
  angle_hist  -> Fill( angle );
  energy_angle_hist -> Fill( energy, angle );
  beam_energy_angle_hist -> Fill( energy, angle, pip_ke );

}

void CrossSection::FillTrueXsecHisto( double pip_ke, double energy, double angle ) {
  truth_beam_energy_angle_hist -> Fill( energy, angle, pip_ke );
}

// ............................................................................
void CrossSection::ExtractXsecEnergy( int total_events, TString& out_file, bool true_xsec ) {

  std::cout << "Calculating cross section" << std::endl;

  // Select histogram
  TH3D* xsec_hist;
  if( true_xsec ) xsec_hist = (TH3D*)truth_beam_energy_angle_hist->Clone();
  else xsec_hist = (TH3D*)beam_energy_angle_hist->Clone();

  // X = pi0 KE
  // Y = pi0 angle
  // Z = pi+ interaction KE
  // xsec calculation: https://ir.uiowa.edu/cgi/viewcontent.cgi?article=1518&context=etd (pg 54)

  std::map<std::string, TGraph*> xsec_graphs;

  // N_tgt = 39.9624 / (6.022e23) * (1.3973) = 4.7492e-23 = Ar atomic mass / (Avagadro's number * LAr density)
  double n_tgt = 39.9624 / ( 6.022e23 * 1.3973 );

  // Get the number of bins in energy and angle
  int energy_bins = xsec_hist -> GetNbinsX();
  int angle_bins  = xsec_hist -> GetNbinsY();
  int beam_bins   = xsec_hist -> GetNbinsZ();

  for( size_t bbins = 1; bbins < beam_bins+1; bbins++ ) { // beam KE

    double beam_center = xsec_hist->GetZaxis()->GetBinCenter( bbins );
    double bbin_width  = xsec_hist->GetZaxis()->GetBinWidth( bbins );

    for ( size_t abins = 1; abins < angle_bins + 1; abins++ ) { // angular xsec
      // Get the angular bin center and width
      double angle_center = xsec_hist->GetYaxis()->GetBinCenter( abins );
      double abin_width   = xsec_hist->GetYaxis()->GetBinWidth( abins );

      // Project out the energy so we can get the error bars
      TH1D * h_xerr = xsec_hist -> ProjectionX( "h_xerr", abins, abins, bbins, bbins, "e" );
      std::vector<double> xsec, true_xsec, xsec_yerr, true_xsec_yerr, energy, xsec_xerr;

      for ( size_t ebins = 1; ebins < energy_bins + 1; ebins++ ) { // energy xsec
        // Get the energy bin center and width
        double energy_center = xsec_hist->GetXaxis()->GetBinCenter( ebins );
        double ebin_width = xsec_hist->GetXaxis()->GetBinWidth( ebins );
        double Ni = xsec_hist->GetBinContent( ebins, abins, bbins );

        // xsec calculation
        double xsec_calc = (( Ni * n_tgt ) / ( total_events * ebin_width * abin_width * bbin_width )) * 1.e27;
        xsec.emplace_back( xsec_calc ); // [milli-barn (mb)]
        xsec_yerr.emplace_back(( xsec_calc / Ni ) * h_xerr->GetBinError( ebins ));

        // True xsec
        true_xsec.emplace_back( GetDDXsec( energy_center, angle_center, beam_center ) );
        true_xsec_yerr.emplace_back(0.);

        energy.emplace_back( energy_center );
        xsec_xerr.emplace_back( ebin_width / 2. );

        std::cout << "Ebin width " << ebin_width << " Abin width " << abin_width << " Ni " << Ni
                  << " Energy " << energy_center << " Angle " << angle_center << " Xsec " << xsec_calc << std::endl;

      } //energy
      // Get TGraph
      std::stringstream gr_name;
      gr_name << "beam_" << (int)beam_center << "_angle_" << (int)( TMath::ACos( angle_center ) * TMath::RadToDeg() );
      xsec_graphs[gr_name.str()] = SinglePlotXsecEnergy( xsec, angle_center, energy, beam_center, xsec_xerr, xsec_yerr, false );

      std::stringstream true_gr_name;
      true_gr_name << "true_beam_" << (int)beam_center << "_angle_" << (int)( TMath::ACos( angle_center ) * TMath::RadToDeg() );
      xsec_graphs[true_gr_name.str()] = SinglePlotXsecEnergy( true_xsec, angle_center, energy, beam_center, xsec_xerr, true_xsec_yerr, true );
    } //angle
  }

  // Write all TGraphs to file
  auto ofile = std::make_unique<TFile>(out_file, "recreate");
  for( auto& graph : xsec_graphs ) graph.second->Write( graph.first.c_str() );
  beam_energy_hist -> Write();
  energy_hist -> Write();
  angle_hist -> Write();
  energy_angle_hist -> Write();
  beam_energy_angle_hist -> Write();
  truth_beam_energy_angle_hist -> Write();
  ofile -> Close();

}

// ............................................................................
TGraph* CrossSection::SinglePlotXsecEnergy( std::vector<double> &xsec, double angle, std::vector<double> &energy, double beam_energy,
                                      std::vector<double> &xerr, std::vector<double> &yerr, bool true_xsec ) {

  std::cout << "Writing Xsec to file" << std::endl;

  auto xsec_graph = new TGraphErrors( energy.size(), energy.data(), xsec.data(), xerr.data(), yerr.data());
  xsec_graph->SetLineWidth( 1 ); xsec_graph->SetMarkerStyle( 8 ); xsec_graph->SetMarkerSize( 0.5 );
  if( true_xsec ) xsec_graph->SetLineColor( 64 );
  else xsec_graph->SetLineColor( 46 );

  std::stringstream title;
  title << "1 GeV/c #pi^{+} CEX Cross-section (Angle = " << (int)( TMath::ACos( angle ) * TMath::RadToDeg() ) << "#circ)"
        << " (T_{#pi^{+}} = " << (int)beam_energy << ")";
  xsec_graph->SetTitle( title.str().c_str() );
  xsec_graph->GetXaxis()->SetTitle( "T_{#pi^{0}} [MeV/c]" );
  xsec_graph->GetYaxis()->SetTitle( "#frac{d^{2}#sigma}{dT_{#pi^{0}}d#Omega_{#pi^{0}}} [mb/MeV/c#upointsr]" );

  return xsec_graph;

}

// ............................................................................
void CrossSection::ExtractXsecAngle( int total_events, TString& out_file, bool true_xsec ) {

  std::cout << "Calculating cross section" << std::endl;

  // Select histogram
  TH3D* xsec_hist;
  if( true_xsec ) xsec_hist = (TH3D*)truth_beam_energy_angle_hist->Clone();
  else xsec_hist = (TH3D*)beam_energy_angle_hist->Clone();

  // X = pi0 KE
  // Y = pi0 angle
  // Z = pi+ interaction KE
  // xsec calculation: https://ir.uiowa.edu/cgi/viewcontent.cgi?article=1518&context=etd (pg 54)

  std::map<std::string, TGraph*> xsec_graphs;

  // N_tgt = 39.9624 / (6.022e23) * (1.3973) = 4.7492e-23 = Ar atomic mass / (Avagadro's number * LAr density)
  double n_tgt = 39.9624 / ( 6.022e23 * 1.3973 );

  // Get the number of bins in energy and angle
  int energy_bins = xsec_hist -> GetNbinsX();
  int angle_bins  = xsec_hist -> GetNbinsY();
  int beam_bins   = xsec_hist -> GetNbinsZ();

  for( size_t bbins = 1; bbins < beam_bins+1; bbins++ ) { // beam KE

    double beam_center = xsec_hist->GetZaxis()->GetBinCenter( bbins );
    double bbin_width = xsec_hist->GetZaxis()->GetBinWidth( bbins );

    for ( size_t ebins = 1; ebins < energy_bins + 1; ebins++ ) { // energy xsec
      // Get the energy bin center and width
      double energy_center = xsec_hist->GetXaxis()->GetBinCenter( ebins );
      double ebin_width = xsec_hist->GetXaxis()->GetBinWidth( ebins );

      // Project out the angle so we can get the error bars
      TH1D * h_xerr = xsec_hist->ProjectionY( "h_xerr", ebins, ebins, bbins, bbins, "e" );

      std::vector<double> xsec, true_xsec, xsec_yerr, true_xsec_yerr, angle, xsec_xerr;

      for ( size_t abins = 1; abins < angle_bins + 1; abins++ ) { // angular xsec
        // Get the angle bin center and width
        double angle_center = xsec_hist->GetYaxis()->GetBinCenter( abins );
        double abin_width = xsec_hist->GetYaxis()->GetBinWidth( abins );
        double Ni = xsec_hist->GetBinContent( ebins, abins, bbins );

        // xsec calculation
        double xsec_calc = (( Ni * n_tgt ) / ( total_events * ebin_width * abin_width * bbin_width )) * 1.e27; // [mb]
        xsec.emplace_back( xsec_calc );
        xsec_yerr.emplace_back(( xsec_calc / Ni ) * h_xerr->GetBinError( abins ));

        // True xsec
        true_xsec.emplace_back( GetDDXsec( energy_center, angle_center, beam_center ) );
        true_xsec_yerr.emplace_back(0.);

        angle.emplace_back( angle_center );
        xsec_xerr.emplace_back( abin_width / 2. );

        std::cout << "Ebin width " << ebin_width << " Abin width " << abin_width << " Ni " << Ni
                  << " Energy " << energy_center << " Angle " << angle_center << " Xsec " << xsec_calc << std::endl;
      }
      std::stringstream gr_name;
      gr_name << "beam_" << (int)beam_center << "_energy_" << (int)energy_center;
      xsec_graphs[gr_name.str()] = SinglePlotXsecAngle( xsec, energy_center, angle, beam_center, xsec_xerr, xsec_yerr, false );

      std::stringstream true_gr_name;
      true_gr_name << "true_beam_" << (int)beam_center << "_angle_" << (int)energy_center;
      xsec_graphs[true_gr_name.str()] = SinglePlotXsecAngle( true_xsec, energy_center, angle, beam_center, xsec_xerr, true_xsec_yerr, true );
      }
    }

  // Write all TGraphs to file
  auto ofile = std::make_unique<TFile>(out_file,"recreate");
  for( auto& graph : xsec_graphs ) graph.second->Write( graph.first.c_str() );
  beam_energy_hist -> Write();
  energy_hist -> Write();
  angle_hist -> Write();
  energy_angle_hist -> Write();
  beam_energy_angle_hist -> Write();
  truth_beam_energy_angle_hist -> Write();
  ofile -> Close();

}

// ............................................................................
TGraph* CrossSection::SinglePlotXsecAngle( std::vector<double> &xsec, double energy, std::vector<double> &angle, double beam_energy,
                                            std::vector<double> &xerr, std::vector<double> &yerr, bool true_xsec ) {

  std::cout << "Writing Xsec to file" << std::endl;

  auto xsec_graph = new TGraphErrors( angle.size(), angle.data(), xsec.data(), xerr.data(), yerr.data());
  xsec_graph->SetLineWidth( 1 ); xsec_graph->SetMarkerStyle( 8 ); xsec_graph->SetMarkerSize( 0.5 );
  if( true_xsec ) xsec_graph->SetLineColor( 64 );
  else xsec_graph->SetLineColor( 46 );

  std::stringstream title;
  title << "#pi^{+} CEX Cross-section (T_{#pi^{0}} = " << (int)energy << " [MeV/c])"
        << " (T_{#pi^{+}} = " << beam_energy << ")";
  xsec_graph->SetTitle( title.str().c_str() );
  xsec_graph->GetXaxis()->SetTitle( "cos#theta_{#pi^{0}}" );
  xsec_graph->GetYaxis()->SetTitle( "#frac{d^{2}#sigma}{dT_{#pi^{0}}d#Omega_{#pi^{0}}} [mb/MeV/c#upointsr]" );

  return xsec_graph;

}

// ............................................................................
double CrossSection::GetDDXsec( double energy, double angle, double beam ) {

  int global_bin = geant_dd_xsec -> FindBin( energy, angle );
  return geant_dd_xsec -> GetBinContent( global_bin );

}