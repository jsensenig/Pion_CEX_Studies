//
// Created by Jon Sensenig on 12/5/20.
//

#include "Histograms.hpp"


Histograms::Histograms()
{ ; }

Histograms::~Histograms()
{ ; }

bool Histograms::ConfigureHistos( std::string config_file ) {

  json conf = utils::LoadConfig( config_file );

  if ( conf == nullptr ) return false; // Failed to load config

  json hists_1d = conf.at("1d_hists");

  for ( auto & h : hists_1d ) {

    auto name = h.at("name").get<std::string>();
    th1_hists[name] = Create1dHist( h );

  }

  json hists_2d = conf.at("2d_hists");

  for ( auto & h : hists_2d ) {

    auto name = h.at("name").get<std::string>();
    th2_hists[name] = Create2dHist( h );

  }

  return true;

}

std::unique_ptr<TH1> Histograms::Create1dHist( json & c ) {

  auto name = c.at( "name" ).get<std::string>();
  auto type = c.at( "type" ).get<std::string>();
  auto axes = c.at( "axes" ).get<std::string>();

  int bins = c.at( "bins" ).get<int>();
  auto u_lim = c.at( "u_lim" ).get<double>();
  auto l_lim = c.at( "l_lim" ).get<double>();

  if ( type == "TH1I" ) return std::make_unique<TH1I>( name.c_str(), axes.c_str(), bins, l_lim, u_lim );
  else if ( type == "TH1D" ) return std::make_unique<TH1D>( name.c_str(), axes.c_str(), bins, l_lim, u_lim );
  else if ( type == "TH1F" ) return std::make_unique<TH1F>( name.c_str(), axes.c_str(), bins, l_lim, u_lim );
  else return nullptr;

}

std::unique_ptr<TH2> Histograms::Create2dHist( json & c ) {

  auto name = c.at( "name" ).get<std::string>();
  auto type = c.at( "type" ).get<std::string>();
  auto axes = c.at( "axes" ).get<std::string>();

  auto xbins = c.at( "xbins" ).get<int>();
  auto xu_lim = c.at( "xu_lim" ).get<double>();
  auto xl_lim = c.at( "xl_lim" ).get<double>();
  auto ybins = c.at( "ybins" ).get<int>();
  auto yu_lim = c.at( "yu_lim" ).get<double>();
  auto yl_lim = c.at( "yl_lim" ).get<double>();

  if ( type == "TH2I" )
    return std::make_unique<TH2I>( name.c_str(), axes.c_str(), xbins, xl_lim, xu_lim, ybins, yl_lim, yu_lim );
  else if ( type == "TH2D" )
    return std::make_unique<TH2D>( name.c_str(), axes.c_str(), xbins, xl_lim, xu_lim, ybins, yl_lim, yu_lim );
  else if ( type == "TH2F" )
    return std::make_unique<TH2F>( name.c_str(), axes.c_str(), xbins, xl_lim, xu_lim, ybins, yl_lim, yu_lim );
  else return nullptr;

}

void Histograms::WriteHistos( TString & out_file ) {

  //if ( ! OpenFile( out_file ) ) return;
  auto ofile = std::make_unique<TFile>(out_file, "recreate");

  for ( auto & h : th1_hists ) {
    h.second -> Write();
  }

  for ( auto & h : th2_hists ) {
    h.second -> Write();
  }

  ofile -> Close();

}

bool Histograms::OpenFile( TString & out_file ) {

  ofile = std::make_unique<TFile>(out_file,"recreate");

  if ( !ofile -> IsOpen() ) {
    std::cout << " Failed to open output file " << out_file << std::endl;
    return false;
  }
  return true;
}
