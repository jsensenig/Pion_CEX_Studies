//
// Created by Jon Sensenig on 12/5/20.
//

#ifndef PION_QE_HISTOGRAMS_HPP
#define PION_QE_HISTOGRAMS_HPP

#include <map>
#include <iostream>
#include "TH1I.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TFile.h"
#include "Utilities.hpp"


class Histograms {

public:

  Histograms();
  virtual ~Histograms();

  bool ConfigureHistos( std::string config_file );
  std::unique_ptr<TH1> Create1dHist( json & c );
  std::unique_ptr<TH2> Create2dHist( json & c );
  void WriteHistos( TString & out_file );

  // Hold all the histograms in these maps
  std::map< std::string, std::unique_ptr<TH1> > th1_hists;
  std::map< std::string, std::unique_ptr<TH2> > th2_hists;

protected:

  bool OpenFile( TString & out_file );

private:

  std::unique_ptr<TFile> ofile;

};


#endif //PION_QE_HISTOGRAMS_HPP
