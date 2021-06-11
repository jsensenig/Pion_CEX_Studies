//
// Created by Jon Sensenig on 12/3/20.
//

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "../contrib/json.hpp"

#include "TString.h"

#include <iostream>
#include <fstream>
#include <vector>

using json = nlohmann::json;

namespace utils {

  namespace pdg {

    const int kPdgPositron = -11;
    const int kPdgElectron = 11;
    const int kPdgMuon = 13;
    const int kPdgAntiMuon = -13;
    const int kPdgGamma = 22;
    const int kPdgPiP = 211; // pi+
    const int kPdgPiM = -211; // pi-
    const int kPdgPi0 = 111; // pi0
    const int kPdgKP = 321; // K+
    const int kPdgKM = -321; // K-
    const int kPdgK0 = 311; // K0
    const int kPdgProton = 2212; //
    const int kPdgAntiProton = -2212; //
    const int kPdgNeutron = 2112; //
    const int kPdgAntiNeutron = -2112; //

    static std::string pdg2string( int p ) {

      switch ( p ) {
        case kPdgMuon     :
          return "Muon";
        case kPdgAntiMuon :
          return "AntiMuon";
        case kPdgPiP      :
          return "Pion";
        case kPdgPi0      :
          return "Pi0";
        case kPdgKP       :
          return "Kaon";
        case kPdgProton   :
          return "Proton";
        case kPdgNeutron  :
          return "Neutron";
        case kPdgPositron :
          return "Positron";
        case kPdgElectron :
          return "Electron";
        case kPdgGamma :
          return "Gamma";
        default           :
          return "Unknown";
      }

    }

    // Mass in units MeV/c
    static double pdg2mass( int p ) {

      switch ( p ) {
        case kPdgMuon     :
          return 105.658;
        case kPdgAntiMuon :
          return 105.658;
        case kPdgPiP      :
          return 139.570;
        case kPdgPi0      :
          return 134.977;
        case kPdgKP       :
          return 493.677;
        case kPdgProton   :
          return 938.272;
        case kPdgNeutron  :
          return 939.565;
        case kPdgPositron :
          return 0.511;
        case kPdgElectron :
          return 0.511;
        case kPdgGamma :
          return 0.0;
        default           :
          return -1.0;
      }

    }

  }


//@brief Utility to load a list of files to process from a file.
//
//----------------------------------------------------------------
  static std::vector<TString> LoadFileList( const std::string & file_list ) {

    std::vector<TString> file_vec;
    std::string line;

    std::ifstream file( file_list );

    if ( file.is_open()) {
      while ( getline( file, line )) {
        std::cout << "Loading file: " << line << std::endl;
        file_vec.emplace_back( line );
      }
      file.close();
    } else std::cout << "Unable to open file " << file_list << std::endl;

    return file_vec;

  }

  static json LoadConfig( const std::string & file ) {

    std::ifstream injson( file );

    if ( !injson.is_open()) {

      std::cout << "Failed to open config file " << file << std::endl;
      return nullptr;

    }

    json parsed;
    try {

      parsed = json::parse( injson );

    } catch ( json::exception & ex ) {

      std::cout << " Parse error in file " << file << " " << ex.what() << std::endl;
      return nullptr;

    }

    return parsed;

  }

  //------------------------------------------------------------
  // Rotate vector v into dir reference frame then return the theta angle in that frame
//  static double ThetaAngle( const Vector_t & dir, Vector_t v ) {
//
//    //v.RotateUz( dir.Unit() ) ;   FIXME
//    //ROOT::Math::AxisAngle r(dir, 0);
//    //return (r * v).Theta();
//    return v.Theta();
//
//  }

  //-------------------------------------------------------------
  // Rotate vector v into dir reference frame then return the phi angle in that frame
//  static double PhiAngle( const Vector_t & dir, Vector_t v ) {
//
//    //v.RotateUz( dir.Unit() ) ;  FIXME
//    //ROOT::Math::AxisAngle r(dir, 0);
//    //return (r * v).Phi();
//    return v.Phi();
//
//  }

  //-------------------------------------------------------------
//  // Calculate kinetic energy for a given momentum and mass
//  static double CalculateKE( const Vector_t & p, const double m ) {
//
//    return sqrt( p.Mag2() + pow( m, 2 )) - m;
//
//  }

  //-------------------------------------------------------------
  // Calculate kinetic energy for a given momentum and mass
  static double CalculateKE( const double p, const double m ) {

    return sqrt( pow( p, 2 ) + pow( m, 2 ) ) - m;

  }

  //-------------------------------------------------------------
  // Calculate kinetic energy for a given momentum and mass
  static double CalculateE( const double p, const double m ) {
                                                                  
    return sqrt( pow( p, 2 ) + pow( m, 2 ) );
                                                                  
  }

  // Template function to count the occurrence of a specified TTree array elements
  // ----------------------------------------------------------------
  template<typename T>
  int Count( std::vector<T> & vec, T a ) {

    return std::count( vec.begin(), vec.end(), a );

  }

  // Template function to find an element in a TTree array
  //----------------------------------------------------------------
//  template<typename T>
//  T Find( TTreeReaderArray<T> & arr, T a ) {
//
//    auto it = std::find( arr.begin(), arr.end(), a );
//
//    if ( it != arr.end()) return *it;
//    else return -1;
//
//  }

// Template function to find and return an element index in a TTree array
//----------------------------------------------------------------
  template<typename T>
  int FindIndex( std::vector<T> & vec, T a ) {

    auto it = std::find( vec.begin(), vec.end(), a );

    if ( it != vec.end()) return it - vec.begin();
    else return -1;

  }

// Template function to find and return an element index in a TTree array
//----------------------------------------------------------------
  static double dR( double xbeam, double ybeam, double zbeam, double xdaughter, double ydaughter, double zdaughter ) {

    return sqrt( (xbeam - xdaughter)*(xbeam - xdaughter) + 
                 (ybeam - ydaughter)*(ybeam - ydaughter) +
                 (zbeam - zdaughter)*(zbeam - zdaughter) );
  }

  static double Distance( double x, double y, double z ) {

    return sqrt( pow(x,2) + pow(y,2) + pow(z,2) );

  }

}

#endif //UTILITIES_HPP
