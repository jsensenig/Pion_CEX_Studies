//
// Created by Jon Sensenig on 2/15/21.
//

#ifndef PION_QE_CHI2PIDALG_HPP
#define PION_QE_CHI2PIDALG_HPP

#include "Utilities.hpp"
#include "TProfile.h"

#include <string>
#include <vector>

class TProfile;

namespace pid {

  class Chi2PIDAlg {

   public:

     Chi2PIDAlg();
     ~Chi2PIDAlg();

//    std::pair<double, int> Chi2PID( const std::vector<double> & track_dedx, const std::vector<double> & range, int pdg);
    double Chi2PID( const std::vector<double> & track_dedx, const std::vector<double> & range, int pdg );
    TProfile * SelectProfile( int pdg );
    double PidaFit( const std::vector<double> & track_dedx, const std::vector<double> & range );
    void Config();

   private:

     std::string fTemplateFile = "../../contrib/dEdxrestemplates.root";

     TProfile * dedx_range_pro;
     TProfile * dedx_range_ka;
     TProfile * dedx_range_pi;
     TProfile * dedx_range_mu;

     std::string _config_file = "../../utilities/Chi2PIDAlg.json";

    double MinResRange = 0;
    double MaxResRange = 30;
    double MaxPIDAValue = 50;
    double ExponentConstant = 0.42;

   };
}


#endif //PION_QE_CHI2PIDALG_HPP
