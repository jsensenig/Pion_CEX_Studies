//
// Created by Jon Sensenig on 4/23/21.
//

#ifndef PI0_MC_STUDY_H
#define PI0_MC_STUDY_H

#include "../utilities/Histograms.hpp"
#include "../utilities/CrossSection.h"

/// Pi0 Data Structure
struct pi0_proxy {
  std::pair<int,int> id;            // Gamma ID
  std::pair<double,double> energy;  // Gamma E [MeV]
  std::pair<double,double> angle;   // Gamma angle [rad]
  double open_angle;                // Gamma open angle [rad]
  double pip_energy;                // Incident pi+ kinetic energy

  void reset() { // Reset all variables
    id         = {0, 0};
    energy     = {0., 0.};
    angle      = {0., 0.};
    open_angle = 0.;
    pip_energy = 0.;
  }
};

bool parseArgs(int argc, char ** argv);

void run_pi0_mc_xsec( const std::string& in_file, const std::string& in_tree, Histograms &hists, bool truth_xsec );

std::map<int, double> daughter_pi0_energy( Histograms& hists );

double pip_interaction_ke();

double pi0_cos_angle( pi0_proxy& pi0 );

double pi0_mom( pi0_proxy& pi0 );

double pi0_energy( pi0_proxy& pi0 );

double pi0_kinetic_energy( pi0_proxy& pi0 );

void plot_single_pi0( pi0_proxy& pi0, Histograms& hists, size_t idx );

void plot_all_pi0( pi0_proxy& pi0, double pi0_energy, Histograms& hists, size_t idx );

double open_angle( double px1, double py1, double pz1, double px2, double py2, double pz2 );

void clean_pointers();

/// Input args
int erange_arg;
bool truth_xsec_arg;
std::string output_file_override;

int true_cex_count;
int true_inc_piplus_count;
CrossSection cs;
Histograms hists;

/// Data Variables
int true_daughter_nPi0;
int true_daughter_nNeutron;
int true_daughter_nProton;
int true_daughter_nPiMinus;
int true_daughter_nPiPlus;
int true_beam_PDG;
double true_beam_startP, true_beam_endP;
double true_beam_endX, true_beam_endY, true_beam_endZ;
double true_beam_endPx, true_beam_endPy, true_beam_endPz;

auto true_beam_daughter_PDG = new std::vector<int>;
auto true_beam_daughter_ID = new std::vector<int>;
auto true_beam_daughter_startP = new std::vector<double>;
auto true_beam_daughter_startPx = new std::vector<double>;
auto true_beam_daughter_startPy = new std::vector<double>;
auto true_beam_daughter_startPz = new std::vector<double>;
auto true_beam_Pi0_decay_PDG = new std::vector<int>;
auto true_beam_Pi0_decay_parID = new std::vector<int>;
auto true_beam_Pi0_decay_startP = new std::vector<double>;
auto true_beam_Pi0_decay_startPx = new std::vector<double>;
auto true_beam_Pi0_decay_startPy = new std::vector<double>;
auto true_beam_Pi0_decay_startPz = new std::vector<double>;
auto true_beam_Pi0_decay_startX = new std::vector<double>;
auto true_beam_Pi0_decay_startY = new std::vector<double>;
auto true_beam_Pi0_decay_startZ = new std::vector<double>;
auto true_beam_Pi0_decay_len = new std::vector<double>;

auto true_beam_endProcess = new std::string;

double reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ;
double beam_inst_P;
std::vector<double> *reco_beam_calibrated_dEdX_SCE = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_energy = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_len = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_dirX = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_dirY = new std::vector<double>;
std::vector<double> *reco_daughter_allShower_dirZ = new std::vector<double>;

#endif //PI0_MC_STUDY_H

