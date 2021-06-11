cmake_minimum_required(VERSION 3.15)
project(pi0_mc_xsec_variables)

set( CMAKE_CXX_STANDARD 17 )

find_package( ROOT 6.22 CONFIG REQUIRED )

include_directories(. ../contrib ../utilities)

add_executable( pi0_mc_xsec_study
        CMakeLists.txt
        pi0_mc_xsec.cpp
        pi0_mc_xsec.h
        ../utilities/Utilities.hpp
        ../utilities/Histograms.cpp
        ../utilities/Histograms.hpp
        ../utilities/CrossSection.cpp
        ../utilities/CrossSection.h)



target_link_libraries(pi0_mc_xsec_study
        PUBLIC ROOT::Core
        ROOT::Physics
        ROOT::Tree
        ROOT::Hist
        ROOT::TreePlayer
        )
