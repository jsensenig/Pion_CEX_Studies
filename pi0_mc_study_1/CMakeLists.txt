cmake_minimum_required( VERSION 3.15 )
project(pi0_study)

set( CMAKE_CXX_STANDARD 17 )

find_package( ROOT 6.22 CONFIG REQUIRED )

include_directories(. ../contrib ../utilities)

add_executable( pi0_mc_study
        CMakeLists.txt
        pi0_mc_study.cpp
        pi0_mc_study.h
        ../utilities/Utilities.hpp
        ../utilities/Histograms.cpp
        ../utilities/Histograms.hpp)


target_link_libraries(pi0_mc_study
        PUBLIC ROOT::Core
        ROOT::Physics
        ROOT::Tree
        ROOT::Hist
        ROOT::TreePlayer
        )
