cmake_minimum_required(VERSION 3.15)
project(pi0_reco_xsec_proj)

set( CMAKE_CXX_STANDARD 17 )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

find_package( ROOT 6.22 CONFIG REQUIRED )

include_directories(. ../contrib ../utilities datatypes)

# Add the custom ROOT data dictionary
add_subdirectory( datatypes )

add_executable( pi0_reco_xsec
        CMakeLists.txt
        pi0_reco_xsec.cpp
        pi0_reco_xsec.h
        EventSelection.cpp
        EventSelection.h
        selection.h
        ../utilities/CrossSection.cpp
        ../utilities/CrossSection.h
        ../utilities/Utilities.hpp
        ../utilities/Histograms.cpp
        ../utilities/Histograms.hpp
        ../utilities/Chi2PIDAlg.cxx
        ../utilities/Chi2PIDAlg.hpp)

add_library( libVecs SHARED IMPORTED )
set_target_properties( libVecs PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/build/datatypes/libVecs.dylib )

target_link_libraries(pi0_reco_xsec
        PUBLIC Vecs
        libVecs
        ROOT::Core
        ROOT::Physics
        ROOT::Tree
        ROOT::Hist
        ROOT::TreePlayer
        )
