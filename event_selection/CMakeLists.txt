cmake_minimum_required( VERSION 3.15 )
project(event_selection)

set( CMAKE_CXX_STANDARD 17 )

find_package( ROOT 6.22 CONFIG REQUIRED )

include_directories(. ../contrib ../utilities)

add_executable( evt_selection
        CMakeLists.txt
        event_selection.cpp
        event_selection.h
        ../utilities/Utilities.hpp
        ../utilities/Histograms.cpp
        ../utilities/Histograms.hpp)


target_link_libraries(evt_selection
        PUBLIC ROOT::Core
        ROOT::Physics
        ROOT::Tree
        ROOT::Hist
        ROOT::TreePlayer
        )
