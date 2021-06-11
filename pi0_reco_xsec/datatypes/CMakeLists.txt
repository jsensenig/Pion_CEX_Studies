project(custom_dicts)


include_directories(.)

#Generate a ROOT dictionary.
ROOT_GENERATE_DICTIONARY(VecsDict
        vector_vector.h
        LINKDEF LinkDef.h
        MODULE Vecs)

# Create the Vecs library with the ROOT dictionary.
# Must be a shared library since ROOT will load it dynamically at runtime
add_library(Vecs SHARED
            VecsDict.cxx)

target_link_libraries(Vecs
        PUBLIC ROOT::Core
        ROOT::Tree
        ROOT::TreePlayer
        )
