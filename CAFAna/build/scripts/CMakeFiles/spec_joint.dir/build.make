# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cvmfs/dune.opensciencegrid.org/products/dune/cmake/v3_12_2/Linux64bit+3.10-2.17/bin/cmake

# The command to remove a file.
RM = /cvmfs/dune.opensciencegrid.org/products/dune/cmake/v3_12_2/Linux64bit+3.10-2.17/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build

# Include any dependencies generated for this target.
include scripts/CMakeFiles/spec_joint.dir/depend.make

# Include the progress variables for this target.
include scripts/CMakeFiles/spec_joint.dir/progress.make

# Include the compile flags for this target's objects.
include scripts/CMakeFiles/spec_joint.dir/flags.make

scripts/CMakeFiles/spec_joint.dir/spec_joint.C.o: scripts/CMakeFiles/spec_joint.dir/flags.make
scripts/CMakeFiles/spec_joint.dir/spec_joint.C.o: ../scripts/spec_joint.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object scripts/CMakeFiles/spec_joint.dir/spec_joint.C.o"
	cd /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts && /cvmfs/larsoft.opensciencegrid.org/products/gcc/v6_4_0/Linux64bit+3.10-2.17/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/spec_joint.dir/spec_joint.C.o -c /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/spec_joint.C

scripts/CMakeFiles/spec_joint.dir/spec_joint.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/spec_joint.dir/spec_joint.C.i"
	cd /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts && /cvmfs/larsoft.opensciencegrid.org/products/gcc/v6_4_0/Linux64bit+3.10-2.17/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/spec_joint.C > CMakeFiles/spec_joint.dir/spec_joint.C.i

scripts/CMakeFiles/spec_joint.dir/spec_joint.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/spec_joint.dir/spec_joint.C.s"
	cd /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts && /cvmfs/larsoft.opensciencegrid.org/products/gcc/v6_4_0/Linux64bit+3.10-2.17/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/spec_joint.C -o CMakeFiles/spec_joint.dir/spec_joint.C.s

# Object files for target spec_joint
spec_joint_OBJECTS = \
"CMakeFiles/spec_joint.dir/spec_joint.C.o"

# External object files for target spec_joint
spec_joint_EXTERNAL_OBJECTS =

scripts/spec_joint: scripts/CMakeFiles/spec_joint.dir/spec_joint.C.o
scripts/spec_joint: scripts/CMakeFiles/spec_joint.dir/build.make
scripts/spec_joint: Analysis/libCAFAnaAnalysis.so
scripts/spec_joint: Fit/libCAFAnaFit.so
scripts/spec_joint: Decomp/libCAFAnaDecomp.so
scripts/spec_joint: Prediction/libCAFAnaPrediction.so
scripts/spec_joint: Core/libCAFAnaCore.so
scripts/spec_joint: Experiment/libCAFAnaExperiment.so
scripts/spec_joint: Systs/libCAFAnaSysts.so
scripts/spec_joint: PRISM/libCAFAnaPRISM.so
scripts/spec_joint: Cuts/libCAFAnaCuts.so
scripts/spec_joint: Extrap/libCAFAnaExtrap.so
scripts/spec_joint: Vars/libCAFAnaVars.so
scripts/spec_joint: libStandardRecord.so
scripts/spec_joint: libOscLibFunc.so
scripts/spec_joint: libUtilitiesFunc.so
scripts/spec_joint: scripts/CMakeFiles/spec_joint.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable spec_joint"
	cd /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spec_joint.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
scripts/CMakeFiles/spec_joint.dir/build: scripts/spec_joint

.PHONY : scripts/CMakeFiles/spec_joint.dir/build

scripts/CMakeFiles/spec_joint.dir/clean:
	cd /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts && $(CMAKE_COMMAND) -P CMakeFiles/spec_joint.dir/cmake_clean.cmake
.PHONY : scripts/CMakeFiles/spec_joint.dir/clean

scripts/CMakeFiles/spec_joint.dir/depend:
	cd /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/scripts/CMakeFiles/spec_joint.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : scripts/CMakeFiles/spec_joint.dir/depend

