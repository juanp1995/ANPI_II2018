# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/juanp1995/ANPI_II2018/Proyecto 1/code"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build"

# Include any dependencies generated for this target.
include benchmarks/CMakeFiles/proyecto1_benchmarks.dir/depend.make

# Include the progress variables for this target.
include benchmarks/CMakeFiles/proyecto1_benchmarks.dir/progress.make

# Include the compile flags for this target's objects.
include benchmarks/CMakeFiles/proyecto1_benchmarks.dir/flags.make

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/flags.make
benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o: ../benchmarks/benchmarkMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o -c "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks/benchmarkMain.cpp"

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.i"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks/benchmarkMain.cpp" > CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.i

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.s"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks/benchmarkMain.cpp" -o CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.s

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.requires:

.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.requires

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.provides: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.requires
	$(MAKE) -f benchmarks/CMakeFiles/proyecto1_benchmarks.dir/build.make benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.provides.build
.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.provides

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.provides.build: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o


benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/flags.make
benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o: ../benchmarks/benchmarkMatrixAdd.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o -c "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks/benchmarkMatrixAdd.cpp"

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.i"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks/benchmarkMatrixAdd.cpp" > CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.i

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.s"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks/benchmarkMatrixAdd.cpp" -o CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.s

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.requires:

.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.requires

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.provides: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.requires
	$(MAKE) -f benchmarks/CMakeFiles/proyecto1_benchmarks.dir/build.make benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.provides.build
.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.provides

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.provides.build: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o


# Object files for target proyecto1_benchmarks
proyecto1_benchmarks_OBJECTS = \
"CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o" \
"CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o"

# External object files for target proyecto1_benchmarks
proyecto1_benchmarks_EXTERNAL_OBJECTS =

bin/proyecto1_benchmarks: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o
bin/proyecto1_benchmarks: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o
bin/proyecto1_benchmarks: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/build.make
bin/proyecto1_benchmarks: src/libanpi.a
bin/proyecto1_benchmarks: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
bin/proyecto1_benchmarks: /usr/lib/x86_64-linux-gnu/libboost_system.so
bin/proyecto1_benchmarks: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
bin/proyecto1_benchmarks: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../bin/proyecto1_benchmarks"
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/proyecto1_benchmarks.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
benchmarks/CMakeFiles/proyecto1_benchmarks.dir/build: bin/proyecto1_benchmarks

.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/build

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/requires: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMain.cpp.o.requires
benchmarks/CMakeFiles/proyecto1_benchmarks.dir/requires: benchmarks/CMakeFiles/proyecto1_benchmarks.dir/benchmarkMatrixAdd.cpp.o.requires

.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/requires

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/clean:
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" && $(CMAKE_COMMAND) -P CMakeFiles/proyecto1_benchmarks.dir/cmake_clean.cmake
.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/clean

benchmarks/CMakeFiles/proyecto1_benchmarks.dir/depend:
	cd "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/juanp1995/ANPI_II2018/Proyecto 1/code" "/home/juanp1995/ANPI_II2018/Proyecto 1/code/benchmarks" "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build" "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks" "/home/juanp1995/ANPI_II2018/Proyecto 1/code/build/benchmarks/CMakeFiles/proyecto1_benchmarks.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : benchmarks/CMakeFiles/proyecto1_benchmarks.dir/depend

