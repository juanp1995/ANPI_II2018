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
CMAKE_SOURCE_DIR = /home/juanp1995/Proyecto_3/Proyecto3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/juanp1995/Proyecto_3/Proyecto3/build

# Include any dependencies generated for this target.
include src/CMakeFiles/proyecto3.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/proyecto3.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/proyecto3.dir/flags.make

src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o: src/CMakeFiles/proyecto3.dir/flags.make
src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o: ../src/Proyecto3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/juanp1995/Proyecto_3/Proyecto3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o"
	cd /home/juanp1995/Proyecto_3/Proyecto3/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/proyecto3.dir/Proyecto3.cpp.o -c /home/juanp1995/Proyecto_3/Proyecto3/src/Proyecto3.cpp

src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proyecto3.dir/Proyecto3.cpp.i"
	cd /home/juanp1995/Proyecto_3/Proyecto3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/juanp1995/Proyecto_3/Proyecto3/src/Proyecto3.cpp > CMakeFiles/proyecto3.dir/Proyecto3.cpp.i

src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proyecto3.dir/Proyecto3.cpp.s"
	cd /home/juanp1995/Proyecto_3/Proyecto3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/juanp1995/Proyecto_3/Proyecto3/src/Proyecto3.cpp -o CMakeFiles/proyecto3.dir/Proyecto3.cpp.s

src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.requires:

.PHONY : src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.requires

src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.provides: src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/proyecto3.dir/build.make src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.provides.build
.PHONY : src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.provides

src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.provides.build: src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o


# Object files for target proyecto3
proyecto3_OBJECTS = \
"CMakeFiles/proyecto3.dir/Proyecto3.cpp.o"

# External object files for target proyecto3
proyecto3_EXTERNAL_OBJECTS =

bin/proyecto3: src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o
bin/proyecto3: src/CMakeFiles/proyecto3.dir/build.make
bin/proyecto3: src/libanpi.a
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libboost_system.so
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libboost_program_options.so
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libGL.so
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libGLU.so
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libglut.so
bin/proyecto3: /usr/lib/x86_64-linux-gnu/libXi.so
bin/proyecto3: src/CMakeFiles/proyecto3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/juanp1995/Proyecto_3/Proyecto3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/proyecto3"
	cd /home/juanp1995/Proyecto_3/Proyecto3/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/proyecto3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/proyecto3.dir/build: bin/proyecto3

.PHONY : src/CMakeFiles/proyecto3.dir/build

src/CMakeFiles/proyecto3.dir/requires: src/CMakeFiles/proyecto3.dir/Proyecto3.cpp.o.requires

.PHONY : src/CMakeFiles/proyecto3.dir/requires

src/CMakeFiles/proyecto3.dir/clean:
	cd /home/juanp1995/Proyecto_3/Proyecto3/build/src && $(CMAKE_COMMAND) -P CMakeFiles/proyecto3.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/proyecto3.dir/clean

src/CMakeFiles/proyecto3.dir/depend:
	cd /home/juanp1995/Proyecto_3/Proyecto3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/juanp1995/Proyecto_3/Proyecto3 /home/juanp1995/Proyecto_3/Proyecto3/src /home/juanp1995/Proyecto_3/Proyecto3/build /home/juanp1995/Proyecto_3/Proyecto3/build/src /home/juanp1995/Proyecto_3/Proyecto3/build/src/CMakeFiles/proyecto3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/proyecto3.dir/depend

