# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zenop/astrotools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zenop/astrotools/_build

# Include any dependencies generated for this target.
include tools/zeno/CMakeFiles/propAidaValidate.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tools/zeno/CMakeFiles/propAidaValidate.dir/compiler_depend.make

# Include the progress variables for this target.
include tools/zeno/CMakeFiles/propAidaValidate.dir/progress.make

# Include the compile flags for this target's objects.
include tools/zeno/CMakeFiles/propAidaValidate.dir/flags.make

tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/flags.make
tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o: ../tools/zeno/propAidaValidate.cpp
tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zenop/astrotools/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o -MF CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o.d -o CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o -c /home/zenop/astrotools/tools/zeno/propAidaValidate.cpp

tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.i"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zenop/astrotools/tools/zeno/propAidaValidate.cpp > CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.i

tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.s"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zenop/astrotools/tools/zeno/propAidaValidate.cpp -o CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.s

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/flags.make
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o: ../modules/cfg/Config.cpp
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zenop/astrotools/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o -MF CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o.d -o CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o -c /home/zenop/astrotools/modules/cfg/Config.cpp

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.i"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zenop/astrotools/modules/cfg/Config.cpp > CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.i

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.s"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zenop/astrotools/modules/cfg/Config.cpp -o CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.s

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/flags.make
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o: ../modules/cfg/Environment.cpp
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zenop/astrotools/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o -MF CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o.d -o CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o -c /home/zenop/astrotools/modules/cfg/Environment.cpp

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.i"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zenop/astrotools/modules/cfg/Environment.cpp > CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.i

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.s"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zenop/astrotools/modules/cfg/Environment.cpp -o CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.s

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/flags.make
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o: ../modules/geco/JsonParser.cpp
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zenop/astrotools/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o -MF CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o.d -o CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o -c /home/zenop/astrotools/modules/geco/JsonParser.cpp

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.i"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zenop/astrotools/modules/geco/JsonParser.cpp > CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.i

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.s"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zenop/astrotools/modules/geco/JsonParser.cpp -o CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.s

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/flags.make
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o: ../modules/geco/MapEnum.cpp
tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o: tools/zeno/CMakeFiles/propAidaValidate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zenop/astrotools/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o -MF CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o.d -o CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o -c /home/zenop/astrotools/modules/geco/MapEnum.cpp

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.i"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zenop/astrotools/modules/geco/MapEnum.cpp > CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.i

tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.s"
	cd /home/zenop/astrotools/_build/tools/zeno && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zenop/astrotools/modules/geco/MapEnum.cpp -o CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.s

# Object files for target propAidaValidate
propAidaValidate_OBJECTS = \
"CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o" \
"CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o" \
"CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o" \
"CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o" \
"CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o"

# External object files for target propAidaValidate
propAidaValidate_EXTERNAL_OBJECTS =

tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/propAidaValidate.cpp.o
tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Config.cpp.o
tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/cfg/Environment.cpp.o
tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/JsonParser.cpp.o
tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/__/__/modules/geco/MapEnum.cpp.o
tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/build.make
tools/zeno/propAidaValidate: /usr/local/lib/libdace.so
tools/zeno/propAidaValidate: modules/atmos/libatmos.so.0.1.0
tools/zeno/propAidaValidate: modules/cfg/libcfg.so.0.1.0
tools/zeno/propAidaValidate: modules/astro/libastro.so.0.1.0
tools/zeno/propAidaValidate: /usr/local/lib/libdace.so
tools/zeno/propAidaValidate: modules/geco/libgeco.so.0.1.0
tools/zeno/propAidaValidate: /usr/local/lib/cspice.a
tools/zeno/propAidaValidate: /usr/lib/x86_64-linux-gnu/libjsoncpp.so
tools/zeno/propAidaValidate: tools/zeno/CMakeFiles/propAidaValidate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zenop/astrotools/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable propAidaValidate"
	cd /home/zenop/astrotools/_build/tools/zeno && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/propAidaValidate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/zeno/CMakeFiles/propAidaValidate.dir/build: tools/zeno/propAidaValidate
.PHONY : tools/zeno/CMakeFiles/propAidaValidate.dir/build

tools/zeno/CMakeFiles/propAidaValidate.dir/clean:
	cd /home/zenop/astrotools/_build/tools/zeno && $(CMAKE_COMMAND) -P CMakeFiles/propAidaValidate.dir/cmake_clean.cmake
.PHONY : tools/zeno/CMakeFiles/propAidaValidate.dir/clean

tools/zeno/CMakeFiles/propAidaValidate.dir/depend:
	cd /home/zenop/astrotools/_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zenop/astrotools /home/zenop/astrotools/tools/zeno /home/zenop/astrotools/_build /home/zenop/astrotools/_build/tools/zeno /home/zenop/astrotools/_build/tools/zeno/CMakeFiles/propAidaValidate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/zeno/CMakeFiles/propAidaValidate.dir/depend

