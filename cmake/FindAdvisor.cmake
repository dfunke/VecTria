
# Copyright (c) 2013-2014 Stefan.Eilemann@epfl.ch
# finds the ittnotify API

#  Advisor_FOUND - system has Advisor
#  Advisor_INCLUDE_DIRS - the Advisor include directories
#  Advisor_LIBRARIES - link these to use Advisor

include(LibFindMacros)

find_program(Advisor_EXECUTABLE advixe-cl
  HINTS $ENV{ADVISOR_2019_DIR}
  PATH_SUFFIXES bin64
)

if(NOT Advisor_EXECUTABLE)
  return()
endif()

get_filename_component(Advisor_DIR ${Advisor_EXECUTABLE} PATH)
set(Advisor_DIR "${Advisor_DIR}/..")

# Include dir
find_path(Advisor_INCLUDE_DIR
  NAMES advisor-annotate.h
  PATHS ${Advisor_DIR}
  PATH_SUFFIXES include
)

# Finally the library itself
find_library(Advisor_LIBRARY
  NAMES advisor
  PATHS ${Advisor_DIR}
  PATH_SUFFIXES lib64
)

if(UNIX)
  list(APPEND Advisor_LIBRARY dl)
endif()

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(Advisor_PROCESS_INCLUDES Advisor_INCLUDE_DIR)
set(Advisor_PROCESS_LIBS Advisor_LIBRARY)
libfind_process(Advisor)