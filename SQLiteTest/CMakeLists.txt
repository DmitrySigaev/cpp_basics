# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)

if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
message("CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
endif()


project(test)
message("VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")

message("find_package SQLITE3")
find_package(sqlite3)
message("SQLITE3_FOUND=${SQLITE3_FOUND}")
message("SQLITE3_INCLUDE_DIRS=${SQLITE3_INCLUDE_DIRS}")
message("SQLITE3_LIBRARIES=${SQLITE3_LIBRARIES}")
message("SQLITE3_INCLUDE_DIR=${SQLITE3_INCLUDE_DIR}")
message("SQLITE3_LIBRARY=${SQLITE3_LIBRARY}")



# Look for the header file.
message("Look for the header file SQLITE3.h")

FIND_PATH(SQLITE3_INCLUDE_DIR NAMES sqlite3.h)

# Look for the library.
FIND_LIBRARY(SQLITE3_LIBRARY NAMES sqlite3)
          
message("SQLITE3_FOUND=${SQLITE3_FOUND}")
message("SQLITE3_INCLUDE_DIR=${SQLITE3_INCLUDE_DIR}")
message("SQLITE3_LIBRARY=${SQLITE3_LIBRARY}")

# Look for the header file.
FIND_PATH(SQLITE3_INCLUDE_DIR NAMES sqlite3.h)

# Look for the library.
FIND_LIBRARY(SQLITE3_LIBRARY NAMES sqlite3)

message("SQLITE3_FOUND=${SQLITE3_FOUND}")
message("SQLITE3_INCLUDE_DIR=${SQLITE3_INCLUDE_DIR}")
message("SQLITE3_LIBRARY=${SQLITE3_LIBRARY}")


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SQLITE3 DEFAULT_MSG SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR)

message("SQLITE3_FOUND=${SQLITE3_FOUND}")
message("SQLITE3_INCLUDE_DIR=${SQLITE3_INCLUDE_DIR}")
message("SQLITE3_LIBRARY=${SQLITE3_LIBRARY}")

add_executable(main main.cpp)

if (SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIR})
  target_link_libraries (main sqlite3)
endif (SQLITE3_FOUND)         
