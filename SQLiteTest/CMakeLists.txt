# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(test)
message("VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")

find_package (SQLite3 REQUIRED)
message("SQLITE3_FOUND=${SQLITE3_FOUND}")
message("SQLITE3_INCLUDE_DIRS=${SQLITE3_INCLUDE_DIRS}")
message("SQLITE3_LIBRARIES=${SQLITE3_LIBRARIES}")

# Look for the header file.
FIND_PATH(SQLITE3_INCLUDE_DIR NAMES sqlite3.h)

# Look for the library.
FIND_LIBRARY(SQLITE3_LIBRARY NAMES sqlite)

message("SQLITE3_FOUND=${SQLITE3_FOUND}")
message("SQLITE3_INCLUDE_DIRS=${SQLITE3_INCLUDE_DIRS}")
message("SQLITE3_LIBRARIES=${SQLITE3_LIBRARIES}")


if (SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIRS})
  target_link_libraries (main ${SQLITE3_LIBRARIES})
endif (SQLITE3_FOUND)
