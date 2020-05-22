# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(test)
message("VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")

find_package (SQLite3 CONFIG REQUIRED)
if (SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIRS})
  target_link_libraries (tutorial ${SQLITE3_LIBRARIES})
endif (SQLITE3_FOUND)
