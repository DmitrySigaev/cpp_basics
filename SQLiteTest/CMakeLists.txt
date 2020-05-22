# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(test)
message("VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")

find_package(Sqlite3 REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main sqlite3)
