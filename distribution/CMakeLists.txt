cmake_minimum_required(VERSION 3.0)
project(Distribution)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

set(CMAKE_CXX_STANDARD 11)

# Tell cmake where to find our modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(HepMC REQUIRED)
if(HepMC_FOUND)
    include_directories(${HepMC_INCLUDE_DIRS})
    link_directories(${HepMC_LIBRARY_DIRS})
endif()

find_package(Fastjet REQUIRED)
if(Fastjet_FOUND)
    include_directories(${Fastjet_INCLUDE_DIR})
    link_directories(${Fastjet_LIBRARIES})
endif()

# find_package(MPI REQUIRED)
# # add this line only when you are using openmpi which has a different c++ bindings
# add_definitions(-DOMPI_SKIP_MPICXX)

# include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# 添加头文件目录
include_directories("include")

# 添加编译选项
add_compile_options("-O3")

add_executable(dist src/distribution.cpp src/histo.cpp)
target_link_libraries(dist HepMC m fastjet)