cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17) 

project(interval_branch_and_bound)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE RELEASE)
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-fmax-errors=1 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(BUILD_GLO_OPT ON CACHE BOOL "whether to build dco boost interval examples")

add_subdirectory(BnBCore)
add_subdirectory(global_optimization)
