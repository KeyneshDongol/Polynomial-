

set(EIGEN_DIR /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

cmake_minimum_required(VERSION 3.10)

project(Polynomial)

include_directories(${EIGEN_DIR})

add_executable(Polynomial main.cpp)
#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(Polynomial PUBLIC OpenMP::OpenMP_CXX)
#endif()

#target_link_libraries(Polynomial /opt/homebrew/Cellar/llvm/16.0.6/lib/ -lomp)
