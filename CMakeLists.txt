cmake_minimum_required(VERSION 4.1)
project(FluidSolver)

set(CMAKE_CXX_STANDARD 20)

add_executable(FluidSolver main.cpp
        src/physics.cpp
        include/constants.h
        include/physics.h
        include/numerical.h
        include/analytical.h
        src/numerical.cpp
        src/analytical.cpp)
