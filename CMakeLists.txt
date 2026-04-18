cmake_minimum_required(VERSION 4.1)
project(FluidSolver)

set(CMAKE_CXX_STANDARD 20)

add_executable(FluidSolver main.cpp
        src/physics.cpp
        src/numerical.cpp
        src/analytical.cpp
        src/problems.cpp
        include/constants.h
        include/physics.h
        include/numerical.h
        include/analytical.h
        include/problems.h)
