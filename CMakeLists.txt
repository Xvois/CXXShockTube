cmake_minimum_required(VERSION 3.10)
project(FluidSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source and include directories
set(SOURCES
    src/physics.cpp
    src/numerical.cpp
    src/analytical.cpp
    src/problems.cpp
    src/solver.cpp
    main.cpp
)

set(HEADERS
    include/constants.h
    include/physics.h
    include/numerical.h
    include/analytical.h
    include/problems.h
    include/solver.h
)

add_executable(FluidSolver ${SOURCES} ${HEADERS})
target_include_directories(FluidSolver PUBLIC include)
