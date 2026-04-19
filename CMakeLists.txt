cmake_minimum_required(VERSION 3.10)
project(FluidSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP support (built into GCC, needs libomp on macOS)
option(USE_OPENMP "Enable OpenMP parallelisation" ON)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
        set(OPENMP_AVAILABLE ON)
    else()
        message(STATUS "OpenMP not found (macOS clang does not support -fopenmp natively)")
        message(STATUS "Install via: brew install llvm && brew install libomp")
        message(STATUS "Falling back to sequential execution (no parallel overhead)")
        set(OPENMP_AVAILABLE OFF)
    endif()
else()
    set(OPENMP_AVAILABLE OFF)
endif()

add_definitions(-DOPENMP_AVAILABLE=${OPENMP_AVAILABLE})

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
if(OPENMP_AVAILABLE)
    target_link_libraries(FluidSolver PRIVATE OpenMP::OpenMP_CXX)
endif()