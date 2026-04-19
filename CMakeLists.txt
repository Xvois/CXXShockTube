cmake_minimum_required(VERSION 3.10)
project(FluidSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Parallel execution option (requires GCC >= 11 or Clang >= 14 for std::execution)
option(USE_PARALLEL "Enable C++20 parallel execution" ON)

if(USE_PARALLEL)
    # GCC: -fopenmp-simd enables parallel execution without OpenMP runtime dependency
    # Clang: uses libc++ parallel algorithms
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # Clang uses libc++ parallelism; may need -fopenmp-simd for full support
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd")
    endif()
endif()

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

# Link TBB for std::execution::par_unseq support
find_package(TBB REQUIRED)
target_link_libraries(FluidSolver PRIVATE tbb)
