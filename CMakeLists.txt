cmake_minimum_required(VERSION 3.14)
project(elda VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Main Library Target
add_library(elda
    src/matrix.cpp
    src/transforms.cpp
    src/vector_utils.cpp
)

# Specify include directories
target_include_directories(elda PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Main Executable Target
add_executable(main main.cpp)
target_link_libraries(main PRIVATE elda)
