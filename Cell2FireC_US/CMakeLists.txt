cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(fuelmodelBP)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost COMPONENTS filesystem REQUIRED) # Needed for M1 arm64 build
include_directories(${Boost_INCLUDE_DIRS})

add_compile_options(-std=c++11)


##### TORCH
set(CMAKE_PREFIX_PATH /home/minho/fires/cell2fireML/libtorch ..)
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

set(SOURCES Cell2Fire.cpp CellsFBP.cpp SpottingFBP.cpp ReadCSV.cpp ReadArgs.cpp WriteCSV.cpp Ellipse.cpp fuelmodelBP.cpp)
add_executable(Cell2Fire ${SOURCES})

target_link_libraries(Cell2Fire "${TORCH_LIBRARIES}")
set_property(TARGET Cell2Fire PROPERTY CXX_STANDARD 14)

##### Link libraries to executable
set(LIBS -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -lm -lpthread -ldl)
target_link_libraries(Cell2Fire ${LIBS})