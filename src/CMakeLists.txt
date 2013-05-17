CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(Sequential_Experimental_Design)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

ADD_DEFINITIONS(-Wall -O3 -g -pg -llapack -lblas -lpthread -lgsl -lgslcblas -lm)

SET(SRC_SOURCE_FILES main.cpp fileRead.cpp inputParams.cpp
linearArch.cpp stochasticSearch.cpp structDef.cpp tools.cpp
trajectories.cpp userFunctions.cpp)
ADD_EXECUTABLE(main.x ${SRC_SOURCE_FILES})

# Required packages for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)

include_directories(${MPI_INCLUDE_PATH} ${MPI_INCLUDE_PATH} ${LAPACK_INCLUDE_PATH} ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
target_link_libraries(main.x ${MPI_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
