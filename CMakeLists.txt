cmake_minimum_required(VERSION 2.6)
PROJECT(zombie_coldatoms)

find_package(MPI)
Find_Package(SDL)

#SET(CMAKE_CXX_FLAGS_RELEASE "-fopenmp -O3 -march=corei7 -mtune=corei7 -msse4 -msse4.1 -msse4.2 -mpopcnt -maes -mpclmul -mcx16 -msahf -pipe -fipa-pta -ftree-loop-linear -floop-interchange -floop-block -fgcse-after-reload -ftree-vectorize -fivopts -funswitch-loops -ftree-loop-im -ftree-loop-ivcanon -ftracer -fweb -std=c++0x")
#SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=core2 -mtune=core2 -msse4 -msse4.1 -pipe")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=corei7 -mtune=corei7 -Wconversion -pipe -mfpmath=sse -fopenmp -fgcse-after-reload -fpredictive-commoning -fipa-cp-clone -funswitch-loops -Wunsafe-loop-optimizations")
SET(CMAKE_C_FLAGS_RELEASE "-O2 -march=corei7 -mtune=corei7 -Wconversion -mfpmath=sse -fopenmp -Wunsafe-loop-optimizations")
SET(CMAKE_CXX_FLAGS_DEBUG "-march=corei7 -gstabs+ -Wall -Wabi -Wconversion")
SET(CMAKE_C_FLAGS_DEBUG "-march=core2 -O3 -mmmx -finline-functions -Wconversion")
SET(CMAKE_CXX_FLAGS_PROFILED "-O2 -march=corei7 -pg -Wconversion")
SET(CMAKE_C_FLAGS_PROFILED "-O3 -march=corei7 -pg -Wconversion")

include_directories(common/ MPI_Zombie_ColdAtoms SDL_Zombie_ColdAtoms)

Set(MPI_Zombie_src ./MPI_Zombie_ColdAtoms/zombie_coldatoms.cpp)

Set(SDL_Zombie_src ./SDL_Zombie_ColdAtoms/zombie_coldatoms.cpp ./SDL_Zombie_ColdAtoms/stripped_SDL_Net.c)

ADD_EXECUTABLE(MPI_zombie_coldatoms ${MPI_Zombie_src})

ADD_EXECUTABLE(SDL_zombie_coldatoms ${SDL_Zombie_src})

SET(MPI_Zombie_Defines "-DUSE_MPI -I.")

set(tempflags "")

#let's build the proper include file string from MPI_INCLUDE_PATH
foreach(elem ${MPI_INCLUDE_PATH})
set(var "-I${elem}")
set(tempflags "${tempflags} ${var}")
endforeach(elem)

set(MPI_files_flags "-DUSE_MPI ${tempflags}")

#let's retrieve the mpi library string...

set(templinkflags "")
foreach(elem ${MPI_EXTRA_LIBRARY})
get_filename_component(fln ${elem} NAME_WE)#strip path
string(REPLACE "lib" "" fln ${fln}) # strip the lib prefix of the libraries
set(templinkflags "${templinkflags} -l${fln}")
endforeach(elem)


# set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
# include_directories(MPI_INCLUDE_PATH)
# target_link_libraries(mytest ${MPI_LIBRARIES}) 

IF(CMAKE_COMPILER_IS_GNUCXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

set_source_files_properties(${MPI_Zombie_src} PROPERTIES COMPILE_FLAGS "${MPI_files_flags}")

set_target_properties(MPI_zombie_coldatoms PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/MPI_Zombie_ColdAtoms COMPILE_FLAGS "${MPI_Zombie_Defines}" LINK_FLAGS "${MPI_LINK_FLAGS} ${templinkflags}")

set_target_properties(SDL_zombie_coldatoms PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SDL_Zombie_ColdAtoms COMPILE_FLAGS " ${COMPILE_FLAGS} -DUSE_SDL -I." )

target_link_libraries(SDL_zombie_coldatoms ${SDL_LIBRARY})
target_link_libraries(MPI_zombie_coldatoms ${MPI_LIBRARIES}) 