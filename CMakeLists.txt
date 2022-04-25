cmake_minimum_required(VERSION 3.10)
project(MaxInf_BRGSTkNN_2020)

set(CMAKE_CXX_STANDARD 11)



# 修改clion 内存
MATH(EXPR stack_size "500*1024*1024")
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-stack,${stack_size}")


FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB TENINDRA_SRC
        "./Tenindra/command/*.h"
        "./Tenindra/command/*.cpp"
        "./Tenindra/utility/*.h"
        "./Tenindra/utility/*.cpp"
        "./Tenindra/queue/*.h"
        "./Tenindra/queue/*.cpp"
        "./Tenindra/processing/*.h"
        "./Tenindra/processing/*.cpp"
        "./Tenindra/tuple/*.h"
        "./Tenindra/tuple/*.cpp"
        "./Tenindra/external/*.h"
        "./Tenindra/external/*.cpp"
        "./Tenindra/external/command/*.h"
        "./Tenindra/external/command/*.cpp"
        "./Tenindra/external/datastr/*.h"
        "./Tenindra/external/datastr/*.cpp"
        "./Tenindra/external/datastr/graph/*.h"
        "./Tenindra/external/datastr/graph/*.cpp"
        "./Tenindra/external/datastr/pqueue/*.h"
        "./Tenindra/external/datastr/pqueue/*.cpp"
        "./Tenindra/external/io/*.h"
        "./Tenindra/external/io/*.cpp"
        "./Tenindra/external/processing/*.h"
        "./Tenindra/external/processing/*.cpp"
        "./Tenindra/external/stats/*.h"
        "./Tenindra/external/stats/*.cpp"
        )

add_executable(MaxInf_BRGSTkNN_2020   process.h
        ./SFMT/dSFMT/dSFMT.c
        InfluenceModels.cc
        ./SocialInfluenceModule/CELF_PLUS/common.cc
        ./SocialInfluenceModule/CELF_PLUS/MC.cc
        ./SocialInfluenceModule/CELF_PLUS/anyoption.cc
        ./SocialInfluenceModule/PMC/pmc.cpp
        imm.h OPIMC.h
        main.cpp
        ${TENINDRA_SRC})

target_link_libraries(MaxInf_BRGSTkNN_2020 boost_serialization /usr/local/lib/libmetis.a)