cmake_minimum_required(VERSION 3.13)
project(Lab1 C)

set(CMAKE_C_STANDARD 99)

add_executable(Lab1 main.c def.h
        mtwister.h mtwister.c
        matrixutils.c matrixutils.h
        gaussjordan.c gaussjordan.h
        gauss.c gauss.h lup.c lup.h)
