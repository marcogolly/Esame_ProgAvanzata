cmake_minimum_required(VERSION 3.16)

project(ESAME_MATRIX_CLASS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# cerco la libreria OpenCL e gli header corrispondenti
find_package(OpenCL REQUIRED)  # per usare <CL/cl.h>

# includo gli header che verranno usati per la compilazione
INCLUDE_DIRECTORIES(${OpenCL_INCLUDE_DIRS})

LINK_DIRECTORIES(${OpenCL_LIBRARY})


add_executable(matrici matrici.cpp)


add_executable(main main.cpp)
# linko la libreria all'eseguibile che sto compilando
target_link_libraries(main ${OpenCL_LIBRARY})