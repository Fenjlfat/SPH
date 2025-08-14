cmake_minimum_required(VERSION 3.10)
project(SPH)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Явное указание всех исходников (лучше чем GLOB)
set(SOURCES 
    source/SPH.cpp
    source/SPOUT.cpp
    source/MOVE.cpp
    source/INIT.cpp
    source/EOS.cpp
    source/NLIST.cpp
    source/DISLOC.cpp
    # другие .cpp файлы
)

# Создание исполняемого файла
add_executable(SPH ${SOURCES})

# Жёсткое указание абсолютного пути к заголовкам
target_include_directories(SPH PRIVATE 
    ${CMAKE_SOURCE_DIR}/headers
    # другие пути при необходимости
)
set(CMAKE_BUILD_TYPE Debug)