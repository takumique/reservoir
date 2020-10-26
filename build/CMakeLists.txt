cmake_minimum_required(VERSION 3.10)

add_subdirectory(../app ${CMAKE_BINARY_DIR}/app)
add_subdirectory(../tool/data-gen ${CMAKE_BINARY_DIR}/tool/data-gen)
add_subdirectory(../tool/weights-gen ${CMAKE_BINARY_DIR}/tool/weights-gen)
