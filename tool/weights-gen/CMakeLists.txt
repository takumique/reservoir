cmake_minimum_required(VERSION 3.10)
project(weights-gen)

add_executable(${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/../../app/reservoir.c
  ${PROJECT_SOURCE_DIR}/../../app/generic/mat.c
  ${PROJECT_SOURCE_DIR}/main.c
)

target_include_directories(${PROJECT_NAME} PUBLIC
  ${PROJECT_SOURCE_DIR}/../../app
  ${PROJECT_SOURCE_DIR}/../../app/generic
)

target_compile_options(${PROJECT_NAME} PUBLIC
)

target_compile_features(${PROJECT_NAME} PUBLIC
  c_std_99
)

target_compile_definitions(${PROJECT_NAME} PUBLIC
  PRECISION_F32
#  PRECISION_F64
)

target_link_libraries(${PROJECT_NAME}
  m
)
