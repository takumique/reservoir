cmake_minimum_required(VERSION 3.10)
project(data-gen)

add_executable(${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/main.c
)

target_include_directories(${PROJECT_NAME} PUBLIC
)

target_compile_options(${PROJECT_NAME} PUBLIC
)

target_compile_features(${PROJECT_NAME} PUBLIC
  c_std_99
)

target_compile_definitions(${PROJECT_NAME} PUBLIC
)

target_link_libraries(${PROJECT_NAME}
)
