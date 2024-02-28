project(empty)

add_executable(empty empty.cpp)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
