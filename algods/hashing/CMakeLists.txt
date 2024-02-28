project(hashing)

add_executable(hashing
    test/test-hashmap.cpp
    test/test.cpp
    hashmap.h)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
