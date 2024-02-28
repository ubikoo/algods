project(arrays)

add_executable(arrays
    test/test-array.cpp
    test/test-deque.cpp
    test/test.cpp
    array.h
    deque.h)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
