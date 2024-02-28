project(linked)

add_executable(linked
    test/test-bag.cpp
    test/test-list.cpp
    test/test-queue.cpp
    test/test-stack.cpp
    test/test.cpp
    bag.h
    doubly-node.h
    list.h
    queue.h
    singly-node.h
    stack.h)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
