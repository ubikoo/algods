project(searching)

add_executable(searching
    test/test-bstree.cpp
    test/test-kdtree.cpp
    test/test-kstree.cpp
    test/test-rbtree.cpp
    test/test.cpp
    test/bstree.h
    test/kdtree-dynamic.h
    test/kdtree-static.h
    test/rbtree.h
    bstree.h
    kdtree-dynamic.h
    kdtree-static.h
    rbtree.h)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
