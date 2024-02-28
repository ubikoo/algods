project(unionfind)

add_executable(unionfind
    test/test-unionfind.cpp
    test/test.cpp
    test/test-unionfind.h
    unionfind-templ.h
    unionfind.h)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
