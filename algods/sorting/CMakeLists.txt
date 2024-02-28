project(sorting)

add_executable(sorting
    test/test-priority-queue.cpp
    test/test-sort.cpp
    test/test.cpp
    test/test-insertsort.h
    test/test-mergesort.h
    test/test-priority-queue.h
    test/test-quicksort.h
    test/test-selectsort.h
    test/test-shellsort.h
    test/test-sort.h
    insertsort.h
    is-ordered.h
    merge.h
    mergesort.h
    partition.h
    priority-queue.h
    quicksort.h
    selectsort.h
    shellsort.h
    sort.h)

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
