project(numeric)

add_executable(numeric
    test/test-array1d.cpp
    test/test-array2d.cpp
    test/test-array3d.cpp
    test/test-histogram1d.cpp
    test/test-histogram2d.cpp
    test/test-histogram3d.cpp
    test/test-linalg.cpp
    test/test-ode-kepler.cpp
    test/test-ode-spring.cpp
    test/test.cpp
    test/test-array1d.h
    test/test-array2d.h
    test/test-array3d.h
    test/test-linalg.h
    array.h
    array1d.h
    array2d.h
    array3d.h
    histogram.h
    histogram1d.h
    histogram2d.h
    histogram3d.h
    linalg-arithmetic.h
    linalg-gauss.h
    linalg-jacobi.h
    linalg.h
    ode-euler.h
    ode-gauss.h
    ode-gear.h
    ode-kutta.h
    ode.h
    sde-ito.h
    sde-strat.h
    sde-wiener.h
    sde.h)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})

# target_link_libraries(${PROJECT_NAME} PRIVATE coremath)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
