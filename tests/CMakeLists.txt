include(GoogleTest)

find_package(GTest REQUIRED)
find_package(benchmark REQUIRED)

macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${TESTNAME} bioutils ${GTest_LIBRARIES})
    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

macro(package_add_bench TESTNAME)
    set(TARGET_NAME ${TESTNAME})
    add_executable(${TARGET_NAME} ${ARGN})
    set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME ${TESTNAME})
    set_target_properties(${TARGET_NAME} PROPERTIES FOLDER tests)
    target_link_libraries(${TARGET_NAME} bioutils ${benchmark_LIBRARIES})
    #add_test(NAME ${TARGET_NAME} COMMAND ${TARGET_NAME})
endmacro()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${GTest_INCLUDE_DIRS}
    ${benchmark_INCLUDE_DIRS}
)

package_add_test(TestPattern test-pattern.cpp)
package_add_bench(BenchPattern bench-pattern.cpp)

