aux_source_directory(. FIND_FILES)
add_executable(find ${FIND_FILES})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${argparse_INCLUDE_DIRS}
    ${CLI11_INCLUDE_DIRS}
)

link_directories(
    ${CMAKE_BINARY_DIR}/lib
)

target_link_libraries(find bioutils ${argparse_LIBRARIES} ${CLI11_LIBRARIES})
set_target_properties(find PROPERTIES OUTPUT_NAME "${EXE_PREFIX}find")