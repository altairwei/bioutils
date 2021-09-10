aux_source_directory(. SEQ_FILES)
add_executable(seq ${SEQ_FILES})

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/lib
    ${CLI11_INCLUDE_DIRS}
)

link_directories(
    ${CMAKE_BINARY_DIR}/lib
)

target_link_libraries(seq bioutils ${CLI11_LIBRARIES})
set_target_properties(seq PROPERTIES OUTPUT_NAME "${EXE_PREFIX}seq")