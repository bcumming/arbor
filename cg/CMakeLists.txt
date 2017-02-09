set(CG_SOURCE
    cg.cpp
)

add_executable(cg ${CG_SOURCE})

target_link_libraries(cg LINK_PUBLIC nestmc gtest)
target_link_libraries(cg LINK_PUBLIC ${EXTERNAL_LIBRARIES})

set_target_properties(cg
   PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cg"
)

