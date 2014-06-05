
# add_subdirectory(colmajor)
# add_subdirectory(rowmajor)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_C_COMPILER g++)

#
# Use sse, sse2 in gcc
# Reference: http://stackoverflow.com/questions/5217812/c-compilation-issue-with-emmintrin-h-on-linux-gcc
#
add_definitions("-msse")
add_definitions("-msse2")

add_definitions("-D__MMX__")
add_definitions("-D__SSE__")
add_definitions("-D__SSE2__")
add_definitions("-D__SSE3__")

# set(SRC_FAST_MATMULT_LIST main.cpp fast_matmult.cpp cblas_dgemm.cpp aligned_malloc.cpp lang_def.cpp huge_tlb.cpp get_char.cpp ${SRC_FAST_MATMULT_COLMAJOR_LIST} ${SRC_FAST_MATMULT_ROWMAJOR_LIST})

file(GLOB_RECURSE SRC_FAST_MATMULT_LIST "${PROJECT_SOURCE_DIR}/src/fast_matmult/*.cpp")

# aux_source_directory(. SRC_FAST_MATMULT_LIST)

# aux_source_directory(./colmajor SRC_FAST_MATMULT_COLMAJOR_LIST)
# aux_source_directory(./rowmajor SRC_FAST_MATMULT_ROWMAJOR_LIST)

include_directories(${PROJECT_SOURCE_DIR}/include)
# add_library(util STATIC ${SRC_LIST})
# add_library(libhello SHARED ${LIB_SRC})
# link_directories(${PROJECT_SOURCE_DIR}/lib)

message(STATUS "This is BINARY dir " ${PROJECT_BINARY_DIR})
message(STATUS "This is SOURCE dir " ${PROJECT_SOURCE_DIR})

message(STATUS "This is SOURCE colmajor dir " ${SRC_FAST_MATMULT_COLMAJOR_LIST})
message(STATUS "This is SOURCE rowmajor dir " ${SRC_FAST_MATMULT_ROWMAJOR_LIST})

message(STATUS "This is SOURCE fast_matmult dir " ${SRC_FAST_MATMULT_LIST})

add_executable(fast_matmult ${SRC_FAST_MATMULT_LIST})
# add_executable(fast_matmult ${SRC_FAST_MATMULT_LIST} ${SRC_FAST_MATMULT_COLMAJOR_LIST} ${SRC_FAST_MATMULT_ROWMAJOR_LIST})

# target_link_libraries(hello util)
# set_target_properties(libhello PROPERTIES OUTPUT_NAME "hello")
