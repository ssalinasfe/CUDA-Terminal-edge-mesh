cmake_minimum_required(VERSION 3.0.0)
project(DelaunayPolyGenerator VERSION 0.1.0)
project(cmake_cuda LANGUAGES CUDA)
set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch")


add_subdirectory(detri2)
add_executable(DelaunayPolyGenerator main.cu triang_edge.cu polygon.cu io.cu BET_elimination.cu)

set_target_properties(DelaunayPolyGenerator PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(DelaunayPolyGenerator detri2)



