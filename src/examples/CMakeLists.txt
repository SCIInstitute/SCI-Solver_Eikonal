project(Example1)

include_directories(
	"${CMAKE_SOURCE_DIR}/core/include"
)

CUDA_ADD_EXECUTABLE(Example1 example1.cu)
TARGET_LINK_LIBRARIES(Example1 Eikonal-3D_CORE)


add_custom_command(TARGET Example1 PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/../example_data ${CMAKE_BINARY_DIR}/example_data)
