#CMAKE SETTING CUDA_NVCC_FLAGS
#############################
# Check for GPUs present and their compute capability
# based on http://stackoverflow.com/questions/2285185/easiest-way-to-test-for-existence-of-cuda-capable-gpu-from-cmake/2297877#2297877 (Christopher Bruns)
if(CUDA_FOUND AND NOT CUDA_HAVE_GPU)
    message(STATUS "Checking CUDA compute capability.")
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR} 
        ${CMAKE_SOURCE_DIR}/cuda_compute_capability.c
        CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${CUDA_TOOLKIT_INCLUDE}
        -DLINK_LIBRARIES:STRING=${CUDA_CUDART_LIBRARY}
        COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE RUN_OUTPUT_VAR)
    # COMPILE_RESULT_VAR is TRUE when compile succeeds
    # RUN_RESULT_VAR is zero when a GPU is found
    if(COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR)
        set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
        set(CUDA_COMPUTE_CAPABILITY ${RUN_OUTPUT_VAR} CACHE STRING "Compute capability of CUDA-capable GPU present")
        set(CUDA_GENERATE_CODE "arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY}" CACHE STRING "Which GPU architectures to generate code for (each arch/code pair will be passed as --generate-code option to nvcc, separate multiple pairs by ;)")
        mark_as_advanced(CUDA_COMPUTE_CAPABILITY CUDA_GENERATE_CODE)
        LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21" "-arch=sm_30" "-arch=sm_31")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-arch=sm_${CUDA_COMPUTE_CAPABILITY}"  CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
	message(STATUS "Automatic Compute Capability Set: ${CUDA_COMPUTE_CAPABILITY}")
    else()
        message(STATUS "You must set your compute capability manually")
        set(CUDA_HAVE_GPU FALSE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    endif()
endif()

IF (NOT CUDA_HAVE_GPU)
    SET(CUDA_USE_COMPUTE_CAPABILITY_10 false CACHE BOOL "Use Compute Capability 1.0")
    SET(CUDA_USE_COMPUTE_CAPABILITY_11 false CACHE BOOL "Use Compute Capability 1.1")
    SET(CUDA_USE_COMPUTE_CAPABILITY_12 false CACHE BOOL "Use Compute Capability 1.2")  
    SET(CUDA_USE_COMPUTE_CAPABILITY_13 false CACHE BOOL "Use Compute Capability 1.3")
    SET(CUDA_USE_COMPUTE_CAPABILITY_20 true CACHE BOOL "Use Compute Capability 2.0")
    SET(CUDA_USE_COMPUTE_CAPABILITY_21 false CACHE BOOL "Use Compute Capability 2.1")
    
    IF(CUDA_USE_COMPUTE_CAPABILITY_10)
      LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21")
      LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_10")
      SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
      set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    ELSEIF(CUDA_USE_COMPUTE_CAPABILITY_11)
      LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21")
      LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_11")
      SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
      set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    ELSEIF(CUDA_USE_COMPUTE_CAPABILITY_12)
      LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21")
      LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_12")
      SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
      set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    ELSEIF(CUDA_USE_COMPUTE_CAPABILITY_13)
      LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21")
      LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_13")
      SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
      set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    ELSEIF(CUDA_USE_COMPUTE_CAPABILITY_20)
      LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21")
      LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_20")
      SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
      set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    ELSEIF(CUDA_USE_COMPUTE_CAPABILITY_21)
      LIST(REMOVE_ITEM CUDA_NVCC_FLAGS "-arch=sm_10" "-arch=sm_11" "-arch=sm_12" "-arch=sm_13" "-arch=sm_20" "-arch=sm_21")
      LIST(APPEND CUDA_NVCC_FLAGS "-arch=sm_21")
      SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} CACHE LIST "Semi-colon delimit multiple arguments." FORCE)
      set(CUDA_HAVE_GPU TRUE CACHE INTERNAL BOOL "Whether CUDA-capable GPU is present")
    ENDIF()
ENDIF ()

IF(NOT CUDA_HAVE_GPU)
	message( FATAL_ERROR "You must set your compute capability" )
ENDIF()
#END CMAKE SETTING CUDA_NVCC_FLAGS
###################

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${METIS_PATH}/Lib
)

set(CORE_CUDA_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/cuda/CUDA_kernels.cu"
  "${CMAKE_CURRENT_SOURCE_DIR}/cuda/meshFIM.cu"
	"${CMAKE_CURRENT_SOURCE_DIR}/cuda/predicates.cu"
	"${CMAKE_CURRENT_SOURCE_DIR}/cuda/tetgen.cu"
	"${CMAKE_CURRENT_SOURCE_DIR}/cuda/tetmesh.cu"
)

set(CORE_CUDA_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/CUDADefines.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/meshFIM.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/mycutil.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/tetgen.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/tetmesh.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/Vec.h"
)

CUDA_ADD_LIBRARY(Eikonal-3D_CORE STATIC ${CORE_CUDA_SOURCES} ${CORE_CUDA_HEADERS})
ADD_DEPENDENCIES(Eikonal-3D_CORE metis)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  TARGET_LINK_LIBRARIES(Eikonal-3D_CORE -ldl ${METIS_LIBRARY})
else()
  TARGET_LINK_LIBRARIES(Eikonal-3D_CORE ${METIS_LIBRARY})
endif()
