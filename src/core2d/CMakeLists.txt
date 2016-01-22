include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${METIS_PATH}/Lib
  )

FILE(GLOB src cuda/*)
FILE(GLOB hdr include/*)

CUDA_ADD_LIBRARY(Eikonal-2D_CORE STATIC ${src} ${hdr})
ADD_DEPENDENCIES(Eikonal-2D_CORE metis cusp)
TARGET_LINK_LIBRARIES(Eikonal-2D_CORE ${METIS_LIBRARY})
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  TARGET_LINK_LIBRARIES(Eikonal-2D_CORE -ldl)
endif()
