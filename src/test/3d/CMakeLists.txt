########################################################################

# Definitions for Tests
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/bin)
add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../test_data/")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../core2d/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../core3d/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../..)
#test macro
function(NEWTEST name)
  add_executable(${name} ${name}.cc)
  if (NOT WIN32)
    target_link_libraries(${name} ${GTEST_LIBRARY} EIKONAL_CORE -lpthread)
  else ()
    target_link_libraries(${name} ${GTEST_LIBRARY} EIKONAL_CORE )
  endif()
  add_test(${name} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name}) 
  add_dependencies(${name} gtest EIKONAL_CORE)
endfunction()
#the tests
newtest(3D-sanity)
newtest(3D-small)
newtest(3D-basic)
newtest(3D-converge)
newtest(3D-improve)
