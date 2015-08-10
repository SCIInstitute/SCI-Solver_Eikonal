########################################################################

# Definitions for Tests
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/bin)
add_definitions(-DTEST_BINARY_DIR="${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cli")
add_definitions(-DTEST_DATA_DIR="${CMAKE_SOURCE_DIR}/test/test_data/")
include_directories(${CMAKE_SOURCE_DIR}/core2d/include)
#test macro
function(NEWTEST name)
  add_executable(${name} ${name}.cc)
  if (NOT WIN32)
    target_link_libraries(${name} ${GTEST_LIBRARY} Eikonal-2D_CORE -lpthread)
  else ()
    target_link_libraries(${name} ${GTEST_LIBRARY} Eikonal-2D_CORE )
  endif()
  add_test(${name} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name}) 
  add_dependencies(${name} gtest Eikonal-2D_CORE)
endfunction()
#the tests
newtest(2D-sanity)
newtest(2D-basic)
newtest(2D-converge)
newtest(2D-improve)
