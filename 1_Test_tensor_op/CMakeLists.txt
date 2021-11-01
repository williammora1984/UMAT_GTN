cmake_minimum_required(VERSION 3.12)

project (PFUNIT_IMPLEMENTATION_UMAT_GTN_DAMAGE_MODEL_WITH_A_BACK_EULER_INTEGRATION
  VERSION 1.0.0
  LANGUAGES Fortran) #This line can not be changed

find_package(/home/william/pFUnit/build/PFUNIT REQUIRED)
enable_testing()

#FILE()
#FILE (COPY simple.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# system under test

add_library (sut
  tensor_ope_module.for 
  /usr/lib/x86_64-linux-gnu/liblapack.a #Name of the file to test
  )
target_include_directories(sut PUBLIC ${CMAKE_CURRENT_BINARY_DIR})


# tests
set (test_srcs test_tensor_mod.pf) #Name of the file .pf
add_pfunit_ctest (my_test_tensor_op1  #Name of the created executable
  TEST_SOURCES ${test_srcs}
  LINK_LIBRARIES sut
  )

#add_pfunit_ctest (my_test_tensor_llpack  #Name of the created executable
#  TEST_SOURCES ${test_srcs}
#   OTHER_SOURCES a.F90 b.F90 ...
#  a.F90 b.F90 ... 
#  LINK_LIBRARIES  /usr/lib/x86_64-linux-gnu/liblapack.a
#  )

#add_pfunit_ctest(my_test_tensor_llpack
#  TEST_SOURCES ${test_srcs}
#  a.F90 b.F90 ... 
#  LINK_LIBRARIES lapack
#  ) 
#  target_link_directories(my_test_tensor_llpack PRIVATE /usr/lib/x86_64-linux-gnu)
