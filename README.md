# UMAT_GTN

This is a routine writen in fortran77 to use as UMAT subroutine in Fortran.
 
22.09.2021
First commit UMAT GTN  includes

tensor_ope_module.for version 0.9.1 : Includes functions and routines to perform the tensor operations required in a material law, without finish the unit test 


6.10.2021
Implementation of elastic preditor

Updates

tesor_ope_module.for version 0.9.2

Implementation of elastic preditor contained in the file

material_law_GTN.for version 0.9.1

In order to test qualitatively the evaluation of the elastic step the following files have been added 

test_isochoric_tension.for version 0.9.1

figures_test_MatLaw.py version 0.9.1

Addtionally the UMAT subroutine, which call the material subroutine, has been upload.
 
UMAT_GTN_W.for version 0.9.1

To compile material routine and see the result of the qualitative test, use the following lines in the Comand line. This will generate some cvs and png files in the folder where the code files are saved.

gfortran tensor_ope_module.for material_law_GTN.for test_isochoric_tension.for -llapack

./a.out

figures_test_MatLaw.py

18.10.2021
Upload of unit tests and first version of input files to meake the UMAT integration test

The unit test were done using the package PFunit, to run the test without installation of PFunit use the executable file. See more information on the documentation and readme file of this folder. 

The Unit test folder /1_Test_tensor_op contains 
my_test_tensor_op1 :executable file

test_tensor_mod.pf: file to specify the testing

build_with_cmake_and_run.x

CMakeLists.txt

tensor_ope_module.for: module tested

generate_values_test.py: file that generate some of the referece values.

data_TT01.csv,...,data_TT09.csv: 9 files with data to test the different tensor operations
 
The integration test folder /2_Test_Abaqus_input_files contains that defines the geometries, type of elements, restriction steps to make the tension test and hidrostatic in single elements.  In this commit are uploaded the following files

test_ht_vol.inp, test_uni_tens_pla_strain.inp, test_uni_tens_pla_strain_rot90.inp, test_uni_tens_pla_stress.inp

19.10.2021 
Implementation of plastic return and plain stress

Update of the material routine with implementation of the plastic return, plain stress and algorithmic tangent stiffness. All the related files that had modification also are updated.

Uplaoad of the folder 3_UMAT/ This folder contains the files that must be compiled for Abaqus and therefore has some little differences with respect to the files that are compiled only with gfortran.