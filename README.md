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

 
