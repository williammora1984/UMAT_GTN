!@test
subroutine test_tensor_operations()
   use tensor_operations   
   use funit
   
   implicit none
   
   double precision, dimension(2):: vector2
   integer:: i, j,stat     !iterators
   double precision:: a, b, c     !constants
   integer:: n1, n2  !dimension matrix
   double precision, dimension(:,:), allocatable:: m1, m1_T !2D test tensor alloatable
   double precision, dimension(3,3):: m2,m2p,m2pp, T_2nd_3_1, T_2nd_3_2, Re_2nd_3 !2D test tensor 3x3
   double precision, dimension(6,6):: T_2nd_6, T_2nd_6p
   double precision, dimension(9,9)::T_2nd_9_1, Re_2nd_9,Re_2nd_9_2
   double precision, dimension(3,3,3,3)::T4th, T4thp, T4thpp, T4thppp  
   double precision, dimension(3,3):: I_3D, Inot
   character(len=13) :: file_name
   double precision, dimension(3):: vector3
   double precision, dimension(5):: vector5
   double precision, dimension(12):: vector12
   double precision :: Re_0th
   double precision, dimension(3,3):: import1
   double precision, dimension(11,38):: import2
   double precision, dimension(2,162):: import3
   double precision, dimension(2,19):: import4
   double precision, dimension(1,99):: import5
   double precision, dimension(3,99):: import6
   double precision, dimension(5,117):: import7
   double precision, dimension(7,117):: import8

   a=0.0
   b=2.0
   c=3.0
!-----------------------------------------------------------------------   
!                         test of linespace
!-----------------------------------------------------------------------  
   call linespace (a, b, vector2,0)
#line 39 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual((/1.0, 2.0/), vector2, 'linespace(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   call linespace (a, b, vector3,1)
#line 41 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual((/0.0, 1.0, 2.0/), vector3, 'linespace(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   call linespace (a, c, vector12,0)
#line 43 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual((/0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00/), vector12, 'linespace(2)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   call linespace (a, b, vector2,0)
!  @assertEqual((/1.0, 2.0/), vector2, 'Linearspace(1)')

!-----------------------------------------------------------------------
!                        test of Trace
!-----------------------------------------------------------------------
!Import the data for the test of Trace, diadic product, inverse.
   n1=size(import2,1)
   n2=size(import2,2)
   file_name="data_TT02.csv"
   call read2 (import2, n1,n2, file_name)

   do i=1,n1
     print*, "Test trace", i
     m2=reshape ((import2(i,1:9)),shape(m2), order=(/2,1/))
#line 59 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(import2(i,19),Trace (m2), 'Trace(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 59) )
  if (anyExceptions()) return
#line 60 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do 

!----------------------------------------------------------------------- 
!                      test of Ident1
!-----------------------------------------------------------------------
   n1=2; n2=2
   allocate (m1(n1,n2)); allocate (m1_T(n1,n2))
   call Ident1(m1,n1)
   m1_T=reshape((/1.0, 0.0, 0.0, 1.0/), shape(m1_T), order=(/2,1/) )
#line 69 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(m1_T, m1, 'Ident1(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   deallocate (m1); deallocate (m1_T)

   n1=3; n2=3
   allocate (m1(n1,n2)); allocate (m1_T(n1,n2))
   call Ident1(m1,n1)
   m1_T=reshape((/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /), shape(m1_T), order=(/2,1/) )
#line 76 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(m1_T, m1, 'Ident1(2)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   deallocate (m1); deallocate (m1_T)
   
   do i=1,2
   n1=6; n2=6
   allocate (m1(n1,n2)); allocate (m1_T(n1,n2))
   call Ident1(m1,n1)
   m1_T=reshape((/1.0, 0.0, 0.0, 0.0, 0.0, 0.0,&  
                  0.0, 1.0, 0.0, 0.0, 0.0, 0.0,&  
                  0.0, 0.0, 1.0, 0.0, 0.0, 0.0,& 
                  0.0, 0.0, 0.0, 1.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0, 0.0, 1.0, 0.0,&
                  0.0, 0.0, 0.0, 0.0, 0.0, 1.0/),&
                  shape(m1_T), order=(/2,1/) )
#line 90 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(m1_T, m1, 'Ident1(3)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   deallocate (m1); deallocate (m1_T)
   end do

!-----------------------------------------------------------------------
!                    test of inv_T_2D
!-----------------------------------------------------------------------
   n1=size(import2,1)
   do i=1,n1
      print*, "Test inv_T_2D", i
      m2=reshape ((import2(i,1:9)),shape(m2), order=(/2,1/))
      m2p=reshape ((import2(i,21:29)),shape(m2), order=(/2,1/))
!     @assertEqual(m2p,inv_T_2D (m2), 'Inverse(1)')
   end do

!-----------------------------------------------------------------------
!                         test of I4sym
!-----------------------------------------------------------------------
!Import the data for the test of I4sym, P4sym and their inverses
!data_TT03.csv=Test_4th_In_1.csv
   n1=size(import3,1)
   n2=size(import3,2)
   file_name="data_TT03.csv"
   call read2 (import3, n1,n2, file_name)
   call I4sym(T4thp)  

   print*, "Test I4sym"
   T4th=reshape ((import3(1,1:81)),shape(T4th), order=(/4,3,2,1/))
#line 118 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(T4th,T4thp, 'I4sym(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
!-----------------------------------------------------------------------
!                        test of P4sym
!-----------------------------------------------------------------------
   call P4sym(T4thp)  
   !T4th=reshape ((import3(2,1:81)),shape(T4th), order=(/4,3,2,1/))
   T4th=reshape ((import3(2,1:81)),shape(T4th), order=(/1,2,3,4/))
   print*, "Test P4sym"
#line 126 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(T4th,T4thp,'P4sym(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"

!-----------------------------------------------------------------------
!                         test of invT4
!-----------------------------------------------------------------------
   n1=size(import3,1)
   do i=1,n1
      print*, "Test inv_T4", i
      T4th=reshape ((import3(i,1:81)),shape(T4th), order=(/4,3,2,1/))
      T4thp=reshape ((import3(i,82:162)),shape(T4thp), order=(/4,3,2,1/))
!      call invT4(T4th,T4thpp)
!      @assertEqual(T4thp,T4thppp, 'invT4(1)')
   end do

!-----------------------------------------------------------------------
!                test of contrac_2nd_2nd
!-----------------------------------------------------------------------
   n1=size(import2,1)
   do i=1,n1
      print*, "Test contrac_2nd_2nd", i
      m2=reshape ((import2(i,1:9)),shape(m2), order=(/2,1/))
      m2p=reshape ((import2(i,10:18)),shape(m2p), order=(/2,1/))
#line 148 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(import2(i,20),contrac_2nd_2nd(m2,m2p),tolerance=1.e-10, &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do

!Second data set for test
   !data_TT04.csv=Test_Cont_ML.csv !ML =Matlab
   n1=size(import4,1)
   n2=size(import4,2)
   file_name="data_TT04.csv"
   call read2 (import4, n1,n2, file_name)
   do i=1,n1
      print*, "Test contrac_2nd_2nd", i
      m2=reshape ((import4(i,1:9)),shape(m2), order=(/2,1/))
      m2p=reshape ((import4(i,10:18)),shape(m2p), order=(/2,1/))
#line 161 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(import4(i,19),contrac_2nd_2nd(m2,m2p), 'contrac_2nd_2nd(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do
!-----------------------------------------------------------------------
!             test of diadic_prod_T2_T2
!-----------------------------------------------------------------------
!Import the data for test of diadic_prod_T2_T2
!data_TT05.csv=Test_Diad_ML.csv !ML =Matlab
   n1=size(import5,1)
   n2=size(import5,2)
   file_name="data_TT05.csv"
   call read2 (import5, n1,n2, file_name)
   do i=1,n1
      print*, "Test diadic_prod_2nd_2nd", i
      m2=reshape ((import5(i,1:9)),shape(m2), order=(/2,1/))
      m2p=reshape ((import5(i,10:18)),shape(m2p), order=(/2,1/))
      T4th=reshape ((import5(i,19:99)),shape(T4th), order=(/4,3,2,1/))
#line 177 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(T4th,diadic_prod_T2_T2(m2,m2p), 'diad(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do
!-----------------------------------------------------------------------
!                  test of contrac_4th_2nd
!-----------------------------------------------------------------------
!data_TT06.csv=Test_Cont_4_2.csv
   n1=size(import6,1)
   n2=size(import6,2)
   file_name="data_TT06.csv"
   call read2 (import6, n1,n2, file_name)
   do i=1,n1
      print*, "Test contrac_4th_2nd", i
      T4th=reshape ((import6(i,1:81)),shape(T4th), order=(/1,2,3,4/))
      m2=reshape ((import6(i,82:90)),shape(m2), order=(/1,2/))
      m2p=reshape ((import6(i,91:99)),shape(m2p), order=(/1,2/))
      m2pp=contrac_4th_2nd(T4th,m2)
#line 193 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(m2p, m2pp,'contrac_4th_2nd(1)', &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
      call emp_2D(m2pp) 
   end do
!-----------------------------------------------------------------------
!                  test of voig_2_T4th
!-----------------------------------------------------------------------
!data_TT07.csv=Test_kmt1.csv   !Matrix with integers
   n1=size(import7,1)
   n2=size(import7,2)
   file_name="data_TT07.csv"
   call read2 (import7, n1,n2, file_name)
   do i=1,n1
      print*, "Test voig_2_T4th", i
      T_2nd_6=reshape ((import7(i,1:36)),shape(T_2nd_6), order=(/1,2/))
      call voig_2_T4th(T_2nd_6,T4th)
      T4thp=reshape ((import7(i,37:117)),shape(T4th), order=(/1,2,3,4/))
#line 209 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(T4thp,T4th,tolerance=0.5e-6, &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do

!data_TT08.csv=Test_kmt2.csv   !Matrix with reals and representation of 4th order tensors
   n1=size(import8,1)
   n2=size(import8,2)
   file_name="data_TT08.csv"
   call read2 (import8, n1,n2, file_name)
   do i=1,n1
      print*, "test_voig_2_T4th_2", i 
      T_2nd_6=reshape ((import8(i,1:36)),shape(T_2nd_6), order=(/1,2/))
      call voig_2_T4th(T_2nd_6,T4th)
      T4thp=reshape ((import8(i,37:117)),shape(T4th), order=(/1,2,3,4/))
#line 222 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(T4thp,T4th,tolerance=0.5e-6, &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do 
!-----------------------------------------------------------------------
!               test of T4th_2_Voig
!-----------------------------------------------------------------------
   n1=size(import7,1)
   do i=1,n1
      print*, "Test T4th_2_voig", i
      T_2nd_6=reshape ((import7(i,1:36)),shape(T_2nd_6), order=(/1,2/))
      T4th=reshape ((import7(i,37:117)),shape(T4th), order=(/1,2,3,4/))
      call T4th_2_Voig(T4th, T_2nd_6p)
#line 233 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
  call assertEqual(T_2nd_6,T_2nd_6p,tolerance=0.6e-6, &
 & location=SourceLocation( &
 & 'test_tensor_mod.pf', &
 & 233) )
  if (anyExceptions()) return
#line 234 "/home/william/Code_Sem_IV/PPP_II/Code_PPP/1_Test_tensor_op/test_tensor_mod.pf"
   end do




end subroutine test_tensor_operations



module Wraptest_tensor_mod
   use FUnit
   implicit none
   private

contains


end module Wraptest_tensor_mod

function test_tensor_mod_suite() result(suite)
   use FUnit
   use Wraptest_tensor_mod
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   external test_tensor_operations


   suite = TestSuite('test_tensor_mod_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_tensor_operations', &
      test_tensor_operations))
   call suite%addTest(t)


end function test_tensor_mod_suite

