      program test_tensor
        use tensor_operations
        implicit none
        
        integer:: i, j     !iterators
        double precision:: a, b, c     !constants
        integer:: n1, n2  !dimension matrix
        double precision, dimension(:,:), allocatable:: m1 !2D test tensor alloatable
        double precision, dimension(3,3):: m2, T_2nd_3_1, T_2nd_3_2, 
     $                                      Re_2nd_3 !2D test tensor 3x3
        double precision, dimension(6,6):: T_2nd_6_1, T_2nd_6_2 !2D test tensor 6x6
        double precision, dimension(9,9)::T_2nd_9_1, Re_2nd_9,Re_2nd_9_2
        double precision, dimension(3,3,3,3)::I4sym_sub, P4sym_o
     $                                       ,P4sym_Inv , Re_4th_3
        double precision, dimension(3,3):: I_3D, Inot
        double precision, dimension(5):: vector5
        double precision, dimension(12):: vector12
        double precision :: Re_0th
c Variables loading
        integer:: ltype, type_linespace
        integer,dimension(5) :: posi_last_time =(/2,3,3,4,7/)
        integer::n_step
        double precision :: dt,mu, n_ampl
        double precision, dimension(7) :: t
        double precision, dimension(3) ::lam
        double precision, dimension(:), allocatable:: time
        double precision, dimension(:), allocatable:: e11 !11 strain

        n1=3
        n2=3
        
        a=0.0
        b=5.0
        c=3.0
        call linespace (a, b, vector5, 1)
        print*, "linespace", vector5
        call linespace (a, c, vector12, 1)
        print*, "linespace", vector12
        allocate (m1(n1,n2))
        !arraySize = size(m1)
        call Ident1(M1,n1)
        !call print_matrix(M1,n1)
        m2=m1
        deallocate (m1)
        m2=m2+m2
        !call print_matrix(M2,n1)
        
        !call print_matrix(M2,n1)
        call I4sym(I4sym_sub)
        !call print_4D_tensor(I4sym_sub,3)
        !print*,I4sym_sub
        call P4sym(P4sym_o)
        !call print_4D_tensor(P4sym,3)
        !print*, P4sym

        call invT4(P4sym_o,P4sym_Inv) !Using a declaration as routine
        !P4sym_Inv= invT4(P4sym) !Using a declaration as function
        print*, "Inverse of I4sym"
        call print_4D_tensor(P4sym_Inv,3)
        a=1.00
        b=0.0
        c=2.0

        T_2nd_3_1=reshape((
     &        / a,   b,   b,
     &          b,   a,   b,
     &          b,   b,   a /),
     &      shape(T_2nd_3_1), order=(/2,1/) )
       print*, "inverse matrix 1"
       Re_2nd_3=inv_T_2D(T_2nd_3_1)
       call print_matrix(Re_2nd_3,3,3)   

        a=1.0
        b=2.0
        c=3.0

        T_2nd_3_1=reshape((
     &        / a,   b,   c,
     &          c,   b,   a,
     &          b,   a,   c /),
     &      shape(T_2nd_3_1), order=(/2,1/) )
        print*, "inverse matrix 2"
        Re_2nd_3=inv_T_2D(T_2nd_3_1)
        call print_matrix(Re_2nd_3,3,3)

        a=1.0
        b=2.0
        c=3.0 
        T_2nd_3_2=reshape((
     &        / a,   b,   c,
     &          b,   a,   b,
     &          c,   b,   a /),
     &      shape(T_2nd_3_1), order=(/2,1/) )
  
        a=2.0
        b=0.0
        c=2**0.5
  
        T_2nd_6_1=reshape((
     &        / a,   a,   a, c, c, c,
     &          a,   a,   a, c, c, c,
     &          a,   a,   a, c, c, c,
     &          c,   c,   c, a, a, a,
     &          c,   c,   c, a, a, a,
     &          c,   c,   c, a, a, a /),
     &      shape(T_2nd_6_1), order=(/2,1/) )
  
        a=2.0
        b=0.0
        c=2**0.5
        T_2nd_6_2=reshape((
     &        / a,   a,   a, c, c, c,
     &          b,   a,   b, c, b, c,
     &          a,   b,   a, b, c, b,
     &          b,   c,   b, a, b, a,
     &          c,   b,   c, b, a, b,
     &          c,   c,   c, a, a, a /),
     &      shape(T_2nd_6_1), order=(/2,1/) )

      T_2nd_9_1=reshape((                                                     
     &       /1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6, 
     &        1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4, 
     &        1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2, 
     &        1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0, 
     &        1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8, 
     &        2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6, 
     &        2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4, 
     &        2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2, 
     &        2.6,  2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0 /),
     &        shape(T_2nd_9_1,1), order=(/2,1/))

        !print*, "inverse matrix 9"
        Re_2nd_9=inv_T_2D(T_2nd_9_1)
        !print*, inv_T_2D(T_2nd_3_1)
        !call print_matrix(Re_2nd_9,9)
c Example to test 9x9 matrix.
c     https://www.ibm.com/docs/en/essl/6.2?topic=subroutines-sgetrf-dgetrf-cgetrf-zgetrf-general-matrix-factorization

        !call contrac_2nd_2nd(T_2nd_9_1,Re_2nd_9,Re_0th) !Using a declaration as routine
        Re_0th=contrac_2nd_2nd(T_2nd_9_1,Re_2nd_9) !Using a declaration as function
        print*, "contraction 2nd order tensors"
        print*, Re_0th
        n2=size(Re_4th_3,4)
        print*, "n2 outside", n2
        !call diadic_prod_T2_T2(T_2nd_3_1,T_2nd_3_1,Re_4th_3) !Using a declaration as routine
        Re_4th_3=diadic_prod_T2_T2(T_2nd_3_1,T_2nd_3_1) !Using a declaration as function
        print*, "Diadic Product two 2nd order 3x3"
        call print_4D_tensor(Re_4th_3,3)
!
        !call contrac_4th_2nd(Re_4th_3, T_2nd_3_1, Re_2nd_3) !Using a declaration as routine
        Re_2nd_3= contrac_4th_2nd(Re_4th_3, T_2nd_3_1) !Using a declaration as function
        print*, "Contraction 4th order and 2nd order 3x3 tensors"
        call print_matrix(Re_2nd_3,3,3)

c Test loading
         dt=0.2
        !te_n1=0.0
        ltype=1; !4
        if (ltype==1) then
            t(1:2)=(/0.0, 1.0/)
            lam(1:2)=(/0.0, 10.0/)
        elseif (ltype==2) then
            t(1:3)=(/0.0, 1.0, 2.0/)
            lam(1:2)=(/0.0, 1.0/)
        elseif (ltype==3) then
            t(1:3)=(/0.0, 1.0, 2.0/)
            lam=(/0.0, 10.0, 2.0/)
        elseif (ltype==4) then
            t(1:4)=(/0.0, 1.0, 3.0, 4.0/)
            lam=(/0.0,  10.0,  -10.0 /)
        elseif (ltype==5) then
            t=(/0.0, 2.5, 7.5, 12.5, 17.5, 22.5, 25.0/)
            lam=(/0.0, 0.005, -0.005/)
        end if

c        --------------------------------------------------------------------------

c        prescribed load/time step
        !dt=0.1;
c       !start and end-time of loading, time-scale, no. of steps

        
        n_step=nint((t(posi_last_time(ltype))-t(1))/dt)
        type_linespace=1
        print*, "steps", n_step
         
        allocate (time(n_step+1)) 
        call linespace(t(1),t(posi_last_time(ltype)),time,
     &      type_linespace)
        print*,time
        
c        Compute the strain load ramp
        allocate (e11(n_step+1))
c    
          call loading2 (ltype,posi_last_time , dt, t ,  
     $     lam, e11)
         print*,"loading 2", e11

      end program test_tensor
        





