c======================================================================
c                    Module Tensor operations
c=======================================================================
c  This module contains all the operations required to do tensor operations, 
c  map into equivalent notations of tensor. print in a sorted way in the
c termianl, and auxiliary functions used in testing of the material routine.  
c-----------------------------------------------------------------------
c Version 0.9.2
c Oct 2021
c-----------------------------------------------------------------------
      module tensor_operations
       INCLUDE 'ABA_PARAM.INC'        
        implicit none
      contains

c=======================================================================
c                    subroutine linespace
c=======================================================================
c=======================================================================  
c subroutine linespace(ini, end, vector): returns  a  1st order tensor with 
C even spaced double precision numbers according to the size of vector 
c-----------------------------------------------------------------------
c Inputs
c   ini: initial value 
c   end: end value
c   vector: 1st order tensor with the right size to store the result
c   type: select if incluedes lower end
c       type= 0: not include loweer end
c       type= 1: includes lower end
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=======================================================================      
      subroutine linespace (ini, end, vector, type)
        double precision ::ini, end
        double precision, dimension (:) :: vector
        integer :: i, type

        if (type==1) then
          do i=1, (size(vector,1))
            vector(i)=ini+(end-ini)*(i-1)/(size(vector,1)-1)
          end do
        else
          do i=1, (size(vector,1))
            vector(i)=ini+(end-ini)*(i)/(size(vector,1))
          end do
        end if
      end subroutine linespace  

c=======================================================================
c                    funcion Trace
c=======================================================================
c=======================================================================  
c Returns the trace of a  2nd order tensor
c------------------------------------------------------------------]
c Inputs
c   array: array of size n1xn1
c-----------------------------------------------------------------
c Output
c   tr_A: trace
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c==================================================================      
      function Trace(A) result(Tr_A)
        implicit none
        double precision, dimension(:,:):: A
        double precision :: Tr_A
        integer n1
        Tr_A=0
        do n1=1,size (A,1)
          Tr_A=Tr_A+A(n1,n1)
        end do
      end function Trace

c=======================================================================
c                    subroutine Ident1
c=======================================================================
c==================================================================  
c subroutine Ident1(array,n1): returns a 2nd order identity  tensor
c------------------------------------------------------------------]
c Inputs
c   array: array of size n1xn1
c   n1 :   dimension of the array
c-----------------------------------------------------------------
c Output
c   array: identy tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c==================================================================      

      subroutine Ident1(array,n1)
        implicit none
        integer:: i, j     !iterators
        integer:: n1  !dimension square tensor
        double precision :: array(n1,n1)

        do i=1,n1
            do j=1,n1
                array(j,i)=0.0
            end do
        end do

        do i=1,n1
            array(i,i)=1.0
        end do
      end subroutine Ident1

c======================================================================
c                    function inv_T_2D
c=======================================================================
c==================================================================  
c Returns the inverse of a 2nd order tensor using the LU decomposition.  
c Depends on the LAPACK Fortran package
c http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
c------------------------------------------------------------------]
cInputs
c    T_2D: 6x6 2nd order tensor in voigth notation
c    T_1D: vector 6 size in voigth notation
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      function inv_T_2D(A) result(Ainv)
        external:: DGETRF
        external:: DGETRI 
        double precision, dimension(:,:), intent(in) :: A
        double precision, dimension(size(A,1),size(A,2)) :: Ainv
        double precision, dimension(size(A,1),size(A,2)) :: Ainv_temp
        double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
       
         ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv_temp =dble(A)
        n = size(A,1)
       
         ! DGETRF computes an LU factorization of a general M-by-N matrix A
         ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv_temp, n, ipiv, info)
        !Error mesage genetrated in the LU decomposition, see the possible
        !option of values for info in LAPACK documentation
        if (info /= 0) then
           stop 'function inv_T_2D- Matrix is numerically singular!'
        end if
       
         ! DGETRI computes the inverse of a matrix using the LU factorization
        call DGETRI(n, Ainv_temp, n, ipiv, work, n, info)
       !Error mesage genetrated in the LU decomposition, see the possible
        !option of values for info in LAPACK documentation 
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if

        Ainv=Ainv_temp
      end function inv_T_2D

c=======================================================================  
c  Returns the 4th order unit tensor Iikjl
c-----------------------------------------------------------------------
cInputs
c    I4dikdjl_o: 4th order tensor
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=====================================================================

      subroutine I4dikdjl(I4dikdjl_o)
        implicit none
         
         integer:: i, j, k , l     !iterators
         integer:: n1              !matrix dimension 
         double precision, dimension(3,3,3,3)::I4dikdjl_o
         double precision, dimension(3,3):: I_3D !Auxiliary 2D order tensors
         
         n1=3
         call Ident1(I_3D,n1)
      
         do i=1,n1
             do j=1,n1
                 do k=1,n1
                     do l=1,n1
                      I4dikdjl_o(i,j,k,l)=I_3D(i,k)*
     &                   I _3D(j,l)
                     end do
                 end do
             end do
         end do
      end subroutine I4dikdjl

c=======================================================================  
c  Returns the 4th order unit tensor Iiljk
c-----------------------------------------------------------------------
cInputs
c    I4dildjk_o: 4th order tensor
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=====================================================================

      subroutine I4dildjk(I4dildjk_o)
        implicit none
         
         integer:: i, j, k , l     !iterators
         integer:: n1              !matrix dimension 
         double precision, dimension(3,3,3,3)::I4dildjk_o
         double precision, dimension(3,3):: I_3D !Auxiliary 2D order tensors
         
         n1=3
         call Ident1(I_3D,n1)
      
         do i=1,n1
             do j=1,n1
                 do k=1,n1
                     do l=1,n1
                      I4dildjk_o(i,j,k,l)=I_3D(i,l)*
     &                   I _3D(j,k)
                     end do
                 end do
             end do
         end do
      end subroutine I4dildjk

c======================================================================
c                    Subroutine I4sym
c=======================================================================
c==================================================================  
c  Returns the 4th order symmetric identity tensor
c-------------------------------------------------------------------
cInputs
c    I4sym_o: 4th order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=====================================================================

      subroutine I4sym(I4sym_o)
       implicit none
        
        integer:: i, j, k , l     !iterators
        integer:: n1              !matrix dimension 
        double precision, dimension(3,3,3,3)::I4sym_o
        double precision, dimension(3,3):: I_3D !Auxiliary 2D order tensors
        
        n1=3
        call Ident1(I_3D,n1)

        do i=1,n1
            do j=1,n1
                do k=1,n1
                    do l=1,n1
                        I4sym_o(i,j,k,l)=(1.0/2.0)*(I_3D(i,k)*
     &                   I_3D(j,l)+I_3D(i,l)*I_3D(j,k))
                    end do
                end do
            end do
        end do
      end subroutine I4sym


c======================================================================+
c                    Subroutine I4dikdjl
c=======================================================================

c======================================================================
c                    subroutine P4sym
c=======================================================================
c==================================================================  
c  Returns 4th order deviatoric projection tensor
c------------------------------------------------------------------]
cInputs
c    P4sym: 4th order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      subroutine P4sym(P4sym_o)

       implicit none
       integer:: i, j, k , l     !iterators
       integer:: n1              !matrix dimension 
       double precision, dimension(3,3,3,3)::I4sym_o, P4sym_o !4th order tensors
       double precision, dimension(3,3):: I_3D            !Auxiliary 2D order tensors

       n1=3
       call Ident1(I_3D,n1)
       
       call I4sym(I4sym_o)    

          do i = 1,3
              do j = 1,3
                  do k = 1,3
                      do l = 1,3
                        P4sym_o(i,j,k,l) =0  
                        P4sym_o(i,j,k,l) =I4sym_o(i,j,k,l)-1.0/3.0*         
     &                     I_3D(i,j)*I_3D(k,l) 
                      end do
                    end do
                  end do
                end do  
      end subroutine P4sym

c======================================================================
c                    Subroutine invT4
c=======================================================================
c==================================================================  
c  inverts a 4th order tensor
c------------------------------------------------------------------]
cInputs
c    tensor:    4th order tensor to invert
c    tensor_In: inverted tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      subroutine invT4(tensor,tensor_Inv)

       implicit none
       integer:: i, j, k , l     !iterators
       integer:: n1               !matrix dimension 
       double precision, dimension(3,3,3,3)::tensor, tensor_Inv !4th order tensors 
       double precision, dimension(6,6)::A6x6, A6x6_inv         !Auxiliary 2D order

        call T4th_2_Voig(tensor, A6x6)
        A6x6_inv=inv_T_2D(dble(A6x6))
        call Voig_2_T4th(A6x6_inv,tensor_Inv)

      end subroutine invT4

c=======================================================================
c                    Function contrac_2nd_2nd
c=======================================================================
c=======================================================================  
c Returns the contraction between two second order tensors
c-----------------------------------------------------------------------
cInputs
c    T_2D_1: 2nd order tensor
c    T_2D_2: 2nd order tensor
c------------------------------------------------------------------------
cOutput
c    Contr: result of the contraction, 0 order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      function contrac_2nd_2nd(T_2D_1, T_2D_2) result(Contr)

       implicit none 
       double precision, dimension(:,:) :: T_2D_1, T_2D_2 !2nd order tensor
       double precision :: Contr !Result of the contraction
       integer:: i, j !Iterators
       integer:: n    !dimension tensor. 

        Contr =0.0
        n=size(T_2D_1,1)

        do i=1,n
          do j=1,n
            Contr =Contr+ T_2D_1(j,i)*T_2D_2(j,i)
          end do
        end do

      end function contrac_2nd_2nd

c======================================================================
c                  Fucntion diadic_prod_T2_T2
c=======================================================================
c==================================================================  
c Returns the diadic product between two second order tensors
c------------------------------------------------------------------]
cInputs
c    T_2D_1: 2nd order tensor
c    T_2D_2: 2nd order tensor
c    Diad: result of the diadic product, 4th order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      function diadic_prod_T2_T2(T_2D_1,T_2D_2) result(Diad)  
      
        double precision, dimension(:,:) :: T_2D_1, T_2D_2  !2nd order tensor
        double precision, dimension(size(T_2D_1,1),size(T_2D_1,1),
     &                 size(T_2D_1,1),size(T_2D_1,1)):: Diad !4th order tensor   
        integer:: i, j, k, l !Iterators
        integer:: n,n2       !size dimension tensor

        n=size(T_2D_1,1)
        n2=size(Diad,4)  

        do i=1,n
            do j=1,n
                do k=1,n
                    do l=1,n
                      Diad(j,i,l,k) = T_2D_1(j,i)*T_2D_2(l,k)
                    end do
                  end do
                end do
              end do

      end function diadic_prod_T2_T2

c=======================================================================
c                    function contrac_4th_2nd
c=======================================================================
c=======================================================================
c Returns the contraction between a 4th order tensor a a 2nd order tensor
c-----------------------------------------------------------------------
cInputs
c    T_4th: 4th order tensor
c    T_2D: 2nd order tensor
c-----------------------------------------------------------------------
cOutput
c    Contr: result of the contraction, a 2nd order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=======================================================================
      function contrac_4th_2nd(T_4th,T_2D) result(Contr)

        double precision, dimension(:,:,:,:):: T_4th                      !2nd order tensor
        double precision, dimension(size(T_4th,3),size(T_4th,4)) :: T_2D
     $   , Contr                                                          !4th order tensor
        integer:: i, j, k, l 
        integer:: n  !size dimension tensor
    
        n=size(T_4th,1)
        do i=1,3
            do j=1,3
                do k=1,3
                    do l=1,3
                      Contr(i,j)=0  
                      Contr(j,i) =T_4th(j,i,l,k)*T_2D(l,k)
                    end do
                  end do
                end do
              end do

      end function contrac_4th_2nd

c======================================================================
c                    subroutine Voig_2_T4th
c=======================================================================
c==================================================================  
c  maps a 2D (6x6) tensor in voigth notation into a 4th order tensor 
c------------------------------------------------------------------]
cInputs
c    tensor: 4th order tensor to map
c    T_2D: Voigth notation (6x6) 2nd order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================

      subroutine voig_2_T4th(T_2D, tensor_4th)
       implicit none
       
       integer:: i, j, k , l                            !iterators
       integer:: n1,n2,elements                         !dimension matrix
       integer, dimension(3,3)::nn                      !Auxiliar iterators
       double precision, dimension(6,6)::T_2D           !2nd order tensor 
       double precision, dimension(6,6):: Aux,temp      !2nd order auxiliary tensor
       double precision, dimension(3,3,3,3)::tensor_4th !4th order
       double precision:: sq2, tem_n1, tem_n2           !Auxiliary scalars
       
       sq2 =2.00**0.5
       tem_n1=1.0
       tem_n2=2.0

       Aux=reshape((
     &        / tem_n1,   tem_n1,   tem_n1, sq2,    sq2,    sq2,
     &          tem_n1,   tem_n1,   tem_n1, sq2,    sq2,    sq2,
     &          tem_n1,   tem_n1,   tem_n1, sq2,    sq2,    sq2,
     &          sq2,      sq2,      sq2,    tem_n2, tem_n2, tem_n2,
     &          sq2,      sq2,      sq2,    tem_n2, tem_n2, tem_n2,
     &          sq2,      sq2,      sq2,    tem_n2, tem_n2, tem_n2 /),
     &      shape(Aux), order=(/2,1/) )
 
        nn=reshape(( /1, 4, 6,
     &                4, 2, 5,
     &                6, 5, 3/),shape(nn), order =(/2,1/)) 
  
        do i=1,6
          do j=1,6
            temp(j,i)=0
            temp(j,i)=T_2D(j,i)/Aux(j,i)
          end do
        end do
        
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                tensor_4th(i,j,l,k)=temp(nn(i,j),nn(l,k))
              end do
            end do
          end do
        end do
      end subroutine voig_2_T4th

c======================================================================
c                    subroutine T4th_2_Voig
c=======================================================================
c==================================================================  
c  maps a 4th order tensor into a 2th order tensor according 
cthe voight notation
c------------------------------------------------------------------]
cInputs
c    tensor: 4th order tensor to map
c    Vn_2D: kelvin notation 6x6 2nd order tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================

      subroutine T4th_2_Voig(tensor_4th, Vn_2D)
        implicit none
        integer:: i, j, k , l                              !iterators
        integer:: n1,n2,elements                           !dimension matrix
        double precision, dimension(3,3,3,3)::tensor_4th   !4nd order tensor
        double precision, dimension(6,6)::Vn_2D            !2nd order tensor 
        double precision:: sq2                             !Auxiliary scalar

        sq2 =2.0**0.5
        elements=size(tensor_4th)
        n1=size(tensor_4th,1)
        n2=size(tensor_4th,2)
          if (elements .eq. 81) then         
           Vn_2D=reshape( (  
     &               /tensor_4th(1,1,1,1),     tensor_4th(1,1,2,2),     
     &                tensor_4th(1,1,3,3), sq2*tensor_4th(1,1,1,2), 
     &            sq2*tensor_4th(1,1,2,3), sq2*tensor_4th(1,1,1,3),
     &                tensor_4th(2,2,1,1),     tensor_4th(2,2,2,2),     
     &                tensor_4th(2,2,3,3), sq2*tensor_4th(2,2,1,2), 
     &            sq2*tensor_4th(2,2,2,3), sq2*tensor_4th(2,2,1,3),
     &                tensor_4th(3,3,1,1),     tensor_4th(3,3,2,2),
     &                tensor_4th(3,3,3,3), sq2*tensor_4th(3,3,1,2), 
     &            sq2*tensor_4th(3,3,2,3), sq2*tensor_4th(3,3,1,3),
     &            sq2*tensor_4th(1,2,1,1), sq2*tensor_4th(1,2,2,2), 
     &            sq2*tensor_4th(1,2,3,3),   2*tensor_4th(1,2,1,2),
     &              2*tensor_4th(1,2,2,3),   2*tensor_4th(1,2,1,3),
     &            sq2*tensor_4th(2,3,1,1), sq2*tensor_4th(2,3,2,2), 
     &            sq2*tensor_4th(2,3,3,3),   2*tensor_4th(2,3,1,2),   
     &              2*tensor_4th(2,3,2,3),   2*tensor_4th(2,3,1,3),
     &            sq2*tensor_4th(1,3,1,1), sq2*tensor_4th(1,3,2,2), 
     &            sq2*tensor_4th(1,3,3,3),   2*tensor_4th(1,3,1,2), 
     &              2*tensor_4th(1,3,2,3),   2*tensor_4th(1,3,1,3)/),
     &          shape(Vn_2D), order=(/2,1/) )

          end if

      end subroutine T4th_2_Voig

c======================================================================
c                    T2nd_2_Voig
c=======================================================================
c==================================================================  
c  maps 2nd order tensor to 1st order tensor according the voight 
cnotation
c------------------------------------------------------------------]
cInputs
c    T_2D: 6x6 2nd order tensor in voigth notation
c    T_1D: vector 6 size in voigth notation
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      subroutine T2nd_2_Voig(T_2D,T_1D)
        
        implicit none
        integer:: i, j                          !iterators
        integer:: n1,n2,elements                !dimension matrix
        double precision, dimension(6,6)::T_2D  !2nd order tensor
        double precision, dimension(6)::T_1D    !1st order tensor
        double precision:: sq2                  !Auxiliar scalar
        
        sq2 =2.0**0.5
        elements=size(T_2D)
        if (elements==9) then

            T_1D=(/T_2D(1,1), T_2D(2,2), T_2D(3,3), sq2*T_2D(1,2), 
     &         sq2*T_2D(2,3), sq2*T_2D(1,3)/)
        end if

      end subroutine T2nd_2_Voig  

c======================================================================
c                    subroutine print_matrix
c=======================================================================
c==================================================================  
c subroutine print in the terminal a 2D tensor of dimension n1 x n2
c------------------------------------------------------------------]
c Inputs
c   array:     array of size n1xn2
c   n1, n2 :   array dimensions 
cc------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c==================================================================      
      subroutine print_matrix(array,n1,n2)
        implicit none
        integer:: i, j     !iterators
        integer:: n1, n2  !dimension matrix
        double precision :: array(n1,n2) !2D test tensors

        do i=1,n1
          write(*,*) (array(i,j), j = lbound(array,2), ubound(array,2))
        end do
      end subroutine print_matrix  
      
c======================================================================
c                    subroutine print_4D_tensor
c=======================================================================
c==================================================================  
c Prints in the terminal a 4th grade tensor of dimension n1 x n1 x n1 x n1
c------------------------------------------------------------------]
cInputs
c    tensor :
c     n1 : size of the inner 2nd grade tensor
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c===================================================================
      subroutine print_4D_tensor(tensor,n1)
        implicit none
        integer:: i, j,k,l                            !iterators
        integer:: n1, n2                              !dimension matrix
        double precision , dimension(3,3,3,3)::tensor !4th order tensor

        print*,"tensor 4th grade"
        do i = 1,3
          do j = 1,3
            print*,"submatrix", j,",",i
            call print_matrix(tensor(j,i,:,:),n1,n1)
          end do
        end do

      end subroutine print_4D_tensor  

c======================================================================
c                 Subroutine loading2
c=======================================================================
c==================================================================  
c Return the ramp of the strain load for a strain driven tension test
c------------------------------------------------------------------
c Inputs
c       ltype: load tyoe
c           1. Simple load
c           2. Load and unload (only compression)
c           3.
c           4. Load and compression (full cycle)
c       posi_last_time: indeces that indicates the position in the vector 
c       "t" of last data for each kind of load
c       dt: time step size
c       t: vector with the range of time for the load, its contents depends 
c          on each load type 
c       lam: vector with the range of strain for the load, its contents depends 
c          on each load type
c       load : vector with the strain load ramp
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c==================================================================

      subroutine loading2(ltype,posi_last_time,dt,t,lam,load) 
        !use tensor_operations
      
        implicit none
        integer:: ltype, type  !
        integer, dimension(5) :: posi_last_time
        integer::n_step_a, n_step_b, n_step_c
        double precision :: dt
        double precision, dimension(7) :: t
        double precision, dimension(3) ::lam
        double precision, dimension(:):: load
        double precision, dimension(:), allocatable :: xa, ya, xb, yb, 
     $                                                 xc, yc
        double precision ::xa1, ya1, xb1, yb1, xc1, yc1 
        double precision ::xa2, ya2, xb2, yb2, xc2, yc2
        print*, "Enter load general"
        
        !ltype==1 correspond to a load with one linear load
        if (ltype==1) then
          type=1
          n_step_a=nint((t(2)-t(1))/dt)
          allocate (xa(n_step_a+1))
          allocate (ya(n_step_a+1))
          
          xa1=t(1); ya1=lam(1)
          xa2=t(2); ya2=lam(2)
          call linespace(ya1,ya2,ya,type)

          load=ya
          deallocate(xa)
          deallocate(ya)    

        !ltype==2 correspond to a semicycle with one linear load and a linear unload
        elseif (ltype==2) then
          type=1
          n_step_a=nint((t(2)-t(1))/dt)
          allocate (xa(n_step_a+1))
          allocate (ya(n_step_a+1))
          call linespace(t(1),t(2),xa,type)
          
          xa1=t(1); ya1=lam(1)
          xa2=t(2); ya2=lam(2)
          call linespace(ya1,ya2,ya,type)
          
          type=0
          n_step_b=nint((t(3)-t(2))/dt)
          allocate (xb(n_step_b))
          allocate (yb(n_step_b))
          call linespace(t(2),t(3),xb, type)
          xb1=t(2); yb1=lam(2)
          xb2=t(3); yb2=lam(1)
          yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1

          deallocate (xa)
          deallocate (xb)
          load=(/ya,yb/)
          deallocate (ya)
          deallocate (yb)   
          print*, "End load 2"

        !ltype==4 correspond to a load with 1 load, 1 unload including
        ! a negative part and a load to the zero point
        elseif (ltype==4) then
          type=1
          n_step_a=nint((t(2)-t(1))/dt)
          allocate (xa(n_step_a+1))
          allocate (ya(n_step_a+1))
          call linespace(t(1),t(2),xa,type)
          
          xa1=t(1); ya1=lam(1)
          xa2=t(2); ya2=lam(2)
          call linespace(ya1,ya2,ya,type)

          type=0
          n_step_b=nint((t(3)-t(2))/dt)
          allocate (xb(n_step_b))
          allocate (yb(n_step_b))
          call linespace(t(2),t(3),xb, type)
          xb1=t(2); yb1=lam(2)
          xb2=t(3); yb2=lam(3)
          yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1
          
          n_step_c=nint((t(4)-t(3))/dt)
          allocate (xc(n_step_c))
          allocate (yc(n_step_c))
          call linespace(t(3),t(4),xc, type)
          xc1=t(3); yc1=lam(3)
          xc2=t(4); yc2=lam(1)
          yc=(yc2-yc1)/(xc2-xc1)*(xc-xc1)+yc1

          deallocate (xa)
          deallocate (xb)
          deallocate (xc)
          load=(/ya, yb, yc/)
          deallocate (ya)
          deallocate (yb)
          deallocate (yc)
          print*, "End load 4"

        end if
      end subroutine loading2

c======================================================================
c                    function partition_1D
c=======================================================================
c==================================================================  
c function partition_1D (T_1D)
c------------------------------------------------------------------]
c Inputs
c   T_1D: 1st order tensor  dimension 1x6
c-----------------------------------------------------------------
c Outputs
c   T_1D_Par: 1st order tensor  dimension 1x5
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c==================================================================
      function partition_1D(T_1D) result(T_1D_Par)

        implicit none
        double precision, dimension (6) ::T_1D
        double precision, dimension (5) ::T_1D_Par
        
        if (size(T_1D)==6) then
            T_1D_Par=T_1D(2:) 
        end if
      end function partition_1D

c=======================================================================
c                 function partition_2D_1D
c=======================================================================
c=======================================================================  
c returns a vecotr notation of a 3x3 simmetric 2nd order tensor without
c  the element 1,1.
c-----------------------------------------------------------------------
c Inputs
c   T_1D : 2 order symmetric tensor  
c-----------------------------------------------------------------------
c Outputs
c   T_1D_Par : 1st order voigth notation with out element 1,1.
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=======================================================================
      function partition_2D_1D(T_2D) result(T_1D_Par)
        implicit none

        double precision, dimension (3,3) ::T_2D    !2nd order tensor
        double precision, dimension (5) ::T_1D_Par  !1st order tensor
        
        if (size(T_2D)==9) then
            T_1D_Par=(/T_2D(2,2), 
     $                 T_2D(3,3), 
     $                 T_2D(1,2), 
     $                 T_2D(2,3), 
     $                 T_2D(1,3)/)
        end if

      end function partition_2D_1D

c======================================================================
c                 function partition_2D_2D
c=======================================================================
c==================================================================  
c returns a 2nd order tensor of dimension 5x5 from a 2nd order tensor of
c dimension 6x6 without include the forst row and colum.
c------------------------------------------------------------------]
c Inputs
c   T_2D: 2nd order tensor 5x5  
c-----------------------------------------------------------------
c Outputs
c   T_2D_Par: 2nd order tensor 6x6
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c==================================================================

      function partition_2D_2D(T_2D) result(T_2D_Par)
        implicit none

        double precision, dimension (6,6) ::T_2D     !2nd order tensor
        double precision, dimension (5,5) ::T_2D_Par !2nd order tensor

        if (size(T_2D)==36) then
          T_2D_Par=reshape(
     $       (/T_2D(2,2), T_2D(2,3), T_2D(2,4), T_2D(2,5), T_2D(2,6),
     $         T_2D(3,2), T_2D(3,3), T_2D(3,4), T_2D(3,5), T_2D(3,6),
     $         T_2D(4,2), T_2D(4,3), T_2D(4,4), T_2D(4,5), T_2D(4,6),
     $         T_2D(5,2), T_2D(5,3), T_2D(5,4), T_2D(5,5), T_2D(5,6),
     $         T_2D(6,2), T_2D(6,3), T_2D(6,4), T_2D(6,5), T_2D(6,6)/),
     $          shape(T_2D_Par),order=(/2,1/))
        end if

      end function partition_2D_2D

c======================================================================
c                 function Kron_d
c=======================================================================
c==================================================================  
c returns the evaluation of the Kronecker delta function given the indices.
c------------------------------------------------------------------]
c Inputs
c   i:index 1
c   j:index 2
c-----------------------------------------------------------------
c Outputs
c   kd_ij: result of the evaluation
c------------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c==================================================================

      function Kron_d(i,j) result(kd_ij)
        implicit none
      
        integer::i, j
        double precision :: kd_ij
        if (i==j) then
          kd_ij=1.0
        else
          kd_ij=0.0
        end if
      end function Kron_d

      end module tensor_operations 


