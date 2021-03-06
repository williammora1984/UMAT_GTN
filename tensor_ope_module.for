c======================================================================
c                    Module Tensor operations
c=======================================================================
c  This module contains all the operations required to do tensor operations, 
c  map into equivalent notations of tensor. print in a sorted way in the
c termianl, and auxiliary functions used in testing of the material routine.  
c-----------------------------------------------------------------------
c Version 0.9.4  Include Inverse without external package
c Oct 2021
c-----------------------------------------------------------------------
      module tensor_operations
c       INCLUDE 'ABA_PARAM.INC'        
        implicit none
      contains

c=======================================================================
c                    subroutine linespace
c=======================================================================
c=======================================================================  
c subroutine linespace(ini, end, vector, type): returns  a  1st order  
C tensor with even spaced double precision numbers according to the size
c of vector  
c-----------------------------------------------------------------------
c Inputs
c   ini: initial value 
c   end: end value
c   vector: 1st order tensor with the right size to store the result
c   type: select if incluedes lower end
c       type= 0: not include lower end
c       type= 1: includes lower end
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=======================================================================      
      subroutine linespace (ini, end, vector, type)
        implicit none
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
c-----------------------------------------------------------------------
c Inputs
c   array: array of size n1xn1
c-----------------------------------------------------------------------
c Output
c   tr_A: trace
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=======================================================================     
      function Trace(A) result(Tr_A)
        implicit none
        double precision, dimension(:,:):: A
        double precision :: Tr_A
        integer n1
        Tr_A=0.0
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

c=======================================================================
c                    function inv_T_2D
c=======================================================================
c=======================================================================  
c Returns the inverse of a 2nd order tensor using the LU decomposition.  
c Depends on the LAPACK Fortran package
c http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
c-----------------------------------------------------------------------
cInputs
c    A: Tensor 2D se nxn
c-----------------------------------------------------------------------
cOutput
c    A: inverted tensor
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=======================================================================
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
c                    inv_T_2Dg
c=======================================================================
c=======================================================================  
c Returns the inverse of a 2nd order tensor according to the Gaus Jordan
c elimination, the code was modified from
c http://computer-programming-forum.com/49-fortran/6083a0ae451dd206.htm
c-----------------------------------------------------------------------
cInputs
c    A: Tensor 2D se nxn
c-----------------------------------------------------------------------
cOutput
c    A: inverted tensor
c-----------------------------------------------------------------------
c coded by: -----
c=======================================================================

      SUBROUTINE inv_T_2Dg(A,n,Ainv) 
      INTEGER :: m,n,NMAX
      double precision :: a(n,n),Ainv(n,n)
      PARAMETER (NMAX=50)
      INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      double precision :: big,dum,pivinv
      Ainv=A
      do 11 j=1,n
        ipiv(j)=0
 11   continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(Ainv(j,k)).ge.big)then
                  big=abs(Ainv(j,k))
                  irow=j
                  icol=k
                endif

              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
 12          continue
          endif
 13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=Ainv(irow,l)
            Ainv(irow,l)=Ainv(icol,l)
            Ainv(icol,l)=dum
 14        continue
c          do 15 l=1,m
c            dum=b(irow,l)
c            b(irow,l)=b(icol,l)
c            b(icol,l)=dum
c 15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (Ainv(icol,icol).eq.0.) pause 'singular matrix in gaussj'

        pivinv=1./Ainv(icol,icol)
        Ainv(icol,icol)=1.
        do 16 l=1,n
          Ainv(icol,l)=Ainv(icol,l)*pivinv
 16      continue
c        do 17 l=1,m
c          b(icol,l)=b(icol,l)*pivinv
c 17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=Ainv(ll,icol)
            Ainv(ll,icol)=0.
            do 18 l=1,n
              Ainv(ll,l)=Ainv(ll,l)-Ainv(icol,l)*dum
 18          continue
c            do 19 l=1,m
c              b(ll,l)=b(ll,l)-b(icol,l)*dum
c 19          continue
          endif
 21      continue
 22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then

          do 23 k=1,n
            dum=Ainv(k,indxr(l))
            Ainv(k,indxr(l))=Ainv(k,indxc(l))
            Ainv(k,indxc(l))=dum
 23        continue
        endif
 24    continue
      return
      END subroutine inv_T_2Dg

c=======================================================================
c                    inv_T_2Dd
c=======================================================================
c=======================================================================  
c Returns the inverse of a 2nd order tensor of maximum size equal to 3x3
c-----------------------------------------------------------------------
cInputs
c    A: Tensor 2D se 2x2 or 3x3
c-----------------------------------------------------------------------
cOutput
c    A: inverted tensor
c    n: size of square matrix
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=======================================================================

      subroutine inv_T_2Dd(A,n,Ainv)

c        INCLUDE 'ABA_PARAM.INC'
        Implicit none
        integer :: n  
        double precision :: A(n,n), Ainv(n,n)
        double precision :: Det
        double precision::Ainv2(2,2), Ainv3(3,3), Cof(3,3), Adj(3,3) 
        
        if (n==2) then
          Det=A(1,1)*A(2,1)-A(1,2)*A(2,1)
            if (Det==0) then
               WRITE(6,1)
 1             FORMAT(//,30X,'***Error - Singular matrix',
     1          'without inverse')
             end if
          Ainv=reshape((/A(2,2), -A(1,2),
     &                   -A(2,1), A(1,1)/), shape(Ainv2), order=(/2,1/))
     &                  /Det
             
        !In this case Dets is the cofactor matrix
        elseif (n==3) then
          Cof(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
          Cof(2,1)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
          Cof(3,1)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
          Cof(1,2)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
          Cof(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
          Cof(3,2)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
          Cof(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
          Cof(2,3)=(A(1,1)*A(2,3)-A(1,3)*A(2,2))
          Cof(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
          
          !Adjugate Matrix
          Adj(1,1)=Cof(1,1)
          Adj(2,1)=Cof(1,2)
          Adj(3,1)=Cof(1,3)
          Adj(1,2)=Cof(2,1)
          Adj(2,2)=Cof(2,2)
          Adj(3,2)=Cof(2,3)
          Adj(1,3)=Cof(3,1)
          Adj(2,3)=Cof(3,2)
          Adj(3,3)=Cof(3,3)
          Det=(A(1,1)*Cof(1,1)+A(2,1)*Cof(2,1)+A(3,1)*Cof(3,1))
          if (Det==0) then
            WRITE(6,1)
          end if          
          Ainv=Adj/Det
        else
         WRITE(6,2)
 2       FORMAT(//,30X,'***Error -This function can be used only with',
     2          '2x2 and 3x3 matrices')
        end if

      
      end subroutine inv_T_2Dd                   

c=======================================================================  
c                   subroutine I4dikdjl
cc=======================================================================
cc==================================================================  
c Returns the 4th order unit tensor I4ikjl 
c-----------------------------------------------------------------------c-----------------------------------------------------------------------
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
                      I4dikdjl_o(i,j,k,l)=I_3D(i,k)*I_3D(j,l)
                     end do
                 end do
             end do
         end do
      end subroutine I4dikdjl

c=======================================================================  
c                   subroutine I4dildjk
c=======================================================================
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
                      I4dildjk_o(i,j,k,l)=I_3D(i,l)*I_3D(j,k)                   
                     end do
                 end do
             end do
         end do
      end subroutine I4dildjk

c=======================================================================  
c                   subroutine I4dildjk
c=======================================================================
c=======================================================================  
c  Returns the 4th order unit tensor Iiljk
c-----------------------------------------------------------------------
cInputs
c    I4dildjk_o: 4th order tensor
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=====================================================================

      subroutine I4dijdkl(I4dijdkl_o)
        implicit none
         
         integer:: i, j, k , l     !iterators
         integer:: n1              !matrix dimension 
         double precision, dimension(3,3,3,3)::I4dijdkl_o
         double precision, dimension(3,3):: I_3D !Auxiliary 2D order tensors
         
         n1=3
         call Ident1(I_3D,n1)
      
         do i=1,n1
             do j=1,n1
                 do k=1,n1
                     do l=1,n1
                      I4dijdkl_o(i,j,k,l)=I_3D(i,j)*I_3D(k,l)
                     end do
                 end do
             end do
         end do
      end subroutine I4dijdkl

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
                        P4sym_o(i,j,k,l) =0.0  
                        P4sym_o(i,j,k,l) =I4sym_o(i,j,k,l)-         
     &                     I_3D(i,j)*I_3D(k,l)/3.000000
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
        call inv_T_2Dg(A6x6, 6,A6x6_inv)
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

        Contr =0.00000000
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
        implicit none
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
        implicit none
        double precision, dimension(3,3,3,3):: T_4th                      !2nd order tensor
        double precision, dimension(size(T_4th,3),size(T_4th,4)) ::T_2D
     $   , Contr                                                          !4th order tensor
        integer:: i, j, k, l 
        integer:: n  !size dimension tensor
    
        n=size(T_4th,1)
        do i=1,3
            do j=1,3
                do k=1,3
                    do l=1,3  
                      !Contr(j,i) =Contr(j,i)+T_4th(j,i,l,k)*T_2D(l,k)
                      Contr(i,j)= Contr(i,j)+T_4th(i,j,k,l)*T_2D(k,l)
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
       
       sq2 =2.00000**0.5
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
        double precision:: sq2,const1                      !Auxiliary scalar

        sq2 =2.000000**0.5
        const1=2.000000
        elements=size(tensor_4th)
        n1=size(tensor_4th,1)
        n2=size(tensor_4th,2)
          if (elements .eq. 81) then         
           Vn_2D=reshape( (  
     &               /tensor_4th(1,1,1,1),         tensor_4th(1,1,2,2),     
     &                tensor_4th(1,1,3,3),     sq2*tensor_4th(1,1,1,2), 
     &            sq2*tensor_4th(1,1,2,3),     sq2*tensor_4th(1,1,1,3),
     &                tensor_4th(2,2,1,1),         tensor_4th(2,2,2,2),     
     &                tensor_4th(2,2,3,3),     sq2*tensor_4th(2,2,1,2), 
     &            sq2*tensor_4th(2,2,2,3),     sq2*tensor_4th(2,2,1,3),
     &                tensor_4th(3,3,1,1),         tensor_4th(3,3,2,2),
     &                tensor_4th(3,3,3,3),     sq2*tensor_4th(3,3,1,2), 
     &            sq2*tensor_4th(3,3,2,3),     sq2*tensor_4th(3,3,1,3),
     &            sq2*tensor_4th(1,2,1,1),     sq2*tensor_4th(1,2,2,2), 
     &            sq2*tensor_4th(1,2,3,3),  const1*tensor_4th(1,2,1,2),
     &              2*tensor_4th(1,2,2,3),  const1*tensor_4th(1,2,1,3),
     &            sq2*tensor_4th(2,3,1,1),     sq2*tensor_4th(2,3,2,2), 
     &            sq2*tensor_4th(2,3,3,3),  const1*tensor_4th(2,3,1,2),   
     &         const1*tensor_4th(2,3,2,3),  const1*tensor_4th(2,3,1,3),
     &            sq2*tensor_4th(1,3,1,1),     sq2*tensor_4th(1,3,2,2), 
     &            sq2*tensor_4th(1,3,3,3),  const1*tensor_4th(1,3,1,2), 
     &         const1*tensor_4th(1,3,2,3), const1*tensor_4th(1,3,1,3)/),
     &          shape(Vn_2D), order=(/2,1/) )

          end if

      end subroutine T4th_2_Voig

c=======================================================================
c                    T2nd_2_Voig
c=======================================================================
c=======================================================================  
c  maps 2nd order tensor to 1st order tensor according the voight 
cnotation
c-----------------------------------------------------------------------
cInputs
c    T_2D: 6x6 2nd order tensor in voigth notation
c    T_1D: vector 6 size in voigth notation
c-----------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c=======================================================================
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

c=======================================================================
c                    4th_2_66_sym
c=======================================================================
c=======================================================================  
c  maps 4th order symmetric tensor to 2nd order tensor 6 x 6
c-----------------------------------------------------------------------
cInputs
c     tensor_4th : 4th order symmetric tensor  
c    M66: 6x6 2nd order tensor in voigth notation
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=======================================================================

      subroutine m4th_2_66_sym (tensor_4th, M66)
         implicit none

         integer, dimension(6):: ii, jj
         integer:: rc, sc
         double precision:: M66(6,6), tensor_4th(3,3,3,3) 
         ii = (/1,2,3,1,2,1/) 
         jj = (/1,2,3,2,3,3/) 
         do rc=1,6
           do sc=1,6
               M66(rc,sc) = tensor_4th(ii(rc),jj(rc),ii(sc),jj(sc))
           end do 
         end do
      end subroutine m4th_2_66_sym

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
          T_2D_Par=T_2D(2:6,2:6)
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

c=======================================================================
c            SUBROUTINES TO MAKE THE TEST
C=======================================================================      
                
c======================================================================
c                 subroutine read2
c=======================================================================
c==================================================================  
c Import the data form a csv file
c------------------------------------------------------------------]
c Inputs
c   import1   : 2D grade tensor to store the readed information
c   n1        : dimesion 1 of the array that contain the information
c   n2        : dimesion 2 of the array that contain the information
c   file_name : information source file name. (must be a csv file)
c------------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c==================================================================

      subroutine read2(import1, n1,n2, file_name)
          implicit none
          integer:: i, j,end,stat,line_no,n1,n2     !iterators
          
          character(len=13) :: file_name
          double precision :: import1(n1,n2)
          n1=size(import1,1)
          n2=size(import1,2)

          open(15, file=file_name,access='sequential',
     &      form="formatted", iostat=stat)
          do i = 1,n1
             read(15,*,iostat=stat) import1(i,:)
          end do
          close (15)

      end subroutine read2

      subroutine read3(import1, n1,n2, file_name)
          implicit none
          integer:: i, j,end,stat,line_no,n1,n2     !iterators
          
          character(len=30) :: file_name
          double precision :: import1(n1,n2)
          n1=size(import1,1)
          n2=size(import1,2)
  
          open(15, file=file_name,access='sequential',
     &      form="formatted", iostat=stat)
          do i = 1,n1
             read(15,*,iostat=stat) import1(i,:)
          end do
          close (15)
  
      end subroutine read3

c======================================================================
c                 subroutine emp_4th
c=======================================================================
c==================================================================  
c Empties a 4th grade tensor
c------------------------------------------------------------------]
c Inputs
c   T_4th   : 4yh grade tensor to empty
c------------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c==================================================================

      subroutine emp_4th(T_4th)
        implicit none
        double precision, dimension(3,3,3,3):: T_4th                      !2nd order tensor
        double precision, dimension(size(T_4th,3),size(T_4th,4)) ::T_2D
     $   , Contr                                                          !4th order tensor
        integer:: i, j, k, l 
        integer:: n  !size dimension tensor
    
        n=size(T_4th,1)
        do i=1,3
            do j=1,3
                do k=1,3
                    do l=1,3  
                      !Contr(j,i) =Contr(j,i)+T_4th(j,i,l,k)*T_2D(l,k)
                      T_4th(i,j,k,l)= 0.000
                    end do
                  end do
                end do
              end do

            end subroutine emp_4th

c======================================================================
c                 subroutine emp_2D
c=======================================================================
c==================================================================  
c Empties a 2nd grade tensor
c------------------------------------------------------------------]
c Inputs
c   T   : 2nd grade tensor to empty
c------------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c==================================================================

      subroutine emp_2D(T)
        implicit none
        double precision, dimension(3,3):: T   !2nd order tensor
        integer:: i, j
        integer:: n  !size dimension tensor
    
        do i=1,3
            do j=1,3
              T(i,j)= 0.000

          end do
        end do

      end subroutine emp_2D

