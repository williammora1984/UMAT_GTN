      program notes
        use tensor_operations
c        use read_data
        implicit none
c        double precision, dimension(11,19):: import1
        character(len=13) :: file_name

        double precision, dimension(3,3):: m2,m2p, m2pp
        double precision, dimension(6,6):: T_2nd_6
        double precision, dimension(9,9)::T_2nd_9_1, Re_2nd_9,Re_2nd_9_2
        double precision, dimension(3,3,3,3)::T4th, T4thp, T4thpp  
        double precision, dimension(11,38):: import2
        double precision, dimension(3,99):: import6
        double precision, dimension(5,117):: import7        
        integer:: n1,n2,i
        
        !data_TT07.csv=Test_kmt1.csv   !Matrix with integers
        n1=size(import7,1)
        n2=size(import7,2)
        file_name="data_TT07.csv"
        call read2 (import7, n1,n2, file_name)
        do i=1,1
           print*, "Test contrac_4th_2nd", i
           T_2nd_6=reshape ((import7(i,1:36)),shape(T_2nd_6), 
     &      order=(/1,2/))
           !call print_matrix (T_2nd_6, 6,6) 
           call voig_2_T4th(T_2nd_6,T4th)
           !call print_4D_tensor(T4th,3)
           T4thp=reshape ((import7(i,37:117)),shape(T4th), 
     &          order=(/1,2,3,4/))
           print*, "Reference"
           !call print_4D_tensor(T4thp,3) 
!           @assertEqual(T4thp,T4th,'voig_2_T4th(1)')
        end do

        n1=size(import6,1)
        n2=size(import6,2)
        file_name="data_TT06.csv"
        call read2 (import6, n1,n2, file_name)
        do i=1,3
           print*, "Test contrac_4th_2nd", i
           T4th=reshape ((import6(i,1:81)),shape(T4th), 
     &          order=(/1,2,3,4/))
           !call print_4D_tensor(T4th,3)
           print*, "input matrix"
           m2=reshape ((import6(i,82:90)),shape(m2), order=(/1,2/))
           call print_matrix (m2, 3,3) 
           m2p=reshape ((import6(i,91:99)),shape(m2p), order=(/1,2/))
           print*, "result reference"
           call print_matrix (m2p, 3,3)
           print*, "result my function"
           m2pp=contrac_4th_2nd(T4th,m2)
           call print_matrix(m2pp,3,3)
           !@assertEqual(m2p,contrac_4th_2nd(T4th,m2), 'contrac_4th_2nd(1)')
           call emp_4th(T4th)
           !call print_4D_tensor(T4th,3)
           call emp_2D(m2pp)  
           call print_matrix(m2pp,3,3)
        end do 
        
      end program

      