Module read 
    !use tensor_operations
    implicit none
    contains
     

    subroutine read2(import1)
        integer:: i, j,end,stat,line_no,n1,n2     !iterators
        double precision, dimension(11,19):: import1
        n1=size(import1,1)
        n2=size(import1,2)   
        open(15, file="diadic2D.csv",access='sequential', form="formatted", iostat=stat)
        do i = 1,n1
            read(15,*,iostat=stat) import1(i,:)
        end do
        close (15)
        print*, import1   
        !call print_matrix (import1, n1,n2)
    end subroutine read1
end module read
