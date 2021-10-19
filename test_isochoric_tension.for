      program test_isoc_ten
        use tensor_operations

        implicit none
        double precision :: te_n1
        double precision :: n_ampl=2.0; !2      ! n_ampl>1
        double precision ::  xE, xnu, sigma_y0, kappa,q_el
        double precision :: mu  !relation for the conversion of load to strain 
        double precision, dimension(13) :: mat_param
        integer :: stat !Variables t print out information
        integer :: ltype, type_linespace !Selector to compute C and linespace respectively
        integer :: ttype !Selector to compute C
        integer :: n_step !number of steps
        integer,dimension(5) :: posi_last_time=(/2,3,3,4,7/)
        integer :: i !iterators
        integer, dimension(3):: options
        double precision, dimension (7) :: t !Time range
        double precision :: dt!step size for strain ramp
        double precision, dimension (3) :: lam ! Strain range
        double precision, dimension(:), allocatable:: time
        double precision, dimension(:), allocatable:: e11 !11 strain
        double precision, dimension (:,:), allocatable:: sdv
        double precision, dimension (21) ::sdvup
        double precision, dimension (:), allocatable::s11, s22, s33 !normal stresses
        double precision, dimension (6):: e6, De6,sigma !Total strain and total stress vectors 
        double precision, dimension(6,6) :: A66
        !double precision, dimension(10):: e11 

c    Define material parameters

        mat_param(1)   = 210.0e3              ! xE
        mat_param(2)   = 0.33                 ! xnu 
        mat_param(3)   = 200.0 !200           ! xsigy0
        mat_param(4)   = 50.0*mat_param(3);   ! 50   ! xH
        mat_param(5)   = 100.0*mat_param(3);  ! 10  ! xh
        mat_param(6)   = 1.5                  ! q1 
        mat_param(7)   = 1.0                  ! q2
        mat_param(8)   = 1.5                  ! q1=q3 Aricle G.Vadillo
        mat_param(9)   = 0.004                ! f_0
        mat_param(10)  = 0.1                  ! f_n
        mat_param(11)  = 0.3                  ! s_n
        mat_param(12)  = 0.2025               ! f_f
        mat_param(13)  = 0.1                  ! E_n

        xE       = mat_param(1)
        xnu      = mat_param(2) 
        sigma_y0 = mat_param(3)
        mu       = xE/(2.0*(1.0+xnu))
        kappa    = xE/(3.0*(1.0-2.0*xnu))
        q_el     =-0.5*(kappa-2.0/3.0*mu)/(kappa+1.0/3.0*mu)

        options(1) = 0 ! ATS selectror 
        options(2) = 0 ! Stress Case
        options(3) = 1 ! Material Law

c     define loading
c     1: linear ramping of load
        ltype=1
        
        te_n1=0.0
        if (ltype==1) then 
            t(1:2)=(/0.0, 10.0/) ! 10
c            %This is like a load which is related with sigma_y0 and an integer n_ampl 
            lam(1:2)=(/te_n1 , (n_ampl*sigma_y0/2/mu)/1.0 /) !3.0 !total Strain considering the value of mu to onvert from load to strain
        elseif (ltype==2) then
            t(1:3)=(/0.0, 5.0, 10.0/)
            lam(1:2)=(/0.0, 1.59155e-3/3.1 /)
        elseif (ltype==3) then
            t(1:3)=(/0.0, 5.0, 10.0/)
            lam=(/0.0, 0.005, 0.001/)
        elseif (ltype==4) then
            t(1:4)=(/0.0, 2.5, 7.5, 10.0/)
            lam=(/te_n1,  n_ampl*sigma_y0/2.0/mu/(1.0-q_el),  -n_ampl*   !
     $           sigma_y0/2.0/mu/(1.0-q_el) /)
        elseif (ltype==5) then
            t=(/0.0, 2.5, 7.5, 12.5, 17.5, 22.5, 25.0/)
            lam=(/0.0, 0.005, -0.005/)
        end if

c       computation of tangent moduli
        
        dt=0.1
c       prescribed load/time step
c       start and end-time of loading, time-scale, no. of steps
        n_step=nint((t(posi_last_time(ltype))-t(1))/dt)
        type_linespace=1

        allocate (time(n_step+1)) 
        call linespace(t(1),t(posi_last_time(ltype)),time,
     &      type_linespace)

c       Compute the strain load ramp
        allocate (e11(n_step+1))
        !call loading1 (ltype,posi_last_time ,dt, t, lam, e11)
        call loading2 (ltype,posi_last_time ,dt, t, lam, e11)
        !print*, "ramp", e11
        

c        initialize internal variables %Here we are going to store all the history, a row for each step
c        The ISV are stored in a row (vector) due to the symmetry
        allocate (sdv(21,n_step+2))
        sdv(20,1)=0.004        
        print*,  "initial sdv", sdv(:,1)
c       initialise quantities for post-processing
c       3 vectors s11, s22, s33 (normal stress). 
c       Here we are going to store all the history of the normal stresses
        allocate (s11 (n_step+2))
        allocate (s22 (n_step+2))
        allocate (s33 (n_step+2))

        ttype = 0 ! 0: analytical 1: numerical
        i=1
        do i=1,n_step 
        !do i=1,2
        
c       Epsilon (strain) in each step is loaded with the value of the ramp 
c       for the respective step
            e6=(/1.0, -0.5, -0.5, 0.0, 0.0, 0.0/)*e11(i);
            De6=(/1.0, -0.5, -0.5, 0.0, 0.0, 0.0/)*(e11(i)-e11(i-1))
            
c        constitutive law: algorithmic stresses and moduli 
        
c        outputs
c        s: sigma_n+1 (stress n+1 as a vector)
c        A: algritmic tangent (4 grade tensor as a 6x6 matrix)
c        sdvup: Accumulated Internal state variables (the updated 3 ISV of hte imput as a vector)
c       
c        Inputs
c        epsilon: ( current total strain in the form of a vector, latter is trasformed to matrix form)
c        sdv: internal state variables (in the form of a vector, contains 3 ISV`(plastic strain,
c        tensorial internal variable alpha, scalar hardening variable alpha) 
c        Latter is trasformed to matrix form
c        ttype: flag for the tangent moduli 
c               ttype = 0, Analitical tangent moduli computation
c               ttype = 1, numerical tangent moduli computation 
         call kGTN (e6,De6,sdv(:,i),options, mat_param, sigma, A66, 
     &        sdvup)
         !print*,"readed",sdv(:,i)       
         !call VM (e6,sdv(:,i),ttype, mat_param, sigma, A66, sdvup)
         !print*,"Algrorithmic tangent"
         !call print_matrix(A66,6,6)
         !print*, "initial internal state variables 2", sdv
         !print*, "iteration ", i," of ", n_step
c        update of internal variables after obtaining convergence
         !print*,"Resul sdv",sdvup
         sdv(:,i+1) = sdvup
         !print*, "sdv"
         !call print_matrix(sdv, 13, n_step+2)
c        store quantities for post-processing.
c        We are only interested in the normal stresses 
         s11(i+1)=sigma(1) 
         s22(i+1)=sigma(2)
         s33(i+1)=sigma(3);
        end do
        print*, "Exit iterations"
        !call print_matrix(SDV,10)
        !print*, sdv
        !print*, s11
        open(unit=1, iostat=stat, file='data1.csv', 
     $         status='old')
             if (stat == 0) close(1, status='delete')

        ! output data into a file as follow: 
             !e11 = Driving strain  
             !s11 = stress 11
             !s22 = stress 22
             !s33 = stress 33 
             !sdv (1)= plastic strain 11
             !sdv (2)= plastic strain 22
             !sdv (3)= plastic strain 33
             !sdv (20)= porosity
             !sdv (21)= microscopy equivalent plastic strain
        open(1, file = 'data1.csv', status = 'new')  
        do i=1,n_step+1  
           write(1,*) e11(i), s11(i), s22(i), s33 (i),
     $                 sdv(1,i),sdv(2,i),sdv(3,i), sdv(20,i),sdv(21,i)   
        end do  
        close(1)
        deallocate (e11)
        deallocate (sdv)
        deallocate (time) !vector with the steps
      end program test_isoc_ten

