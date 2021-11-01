c=======================================================================  
c Program to call the  GTN material law using different material   
c parameters of VM and GTN to compare its behaviour in a isochoric tension  
c-----------------------------------------------------------------------
c Version 0.9.1
c Oct 2021
c=======================================================================      
      program test_isoc_ten
        use tensor_operations

        implicit none
        double precision :: te_n1
        double precision :: n_ampl=2.0; !2      ! n_ampl>1
        double precision ::  xE, xnu, sigma_y0, kappa,q_el
        double precision :: mu  !relation for the conversion of load to strain 
        double precision, dimension(15) :: mat_param
        integer :: stat !Variables t print out information
        integer :: ltype, type_linespace !Selector to compute C and linespace respectively
        integer :: ttype !Selector to compute C
        integer :: n_step !number of steps
        integer,dimension(5) :: posi_last_time=(/2,3,3,4,7/)
        integer :: i,j !iterators
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
        character:: file_na_ex*27 

c   xE             : Young's modulus    
c   xnu            : Poisson's ratio
c   xsigy0         : initial yield stress
c   xmu            : elastic shear modulus       
c   xk             : bulk modulus
c   kappa          : bulk modulus K
c   ltype          : type of load
c   te_n1          : temporal variable 1
c   posi_last_time : vector to help to define the ramp size
c   t              : vector with time range for all the ramp options
c   lamp           : load range
c   n1             : number of test (rows input file)
c   n_test         : number of test
c   n2             : number of parameters (columnss input file)
c   mat_param      : matrix to store imported parameters for test
c   matL           : auxiliary integer to run different material laws
c   type_linespace : Selector functionality of linespace
c   ttype          : Selector tangent stiffness coputation method 
c   dt             : delta time increment
c   time            :time vector for the ramp 
c   e11            : strain 11 ramp 
c   n_step         : number of steps ramp
c   sdv , sdv_1    : tensors to store the state variables history, 
c                    of material law 0 and 1 respectively, a row for 
c                    each step. The vector form comes already from the 
c                    mterial subroutine
c   sdvup          : updated ISV from the material routine
c   A66            : ATS in 6x6 matrix
c   s11, s22, s33  : vetors to stroe the hystorey of the normal stresses
c   e6             : current total strain tensor defined for an specified time 
c   De6            : standard strain increment between time increment 
c                    defined for this test
c   file_name      : Name of the file to write results,******DON'T CHANGE IT ******* 
c   stat           : Variable for the state of the opened file to write
c   options        : flags for options in the material routine

c    Define material parameters

        mat_param(1)   = 200.0e3              ! xE
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
        mat_param(14)   = 0.25                 !f_c 
        mat_param(15)   = 0.1                  !NU

        xE       = mat_param(1)
        xnu      = mat_param(2) 
        sigma_y0 = mat_param(3)
        mu       = xE/(2.0*(1.0+xnu))
        kappa    = xE/(3.0*(1.0-2.0*xnu))
        q_el     =-0.5*(kappa-2.0/3.0*mu)/(kappa+1.0/3.0*mu)

c     define loading
c     1: linear ramping of load
        ltype=1
        
        te_n1=0.0
        if (ltype==1) then 
            t(1:2)=(/0.0, 10.0/) ! 10
            lam(1:2)=(/te_n1 , (n_ampl*sigma_y0/2/mu)/1.7/) !0.2!3.0-2.9 2D, 2.7-2.6 Plane stress 
        elseif (ltype==2) then
            t(1:3)=(/0.0, 5.0, 10.0/)
            lam(1:2)=(/0.0, 1.59155e-3/2.9 /)

        end if
        
        dt=0.2
c       prescribed load/time step
c       start and end-time of loading, time-scale, no. of steps
        n_step=nint((t(posi_last_time(ltype))-t(1))/dt)
        type_linespace=1

        allocate (time(n_step+1)) 
        call linespace(t(1),t(posi_last_time(ltype)),time,
     &      type_linespace)

c       Compute the strain load ramp
        allocate (e11(n_step+1))
        call loading2 (ltype,posi_last_time ,dt, t, lam, e11)
        print*, "ramp", e11
        

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
        options(1) = 0 ! ATS selectror 
        options(2) = 0 ! Stress Case
        options(3) = 1 ! Material Law
        do i=1,n_step 
        j=i+1
c       Epsilon (strain) in each step is loaded with the value of the ramp 
c       for the respective step
            if (options(2)==0) then            
              e6=(/1.0, -0.5, -0.5, 0.0, 0.0, 0.0/)*e11(j)
              De6=(/1.0, -0.5, -0.5, 0.0, 0.0, 0.0/)*(e11(j)-e11(i))
            elseif (options(2)==1) then
              e6=(/1.0, -0.5, 0.0, 0.0, 0.0, 0.0/)*e11(j)
              De6=(/1.0, -0.5, 0.0, 0.0, 0.0, 0.0/)*(e11(j)-e11(i))
            end if
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

         sdv(:,i+1) = sdvup

c        We are only interested in the normal stresses 
         s11(i+1)=sigma(1) 
         s22(i+1)=sigma(2)
         s33(i+1)=sigma(3);
        end do
        print*, "Exit iterations"

        file_na_ex='5_Isochoric_test/data1.csv'
        open(unit=1, iostat=stat, file=file_na_ex, 
     $         status='old')
             if (stat == 0) close(1, status='delete')

        ! output data into a file as follow: 
             !e11 = Driving strain  
             !s11 = stress 11
             !s22 = stress 22
             !s33 = stress 33 (21 for plane stress)
             !sdv (1)= plastic strain 11
             !sdv (2)= plastic strain 22
             !sdv (3)= plastic strain 33 (21 for plane stress)
             !sdv (20)= porosity
             !sdv (21)= microscopy equivalent plastic strain
        open(1, file = file_na_ex, status = 'new')  
        do i=1,n_step+1  
           write(1,*) e11(i), s11(i), s22(i), s33 (i),
     $                 sdv(1,i),sdv(2,i),sdv(3,i), sdv(20,i),sdv(21,i)   
        end do  
        close(1)
        deallocate (e11)
        deallocate (sdv)
        deallocate (time) 
      end program test_isoc_ten

