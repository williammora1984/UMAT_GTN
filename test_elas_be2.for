c=======================================================================  
c Program to call the  GTN material law using different material   
c parameters of VM and GTN to compare its behaviour in a isochoric tension  
c-----------------------------------------------------------------------
c In order to call this program in the same folder muste be a csv file of
c at least one file  and maximum 5, with 15 parameters, 3 linear elasticity,
c  2 VM plasticity and 10 fort GTN model. See more information in the 
c documentation.
c-----------------------------------------------------------------------
c Version 0.9.1
c Oct 2021
c=======================================================================
        program test_isoc_ten
        use tensor_operations

        implicit none
        double precision :: te_n1
        double precision :: n_ampl=2.0; !2      ! n_ampl>1
        double precision ::  xE, xnu, sigma_y0, kappa
        double precision :: mu   
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
        double precision, dimension (:,:), allocatable:: sdv, sdv_1
        double precision, dimension (21) ::sdvup, sdvup_1
        double precision, dimension (:), allocatable::s11, s22, s33, !normal stresses
     &       s11_1, s22_1, s33_1
        double precision, dimension (6):: e6, De6, sigma, sigma_1 !Total strain and total stress vectors 
        double precision, dimension(6,6) :: A66, A66_1
        double precision:: mat_param(5,15)
        integer:: n1, n2, n_test, matL
        character:: file_name*13, file_name2*30, file_na_ex*27  !******DON'T CHANGE IT ******* 
c       double precision, dimension(10):: e11 

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
c=======================================================================
c                  Define material parameters
c=======================================================================
c-----------Import data with material parameters for test---------------
        n1=size(mat_param,1)
        n2=size(mat_param,2)
        file_name2="4_Test_lin_elas0/dat_LE_02.csv"
        !file_name ="dat_LE_01.csv"
        call read3 (mat_param, n1,n2, file_name2)
        do matL=0,1
        do n_test=1, n1
                print*, "test", n_test
c-------------Compute elastic constant---------------------------------
        xE       = mat_param(n_test,1)
        xnu      = mat_param(n_test,2) 
        sigma_y0 = mat_param(n_test,3)
        mu       = xE/(2.0*(1.0+xnu))
        kappa    = xE/(3.0*(1.0-2.0*xnu))
c=======================================================================
c                       define loading ramp
c=======================================================================
c----------------Select the time (load) amplitud------------------------
c     1: linear ramp of load
c     2: Linear load and unload toero
        ltype=1
    
        te_n1=0.0
        if (ltype==1) then 
            t(1:2)=(/0.0, 10.0/) ! 10
c            %This is like a load which is related with sigma_y0 and an integer n_ampl 
            lam(1:2)=(/te_n1 , (n_ampl*sigma_y0/2/mu)/2.1 /) !3.0 !total Strain considering the value of mu to onvert from load to strain
        elseif (ltype==2) then
            t(1:3)=(/0.0, 5.0, 10.0/)
            lam(1:2)=(/0.0, 1.59155e-3/3.1 /)
        end if

c-------- computation of number of steps and e11 strain vs time--------
        dt=0.1
        n_step=nint((t(posi_last_time(ltype))-t(1))/dt)
        type_linespace=1
c-------------------Define time ramp------------------------------------ 
        allocate (time(n_step+1)) 
        call linespace(t(1),t(posi_last_time(ltype)),time,
     &      type_linespace)

c--------------Define the strain e11 reference load ramp------------------------------
        allocate (e11(n_step+1))
        call loading2 (ltype,posi_last_time ,dt, t, lam, e11)
        !print*, "ramp", e11
c=======================================================================
c       initiate variables for post-processing
c=======================================================================
c-------initiate the tensor sdv according the number of steps-----------  
        allocate (sdv(21,n_step+2)) ; allocate (sdv_1(21,n_step+2))
        sdv(20,1)=mat_param(n_test,9); sdv_1(20,1)=mat_param(n_test,9)        
        print*,  "initial sdv", sdv(:,1)

c----------Initiate vectors  s11, s22, s33 (normal stress).------------- 
        allocate (s11 (n_step+2)); allocate (s11_1 (n_step+2))
        allocate (s22 (n_step+2)); allocate (s22_1 (n_step+2))
        allocate (s33 (n_step+2)); allocate (s33_1 (n_step+2))

c=======================================================================
c   Start loop to compute material stress behavior due to strain ramp 
c=======================================================================        
        ttype = 0 ! 0: analytical 1: numerical
        do i=1,n_step
c----Define current driven strain at the begining of an iteration------- 
            e6=(/1.0, -0.5, -0.5, 0.0, 0.0, 0.0/)*e11(i);
            De6=(/1.0, -0.5, -0.5, 0.0, 0.0, 0.0/)*(e11(i)-e11(i-1))
c---------Call the material law to compute the stress, ATS, and ISV-----            
c        outputs
c        sigma: stress n+1 as a vector
c        A66: algritmic tangent (4 grade tensor as a 6x6 matrix) at n+1
c        sdvup: Accumulated Internal state variables at n+1 
c       
c        Inputs
c        e      : current total strain
c        De6    : strain increment between time increment      
c        sdv(i) : current internal state variables
c        options: flags for options in the material routine
c        mat_param: material parameters according to the test
         options(1) = 0 ! ATS selectror 
         options(2) = 0 ! Stress Case
         options(3) = matL ! Material Law
         call kGTN (e6,De6,sdv(:,i),options, mat_param(n_test,:), sigma,  
     &        A66, sdvup)

c=======================================================================
c store the updated variables to start the next step and for the 
c post-processing
c=======================================================================
         sdv(:,i+1) = sdvup
 
c        We are only interested in the normal stresses 
         s11(i+1)=sigma(1) 
         s22(i+1)=sigma(2)
         s33(i+1)=sigma(3)
     
        end do

c=======================================================================
c        Export the results to a cs file for post-processing
c=======================================================================
c+++++++++++++++++++++++Material law 1++++++++++++++++++++++++++++++++++
        if (matL==0) then
            file_na_ex='4_Test_lin_elas0/dVM1'//char(48+n_test)//'.csv'
c---------------Erase old files with the same name-----------------------
            open(unit=1, iostat=stat, file=file_na_ex, 
     $             status='old')
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
c------------------------Export the data--------------------------------                 
            open(1, file = file_na_ex, status = 'new')  
            do i=1,n_step+1  
               write(1,*) e11(i), s11(i), s22(i), s33 (i),
     &             sdv(1,i),sdv(2,i),sdv(3,i), sdv(20,i),sdv(21,i)
   
            end do  
            close(1)
        elseif (matL==1) then
c+++++++++++++++++++++++Material law 2++++++++++++++++++++++++++++++++++                
            file_na_ex='4_Test_lin_elas0/dGN1'//char(48+n_test)//'.csv'
c---------------Erase ol files with the same name-----------------------
            open(unit=1, iostat=stat, file=file_na_ex, 
     $             status='old')
                 if (stat == 0) close(1, status='delete')
c------------------------Exprot the data--------------------------------                 
            open(1, file = file_na_ex, status = 'new')  
            do i=1,n_step+1  
               write(1,*) e11(i), s11(i), s22(i), s33 (i),
     &             sdv(1,i),sdv(2,i),sdv(3,i), sdv(20,i),sdv(21,i)   
            end do  
            close(1)
        end if
c-----------------Deallocate loop allocatable variables-----------------        
        deallocate (s11); deallocate (s22); deallocate (s33)
        deallocate (s11_1); deallocate (s22_1); deallocate (s33_1)
        deallocate (e11)
        deallocate (sdv); deallocate (sdv_1)
        deallocate (time) !vector with the steps
       end do
       end do
      end program test_isoc_ten



