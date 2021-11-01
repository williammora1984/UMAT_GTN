c=======================================================================
c                    Progtam Test 1D tension poin materila
c=======================================================================
cThis program calls the material routine to perform in a 1D tension 
ctension test in a material point
c=======================================================================  
c-----------------------------------------------------------------------
c Version 0.9.1
c Oct 2021
c=======================================================================

      program drive_1D_tension
        !subroutine drive
c        use tensor_operations
c        
c        double precision, dimension(3,3):: eps !total strain
c        double precision, dimension (3,3,3,3):: C_ve !Elastic tensor
c        double precision, dimension(3,3):: sig_ve !elastic predictor
c        double precision, dimension(3,3):: eps_ve_t !elastic strain at time t 
c        double precision, dimension(3,3):: D_eps !Delta total strain
c        double precision, dimension(3,3):: sig_t1 !stress t+1
c        double precision, dimension(3,3):: D_eps_p !Delta plastic strain
c        sig_ve=contrac_4th_2nd(c_ve,eps_ve_t+D_eps)
c        sig_t1= sig_ve-contrac_4th_2nd(C_ve,D_eps_p)


        use tensor_operations
        
        implicit none
        double precision :: te_n1,te_n2,te_n3, tol, dt
        double precision :: n_ampl=2.0; !2      ! n_ampl>1
        double precision ::  xE, xnu, sigma_y0, Xmu, kappa, q_el 
        double precision, dimension(15) :: mat_param
        integer,dimension(5) :: posi_last_time =(/2,3,3,4,7/)
        integer :: ltype, ttype, type_linespace, matL, options(3)
        integer :: i, maxit,iter,n_step,stat
        double precision, dimension (7) :: t 
        double precision, dimension (3) :: lam 
        double precision, dimension(:), allocatable:: time, e11, eps22, 
     &       eps33, s11, s22, s33
        double precision, dimension (:,:), allocatable:: sdv
        double precision, dimension (21) ::sdvup
        double precision, dimension (6):: e6, e6n, De6, sigma
        double precision, dimension(6,6) :: A66
        double precision, dimension(5,5) :: Abar, Abar_inv
        double precision, dimension(5) :: epsbar
        double precision, dimension(5) :: sbar=(/1.0,1.0,1.0,1.0,1.0/)
        character:: file_na_ex*40 

c   xE             : Young's modulus    
c   xnu            : Poisson's ratio
c   xsigy0         : initial yield stress
c   xmu            : elastic shear modulus       
c   xk             : bulk modulus
c   kappa          : bulk modulus K
c   ltype          : type of load
c   te_n1,te_n2,te_n3 : temporal variable 1
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
c   time           :time vector for the ramp 
c   e11            : strain 11 ramp 
c   n_step         : number of steps ramp
c   sdv , sdv_1    : tensors to store the state variables history, 
c                    of material law 0 and 1 respectively, a row for 
c                    each step. The vector form comes already from the 
c                    mterial subroutine
c   sdvup          : updated ISV from the material routine
c   A66            : ATS in 6x6 matrix
c   sigma          : current total stress in a vector
c   s11, s22, s33  : vetors to store the history of the normal stresses
c   eps33, eps22   ! Lateral strains
c   e6             : current total strain tensor defined for an specified time
c   e6n            : previous current total strain 
c   De6            : standard strain increment between time increment 
c                    defined for this test
c   
c   file_name      : Name of the file to write results,******DON'T CHANGE IT ******* 
c   stat           : Variable for the state of the opened file to write
c   options        : flags for options in the material routine
c   tol            : Tolerance numerical solution
c   epsbar         : lateral strains
c   sbar           : Lateral strains

c=======================================================================
c                  Define material parameters
c=======================================================================
        mat_param(1)   = 210.0e3              ! xE
        mat_param(2)   = 0.33                 ! xnu 
        mat_param(3)   = 200.0 !200           ! xsigy0
        mat_param(4)   = 50.0*mat_param(3);   !50   ! xH
        mat_param(5)   = 100.0*mat_param(3);  !10  ! xh
        mat_param(6)   = 1.5                  !q1 
        mat_param(7)   = 1.0                  !q2
        mat_param(8)   = 1.5                  !q1=q3 Aricle G.Vadillo
        mat_param(9)   = 0.004                !f_0
        mat_param(10)  = 0.1                  !f_n
        mat_param(11)  = 0.3                  !s_n
        mat_param(12)  = 0.2025               !f_f
        mat_param(13)  = 0.1                  !E_n
        mat_param(14)   = 0.25                 !f_c 
        mat_param(15)   = 0.1                  !NU

        xE = mat_param(1)
        xnu = mat_param(2) 
        sigma_y0 = mat_param(3)
        xmu = xE/(2.0*(1.0+xnu))

c       strain amplitude in terms of multiple of normalized yield stress
c       strain_ampl=sigma_y0/2/mu*n_ampl/(1-q_el)
        kappa = xE/(3.0*(1.0-2.0*xnu))
        q_el=-0.5*(kappa-2.0/3.0*xmu)/(kappa+1.0/3.0*xmu)

c=======================================================================
c                       define loading ramp
c=======================================================================
c 1: linear ramping of load
c 2: half cycle with linear load change (half load and half unload)
c 3: linear loading and unloading, start and end point different
c 4: full cycle with linear load change
c 5: two cycles

        te_n1=0.0
        ltype=2; !4
        if (ltype==1) then
            t(1:2)=(/0.0, 10.0/)
            lam(1:2)=(/0.0, 0.005/1.0 /)  !3.0
        elseif (ltype==2) then
            t(1:3)=(/0.0, 5.0, 10.0/)
            lam(1:2)=(/0.0, 1.59155e-3/)
        end if

c-------- computation of number of steps and e11 strain vs time--------
        dt=0.1;        
        n_step=nint((t(posi_last_time(ltype))-t(1))/dt)
        type_linespace=1
        
        allocate (time(n_step+1)) 
        call linespace(t(1),t(posi_last_time(ltype)),time,
     &      type_linespace)
        
c----------------Define the strain e11 load ramp------------------------
        allocate (e11(n_step+1))
        call loading2 (ltype,posi_last_time ,dt, t, lam, e11)

c=======================================================================
c       initiate variables for post-processing
c=======================================================================
c-------initiate the tensor sdv according the number of steps-----------
        allocate (sdv(21,n_step+2))  !initial internal variables
c----------Initiate vectors  s11 (stress 1D tension test).------------- 
         allocate (s11 (n_step+2)) 
c---------------vectors for the normal strains-------------------------
         allocate (eps22(n_step+2)) 
         allocate (eps33(n_step+2))
c=======================================================================
c   Start loop to compute material behavior due to strain ramp 
c======================================================================= 
c--Define tolerance and maximum no. of iterations for Newton iteration--
         tol=1e-10 !1e-10
         maxit=200
         ttype = 0 ! 0: analytical, 1: numerical tangent moduli computation

         i=1
        !print*, "Enter material model 1" 
         do i=1,n_step

c=======================================================================
c      Start loop to compute material stress and strain in a step using 
c                  a newton methud until convergence
c======================================================================= 
            sbar=(/1,1,1,1,1/)
            iter=0
c------------Store strain old values to find Deps-----------------------              
            e6n(1)=e11(i)
            e6n(2:6)=epsbar
            do iter=1,maxit-1

                print*,"Iteration",iter," of maximum",maxit-1,  
     &            "in the step", i, "of", n_step
                if (iter > maxit) then
                  print*, "error"  
c                  error(['No convergence after ', iter,
c     $                   ' global iterations'])
                end if
c               % 1.) total deformation
c------Update the current total strain  composed by 1, 2 ---------------
c-----------1 = e11 prescribed by ramp----------------------------------
c-----------2 = Solution of lateral strains using the ATS---------------
          !      print*, "enter program uniaxial test 1"
                e6(1) = e11(i+1) !Here I assign the second value of the controlled strain 
                e6(2:6) = epsbar !At the first iteration this part of the strain is equal to the previous, then in the second iteration change.
                De6=e6-e6n
c---------Call the material law to compute the stress, ATS, and ISV-----
c                call kGTN (e6,sdv(:,i),ttype, mat_param, sigma, A66,
c     $                      sdvup)
                options(1) = 0 ! ATS selectror 
                options(2) = 0 ! Stress Case
                options(3) = 0 ! Material Law
                call kGTN (e6,De6,sdv(:,i),options, mat_param, sigma,  
     &              A66, sdvup)     

c----Take the stress and ATS corresponding to the lateral variables-----
                sbar=partition_1D(sigma)
                !call print_matrix(A66,6)
                Abar=partition_2D_2D(A66)                
c-------------update of lateral strains----------------------------
                Abar_inv=inv_T_2D(dble(Abar))

                epsbar=epsbar-matmul(Abar_inv,sbar)
                !sdv(:,i)=sdvup 
            end do
c=======================================================================
c       Store the updated variables after convergence to start the  
c                 next step and dor the post-processing
c=======================================================================
            sdv(:,i+1) = sdvup
c        store quantities for post-processing.
            s11(i+1)=sigma(1) 
            eps22(i+1)=e6(2)
            eps33(i+1)=e6(3)
         end do
c=======================================================================
c        Export the results to a cs file for post-processing
c=======================================================================
c---------------Erase old files with the same name----------------------
         file_na_ex='6_1D_Stra_dri/data_str_drive.csv'
          open(unit=1, iostat=stat, file=file_na_ex, 
     $         status='old')
          if (stat == 0) close(1, status='delete')
c      ! output data into a file as follow: 
             !time = time ramp
             !e11 = Driving strain  
             !s11 = stress 11
             !eps22 = strain 22
             !eps33 = strain 33 
             !sdv (1)= plastic strain 11
             !sdv (2)= plastic strain 22
             !sdv (3)= plastic strain 33
             !sdv (20)= porosity
             !sdv (21)= microscopy equivalent plastic strain
c------------------------Export the data--------------------------------
         open(1, file = file_na_ex, status = 'new')  
         do i=1,n_step+1  
            write(1,*) time(i), e11(i), S11(i), eps22(i), eps33 (i),
     &                 sdv(1,i),sdv(2,i),sdv(3,i)
c     &                 sdv(20,i),sdv(21,i)   
         end do  
         close(1)
c-----------------Deallocate loop allocatable variables-----------------          
         deallocate (S11) 
         deallocate (eps22) 
         deallocate (eps33)
         deallocate (e11)
         deallocate (sdv)
         deallocate (time)
         print*, "program end, review the exported files"

      end program drive_1D_tension


