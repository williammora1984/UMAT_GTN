c======================================================================
c    UMAT_for the Gurson, Tvergaard and Needleman (GTN) damage model
c             coded by  Wiliam Fernando Mora Pulido   
c=======================================================================
c  This is a user defined subroutine UMAT for Abaqus, is and inter-
c face that calls the material subroutine kGTN based on the on the 
c algorithm proposed by N.Aravas and Z.L.Zhang
c     
c-----------------------------------------------------------------------
c Version 0.9.2
c coded by: W. Mora Oct 2021
c-----------------------------------------------------------------------

c      NDI         : number of direct components of DDSDDE, DDSDDT, and DRPLDE
c      NSHR        : number of engineering shear components of DDSDDE, DDSDDT, and DRPLDE
c      NTENS       : NDI + NSHR: Size of the stress or strain component array
c      DDSDDE      : Algorithmic tangent stiffness (Jacobian-stiffness matrix) 
c      stran       : strain from previous increment
c      stress      : stress  
c      Dstran      : increment of strain for current iteration
c
c      PROPS(NPROPS): Array with material property data
c      NPROPS       : number of material properties.
c                     3D     Plane stress 
       !ntens  =      6          3
       !nprops =      15         15
       !nstatv =      3          3
        

c======================================================================
c                    Standard declaration UMAT
c=======================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
c======================================================================
c                    Standard variables UMAT
c=======================================================================
       !implicit none
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3) 
C
      DIMENSION EELAS(6),EPLAS(6),FLOW(6)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
c======================================================================
c              Additional variables used in this UMAT file
c=======================================================================
      double precision:: Tr, Ainv(3,3)
      double precision:: sdv(21),epspr(6), epser(6), eps6(6), Deps6(6), 
     &  ATS66(6,6), sig6(6)
      integer:: ttype,stat,options(3)  
      double precision,dimension(21) ::sdvup
c      WRITE(6,"(/'DDSDDE')") 

C -----------------------------------------------------------

c      sdv   : internal state variables input to material routine
c      options(1) : type of tangent stiffnes 
c              0: analytical 
c              1: numerical
c      options(2): stress_case
c              0: 3D 
c              1: Plane stress
c      options(3): Material law
c              0: Von Mises
c              1: GTN

c      sdvup : updated internal state variables from the material routine
c      epspr : rotated platic strain
c      epser : rotated elastic strain
c--------Variables of standard size to manipulte the information-------- 
c------------in the change from 3D to plain stress and interact--------- 
c----------------------with the material routine------------------------
c      eps6  : previous strain 
c      Deps6 : delta strain 
c      sig6       : stress
c      ATS66(6,6) : Algorithmic tangent stiffness
c------------------MATERIAL PARAMETERS----------------------------------

c       props(1)   ! : 210.0e3            ! xE
c       props(2)   ! : 0.33               ! xnu 
c       props(3)   ! : 200.0 !200         ! xsigy0
c       props(4)   ! : 50.0*mat_param(3); !50   ! xH
c       props(5)   ! : 100.0*mat_param(3); !10  ! xh
c       props(6)   ! : 1.5                 !q1 
c       props(7)   ! : 1.0                 !q2
c       props(8)   ! : 1.5                 !q1=q3 Aricle G.Vadillo
c       props(9)   ! : 0.004                !f_0
c       props(10)  ! :  0.1                 !f_n
c       props(11)  ! :  0.3                 !s_n
c       props(12)  ! :  0.2025              !f_F
c       props(13)  ! : 0.1                  !E_n
c       props(14)  ! : 0.1                  !f_c

C-----------------INTERNAL STATE VARIABLES------------------------------
c       statev(1:6)    : plastic strain
c       statev(7:12)   : Back stress
c       statev(13)     : Drag stress, alpha,
c       statev(14:19)  : Elastic strain
c       statev(20)     : void fraction 
c       statev(21)     : Microstrain eps_b  

       options(1) = 0
       options(3) = 0
c ---Assign 3D or plane stress program
       if (NTENS==6) then
          options(2) = 0
       elseif (NTENS==3) then
          options(2) = 1
       else  
         WRITE(6,1)
 1       FORMAT(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',
     1          'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
         
       end if

c Get the updated rotated elastic and plastic strain
       !The ISV vectors does not change with the streess case 
       CALL ROTSIG(STATEV(14), DROT,epser, 2,NDI,NSHR)
       CALL ROTSIG(STATEV(1),DROT,epspr,2,NDI,NSHR)
       

      sdv(1:13)  = statev(1:13)
      sdv(1:6)   = epspr
      sdv(14:19) = epser
    
      if (options(2)==0) then
        eps6=STRAN
        Deps6=DSTRAN  
      elseif (options(2)==1) then
        eps6(1:3)=STRAN
        Deps6(1:3)=DSTRAN
      end if 
 3       FORMAT(//,30X,'***Subroutine von Mises ',
     3          'Test UMAT')
 4       FORMAT(//,30X,'***Subroutine GTN ',
     4          'Test UMAT')
       WRITE(6,*) "Hello subroutine"
       WRITE(7,*) "Hello subroutine4"
       if (options(3)==0) then
          WRITE(6,3)
          WRITE(7,3)
       elseif (options(3)==1) then
          WRITE(6,4)
          WRITE(7,4)
       end if

      call kGTN (eps6+Deps6,Deps6,sdv,options, props, sig6,  
     &         ATS66, sdvup)
      statev=sdvup

      if (options(2)==0) then
        STRESS=sig6
        DDSDDE=ATS66
  
      elseif (options(2)==1) then
        STRESS=sig6(1:3)
        DDSDDE=ATS66(1:3,1:3)
      end if      


      END SUBROUTINE UMAT


c=======================================================================  
c subroutine GTN : Material model according to GURSON, TVERGAARD AND 
c NEEDLEMAN (GTN) DAMAGE MODEL  
c-----------------------------------------------------------------------
c  This subroutine implement the nmerical solution for the GURSON, 
c TVERGAARD AND NEEDLEMAN (GTN) DAMAGE MODEL according to the N.Aravas
c  article.
c-----------------------------------------------------------------------
c Version 0.9.2
c Oct 2021
c-----------------------------------------------------------------------

c Inputs
c   eps6 :epsilon ( current total strain in the form of a vector, latter
c         is trasformed to matrix form, the vector comes from the history
c          of the properties for the material)
c   sdvl : internal state variables (in the form of a vector, contains 3 
c          ISV`(plastic strain, tensorial internal variable alpha, 
c          scalar hardening variable alpha) 
c   Options
c      options(1) : type of tangent stiffnes computation 
c              0: analytical 
c              1: numerical
c      options(2): stress_case
c              0: 3D 
c              1: Plane stress
c      options(3): Material law
c              0: Von Mises
c              1: GTN
c  inputmat: vector with the material's intrinsic properties. 
c-----------------------------------------------------------------
c Outputs

c   sig6:  sigma_n+1 (stress at n+1 as a vector)
c   A66:   algritmic tangent stiffness (ATS) tensor(4 grade tensor as a 
c          6x6 tensor)
c   sdvup: Updated Internal state variables (the updated 3 ISV of the s
c          dvl input vector)
c==================================================================


      subroutine kGTN (eps6,D_eps6,sdvl,options,inputmat,sig6,A66,sdvup)
      !eps,epsn,epsen
      INCLUDE 'ABA_PARAM.INC'       
       !use tensor_operations
    
       !implicit none
c         material parameters
        external plas_corre_GTN
        double precision, dimension(6) :: eps6, D_eps6, epsn6
        double precision, dimension(3,3) :: eps, D_eps, epsn, sig
c             double precision, dimension(:),intent(in) :: sdvl
        double precision, dimension(21) :: sdvl 
        double precision, dimension(21) :: sdvup
        double precision, dimension(15):: inputmat
        double precision :: xE     ! Young's modulus
        double precision :: xnu    ! Poisson's ratio
        double precision :: xsigy0 ! initial yield stress
        double precision :: xHk    ! kinematic hardening modulus
        double precision :: xhi     ! isotropic hardening modulus
        double precision :: xmu    !shear modulus
        double precision :: xk     ! bulk modulus
        double precision :: q1, q2, q3, f_0, f_n, s_n, E_n, f_c, NU !Parameters GTN model
        double precision, dimension (3,3) :: epspn, epsp, epse,epsen   !plastic and elastic strain at t0 and t1
        double precision, dimension (3,3):: Balphan, Balpha ! (B=bold= Tensor)strain-like ISV that thermodynamically conjugates to the kinematic hardening at time 0 (input) and t+1 (output)
        double precision :: alphan, alpha, fn, f, epsp_b, epsp_b_n !(scalar) strain-like ISV that thermodynamically conjugates to the Isotropic hardening at time t (input) and t+1 (output)
        double precision, dimension(6,6) :: A66
        double precision, dimension(6) :: sig6
        double precision, dimension (3,3)::sigtr, dsigtr !deviatoric part of the trial stress tensor
        double precision, dimension (3,3)::dsig !deviatoric part of the stress tensor
        double precision, dimension (3,3)::d_epsp  !<----New
        double precision, dimension(3,3,3,3):: P4sym_r !Fouth order identity tensor and deviatoriser tensor 
        double precision, dimension(3,3,3,3)::C, Cdev, C_e  !Elasticity tensor
        double precision , dimension (3,3) :: ntr, nxi_a !flow direction from trial state
        double precision :: nxitr !norm of xitr (norm of relative stress tensor)
        double precision :: phitr,q_tr,p_tr !elastic predictor (or trial stress) 
        double precision :: Beta1, Beta2 !Term for the ATS
        double precision :: gamma !incremental plastic multiplier

c-------------------Auxiliary variables---------------------------------
        double precision :: tol =1e-8
        integer, dimension (6) :: ii = (/1,2,3,1,2,1/) !Auxiliar iterators
        integer, dimension (6) :: jj = (/1,2,3,2,3,3/) !Auxiliar iterators
        integer, dimension (3) ::ii_PS = (/1,2,2/)
        integer, dimension (3) ::jj_PS = (/1,2,1/)
        integer ::i, j, ATStype, Mat_L, str_case, rec,pru !Iterators and selectors
        integer, dimension (3) :: options
        double precision, dimension (3,3):: xid !Identity
        real:: e, t1(2), t2(2)
        double precision:: time_su (5,2)
c  Especific variables according to the article N. ARAVAS, 1987
       double precision :: q, p, sig_0, d_sig_0_d_epsp_b
       double precision, dimension (3,3) :: D_eps3
c        s_ve_33
       double precision ::D_eps_p, D_eps_q  !p: hydrotatyx stress and equivalent stress
       double precision:: dg_dp, dg_dq, d2g_d2p, d2g_d2q, d2g_dpdq 

c=======================================================================
c                   Retrieve material parameters
c=======================================================================
        xE         = inputmat(1)  !Young's modulus    
        xnu        = inputmat(2)  !Poison    
        xsigy0     = inputmat(3)  !Yield point     
        xHk        = inputmat(4) !Kinematic hardeniing parameter      
        xhi        = inputmat(5) !Isotropic hardeniing parameter 
        q1         = inputmat(6)
        q2         = inputmat(7)
        q3         = inputmat(8)
        f_0        = inputmat(9)
        f_n        = inputmat(10)
        s_n        = inputmat(11)
        E_n        = inputmat(13)
                !

        xmu  = xE/(2.0*(1.0+xnu))    !mu or G: Elastic Shear moduli 
        xk   = xE/(3.0*(1.0-2.0*xnu))  !k: Compression (Bulk) Moduli


c=======================================================================
c             Retrieve option parameters
c=======================================================================
        ATStype=options(1) !type of tangent stiffnes computation
        str_case=options(2) !stress_case
        Mat_L=options(3) !Material law

c=======================================================================
c      Retrieve strain drive of the current step and last step 
c=======================================================================
        if (str_case==0) then
c=======================================================================
c                            3D 
c=======================================================================
c-----Current strain drive input to compute stresses, ATS and ISV-------
          eps = reshape((/eps6(1),      eps6(4)/2.0,  eps6(6)/2.0,
     $                    eps6(4)/2.0,  eps6(2),      eps6(5)/2.0,
     $                    eps6(6)/2.0,  eps6(5)/2.0,  eps6(3)/),
     $                    shape(eps), order=(/2,1/))
c---------------Strain drive input of the previous step-----------------
c    <----------------------------------------------------------------------Change to implement and review               
          epsn6=eps6-D_eps6
          epsn = reshape((/epsn6(1),      epsn6(4)/2.0,  epsn6(6)/2.0,
     $                     epsn6(4)/2.0,  epsn6(2),      epsn6(5)/2.0,
     $                     epsn6(6)/2.0,  epsn6(5)/2.0,  epsn6(3)/),
     $                     shape(eps), order=(/2,1/))
          D_eps=eps-epsn
c    <----------------------------------------------------------------------Change until here to implement and review
c--------------------plastic strain at t0------------------------------- 
        !print*,"before copy in function",sdvl
        epspn   = reshape((/sdvl(1),      sdvl(4)/2.0,  sdvl(6)/2.0,
     $                      sdvl(4)/2.0,  sdvl(2),      sdvl(5)/2.0,
     $                      sdvl(6)/2.0,  sdvl(5)/2.0,  sdvl(3)/),
     $                      shape(epspn), order=(/2,1/))

c-------------------------Elastic strain--------------------------------        
        epsen  = reshape((/sdvl(1+13),    sdvl(4+13)/2.0,sdvl(6+13)/2.0,
     $                     sdvl(4+13)/2.0,sdvl(2+13),    sdvl(5+13)/2.0,
     $                     sdvl(6+13)/2.0,sdvl(5+13)/2.0,sdvl(3+13)/),
     $                     shape(epspn), order=(/2,1/))

        elseif (str_case==1) then
c=======================================================================
c                      PLANE STRESS 
c=======================================================================
c-----Current strain drive input to compute stresses, ATS and ISV-------
          eps(1,1) = eps6(1)
          eps(2,2) = eps6(2)
          eps(1,2) = eps6(3)/2.0
          eps(2,1) = eps6(3)/2.0

c---------------Strain drive input of the previous step-----------------
c    <----------------------------------------------------------------------Change to implement and review               
          epsn6=eps6-D_eps6
          epsn(1,1) = epsn6(1)
          epsn(2,2) = epsn6(2)
          epsn(1,2) = epsn6(3)/2.0
          epsn(2,1) = epsn6(3)/2.0
          epsn(1,3) = 0.0; epsn(2,3) = 0.0; epsn(3,1) = 0.0 
          epsn(3,2) = 0.0; epsn(3,3) = 0.0
c          epsn(3,1:3,3)=(/0.0, 0.0 ,0.0/)    

          D_eps=eps-epsn
c--------------------plastic strain at t0-------------------------------
          epspn(1,1) = sdvl(1)
          epspn(2,2) = sdvl(2)
          epspn(1,2) = sdvl(3)/2.0
          epspn(2,1) = sdvl(3)/2.0
          epspn(1,3) = 0.0; epspn(2,3) = 0.0; epspn(3,1) = 0.0 
          epspn(3,2) = 0.0; epspn(3,3) = 0.0
c-------------------------Elastic strain--------------------------------
          epsen(1,1) = sdvl(1+13)
          epsen(2,2) = sdvl(2+13)
          epsen(1,2) = sdvl(3+13)/2.0
          epsen(2,1) = sdvl(3+13)/2.0
          epsen(1,3) = 0.0; epsen(2,3) = 0.0; epsen(3,1) = 0.0 
          epsen(3,2) = 0.0; epsen(3,3) = 0.0    
        end if
c=======================================================================
c            Retrieve Internal State variables (last step)
c=======================================================================
c --------strain-like ISV that thermodynamically conjugates-------------
c---------------to the kinematic hardening at time t--------------------
        Balphan = reshape((/sdvl(1+6),     sdvl(4+6)/2,   sdvl(6+6)/2.0,
     $                      sdvl(4+6)/2.0, sdvl(2+6),     sdvl(5+6)/2.0,
     $                      sdvl(6+6)/2.0, sdvl(5+6)/2.0, sdvl(3+6)/),
     $                      shape(Balphan), order=(/2,1/))
c--------strain-like ISV that thermodynamically conjugates--------------
c---------------- to the sotropic hardening at time t-------------------
        alphan = sdvl(13)
c-----------------------Void volume fraction----------------------------        
        fn  = sdvl(20)
c--------------microscopic equivalent plastic strain--------------------
        epsp_b_n= sdvl(21)
c=======================================================================
c     Start computations 
c=======================================================================

c------------Fill the identity and deviatoriser tensor------------------               
        call Ident1(xid,3)
        call P4sym(P4sym_r)
        call Iso_hard1(epsp_b_n, sig_0, d_sig_0_d_epsp_b, inputmat)
        
        pru=0
c=======================================================================
c                      Elastic prediction
c=======================================================================
       call dtime(t1,e)         !  Startup etime - do not use result

        if (Mat_L==0) then
              call trial_step (eps,epspn,xmu,xHk,Balphan,alphan,xhi,
     &                   xsigy0, phitr, dsigtr, ntr, nxitr)       
        else
              call trial_step_GTN (eps,epsn,epsen,xmu,xk,phitr, sig_0,
     &               ntr,q1,q2,q3,fn,dsigtr,sigtr, q_tr,p_tr)
        end if
        call dtime(t1,e)
        time_su(1,:)=t1
 
c      Evalate if it is required a elastic or elastic-plastic step
        if (phitr<tol) then  !Tol is a number very close to zero
            !print*, "elastic"
            sig=sigtr
            dsig=dsigtr
            Beta1=1.0
            Beta2=1.0
            Balpha=Balphan
            alpha=alphan
            epsp=epspn
            epse=eps
            !epse=epsen+(eps-epsn) !<<<-----------Change to implement
            f=fn
            epsp_b=epsp_b_n
            if (ATStype==0) then
              Cdev=2.0*xmu*P4sym_r      
            endif
            if (Mat_L==0) then  !I change this part here to avoid use more complicate conditionals if I put this below
c-------------------Update stress--------------------------------------- 
c-----------Addition of spherical contribution to stress----------------
              call Trace(eps,3,tra_tem)
              sig=dsig+xk*tra_tem*xid
              !print*, "sig", sig
              !print*, "eps", eps
              epse=eps-epsp
              call diadic_prod_T2_T2(xid,xid,3,3,dia_tem)
              C=Cdev+xk*dia_tem
            end if
            if (Mat_L==1) then
c-------------------Update stress--------------------------------------- 
c-----------Addition of spherical contribution to stress----------------
              epse=eps-epsp
              call diadic_prod_T2_T2(xid,xid,3,3,dia_tem)
              C=Cdev+xk*dia_tem
              call m4th_2_66_sym (C, A66)           
            end if                          

c=======================================================================
c                      plastic correction
c=======================================================================
        call dtime(t2,e)
        else 
c           print '(A)'   
c           print*, "enter to material routine plastic return", phitr
c          Mat_L=0 : VM, Mat_L=1 : GTN
           if (Mat_L==0) then
         
             call plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
     &            nxitr, epsp, epspn, Balpha, Balphan, alpha, alphan,
     &             beta1,beta2 )
             if (ATStype==0) then
                call diadic_prod_T2_T2(ntr,ntr,3,3,dia_tem)
                Cdev=2*xmu*beta1*P4sym_r-2*xmu*beta2*dia_tem
c------------   -------Update stress--------------------------------------- 
c-----------A   ddition of spherical contribution to stress----------------
                 call Trace(eps,3,tra_tem) 
                 sig=dsig+xk*tra_tem*xid
                 epse=eps-epsp
c-------------------Compute the ATS for VMs-----------------------------           
                 call diadic_prod_T2_T2(xid,xid,3,3,dia_tem)
                 C=Cdev+xk*dia_tem
             end if                                 
        
          else
              rec=1
             call plas_corre_GTN (rec, q_tr,p_tr, sigtr,inputmat,fn, 
     &       dsigtr, phitr, epsp_b,eps, D_eps, epsp, epse, Balpha, 
     &       alpha,sig_0,sig,dsig,f, ATStype, str_case, epspn, epsp_b_n,
     &       A66, plas_corre_GTN)      
           end if
        end if
        call dtime(t2,e)
        time_su(2,:)=t2

c=======================================================================         
c      Retrieve condensed form of tensors to export information
c=======================================================================
          if (ATStype==0 .and. Mat_L==0 ) then
c-----------------restore stiffness tensor as matrix-------------------- 
            call m4th_2_66_sym (C, A66)
          end if
c=======================================================================         
c                              3D
c=======================================================================
c-----------------------stress tensor as vector-------------------------
        if (str_case==0) then
           do i=1,6
               sig6(i) = sig(ii(i),jj(i))
           end do
c---------------------- store ISV as vector----------------------------- 
c-------------------------plastic strain--------------------------------  
          sdvup(1:6) =(/epsp(1,1),   epsp(2,2),   epsp(3,3),
     $               2.0*epsp(1,2), 2.0*epsp(2,3), 2.0*epsp(1,3)/)
c---------------------------Back stress---------------------------------        
          sdvup(7:12)=(/Balpha(1,1),   Balpha(2,2),   Balpha(3,3),
     $               2.0*Balpha(1,2), 2.0*Balpha(2,3), 2.0*Balpha(1,3)/)
          sdvup(14:19)=(/epse(1,1),   epse(2,2),   epse(3,3),
     $               2.0*epse(1,2), 2.0*epse(2,3), 2.0*epse(1,3)/)

        elseif (str_case==1) then
c=======================================================================         
c                         Plane stress
c=======================================================================
          do i=1,3
                 sig6(i) = sig(ii_PS(i),jj_PS(i))
          end do
c---------------------- store ISV as vector----------------------------- 
c-------------------------plastic strain--------------------------------
          sdvup(1:3) =(/epsp(1,1),   epsp(2,2),  2.0*epsp(1,2)/)
c---------------------------Back stress---------------------------------        
          sdvup(7:9)=(/Balpha(1,1),   Balpha(2,2), 2.0*Balpha(1,2)/)
          sdvup(14:16)=(/epse(1,1),   epse(2,2),   2.0*epse(1,2)/)

        end if

c--------------------------drag stress----------------------------------
        sdvup(13) = alpha        
        sdvup(20)=fn
        sdvup(21)=epsp_b
      end subroutine kGTN


c==================================================================  
c subroutine trial_step (eps,epspn,xmu,xHk,Balphan,alphan,xhi,
c            xsigy0, phitr, dsigtr, ntr, nxitr)
c------------------------------------------------------------------]
c Inputs
c   eps : stain tensor
c   epspn : plastic strain tensor
c   xmu : shear modulus
c   xHk : kinematic hardening modulus, scalar 
c   Balphan: tensorial internal variable alpha
c   alphan: scalar hardening variable alpha
c   xh :  isotropic hardening modulus, scalar
c   xsigy0: initial yield stress, uniaxial test
c-----------------------------------------------------------------
c Outputs
c   phitr : trial yield function 
c   dsigtr : deviatoric part of the trial stress tensor
c   ntr : flow direction from trial state
c   nxitr : norm of xitr
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c======================================================================= 


      subroutine trial_step (eps,epspn,xmu,xHk,Balphan,alphan,xhi,
     $       xsigy0, phitr, dsigtr, ntr, nxitr)
       INCLUDE 'ABA_PARAM.INC'       
       !use tensor_operations
       
       !implicit none
       double precision, dimension (3,3) :: eps, xid,deps, dsigtr, 
     $                                Bbetatr,Balphan, epspn 
       double precision:: betatr, xHk, xhi, xmu, alphan, xsigy0 
       double precision, dimension (3,3) :: xitr !trial value of deviatoric stress difference 
       double precision :: nxitr !norm of xitr           
       double precision , dimension (3,3) :: ntr !flow direction from trial state
       double precision :: phitr, tra_tem,cont2_2_tem !trial value of the yield function 
   
       call Ident1(xid,3) 
       !deviatoric part of the strain tensor
       call Trace(eps,3,tra_tem)
       deps = eps - tra_tem*xid/3.0 
       !deviatoric part of the trial stress tensor
       dsigtr = 2.0*xmu*(deps-epspn)
       !trial value of the back-stress tensor
       Bbetatr = xHk*Balphan  
       !trial value of the increase in yield stress (scalar beta: drag stress)
       betatr =xhi*alphan
       !trial value of deviatoric stress difference
       xitr =dsigtr - Bbetatr
       call contrac_2nd_2nd(xitr,xitr,3,cont2_2_tem)
       nxitr = (cont2_2_tem)**0.5
       !flow direction from trial state
       ntr = xitr/nxitr
c-----------------------trial yield function----------------------------  
       phitr = nxitr-(2.0/3.0)**(0.5)*(xsigy0+betatr)

      end subroutine trial_step


c=======================================================================
c subroutine trial_step_GTN (eps,epsn,epsen,xmu,xk,phitr, sig_0,ntr)
c=======================================================================
c=======================================================================
c Inputs
c   eps : stain tensor acording to load
c   epsn : strain tensor last step
c   epsen : elastic strain tensor last step
c   xmu : shear moduli
c   xk : Bulk moduli 
c   phitr : trial yield function 
c   sig_0 : equivalent tensile flow stress
c   ntr : normalized flow direction from trial state (according to article)
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c======================================================================= 
      subroutine trial_step_GTN (eps,epsn,epsen,xmu,xk,phitr, sig_0,
     $     n_a_tr,q1,q2,q3,f,dsig_tr,sig_e_tr, q_tr,p_tr)
        INCLUDE 'ABA_PARAM.INC'
        !use tensor_operations
        !implicit none
        double precision, dimension (3,3) :: eps, epsn, xid,deps,  
     $                                dsigtr,sig_e_tr, sig_test 
        double precision::   xhi, xmu,xk, alphan, sig_0, q1, q2,q3 
        double precision, dimension (3,3) :: xitr, dsig_tr !trial value of deviatoric stress difference            
        double precision , dimension (3,3) :: s_tr, n_a_tr !flow direction from trial state
        double precision , dimension (3,3) :: C_e          !Elasticity tensor in voigth notation due to simmetry
        double precision :: phitr, q_tr, p_tr, f, tra_tem,cont2_2_tem           !trial value of the yield function 
        double precision , dimension (3,3) :: epsen        !elastic strain last step
        double precision , dimension (3,3):: D_eps         !Change strain
        !Variables plain strain
        double precision:: pe
        double precision , dimension (3,3)::D_eps_b !Delta Strain in the plane of stressing
        double precision , dimension (3,3)::dD_eps_b !viatoric of D_eps_b 
        double precision , dimension (3,3)::D_epse !Delta elastic eps 
        double precision , dimension (3,3):: depse  !Deviatoric elastic eps
        double precision , dimension (3,3)::a,ad !e3e3 and its deviatoric  

        call Ident1(xid,3) 
        !Elasticity tensor

        C_e=2.0*xmu*Xid!-(xk-2.0/3.0*xmu)*Xid   !<----Take care that here maybe we require the 4th order identity tensor and deviator           
        !strain Change  in the step
        D_eps=eps-epsn
        !Elastic predictor (trial). All the strain is supoused to be elastic
        sig_e_tr=C_e*(epsen+D_eps)
        !trial values for p, s and q
c           if (str_case==0) then   
           call Trace(sig_e_tr,3,tra_tem)
           p_tr= - 1.0/3.0*tra_tem
           dsig_tr=sig_e_tr-1.0/3.0*tra_tem*xid
           call contrac_2nd_2nd(dsig_tr,dsig_tr,3,cont2_2_tem)
           q_tr=(3.0/2.0*cont2_2_tem)**0.5              
c           end if   
        !norm according to the article
        n_a_tr = 3.0/2.0*dsig_tr/q_tr
        sig_test=-p_tr*xid+(2.0/3.0)*q_tr*n_a_tr
        phitr = (q_tr/sig_0)**2+2.0*q1*f*cosh(-3.0*q2*p_tr/
     &                (2.0*sig_0))-(1.0+q3*f**2)
      end subroutine trial_step_GTN


c=======================================================================  
c subroutine plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
c       epsp, epspn, Balpha, Balphan, alpha, alphan)
c=======================================================================
c
c------------------------------------------------------------------------
c coded by: W. Mora Sep 2021
c======================================================================= 

      subroutine plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
     &  nxitr, epsp, epspn, Balpha, Balphan, alpha, alphan,beta1,beta2)
         INCLUDE 'ABA_PARAM.INC'     
         !implicit none
         
         double precision, dimension (3,3) ::  dsig, dsigtr,ntr,
     &    epsp, epspn, Balpha, Balphan  
         double precision :: phitr, xmu, xHk, alpha, alphan, beta1, 
     &    beta2,gamma,xhi, nxitr

         gamma=phitr/(2.0*xmu+xHk+2.0*xhi/3.0)
         dsig=dsigtr-2.0*xmu*gamma*ntr
         epsp=epspn+(dsigtr-dsig)/(2*xmu)
c         epsp=epspn+gamma*ntr
         Balpha = Balphan+gamma*ntr        ! tensorial internal variable alpha
         alpha = alphan+(2.0/3.0)**(0.5)*gamma           !scalar hardening variable alpha
         beta1=1.0-(phitr/nxitr)*(1.0/(1.0+xHk/(2.0*xmu)+
     $          xhi/(3.0*xmu)))
         beta2=(1.0-phitr/nxitr)*(1.0/(1.0+xHk/(2.0*xmu)+
     $          xhi/(3.0*xmu))) 
      end subroutine plas_corre_VM
     
c=======================================================================  
c subroutine plas_corre_GTN (q, p, sig,inputmat,fn, dsig,PHI, 
c              d_epsp, epsp_b,epsp,epse)
c=======================================================================
c=======================================================================
c Inputs
c   rec: counter to make a recursive call
c   q : hydrostac stress trial
c   p : equivalent stress trial
c   sig : stress trial 
c   inputmat : vector with the material's intrinsic properties.
c   f : void fraction trial
c   dsig : deviatoric stress trial
c   PHI : yield suface trial
c   epsp_b: micorscopic equivalent plastic strain
c   eps,D_eps, epsp, epse: strain, delta strain, plastic strain, elastic strain
c   Balpha, alpha: strass like variables Von Mises
c   sig_0 : equivalent tensile flow stress
c   sig1 : final stress
c    dsig1 : Change of stress
c    f1 : final void volume fraction  
c    ATStype : Selector for type of ATS computtion 
c    str_case : Selector for the stress case
c    epspn : initial plastic strain
c    epsp_b_n : initial micorscopic equivalent plastic strain
c    dumbsub : varible to declare the function as recursive
c------------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c======================================================================= 

      subroutine plas_corre_GTN (rec, q,p, sig,inputmat,f, dsig,phi, 
     &      epsp_b,eps,D_eps, epsp, epse, Balpha, alpha, sig_0, sig1,
     &      dsig1,f1,ATStype, str_case, epspn, epsp_b_n, ATS66,dumbsub)
            INCLUDE 'ABA_PARAM.INC'            
            !use tensor_operations
            !implicit none
            external dumbsub
c  Especific variables according to the article N. ARAVAS, 1987
            !Parameters material
            double precision, dimension (15)::inputmat
            !Material parameters
            double precision :: q, p, sig_0, q1, q2, q3,f,f1, xE,  
     &         xsigy0,  xk, xmu, xnu, q_e, p_tr, q_tr 
            double precision :: f_c, f_s, f_f, df_s_df, k_b, f_su !Parameters GTN model (Needleman part)
            double precision ::A_c 
            
            !Change of parameter between steps
            double precision ::D_eps_p, D_eps_q, D_epsp_b,D_f, D_eps_pn,
     &       D_eps_qn  
            double precision, dimension (3,3) ::D_epsp,D_epspn1, dsig,
     &       dsig_e, D_eps
            !Current value of ISV
            double precision :: epsp_b, epsp_b1, alpha, alphan
            double precision, dimension (3,3) ::sig ,epsp, eps, epsn, 
     &            eps_pp, epse, n_a, sign_l, Balpha, Balphan, sig1,dsig1
            !derivates defined in article of real type
            double precision:: dg_dp, dg_dq, d2g_d2p, d2g_d2q, d2g_dpdq,
     &        d2g_dpdf, d2g_dqdf, d2g_dq_d_epsp,
     &        dDf_dDeps_p, dDepsp_b_dDeps_p, dDf_dDeps_q, dDepsp_b_dq,
     &        dDf_dp, dDepsp_b_dDeps_q, dPHI_dp, dPHI_dq, dPHI_df,
     &        dPHI_d_epsp_b, d_epsp_b_dDeps_p, d_epsp_b_dp, dDepsp_b_dp,
     &        A11, A12, A21, A22, b1, b2, PHI, A11_D, A12_D, A21_D,
     &        A22_D, B11_D, B12_D, B21_D, B22_D, b1_D, b2_D, df_dp, 
     &        d_epsp_dp, d_epsp_b_dq, df_dq, d_epsp_dq, df_d_epsp_b,
     &        d2g_dq_d_epsp_b, d2g_dp_d_epsp, d2g_dp_d_epsp_b,
     &        df_dDeps_p, dH_alpha_dDeps_p, dDf_dq, df_dDeps_q,
     &        dH_alpha_dDeps_q, d_epsp_b_dDeps_q, d_sig0_d_f, 
     &        d_sig0_d_epsp_b, dA_d_epsp_b, dDepsp_b_df, d2g_dqdsig_0,
     &        dDf_d_epsp_b, dDf_df,dDepsp_b_d_epsp_b, NU, 
     &        d_sig_0_d_epsp_b
            !derivates defined in article of tensorial (3,3) type
            double precision, dimension (3,3) ::  dDepsp_dDeps_p, 
     &         dDepsp_dp,d_epsp_dDeps_p, dDepsp_dDeps_q, dDepsp_dq,
     &          d_epsp_dDeps_q, dPHI_d_epsp
            double precision, dimension (3,3,1) :: dou_sum1_P2, 
     &       dou_sum2_P2 
            !Secondary derivatives and variables
            double precision:: c_f_f, c_f_epsp, c_f_epsp_b
     &        , c_epsp_b_epsp , c_epsp_b_f, c_epsp_b_epsp_b
     &        , sum1, sum2, sum3, sum4 , dou_sum1, dou_sum1_P1,
     &        dou_sum2, dou_sum3, dou_sum4, dou_sum5, dou_sum6, 
     &        sum32, sum42,
     &        dou_sum7, dou_sum8 , dD_eps_p,dD_eps_q   
            double precision, dimension (3,3) :: nxi_a,
     &         c_epsp_f, c_epsp_epsp, c_epsp_epsp_b, D_sig   
            !Parameters nucleation
            double precision :: f_0, f_n, s_n, E_n ! E_n: mean value Parameter used in the nucleation equations
            !Auxiliar variables
            double precision, dimension (3,3) ::  xid, Delta_sig66
            double precision :: pi, tol_conv_cor,conv_D_eps_p,
     &          conv_D_eps_q, pert, con, ATS_ave
            double precision, dimension(2,2) :: M,M_inv,ML,MR,dD_eps_pq
            double precision, dimension(2) :: b_v, c_corr
            double precision, dimension(6) :: eps6_pp, sig6_pp, 
     &          Delta_sig6
            double precision, dimension(6,6) ::  ATS66,A66_t,C_e66_tem   
            double precision, dimension(3,3,3,3)::C, Cdev, C_e, C_e_inv,
     &           M_4T, ATS, dn_d_sig,I4dikdjl_r,I4dijdkl_r 
            integer :: iterartion, i,j, rc,sc, ATStype, str_case,rec
            double precision :: Beta1, Beta2, tol
            integer, dimension(6):: ii, jj

c      Variables for plane stress elastoplastic equations  
c           
            double precision ::A23_PS, A31_PS, A32_PS, A33_PS, b3_PS,
     &         A13_PS, D_eps3, dev_e_33, dq_dDeps3, pe, D_eps3n,
     &         conv_d_eps3, ome31, ome32, ome33
            double precision, dimension(3,3) :: M_ps, M_ps_inv,D_eps_b,
     &          dD_eps_b, depse, se, a, da, ATS_PS
            double precision, dimension(3) :: b_ps_v, c_corr_ps
            
            !Auxiliar variables
            double precision::conv_D_eps_p_ps, conv_D_eps_q_ps,D_eps11,
     &           D_eps22
            double precision, dimension (3,3) :: sigtr, dsigtr, epsen, 
     &       epspn, ntr
            double precision:: epsp_b_n, fn, phitr

c   rec       : number of recursive iterations             
c   ATStype   : Type of computation of algorithmic tangent
c   str_case  : Tyoe of analysis 3D or plane stress
c   conv_     : convergence
c   PS        : variable for the plane stress case 
c   D_ or D.. :      : Delta
c   D_        : variable used in computation of ATS lineaation module 
c   d:        : derivative, total or partial depending of the context
c   d2        : second derivative, total or partial depending of the context
c   sum       ;auxiliary variable to store the result of a summatory 
c   dou_sum   :auxiliary variable store the result of a double summatory
c   _P1, _P2  : auxiliary variable to assign the part of a long equation
c   eps       : total strain
c   epse      : elastic strain
c   epsp      : plastic strain
c   epsp_b    : microscopic eqivalent platic strain  
c   ..tr      : trial value
c   ATS       : Algorithmic tangent stiffness 4th order tensor
c   ATS66     : ATS retrieved in 6x6 2nd order tensor 
c   A66_t     : auxiliary ATS66 matrix to run the recursive code
c   p         : hydrotatic stress
c   q         : equivalent stress
c   dsig      : stress deviatoric
c   A_c       : Parmeter for nucleation strain Chu and Needleman 
c   M_4T      : Variable M, 4th order tensor used in solutions ATS
c   M,ML, MR  : matrix representing linear systems of equations (LSE)
c   b         : colum vector independent terms LSE.
c iterartion  : iterator 
c i,j, rc,sc  : iterators
c   f_s       : f* function of the void volume fraction
c   f_su      : ultimate value of f* at ductile rupture
c   k_b       : Relation in the Needleman model
c   f_c       : critical void volume at which coids coalesce
c   f_F       : void volume fraction at final failure of the       

            xE        = inputmat(1)
            xnu       = inputmat(2)  !Poison    
            xsigy0    = inputmat(3)  !Yield point     
            q1        = inputmat(6)
            q2        = inputmat(7)
            q3        = inputmat(8)
            f_0       = inputmat(9)
            f_n       = inputmat(10)
            s_n       = inputmat(11)
            f_f       = inputmat(12)
            E_n       = inputmat(13)
            f_c       = inputmat(14)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c Constants 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
            xmu  = xE/(2.0*(1.0+xnu))    !(G)Shear modulus 
            xk   = xE/(3.0*(1-2.0*xnu))  !(K) Compression Modulus
            call Ident1(xid,3)
            pi=4.D0*DATAN(1.D0)
            tol_conv_cor=1e-11

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c Initial values to start iteration
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            PHI=(q/sig_0)**2+2.0*f*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f**2)
         if (str_case==0) then
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------

c  return direction (only depends on the trial step            )
            dsig_e=sig-Trace(sig)*xid
            q_e=(3.0/2.0*contrac_2nd_2nd(dsig_e,dsig_e))**0.5
            n_a=3.0/2.0*dsig_e/q_e  !<--------This is computed only with the trial values, does not include correction. Last equation of page 1399
            p_tr=p
            q_tr=q  
            sig1=sig

         else
c-----------------------------------------------------------------------
c                       PLAIN STRAIN CASE
c-----------------------------------------------------------------------

            PHI=(q/sig_0)**2+2.0*f*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f**2)

c  return direction (only depends on the trial step            )
            dsig_e=sig-Trace(sig)*xid
            q_e=(3.0/2.0*contrac_2nd_2nd(dsig_e,dsig_e))**0.5
            n_a=3.0/2.0*dsig_e/q_e  !<--------This is computed only with the trial values, does not include correction. Last equation of page 1399
            
            D_eps_p=0.0
            D_eps_q=0.0
            D_eps11=D_eps(1,1)
            D_eps22=D_eps(2,2)     
            D_eps3=D_eps(3,3)
            pe=-xk*(Trace(epse)+D_eps11+D_eps22)
            
            D_eps_b(:2,:2)=D_eps(:2,:2) !In the aricle is called Delta e bar
            dD_eps_b=D_eps_b-1.0/3.0*Trace(D_eps_b)*xid !In the aricle is called deviatoric Delta e bar
            depse=epse-1.0/3.0*Trace(epse)*xid !In the article is called e^e, deviatoric elastic eps
            se=2*xmu*(depse+dD_eps_b) !se is the article name. This is a deviatoric with out name 
            a=reshape((/0,0,0,0,0,0,0,0,1/),shape(a), order=(/2,1/))
            da=a-1.0/3.0*Trace(a)*xid !Deviatoric of e33 
            dsig=se+2.0*xmu*(D_eps3*da-3.0/2.0*D_eps_q/q*dsig)
            p_tr=pe
            q_tr=q  !This q is qe, which comes form the arguments, does not change before this point
            
         end if
            sig1=sig
         
c  Set guess values to enter to the iterative algorithm 
            D_eps_p=0.000
            D_eps_q=0.0001
            D_eps3=0.000
                        
            iterartion=0
            D_eps_pn=D_eps_p
            D_eps_qn=D_eps_q
            D_eps3n=D_eps3
            conv_D_eps_p=1.0
            conv_D_eps_q=3.0

c=======================================================================
c     Volumetric and deviatoric plastic strain increments according to
c      a Newton scheme, and D_eps3 for in case of the plain stress case
c=======================================================================

   10     if ((abs(c_corr(1)) .gt. tol_conv_cor) .or. 
     $        (abs(c_corr(2)).gt. tol_conv_cor ) .and. (PHI>1e-1)) then       

            iterartion=iterartion+1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C           Materal parameter initial values
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          if (str_case==0) then
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------
            !D_epsp_b=0.01
c Update Hydrostatic and equivalent stress

            D_epsp_b=(-p*D_eps_p+q*D_eps_q)/(1.0-f)/sig_0
            epsp_b=epsp_b+D_epsp_b

            p=p_tr+xk*D_eps_p
            q=q_tr-3.0*xmu*D_eps_q

            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b-E_n)/
     &              s_n)**2) !Chu and Needleman parameter A
            d_f=(1.0-f)*D_eps_p+A_c*D_epsp_b
            f=f+d_f  

          else
c-----------------------------------------------------------------------
c                       PLAIN STRAIN CASE
c-----------------------------------------------------------------------
c Update Hydrostatic and equivalent stress

            D_epsp_b=(-p*D_eps_p+q*D_eps_q)/(1.0-f)/sig_0
            epsp_b=epsp_b+D_epsp_b
            
            p=p_tr-xk*(D_eps3-D_eps_p)
            q=-3.0*xmu*D_eps_q+(q_tr**2+6.0*xmu*se(3,3)*D_eps3+  
     &         4.0*xmu**2*D_eps3**2)**0.5   !here q^e is q_tr 

            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b-E_n)/
     &              s_n)**2) !Chu and Needleman parameter A
            d_f=(1.0-f)*D_eps_p+A_c*D_epsp_b
            f=f+d_f  
 
          end if
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Compute the required derivates to compute the constant of the LSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Derivatives to solve system of equations to solve D_eps_p and D_eps_q 
            dg_dq =2.0*(q/sig_0**2)
            dg_dp = -3.0*q1*q2/sig_0*f*sinh(-3.0*q2*p/(2.0*sig_0))
            d2g_dpdq=0.0
            d2g_d2p=(9.0/2.0) *q1*f*(q2/sig_0)**2
     &               *cosh((-3.0*q2*p)/(2.0*sig_0))  !
            d2g_d2q=(2.0/(sig_0**2))

            d_sig0_d_f=(-p*D_eps_p+q*D_eps_q)/((1-f)**2*D_epsp_b)
            d_sig0_d_epsp_b=-(-p*D_eps_p+q*D_eps_q)/((1-f)*D_epsp_b**2)  !Option 1 using the equation where is sig_0, D_epsp_b
            NU=0.1

            d2g_dpdf=-3.0*q1*q2*(p/sig_0*sinh(-3.0*q2*p/(2.0*sig_0))-
     &        f*p/sig_0**2*d_sig0_d_f*sinh(-q2*3.0*p/(2.0*sig_0))+
     &        f*p/sig_0*cosh(-q2*3.0*p/(2.0*sig_0))*q2*3.0*p/
     &        (2.0*sig_0**2)*d_sig0_d_f)

            dA_d_epsp_b=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-0.5*((epsp_b-
     &        E_n)/s_n)**2.0)*(-((epsp_b-E_n)/s_n)*(1/s_n))
            dDepsp_b_df=(-p*D_eps_p+q*D_eps_q)/((1-f)**2*sig_0)
            df_d_epsp_b=(dA_d_epsp_b*D_epsp_b+A_c*dDepsp_b_df)/
     &       (1.0+D_eps_p)
            dDepsp_b_d_epsp_b=dDepsp_b_df*df_d_epsp_b
            d2g_dp_d_epsp_b= -3.0*q1*q2*(df_d_epsp_b*p/sig_0*sinh(-3.0
     &        *q2*p/(2.0*sig_0))-f*p/sig_0**2*d_sig0_d_epsp_b*sinh(-q2
     &        *3.0*p/(2.0*sig_0))+f*p/sig_0*cosh(-q2*3.0*p/(2.0*sig_0))*
     &        q2*3.0*p/(2.0*sig_0**2)*d_sig0_d_epsp_b) 
            
            !d2g_dp_d_epsp_b=d2g_dpdf*df_d_epsp_b
            d2g_dqdsig_0=(-4.0*q/sig_0**3)
            d2g_dqdf=d2g_dqdsig_0*d_sig0_d_f
            !alpha2 = epsp (epsilon plastic)
            !alpha3 = f

            d2g_dq_d_epsp_b= d2g_dqdsig_0*d_sig0_d_epsp_b
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c    !dH_alpha_dDeps_p = df_dDeps_p + d_epsp_dDeps_p + d_epsp_b_dDeps_p
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  !Compute df_dDeps_p

            dDf_dDeps_p=(1.0-f)            
            !dDf_dp=A_c*(-D_eps_p/((1.0-f)*sig_0))  
            dDf_d_epsp_b=dA_d_epsp_b*D_epsp_b+ A_c*dDepsp_b_d_epsp_b
            d_epsp_b_dp=-D_eps_p/((1.0-f)*sig_0)
            dDf_dp=dDf_d_epsp_b*d_epsp_b_dp

            dDepsp_b_dDeps_p=-p/((1.0-f)*sig_0)
            dDepsp_b_dp=-D_eps_p/((1.0-f)*sig_0)

            dDf_df= -D_eps_p  
            c_f_f = 1.0/(Kron_d(1,1)-dDf_df)

            c_f_epsp_b=1.0/(Kron_d(1,3)-dDf_d_epsp_b)
            df_dDeps_p=c_f_f*(dDf_dDeps_p+xk*dDf_dp) 
     &            +c_f_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp) !<---I'm doing the second term 0, since epsp this is not a H_beta  

c  Compute  d_epsp_b_dDeps_p
            c_epsp_b_f =1.0/(Kron_d(3,1)-dDepsp_b_df) 
            !c_epsp_b_epsp =1/(Kron_d(3,2)-(f-1.0)/A_c)
            c_epsp_b_epsp_b=1/(Kron_d(3,3)-dDepsp_b_d_epsp_b)

            d_epsp_b_dDeps_p=c_epsp_b_f*(dDf_dDeps_p+xk*dDf_dp)
     &            +c_epsp_b_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp) !<---I'm doing the second term 0, since epsp this is not a H_beta

c   Finally compute dH_alpha_dDeps_p
            dH_alpha_dDeps_p= df_dDeps_p + d_epsp_b_dDeps_p 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     dH_alpha_dDeps_q = df_dDeps_q + d_epsp_dDeps_q+ d_epsp_b_dDeps_q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Compute df_dDeps_q

            dDf_dDeps_q=0
            d_epsp_b_dq=D_eps_q/((1-f)*sig_0)

            dDf_dq=dDf_d_epsp_b*d_epsp_b_dq

            dDepsp_b_dDeps_q= q/((1.0-f)*sig_0)
            dDepsp_b_dq     = D_eps_q/((1.0-f)*sig_0)
            df_dDeps_q=c_f_f*(dDf_dDeps_q-3.0*xmu*dDf_dq)
     &         +c_f_epsp_b*(dDepsp_b_dDeps_q-3.0*xmu*d_epsp_b_dq)

c Compute d_epsp_b_dDeps_q
            d_epsp_b_dDeps_q=c_epsp_b_f*(dDf_dDeps_q-3.0*xmu*dDf_dq)
     &        +c_epsp_b_epsp_b*(dDepsp_b_dDeps_q-3.0*xmu*d_epsp_b_dq)

c   Finally compute dH_alpha_dDeps_q 
            dH_alpha_dDeps_q = df_dDeps_q + d_epsp_b_dDeps_q
c     &               +d_epsp_dDeps_q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    dPHI_d
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            dPHI_dp =dg_dp !q1*f*sinh(-3.0*q2*p/(2.0*sig_0))
            dPHI_dq=dg_dq  !2*(q/sig_0**2)
c  Compute all  dPHI_dH_alpha = dPHI_df, dPHI_d_epsp 

            dPHI_df = -2*q**2/sig_0**3*d_sig0_d_f + 2.0*q1*(cosh(-3.0*
     &           q2*p/(2.0*sig_0))+f*sinh(-3.0*q2*p/(2.0*sig_0))*
     &           (3.0*q2*p/(2.0*sig_0**2))*d_sig0_d_f)-2.0*q3*f

            dPHI_d_epsp_b= (-2.0*q**2/sig_0**3)*d_sig0_d_epsp_b+2*q1*
     &        (f*sinh(-3.0/2.0*q2*p/sig_0)*(3.0/2.0*q2*p/sig_0**2)*
     &        d_sig0_d_epsp_b)                         ! Option 1 considering f = cte
     
            PHI=(q/sig_0)**2+2.0*f*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f**2) 
  
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
c Consolidation of Elastoplastic equation constants A11, A12, A21, A22, 
c b1, b2 and solution of LSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   
            sum1=d2g_dqdf*df_dDeps_p+d2g_dq_d_epsp_b*d_epsp_b_dDeps_p

            sum2=d2g_dpdf*df_dDeps_p+d2g_dp_d_epsp_b*d_epsp_b_dDeps_p

            A11=dg_dq+D_eps_p*(xk*d2g_dpdq+sum1)+D_eps_q*(xk*d2g_d2p
     &          +sum2)
            sum3=d2g_dqdf*df_dDeps_q+d2g_dq_d_epsp_b*d_epsp_b_dDeps_q
            sum4=d2g_dpdf*df_dDeps_q+d2g_dp_d_epsp_b*d_epsp_b_dDeps_q

            A12= dg_dp+D_eps_p*(-3.0*xmu*d2g_d2q+sum3)
     &             + D_eps_q*(-3*xmu*d2g_dpdq+sum4)  
            sum32=dPHI_df*df_dDeps_p+dPHI_d_epsp_b*d_epsp_b_dDeps_p

            A21=xk*dPHI_dp+sum32
            sum42=dPHI_df*df_dDeps_q+dPHI_d_epsp_b*d_epsp_b_dDeps_q

            A22=-3.0*xmu*dPHI_dq+sum42

            b1=-D_eps_p*dg_dq-D_eps_q*dg_dp
            b2= -PHI

         if (str_case==1) then
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
c Consolidation of Elastoplastic equation constants A13, A23, A31, A32, 
c b3 / Applies only to plane stress
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            A13_PS=D_eps_p*(-xk*d2g_dpdq+dq_dDeps3*d2g_d2q)+
     &        D_eps_q*(-xk*d2g_d2p+dq_dDeps3*d2g_dpdq)            
            A23_PS=-xk*dPHI_dp+dq_dDeps3*dPHI_dq
            A31_PS=xk*(q+3.0*xmu*D_eps_q)
            A32_PS=3.0*xmu*dev_e_33
            A33_PS=dq_dDeps3*p-xk*(q+3.0*xmu*D_eps_q)-(4.0/3.0)*
     &        xmu*q-(dev_e_33+(4.0/3.0)*xmu*D_eps3)*dq_dDeps3
            dq_dDeps3= (3.0*xmu*dev_e_33+4.0*xmu**2*D_eps3)/
     &        (q+3.0*xmu*D_eps_q)
            b3_PS=-(q+3.0*xmu*D_eps_q)*p+(dev_e_33+
     &        (4.0/3.0)*xmu*D_eps3)*q
            
         end if

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      Solution of LSE for the correction values and update of volumetric 
c       and deviatoric plastic strain increments
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (str_case==0) then
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------
c  Solution of LSE for the correction values
            M = reshape((/A11, A12,
     $                    A21, A22/),shape(M), order=(/2,1/))
            b_v=(/b1,b2/))
            call inv_T_2Dd(M, 2,M_inv)
            c_corr = matmul(M_inv, b_v)


c Update volumetric and deviatoric plastic strain increments
            D_eps_p=D_eps_p+c_corr(1) 
            D_eps_q=D_eps_q+c_corr(2)

c Evaluate step convergence  
            conv_D_eps_p=abs(D_eps_p-D_eps_pn)
            conv_D_eps_q=abs(D_eps_p-D_eps_pn)
            D_eps_pn=D_eps_p
            D_eps_qn=D_eps_q
         else
c-----------------------------------------------------------------------
c                       PLAIN STRESS CASE
c-----------------------------------------------------------------------
c  Solution of LSE for the correction values
            
            M_ps = reshape ((/A11,    A12,    A13_PS,
     &                        A21,    A22,    A23_PS,
     &                        A31_PS, A32_PS, A33_PS/),
     &                       shape(M_ps), order=(/2,1/))
            b_ps_v =(/b1, b2, b3_PS/) 
            call inv_T_2Dd(M_ps,3,M_ps_inv)
            c_corr_ps=matmul(M_ps_inv, b_ps_v)      

c Update volumetric and deviatoric plastic strain increments and D_eps3

            D_eps_p=D_eps_p+c_corr_ps(1) 
            D_eps_q=D_eps_q+c_corr_ps(2)
            D_eps3=D_eps_q+c_corr_ps(3)

c Evaluate step convergence  
            conv_D_eps_p=abs(D_eps_p-D_eps_pn)
            conv_D_eps_q=abs(D_eps_p-D_eps_pn)
            conv_D_eps3=abs(D_eps3-D_eps3n)

            D_eps_pn=D_eps_p
            D_eps_qn=D_eps_q
            D_eps3n=D_eps3

         end if
        goto 10
        endif

c Update stress and strain tensors
            
  
            sig=-p*xid+(2.0/3.0)*q*n_a      
            dsig=sig-1.0/3.0*xid*Trace(sig)
            D_epsp=1.0/3.0*D_eps_p*xid+D_eps_q*n_a
            epsp=epsp + D_epsp
            epse=eps-epsp  

            Balpha=Balphan
            alpha= alphan
            D_sig =sig-sig1 !Delta sigma
            sig1=sig
            dsig1=dsig
            f1=f

         if (ATStype==0) then
c=======================================================================
c          Compute Material tangent stiffness "analitical" solution
c=======================================================================
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Additional derivatives 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            df_dp=1.0 
            d_epsp_dp=1.0
            df_dq=1.0
            d_epsp_dq=1.0
            dou_sum1=(D_eps_p*d2g_dqdf+D_eps_q*d2g_dpdf)*(c_f_f*
     &         df_dDeps_p+c_f_epsp_b*d_epsp_b_dDeps_p)+
     &         (D_eps_p*d2g_dq_d_epsp_b+D_eps_q*d2g_dp_d_epsp_b)* 
     &         (c_epsp_b_f*df_dDeps_p+
     &         c_epsp_b_epsp_b*d_epsp_b_dDeps_p)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Consolidation of constants to compute LSE to solve D: A11, A12, 
c  A21, A22, B11, B12, B21, B22 and solution of LSE 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            A11_D=dg_dq+dou_sum1
            dou_sum2=(D_eps_p*d2g_dqdf+D_eps_q*d2g_dpdf)*(c_f_f*
     &         df_dDeps_q+c_f_epsp_b*d_epsp_b_dDeps_q)+
     &         (D_eps_p*d2g_dq_d_epsp_b+D_eps_q*d2g_dp_d_epsp_b)* 
     &         (c_epsp_b_f*df_dDeps_q+c_epsp_b_epsp_b*d_epsp_b_dDeps_q)
            A12_D=dg_dp+dou_sum2
            A21_D=dPHI_df *(c_f_f*df_dDeps_p+
     &         +c_f_epsp_b*d_epsp_b_dDeps_p)+ 
     &         dPHI_d_epsp_b*(c_epsp_b_f*df_dDeps_p
     &         +c_epsp_b_epsp_b*d_epsp_b_dDeps_p)
            A22_D=dPHI_df *(c_f_f*df_dDeps_q+
     &         c_f_epsp_b*d_epsp_b_dDeps_q)+
     &         dPHI_d_epsp_b*(c_epsp_b_f*df_dDeps_q+
     &         c_epsp_b_epsp_b*d_epsp_b_dDeps_q)
           dou_sum3=(d2g_dqdf)*(c_f_f*df_dp+c_f_epsp_b*d_epsp_b_dp)+
     &         (d2g_dq_d_epsp_b)*(c_epsp_b_f*df_dp+
     &          c_epsp_b_epsp_b*d_epsp_b_dp)
           dou_sum4=(d2g_dpdf)*(c_f_f*df_dp+c_f_epsp_b*d_epsp_b_dp)+
     &         (d2g_dp_d_epsp_b)*(c_epsp_b_f*df_dp 
     &         +c_epsp_b_epsp_b*d_epsp_b_dp)
           B11_D=(1.0/3.0)*(D_eps_p*(d2g_dpdq+dou_sum3)+D_eps_q*(d2g_d2p
     &            +dou_sum4))
           dou_sum5=(d2g_dqdf)*(c_f_f*df_dq+c_f_epsp_b*d_epsp_b_dq)+
     &         (d2g_dq_d_epsp_b)*(c_epsp_b_f*df_dq+ 
     &         c_epsp_b_epsp_b*d_epsp_b_dq)
           dou_sum6=(d2g_dpdf)*(c_f_f*df_dq+
     &         +c_f_epsp_b*d_epsp_b_dq)+
     &         (d2g_dp_d_epsp_b)*(c_epsp_b_f*df_dq 
     &         +c_epsp_b_epsp_b*d_epsp_b_dq)
           B12_D=-D_eps_p*(d2g_d2q+dou_sum5)-D_eps_q*(d2g_dpdq+dou_sum6)
           dou_sum7=dPHI_df *(c_f_f*df_dp+
     &         +c_f_epsp_b*d_epsp_b_dp)+ 
     &         dPHI_d_epsp_b*(c_epsp_b_f*df_dp+
     &         c_epsp_b_epsp_b*d_epsp_b_dp)
           B21_D=(1.0/3.0)*(dPHI_dp+dou_sum7)
           dou_sum8=dPHI_df *(c_f_f*df_dq+
     &         +c_f_epsp_b*d_epsp_b_dq)+ 
     &         dPHI_d_epsp_b*(c_epsp_b_f*df_dq
     &         +c_epsp_b_epsp_b*d_epsp_b_dq)
           B22_D=-(dPHI_dq+dou_sum8)
              
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Solution of LSE to find patial(D_eps_p) and partial(D_eps_q)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ML = reshape((/A11_D, A12_D,
     $               A21_D, A22_D/),shape(M), order=(/2,1/))
       MR = reshape((/B11_D, B12_D,
     $               B21_D, B22_D/),shape(M), order=(/2,1/)) 

        M_inv = inv_T_2D(ML)
        dD_eps_pq = matmul(M_inv, MR)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Find ATS using patial(D_eps_p) and partial(D_eps_q)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ML = reshape((/A11_D, A12_D,
     $               A21_D, A22_D/),shape(M), order=(/2,1/))
       MR = reshape((/B11_D, B12_D,
     $               B21_D, B22_D/),shape(M), order=(/2,1/)) 

        M_inv = inv_T_2D(ML)
        dD_eps_pq = matmul(M_inv, MR)
        !print*, "mpn/3", dD_eps_pq(1,2)/3 ,"mql", dD_eps_pq(1,2)
        !print*, "partial delta eps_p, partial delta eps_p", c_corr

      end if

c-----------------------------------------------------------------------
c                       PLAIN STRESS CASE
c-----------------------------------------------------------------------
       ome31=-ATS(3,3,3,1)/ATS(3,3,3,3)  
       ome32=-ATS(3,3,3,1)/ATS(3,3,3,3)  
       ome33=-2*ATS(3,3,3,1)/ATS(3,3,3,3) 
      ATS_PS= reshape((/ ATS(1,1,1,1)+ome31*ATS(1,1,3,3), 
     & ATS(1,1,2,2)+ome32*ATS(1,1,3,3), ATS(1,1,1,2)+ome33*ATS(1,1,3,3),
     & ATS(2,2,1,1)+ome31*ATS(2,2,3,3),
     & ATS(2,2,2,2)+ome32*ATS(2,2,3,3), ATS(2,2,1,2)+ome33*ATS(2,2,3,3),
     & ATS(1,2,1,1)+ome31*ATS(1,2,3,3), 
     & ATS(1,2,2,2)+ome32*ATS(1,2,3,3), ATS(1,2,1,2)+ome33*ATS(1,2,3,3)
     & /),shape(ATS_PS), order=(/2,1/))

      if (str_case==1) then
         ATS66(1:3,1:3)=ATS_PS
      end if


      end subroutine plas_corre_GTN


c=======================================================================  
c                    subroutine Iso_hard1
c=======================================================================
c returns the yield stress and its derivative with respecto 
c to epsp_b  of the full dense matrix according to an Isotropic hardening  
c funtion of the microscopic equvalent plastic stress according to 
c sig_0/sigy0=(sig_0/sigy0+3G*epsp_b/sigy0)^N 
c------------------------------------------------------------------------
c Inputs
c epsp_b Microscopic equivaent plastic strain
c sig_0 : yield stress of the full dense matrix
c d_sig_0_d_epsp_b : derivative d_sig_0 with respect to epsp_b 
c coded by: W. Mora Oct 2021
c======================================================================= 

      subroutine Iso_hard1(epsp_b, sig_0,d_sig_0_d_epsp_b,inputmat)
        INCLUDE 'ABA_PARAM.INC'

        !implicit none
        double precision:: xE, xnu, mu, sigy0, sig_0, epsp_b
     &   , conv,sig_0n, NU, G, D_Y , d_sig_0_d_epsp_b, xsigy0
        double precision, dimension(15):: inputmat
        integer::i

        xE        = inputmat(1)
        xnu       = inputmat(2)  !Poison    
        sigy0     = inputmat(3)  !Yield point     
        NU        = inputmat(15)        

        mu       = xE/(2.0*(1.0+xnu)) 
        conv=1.0
        sig_0=sigy0
        i=0
   15 if (conv .gt. 1e-6) then
          !i=i+1
          !print*, "Iteration", i  
          G=sig_0/sigy0-(sig_0/sigy0+3*mu/sigy0*epsp_b)**NU
          d_sig_0_d_epsp_b=3.0*mu*NU*(sig_0/sigy0+3.0*mu*epsp_b/sigy0)
     &              **(NU-1.0)/(1.0-NU*(sig_0/sigy0+
     &              3.0*mu*epsp_b/sigy0)**(NU-1))
          D_Y=-G/d_sig_0_d_epsp_b    
          sig_0n=sig_0
          sig_0=sig_0+D_Y
          conv=abs(sig_0-sig_0n)
       goto 15
       endif
      end subroutine Iso_hard1

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
      subroutine linespace (ini, end, vector,n, type)
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
        integer :: i, type,n
        double precision ::ini, end
        double precision, dimension (n) :: vector
        
  
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
       Subroutine Trace(A,n,Tr_A)
       !function Trace(A) result(Tr_A)
         INCLUDE 'ABA_PARAM.INC'
         !implicit none
         integer :: n, n1
         double precision, dimension(n,n):: A
         double precision :: Tr_A
         
         Tr_A=0.0
         do n1=1,size (A,1)
           Tr_A=Tr_A+A(n1,n1)
         end do
       !end function Trace
       end Subroutine Trace   

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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
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
c                    Subroutine inv_T_2Dg
c=======================================================================
c=======================================================================  
c Returns the inverse of a 2nd order tensor according to the Gaus Jordan
c elimination, the code was modified from
c http://computer-programming-forum.com/49-fortran/6083a0ae451dd206.htm
c-----------------------------------------------------------------------
cInputs
c    A: Tensor 2D
c    n: size of square matrix 
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
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (Ainv(icol,icol).eq.0.) pause 'singular matrix in gaussj'

        pivinv=1./Ainv(icol,icol)
        Ainv(icol,icol)=1.
        do 16 l=1,n
          Ainv(icol,l)=Ainv(icol,l)*pivinv
 16      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=Ainv(ll,icol)
            Ainv(ll,icol)=0.
            do 18 l=1,n
              Ainv(ll,l)=Ainv(ll,l)-Ainv(icol,l)*dum
 18          continue
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
c                    function inv_T_2Dd
c=======================================================================
c=======================================================================  
c Returns the inverse of a 2nd order tensor of maximum size equal to 3x3
c http://computer-programming-forum.com/49-fortran/6083a0ae451dd206.htm
c-----------------------------------------------------------------------
cInputs
c    A: Tensor 2D se 2x2 or 3x3
c    n: size of square matrix 
c-----------------------------------------------------------------------
cOutput
c    A: inverted tensor
c-----------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c=======================================================================

      subroutine inv_T_2Dd(A,n,Ainv)

        INCLUDE 'ABA_PARAM.INC'
        !Implicit none    
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
         
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
         
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
        
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
        INCLUDE 'ABA_PARAM.INC'
       !implicit none
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
       INCLUDE 'ABA_PARAM.INC'
       !implicit none
       integer:: i, j, k , l     !iterators
       integer:: n1               !matrix dimension 
       double precision, dimension(3,3,3,3)::tensor, tensor_Inv !4th order tensors 
       double precision, dimension(6,6)::A6x6, A6x6_inv         !Auxiliary 2D order
      
        call T4th_2_Voig(tensor, A6x6)
        call inv_T_2Dg(A6x6,6,A6x6_inv)
        call voig_2_T4th(A6x6_inv,tensor_Inv)
      
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
      subroutine contrac_2nd_2nd(T_2D_1, T_2D_2,n,Contr)
      !function contrac_2nd_2nd(T_2D_1, T_2D_2) result(Contr)
       INCLUDE 'ABA_PARAM.INC'
       !implicit none 
       integer:: n    !dimension tensor.
       double precision, dimension(n,n) :: T_2D_1, T_2D_2 !2nd order tensor
       double precision :: Contr !Result of the contraction
       integer:: i, j !Iterators
        
        Contr =0.00000000
        do i=1,n
          do j=1,n
            Contr =Contr+ T_2D_1(j,i)*T_2D_2(j,i)
          end do
        end do
      
      !end function contrac_2nd_2nd
      end subroutine contrac_2nd_2nd
  
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
      subroutine diadic_prod_T2_T2(T_2D_1,T_2D_2, n1,n2,Diad)
      !function diadic_prod_T2_T2(T_2D_1,T_2D_2) result(Diad)  
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
        integer:: n, n1, n2       !size dimension tensor
        double precision, dimension(n1,n2) :: T_2D_1, T_2D_2  !2nd order tensor
        double precision, dimension(size(T_2D_1,1),size(T_2D_1,1),
     &                 size(T_2D_1,1),size(T_2D_1,1)):: Diad !4th order tensor   
        integer:: i, j, k, l !Iterators
        
      
        n=size(T_2D_1,1)
        !n2=size(Diad,4)  
      
        do i=1,n
            do j=1,n
                do k=1,n
                    do l=1,n
                      Diad(j,i,l,k) = T_2D_1(j,i)*T_2D_2(l,k)
                    end do
                  end do
                end do
              end do
      
      !end function diadic_prod_T2_T2
      end subroutine diadic_prod_T2_T2
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
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
       !implicit none
       INCLUDE 'ABA_PARAM.INC'
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
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
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
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

      subroutine loading2(ltype,posi_last_time,dt,t,lam,load,nlo) 
        !use tensor_operations
        INCLUDE 'ABA_PARAM.INC'
        !implicit none
        integer:: ltype, type , nlo !
        integer, dimension(5) :: posi_last_time
        integer::n_step_a, n_step_b, n_step_c
        double precision :: dt
        double precision, dimension(7) :: t
        double precision, dimension(3) ::lam
        double precision, dimension(nlo):: load
        double precision, dimension(:), allocatable :: xa, ya, xb, yb, 
     $                                                 xc, yc
        double precision ::xa1, ya1, xb1, yb1, xc1, yc1 
        double precision ::xa2, ya2, xb2, yb2, xc2, yc2
        
        !ltype==1 correspond to a load with one linear load
        if (ltype==1) then
          type=1
          n_step_a=nint((t(2)-t(1))/dt)
          allocate (xa(n_step_a+1))
          allocate (ya(n_step_a+1))
          
          xa1=t(1); ya1=lam(1)
          xa2=t(2); ya2=lam(2)
          call linespace(ya1,ya2,ya,n_step_a+1,type)
      
          load=ya
          deallocate(xa)
          deallocate(ya)    
      
        !ltype==2 correspond to a semicycle with one linear load and a linear unload
        elseif (ltype==2) then
          type=1
          n_step_a=nint((t(2)-t(1))/dt)
          allocate (xa(n_step_a+1))
          allocate (ya(n_step_a+1))
          call linespace(t(1),t(2),xa,n_step_a+1,type)
          
          xa1=t(1); ya1=lam(1)
          xa2=t(2); ya2=lam(2)
          call linespace(ya1,ya2,ya,n_step_a+1,type)
          
          type=0
          n_step_b=nint((t(3)-t(2))/dt)
          allocate (xb(n_step_b))
          allocate (yb(n_step_b))
          call linespace(t(2),t(3),xb, n_step_b, type)
          xb1=t(2); yb1=lam(2)
          xb2=t(3); yb2=lam(1)
          yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1
      
          deallocate (xa)
          deallocate (xb)
          load=(/ya,yb/)
          deallocate (ya)
          deallocate (yb)   
      
        !ltype==4 correspond to a load with 1 load, 1 unload including
        ! a negative part and a load to the zero point
        elseif (ltype==4) then
          type=1
          n_step_a=nint((t(2)-t(1))/dt)
          allocate (xa(n_step_a+1))
          allocate (ya(n_step_a+1))
          call linespace(t(1),t(2),xa,n_step_a+1,type)
          
          xa1=t(1); ya1=lam(1)
          xa2=t(2); ya2=lam(2)
          call linespace(ya1,ya2,ya,n_step_a+1,type)
      
          type=0
          n_step_b=nint((t(3)-t(2))/dt)
          allocate (xb(n_step_b))
          allocate (yb(n_step_b))
          call linespace(t(2),t(3),xb,n_step_b, type)
          xb1=t(2); yb1=lam(2)
          xb2=t(3); yb2=lam(3)
          yb=(yb2-yb1)/(xb2-xb1)*(xb-xb1)+yb1
          
          n_step_c=nint((t(4)-t(3))/dt)
          allocate (xc(n_step_c))
          allocate (yc(n_step_c))
          call linespace(t(3),t(4),xc,n_step_c, type)
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
      
        end if
      end subroutine loading2
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
        !implicit none
        INCLUDE 'ABA_PARAM.INC'
        integer::i, j
        double precision :: kd_ij
        if (i==j) then
          kd_ij=1.0
        else
          kd_ij=0.0
        end if
      end function Kron_d
    
c      end module tensor_operations 

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
          !implicit none
          INCLUDE 'ABA_PARAM.INC'
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
          !implicit none
          integer:: i, j,end,stat,line_no,n1,n2     !iterators
          
          character(len=30) :: file_name
          double precision :: import1(n1,n2)
          n1=size(import1,1)
          n2=size(import1,2)
      
          open(15, file=file_name,access='sequential',
     &        form="formatted", iostat=stat)
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
        !implicit none
        INCLUDE 'ABA_PARAM.INC'
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
        !implicit none
        INCLUDE 'ABA_PARAM.INC'
        double precision, dimension(3,3):: T   !2nd order tensor
        integer:: i, j
        integer:: n  !size dimension tensor
      
        do i=1,3
            do j=1,3
              T= 0.000
      
          end do
        end do
      
      end subroutine emp_2D
