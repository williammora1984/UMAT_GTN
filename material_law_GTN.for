c=======================================================================  
c subroutine GTN : Material model according to GURSON, TVERGAARD AND 
c NEEDLEMAN (GTN) DAMAGE MODEL  
c-----------------------------------------------------------------------
c  This subroutine implement the nmerical solution for the GURSON, 
c TVERGAARD AND NEEDLEMAN (GTN) DAMAGE MODEL according to the N.Aravas
c  article.
c-----------------------------------------------------------------------
c Version 0.9.3
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
c       INCLUDE 'ABA_PARAM.INC'       
        use tensor_operations

        implicit none
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


cFill the identity and deviatoriser tensor               
        call Ident1(xid,3)
        call P4sym(P4sym_r)
        call Iso_hard1(epsp_b_n, sig_0, d_sig_0_d_epsp_b, inputmat)
c      sig_0=xsigy0         
        
        pru=0

c=======================================================================
c                      Elastic prediction
c=======================================================================
       call dtime(t1,e)         !  Startup etime - do not use result

       !Mat_L=0 : VM, Mat_L=1 : GTN
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
            f=fn
            epsp_b=epsp_b_n
            if (ATStype==0) then
              Cdev=2.0*xmu*P4sym_r      
            endif
            if (Mat_L==0) then 
c-------------------Update stress--------------------------------------- 
c-----------Addition of spherical contribution to stress----------------
              sig=dsig+xk*Trace(eps)*xid
              epse=eps-epsp
              C=Cdev+xk*diadic_prod_T2_T2(xid,xid)
            end if
            if (Mat_L==1) then
c-------------------Update stress--------------------------------------- 
c-----------Addition of spherical contribution to stress----------------
              epse=eps-epsp
              C=Cdev+xk*diadic_prod_T2_T2(xid,xid)
              call m4th_2_66_sym (C, A66)           
            end if                          
c=======================================================================
c                      plastic correction
c=======================================================================
        call dtime(t2,e)
        else 

c          Mat_L=0 : VM, Mat_L=1 : GTN
           if (Mat_L==0) then
              !print*, "plastic return VM"
             call plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
     &            nxitr, epsp, epspn, Balpha, Balphan, alpha, alphan,
     &             beta1,beta2 )
             if (ATStype==0) then
                Cdev=2*xmu*beta1*P4sym_r-2*xmu*beta2*diadic_prod_T2_T2
     $               (ntr,ntr)
c------------   -------Update stress--------------------------------------- 
c-----------A   ddition of spherical contribution to stress----------------
                 sig=dsig+xk*Trace(eps)*xid
                 epse=eps-epsp
                 C=Cdev+xk*diadic_prod_T2_T2(xid,xid)
             end if
  
c------------------Algorithmic Tangen Stifness-------------------------     
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
c   Due to the symmetry the tangent stiffness 4 grade tensor can be 
c   expressed as a 6x6 matrix
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
        sdvup(20)=f
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
c       INCLUDE 'ABA_PARAM.INC'       
       use tensor_operations
       
       implicit none
       double precision, dimension (3,3) :: eps, xid,deps, dsigtr, 
     $                                Bbetatr,Balphan, epspn 
       double precision:: betatr, xHk, xhi, xmu, alphan, xsigy0 
       double precision, dimension (3,3) :: xitr !trial value of deviatoric stress difference 
       double precision :: nxitr !norm of xitr           
       double precision , dimension (3,3) :: ntr !flow direction from trial state
       double precision :: phitr !trial value of the yield function 
   
       call Ident1(xid,3) 
       !deviatoric part of the strain tensor
       deps = eps - Trace(eps)*xid/3.0 
       !deviatoric part of the trial stress tensor
       dsigtr = 2.0*xmu*(deps-epspn)
       
       !trial value of the back-stress tensor
       Bbetatr = xHk*Balphan  
       !trial value of the increase in yield stress (scalar beta: drag stress)
       betatr =xhi*alphan
       !trial value of deviatoric stress difference
       xitr =dsigtr-Bbetatr
       !norm of xitr 
       nxitr = (contrac_2nd_2nd(xitr,xitr))**0.5 
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
c           INCLUDE 'ABA_PARAM.INC'
           use tensor_operations
           implicit none
           double precision, dimension (3,3) :: eps, epsn, xid,deps,  
     $                                dsigtr,sig_e_tr, sig_test 
           double precision::   xhi, xmu,xk, alphan, sig_0, q1, q2,q3 
           double precision, dimension (3,3) :: xitr, dsig_tr !trial value of deviatoric stress difference            
           double precision , dimension (3,3) :: s_tr, n_a_tr !flow direction from trial state
           double precision , dimension (3,3) :: C_e          !Elasticity tensor in voigth notation due to simmetry
           double precision :: phitr, q_tr, p_tr, f           !trial value of the yield function 
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

           C_e=2.0*xmu*Xid
           !Change strain in the step
           D_eps=eps-epsn

           !Elastic predictor (trial). All the strain is supoused to be elastic
           sig_e_tr=C_e*(epsen+D_eps)
           !trial values for p, s and q   
              p_tr= - 1.0/3.0*Trace(sig_e_tr)
              dsig_tr=sig_e_tr-1.0/3.0*Trace(sig_e_tr)*xid
              q_tr=(3.0/2.0*contrac_2nd_2nd(dsig_tr,dsig_tr))**0.5              

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
c       INCLUDE 'ABA_PARAM.INC'     
       implicit none
       
       double precision, dimension (3,3) ::  dsig, dsigtr,ntr,
     &    epsp, epspn, Balpha, Balphan  
       double precision :: phitr, xmu, xHk, alpha, alphan, beta1, beta2,
     &    gamma,xhi, nxitr
       gamma=phitr/(2.0*xmu+xHk+2.0*xhi/3.0)
       dsig=dsigtr-2.0*xmu*gamma*ntr
       epsp=epspn+(dsigtr-dsig)/(2*xmu) !This is the same, you can solve from the avobe equation this expression  

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
c   q : hydrostac stress trial
c   p : equivalent stress trial
c   sig : stress trial 
c   inputmat : vector with the material's intrinsic properties.
c   fn : void fraction trial
c   dsig : deviatoric stress trial
c   PHI : yield suface trial
c-----------------------------------------------------------------------
c Outputs
c   phitr : 
c------------------------------------------------------------------------
c coded by: W. Mora Oct 2021
c======================================================================= 


      subroutine plas_corre_GTN (rec, q,p, sig,inputmat,f, dsig,phi, 
     &      epsp_b,eps,D_eps, epsp, epse, Balpha, alpha, sig_0, sig1,
     &      dsig1,f1,ATStype, str_case, epspn, epsp_b_n, ATS66,dumbsub)
c            INCLUDE 'ABA_PARAM.INC'            
            use tensor_operations
            implicit none
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
            double precision :: epsp_b, epsp_b1, alpha, alphan  ! microscopic equivalent plastic strain
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
            double precision, dimension(6,6) ::  ATS66,A66_t,
     &       C_e66_tem   
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
            tol_conv_cor=0.5e-10

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c Initial values to start iteration
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            f_su=1.0/q1
            k_b=(f_su-f_c)/(f_F-f_c)
            if (f .le. f_c) then
              df_s_df=1
              f_s=f
            else
              df_s_df=k_b
              f_s=f_c+k_b*(f-f_c)
            end if
 
            
            PHI=(q/sig_0)**2+2.0*f*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f**2)
         if (str_case==0) then
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------
 

c  return direction (only depends on the trial step)
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
            D_eps_q=0.0000001
            D_eps3=0.000 
            iterartion=0
            D_eps_pn=D_eps_p
            D_eps_qn=D_eps_q
            D_eps3n=D_eps3
            conv_D_eps_p=1.0
            conv_D_eps_q=3.0
            c_corr(1)=1.0
            c_corr(2)=2.0
            f1=f  
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
c------------------Relation between f and f*----------------------------
            if (f1 .le. f_c) then
              df_s_df=1
              f_s=f1
            else
              df_s_df=k_b
              f_s=f_c+k_b*(f1-f_c)
            end if              
          if (str_case==0) then
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------

c------- Update microscopic equivalent plasic strain and sig_0----------
            D_epsp_b=(-p*D_eps_p+q*D_eps_q)/(1.0-f1)/sig_0
            epsp_b1=epsp_b+D_epsp_b
            if (epsp_b1 < 0.0) then
              epsp_b1=epsp_b_n
            end if
            call Iso_hard1(epsp_b1, sig_0, d_sig_0_d_epsp_b, inputmat)
c----------- Update Hydrostatic and equivalent stress-------------------
            p=p_tr+xk*D_eps_p
            q=q_tr-3.0*xmu*D_eps_q

            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b1-E_n)/
     &              s_n)**2) !Chu and Needleman parameter A
            D_f=(1.0-f1)*D_eps_p+A_c*D_epsp_b
            f1=f+D_f  
          else
c-----------------------------------------------------------------------
c                       PLAIN STRAIN CASE
c-----------------------------------------------------------------------
c Update Hydrostatic and equivalent stress

            D_epsp_b=(-p*D_eps_p+q*D_eps_q)/(1.0-f1)/sig_0
            epsp_b1=epsp_b+D_epsp_b
            call Iso_hard1(epsp_b1, sig_0, d_sig_0_d_epsp_b, inputmat)
            
            p=p_tr-xk*(D_eps3-D_eps_p)
            q=-3.0*xmu*D_eps_q+(q_tr**2+6.0*xmu*se(3,3)*D_eps3+  
     &         4.0*xmu**2*D_eps3**2)**0.5   !here q^e is q_tr 

            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b1-E_n)/
     &              s_n)**2) !Chu and Needleman parameter A
            D_f=(1.0-f1)*D_eps_p+A_c*D_epsp_b
            f1=f+D_f  
 
c     &        "D_epsp", D_epsp,       

          end if
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Compute the required derivates to compute the constant of the LSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Derivatives to solve system of equations to solve D_eps_p and D_eps_q 
            dg_dq =2.0*(q/sig_0**2)
            dg_dp = -3.0*q1*q2/sig_0*f1*sinh(-3.0*q2*p/(2.0*sig_0))
            d2g_dpdq=0.0
            d2g_d2p=(9.0/2.0) *q1*f1*(q2/sig_0)**2
     &               *cosh((-3.0*q2*p)/(2.0*sig_0))  !
            d2g_d2q=(2.0/(sig_0**2))

            d_sig0_d_epsp_b= 3.0*xmu*NU*(sig_0/xsigy0+3.0*xmu*epsp_b1/
     &       xsigy0)**(NU-1)/(1.0-NU*(sig_0/xsigy0+3.0*xmu*epsp_b1/
     &       xsigy0)**(NU-1))   
            d2g_dpdf=-3.0*q1*q2*sinh(-3.0*q2*p/(2.0*sig_0))/sig_0
     &          *df_s_df        

            dA_d_epsp_b=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-0.5*((epsp_b1-
     &        E_n)/s_n)**2.0)*(-((epsp_b1-E_n)/s_n)*(1/s_n))
            dDepsp_b_df=(-p*D_eps_p+q*D_eps_q)/((1-f1)**2*sig_0)
            df_d_epsp_b=(dA_d_epsp_b*D_epsp_b)/(1.0+D_eps_p) 
            dDepsp_b_d_epsp_b=0.0
            d2g_dp_d_epsp_b= -3.0*q1*q2*( -df_d_epsp_b/sig_0*sinh(-q2*
     &        3.0*p/(2.0*sig_0))-f1/sig_0**2*d_sig0_d_epsp_b*sinh(-q2*
     &        3.0*p/(2.0*sig_0))+f1/sig_0*cosh(-q2*3.0*p/(2.0*sig_0))*
     &        q2*3.0*p/(2.0*sig_0**2)*d_sig0_d_epsp_b) 
            
            d2g_dqdsig_0=(-4.0*q/sig_0**3)
            d2g_dqdf= 0 !d2g_dqdsig_0*d_sig0_d_f

            d2g_dq_d_epsp_b= d2g_dqdsig_0*d_sig0_d_epsp_b

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c    !dH_alpha_dDeps_p = df_dDeps_p + d_epsp_dDeps_p + d_epsp_b_dDeps_p
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  !Compute df_dDeps_p
 
            dDf_dDeps_p=(1.0-f1)            

            dDf_d_epsp_b=dA_d_epsp_b*D_epsp_b+ A_c*dDepsp_b_d_epsp_b
            d_epsp_b_dp=-D_eps_p/((1.0-f1)*sig_0)
            dDf_dp=dDf_d_epsp_b*d_epsp_b_dp

            dDepsp_b_dDeps_p=-p/((1.0-f1)*sig_0)
            dDepsp_b_dp=-D_eps_p/((1.0-f1)*sig_0)

            dDf_df= -D_eps_p  
            c_f_f = 1.0/(Kron_d(1,1)-dDf_df)

            c_f_epsp_b=1.0/(Kron_d(1,3)-dDf_d_epsp_b)
            df_dDeps_p=c_f_f*(dDf_dDeps_p+xk*dDf_dp) 
     &            +c_f_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp)  

            c_epsp_b_f =1.0/(Kron_d(3,1)-dDepsp_b_df) 

            c_epsp_b_epsp_b=1/(Kron_d(3,3)-dDepsp_b_d_epsp_b)

            d_epsp_b_dDeps_p=c_epsp_b_f*(dDf_dDeps_p+xk*dDf_dp)
     &            +c_epsp_b_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp) 

            dH_alpha_dDeps_p= df_dDeps_p + d_epsp_b_dDeps_p 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     dH_alpha_dDeps_q = df_dDeps_q + d_epsp_dDeps_q+ d_epsp_b_dDeps_q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Compute df_dDeps_q

            dDf_dDeps_q=0
            d_epsp_b_dq=D_eps_q/((1-f1)*sig_0)
            !dDf_dq=A_c * D_eps_q/(1.0-f)/sig_0
            dDf_dq=dDf_d_epsp_b*d_epsp_b_dq

            dDepsp_b_dDeps_q= q/((1.0-f1)*sig_0)
            dDepsp_b_dq     = D_eps_q/((1.0-f1)*sig_0)

            df_dDeps_q=c_f_f*(dDf_dDeps_q-3.0*xmu*dDf_dq)
     &         +c_f_epsp_b*(dDepsp_b_dDeps_q-3.0*xmu*d_epsp_b_dq)

c Compute d_epsp_b_dDeps_q
            d_epsp_b_dDeps_q=c_epsp_b_f*(dDf_dDeps_q-3.0*xmu*dDf_dq)
     &        +c_epsp_b_epsp_b*(dDepsp_b_dDeps_q-3.0*xmu*d_epsp_b_dq)
c   Finally compute dH_alpha_dDeps_q 
            dH_alpha_dDeps_q = df_dDeps_q + d_epsp_b_dDeps_q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    dPHI_d
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            dPHI_dp =dg_dp !q1*f*sinh(-3.0*q2*p/(2.0*sig_0))
            dPHI_dq=dg_dq  !2*(q/sig_0**2)
c  Compute all  dPHI_dH_alpha = dPHI_df, dPHI_d_epsp 

            dPHI_df = (2*q1-cosh(-3*q2*p/(2*sig_0))-2*q3*f1)*df_s_df

            dPHI_d_epsp_b= (-2.0*q**2/sig_0**3)*d_sig0_d_epsp_b+2*q1*
     &        (df_d_epsp_b *cosh(-3.0/2.0*q2*p/sig_0) + f_s*sinh
     &        (-3.0/2.0*q2*p/sig_0)*(3.0/2.0*q2*p/sig_0**2)*
     &        d_sig0_d_epsp_b)-2.0*f_s*q3*df_d_epsp_b              ! Option 2 considering relation between f  on epsp_b
     
            PHI=(q/sig_0)**2+2.0*f_s*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f_s**2) 

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
c            A12= dg_dp+D_eps_p*(-3.0*xmu*d2g_d2q+sum1)
c     &             + D_eps_p*(-3*xmu*d2g_dpdq+sum2)  !The article has this error, the derivative is not with respect to dDeps_q
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
            b_v=(/b1,b2/)
            M_inv = inv_T_2D(M)
            c_corr = matmul(inv_T_2D(M), b_v)

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
            
            M_ps = reshape ((/A11, A12, A13_PS,
     &                       A21, A22, A23_PS,
     &                       A31_PS, A32_PS, A33_PS/),
     &                       shape(M_ps), order=(/2,1/))
            b_ps_v =(/b1, b2, b3_PS/) 
            M_ps_inv=inv_T_2D(M_ps)
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
             
            sig=-p*xid+(2.0/3.0)*q*n_a
            dsig=sig-1.0/3.0*xid*Trace(sig)
            D_epsp=1.0/3.0*D_eps_p*xid+D_eps_q*n_a
            epsp=epsp + D_epsp
            epse=eps-epsp  

            Balpha=Balphan
            alpha= alphan
            D_sig =sig-sig1 
            sig1=sig
            dsig1=dsig
            epsp_b=epsp_b_n+D_epsp_b
            f1=f+D_f

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
            df_dp=D_eps_p/(D_epsp_b*sig_0) 
            df_dq=D_eps_q/(D_epsp_b*sig_0) 
            call I4dikdjl(I4dikdjl_r)
            call I4dijdkl(I4dijdkl_r) 
       
            dn_d_sig=1.0/q*(3.0*I4dikdjl_r/2.0-diadic_prod_T2_T2(xid, 
     &         xid)/2.0-diadic_prod_T2_T2(n_a, n_a))

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
        !print*, "mpn/3", dD_eps_pq(1,2)/3 ,"mql", dD_eps_pq(1,2)
        !print*, "partial delta eps_p, partial delta eps_p", c_corr

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Find ATS using patial(D_eps_p) and partial(D_eps_q)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
       M_4T=dD_eps_pq(1,1)*diadic_prod_T2_T2(xid, xid)/3.0
     &     +dD_eps_pq(1,2)*diadic_prod_T2_T2(xid, n_a)/3.0
     &     +dD_eps_pq(2,1)*diadic_prod_T2_T2(n_a, xid)
     &     +dD_eps_pq(2,2)*diadic_prod_T2_T2(n_a, n_a)
     &     +D_eps_q*dn_d_sig 

       C_e=2.0*xmu*I4dikdjl_r-(xk-2.0/3.0*xmu)*I4dijdkl_r

       call invT4(C_e,C_e_inv)

       call invT4(M_4T+C_e_inv,ATS)


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
        use tensor_operations

        implicit none
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
         