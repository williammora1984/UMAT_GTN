c==================================================================  
c subroutine GTN : Material model according to GURSON, TVERGAARD AND 
c NEEDLEMAN (GTN) DAMAGE MODEL  
c------------------------------------------------------------------]
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

c   sig6: sigma_n+1 (stress at n+1 as a vector)
c   A66: algritmic tangent stiffness (ATS) tensor(4 grade tensor as a 6x6 matrix)
c   sdvup: Updated Internal state variables (the updated 3 ISV of the sdvl input vector)   
c==================================================================


      subroutine kGTN (eps6,D_eps6,sdvl,options,inputmat,sig6,A66,sdvup)
       !eps,epsn,epsen
       INCLUDE 'ABA_PARAM.INC'       
        use tensor_operations

        implicit none
c         material parameters
        double precision, dimension(6) :: eps6, D_eps6, epsn6
        double precision, dimension(3,3) :: eps, D_eps, epsn, sig
c             double precision, dimension(:),intent(in) :: sdvl
        double precision, dimension(21) :: sdvl 
        double precision, dimension(21) :: sdvup
        double precision, dimension(13):: inputmat
        !double precision, dimension (size(inputmat,1)):: matp != inputmat
        double precision :: xE     ! Young's modulus
        double precision :: xnu    ! Poisson's ratio
        double precision :: xsigy0 ! initial yield stress
        double precision :: xHk    ! kinematic hardening modulus
        double precision :: xhi     ! isotropic hardening modulus
        double precision :: xmu    !shear modulus
        double precision :: xk     ! bulk modulus
        double precision :: q1, q2, q3, f_0, f_n, s_n, f_f, E_n !Parameters GTN model
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

c     General 
        double precision :: tol =1e-8
        integer, dimension (6) :: ii = (/1,2,3,1,2,1/) !Auxiliar iterators
        integer, dimension (6) :: jj = (/1,2,3,2,3,3/) !Auxiliar iterators
        integer ::i, j, ttype, Mat_L, str_case !Iterators and selectors
        integer, dimension (3) :: options
        double precision, dimension (3,3):: xid !Identity
        
c  Especific variables according to the article N. ARAVAS, 1987
       double precision :: q, p, sig_0
       double precision, dimension (3,3) :: D_eps3
c        s_ve_33
       double precision ::D_eps_p, D_eps_q  !p: hydrotatyx stress and equivalent stress
       double precision:: dg_dp, dg_dq, d2g_d2p, d2g_d2q, d2g_dpdq     
        
        !matp = inputmat
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
        f_f        = inputmat(12)
        E_n        = inputmat(13)


        xmu  = xE/(2.0*(1.0+xnu))    !mu or G: Elastic Shear moduli 
        xk   = xE/(3.0*(1.0-2.0*xnu))  !k: Compression (Bulk) Moduli
          
        eps = reshape((/eps6(1),      eps6(4)/2.0,  eps6(6)/2.0,
     $                  eps6(4)/2.0,  eps6(2),      eps6(5)/2.0,
     $                  eps6(6)/2.0,  eps6(5)/2.0,  eps6(3)/),
     $                  shape(eps), order=(/2,1/))
c    <----------------------------------------------------------------------Change to implement and review               
       epsn6=eps6-D_eps6
       epsn = reshape((/epsn6(1),      epsn6(4)/2.0,  epsn6(6)/2.0,
     $                  epsn6(4)/2.0,  epsn6(2),      epsn6(5)/2.0,
     $                  epsn6(6)/2.0,  epsn6(5)/2.0,  epsn6(3)/),
     $                  shape(eps), order=(/2,1/))
       D_eps=eps-epsn
c    <----------------------------------------------------------------------Change until here to implement and review
c plastic strain at t0 
        !print*,"before copy in function",sdvl
        epspn   = reshape((/sdvl(1),      sdvl(4)/2.0,  sdvl(6)/2.0,
     $                      sdvl(4)/2.0,  sdvl(2),      sdvl(5)/2.0,
     $                      sdvl(6)/2.0,  sdvl(5)/2.0,  sdvl(3)/),
     $                      shape(epspn), order=(/2,1/))
c strain-like ISV that thermodynamically conjugates to the kinematic hardening at time t
        Balphan = reshape((/sdvl(1+6),     sdvl(4+6)/2,   sdvl(6+6)/2.0,
     $                      sdvl(4+6)/2.0, sdvl(2+6),     sdvl(5+6)/2.0,
     $                      sdvl(6+6)/2.0, sdvl(5+6)/2.0, sdvl(3+6)/),
     $                      shape(Balphan), order=(/2,1/))
c strain-like ISV that thermodynamically conjugates to the sotropic hardening at time t
        alphan = sdvl(13)
c Elastic strain        
        epsen  = reshape((/sdvl(1+13),    sdvl(4+13)/2.0,sdvl(6+13)/2.0,
     $                     sdvl(4+13)/2.0,sdvl(2+13),    sdvl(5+13)/2.0,
     $                     sdvl(6+13)/2.0,sdvl(5+13)/2.0,sdvl(3+13)/),
     $                     shape(epspn), order=(/2,1/))
c Void volume fraction         
        fn  = sdvl(20)
c microscopic equivalent plastic strain
        epsp_b_n= sdvl(21)
!Fill the identity and deviatoriser tensor               
        call Ident1(xid,3)
        call P4sym(P4sym_r)
        sig_0=200.00


        Mat_L=options(3) ! Mat_L=0 : VM, Mat_L=1 : GTN
        !print*,"Balphan before trial func", Balphan
        if (Mat_L==0) then
              call trial_step (eps,epspn,xmu,xHk,Balphan,alphan,xhi,
     $                   xsigy0, phitr, dsigtr, ntr, nxitr)       
        else
              call trial_step_GTN (eps,epsn,epsen,xmu,xk,phitr, sig_0,
     &               ntr,q1,q2,q3,fn,dsigtr,sigtr, q_tr,p_tr)
        end if
    
        !print*,"Balphan after trial func", Balphan
        !print*,"Balphan after trial func", Balpha
        !print*, "dsigtr", dsigtr   
        print*, "phitr exit", phitr 
c      Evalate if it is required a elastic or elastic-plastic step
        if (phitr<tol) then  !Tol is a number very close to zero
            print*, "elastic"
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
            print*, "into 2"
            if (ttype==0) then
              Cdev=2.0*xmu*P4sym_r
              
            endif
            
        else ! (plastic correction)
           print '(A)'   
           print*, "enter to material routine plastic return", phitr
           !print*, "phitr enter pl", phitr
c          Mat_L=0 : VM, Mat_L=1 : GTN
           if (Mat_L==0) then
              !Radial return algortihm
              
              call plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
     &           epsp, epspn, Balpha, Balphan, alpha, alphan)
              if (ttype==0) then
                Cdev=2.0*xmu*beta1*P4sym_r
     &                  -2.0*xmu*beta2*diadic_prod_T2_T2(ntr,ntr)
              end if
             
           else
              str_case=options(2)
              !print*, "q_tr", q_tr,"p_tr", p_tr, "fn", fn
             call plas_corre_GTN (q_tr,p_tr, sigtr,inputmat,fn, dsigtr,
     &       phitr, epsp_b,eps, D_eps, epsp, epse, Balpha, alpha,sig_0
     &       ,sig,dsig,f,str_case)            
           end if
        end if
        if (Mat_L==0) then
c             Addition of spherical contribution to stress and moduli:
              sig=dsig+xk*Trace(eps)*xid
              !print*, "sig", sig
              !print*, "eps", eps
              epse=eps-epsp
        end if              
c        epse=eps-epsp

c=================================================================         
c         Analitical solution stiffness tensor
c=============================================================

         ttype=options(1)
        if (ttype==0) then
           
           C=Cdev+xk*diadic_prod_T2_T2(xid,xid)

c               restore stiffness tensor as matrix 
c               The values of these vectors (ii and jj) are equal to the declared at the begining of this file
c                maybe there is a case where I need to have different values
          ii = (/1,2,3,1,2,1/) 
          jj = (/1,2,3,2,3,3/)
c               Due to the symmetry the tangent stiffness 4 grade tensor ca be expressed as a 6x6 matrix
           
          do i=1,6
            do j=1,6
                A66(i,j) = C(ii(i),jj(i),ii(j),jj(j))
            end do 
          end do
        end if

c=====================================================================
c Store information in vectors for postprocessing
c=====================================================================
c     store stress tensor as vector
        do i=1,6
            sig6(i) = sig(ii(i),jj(i))
        end do
        
c    store ISV 
c      plastic strain  
        sdvup(1:6) =(/epsp(1,1),   epsp(2,2),   epsp(3,3),
     $               2.0*epsp(1,2), 2.0*epsp(2,3), 2.0*epsp(1,3)/)
c     Back stress        
        sdvup(7:12)=(/Balpha(1,1),   Balpha(2,2),   Balpha(3,3),
     $               2.0*Balpha(1,2), 2.0*Balpha(2,3), 2.0*Balpha(1,3)/)
c     drag stress
        sdvup(13) = alpha
        sdvup(14:19)=(/epse(1,1),   epse(2,2),   epse(3,3),
     $               2.0*epse(1,2), 2.0*epse(2,3), 2.0*epse(1,3)/)
        sdvup(20)=fn
        sdvup(21)=epsp_b
        !print*,"Balphan before leave func 2", sdvup
        !print*, "xk", xk
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
c==================================================================

      subroutine trial_step (eps,epspn,xmu,xHk,Balphan,alphan,xhi,
     $       xsigy0, phitr, dsigtr, ntr, nxitr)
       INCLUDE 'ABA_PARAM.INC'       
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
       deps = eps - 1.0/3.0*Trace(eps)*xid 
       !deviatoric part of the trial stress tensor
       !print*, "xmu", xmu
       dsigtr = 2.0*xmu*(deps-epspn)
       
       !trial value of the back-stress tensor
       !print*,"Balphan from trial func", Balphan
       Bbetatr = xHk*Balphan  !%---------------------------------> Here was the mistake 
       !trial value of the increase in yield stress (scalar beta: drag stress)
       !print*, "alphan from tiral func", alphan
       betatr =xhi*alphan
       !trial value of deviatoric stress difference
       xitr =dsigtr-Bbetatr
       !norm of xitr (use function t2_contr_t2(A,B) contained in ./tensor/ to compute); The norm is elevated at the square root
       nxitr = (contrac_2nd_2nd(xitr,xitr))**0.5 !in the case of these tensor you must use the double contraction
       !flow direction from trial state
       ntr = xitr/nxitr
       !print*, " \n"
       !print*, "xitr", xitr
       !print*, "nxitr", nxitr
       !print*, "xsigy0", xsigy0
       !print*, "betatr", betatr
       !print*, "part2", (2.0/3.0)**(0.5)*(xsigy0+betatr) 
c             trial yield function  
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
c=======================================================================

      subroutine trial_step_GTN (eps,epsn,epsen,xmu,xk,phitr, sig_0,
     $     n_a_tr,q1,q2,q3,f,dsig_tr,sig_e_tr, q_tr,p_tr)
           INCLUDE 'ABA_PARAM.INC'
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

           C_e=2.0*xmu*Xid!-(xk-2.0/3.0*xmu)*Xid   !<----Take care that here maybe we require the 4th order identity tensor and deviator           
           !Change strain in the step
           D_eps=eps-epsn
           !Elastic predictor (trial). All the strain is supoused to be elastic
           sig_e_tr=C_e*(epsen+D_eps)
           !trial values for p, s and q

c           if (str_case==0) then   
              p_tr= - 1.0/3.0*Trace(sig_e_tr)
              dsig_tr=sig_e_tr-1.0/3.0*Trace(sig_e_tr)*xid
              q_tr=(3.0/2.0*contrac_2nd_2nd(dsig_tr,dsig_tr))**0.5              
c           end if   
           !norm according to the article
           n_a_tr = 3.0/2.0*dsig_tr/q_tr
           sig_test=-p_tr*xid+(2.0/3.0)*q_tr*n_a_tr
           !print*, "sig test using p and q", sig_test
c          trial yield functionn_a_tr
           !print*, "n_a_tr before", n_a_tr
           !print*, "dsig_tr before",dsig_tr   

           phitr = (q_tr/sig_0)**2+2.0*q1*f*cosh(-3.0*q2*p_tr/
     &                (2.0*sig_0))-(1.0+q3*f**2)
          !print*, "phitr inside", phitr, "q_tr", q_tr,"p_tr", p_tr, 
c     &            "fn", f, "sig_0 print",sig_0
      end subroutine trial_step_GTN
           
c=======================================================================  
c subroutine plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
c       epsp, epspn, Balpha, Balphan, alpha, alphan)
c=======================================================================


      subroutine plas_corre_VM (phitr, xmu, xHk, xhi,dsig, dsigtr,ntr,
     &  epsp, epspn, Balpha, Balphan, alpha, alphan)
c       INCLUDE 'ABA_PARAM.INC'     
       implicit none
       
       double precision, dimension (3,3) ::  dsig, dsigtr,ntr,
     &    epsp, epspn, Balpha, Balphan  
       double precision :: phitr, xmu, xHk, alpha, alphan, beta1, beta2,
     &    gamma,xhi, nxitr
       gamma=phitr/(2.0*xmu+xHk+2.0*xhi/3.0)
       !print*, "gamma", gamma
       dsig=dsigtr-2.0*xmu*gamma*ntr
       !print*, "dsig", dsig
       epsp=epspn+gamma*ntr
       !print*, "epsp", epsp
       !print*,"Balphan", Balphan
       Balpha = Balphan+gamma*ntr        ! tensorial internal variable alpha
       !print*,"Balpha changed", Balpha
       alpha = alphan+(2.0/3.0)**(0.5)*gamma           !scalar hardening variable alpha
       !print*, "alpha", alpha
       print*, "nxitr", nxitr !<-------------------------Review, according to this moment its values ero since is not computed neiter is an input argument
       beta1=1.0-(phitr/nxitr)*(1.0/(1.0+xHk/(2.0*xmu)+
     $          xhi/(3.0*xmu)))
       !print*, "beta1", beta1
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
c=======================================================================

      subroutine plas_corre_GTN (q,p, sig,inputmat,f, dsig,phi, 
     &      epsp_b,eps,D_eps, epsp, epse, Balpha, alpha, sig_0, sig1,
     &      dsig1,f1, str_case)
            INCLUDE 'ABA_PARAM.INC'            
            use tensor_operations
            implicit none
c  Especific variables according to the article N. ARAVAS, 1987
            !Parameters material
            double precision, dimension (13)::inputmat
            !Material parameters
            double precision :: q, p, sig_0, q1, q2, q3,f,f1, xE,  
     &           xk, xmu, xnu, q_e, p_tr, q_tr 
            double precision ::A_c  ! Parmeter for nucleation strain Chu and Needleman
            
            !Change of parameter between steps
            double precision ::D_eps_p, D_eps_q, D_epsp_b,d_f, D_eps_pn,
     &       D_eps_qn  !p: hydrotatic stress and equivalent stress, D_epsp_b : delta microscopic eqivalent platic strain
            double precision, dimension (3,3) ::D_epsp,D_epspn1, dsig,
     &       dsig_e, D_eps
            !Current value of ISV
            double precision :: epsp_b, alpha, alphan  ! microscopic equivalent plastic strain
            double precision, dimension (3,3) ::sig ,epsp, eps, epse,
     &                 n_a, sign_l, Balpha, Balphan, sig1,dsig1
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
     &        dDf_d_epsp_b, dDf_df,dDepsp_b_d_epsp_b, NU
            !derivates defined in article of tensorial (3,3) type
            double precision, dimension (3,3) ::  dDepsp_dDeps_p, 
     &         dDepsp_dp,d_epsp_dDeps_p
     &         dDepsp_dDeps_q, dDepsp_dq, d_epsp_dDeps_q, dPHI_d_epsp
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
            double precision :: f_0, f_n, s_n, f_f, E_n ! E_n: mean value Parameter used in the nucleation equations
            !Auxiliar variables
            double precision, dimension (3,3) ::  xid
            double precision :: pi, tol_conv_cor,conv_D_eps_p,
     &          conv_D_eps_q
            double precision, dimension(2,2) :: M, M_inv
            double precision, dimension(2) :: b_v, c_corr, dD_eps_pq 
            double precision, dimension(3,3,3,3)::C, Cdev, C_e, C_e_inv,
     &           M_4T, ATS, dn_d_sig, I4sym_r,I4dikdjl_r,I4dildjk_r
            integer :: iterartion, i, str_case  

c      Variables for plane stress elastoplastic equations  
c           
            double precision ::A23_PS, A31_PS, A32_PS, A33_PS, b3_PS,
     &         A13_PS, D_eps3, dev_e_33, dq_dDeps3, pe, D_eps3n,
     &         conv_d_eps3
            double precision, dimension(3,3) :: M_ps, M_ps_inv,D_eps_b,
     &          dD_eps_b, depse, se, a, da
            double precision, dimension(3) :: b_ps_v, c_corr_ps
            
            !Auxiliar variables
            double precision::conv_D_eps_p_ps, conv_D_eps_q_ps,D_eps11,
     &           D_eps22
             !print*, "q", q,"p", p, "poro",f
             !print*, "stress",sig
             !print*, "mat prop",inputmat
             !print*,"delta_stress", dsig
             !print*,phi 
             !print*,"d_epsp",d_epsp
             print*, "delta epsp_b",D_epsp_b, "epsp_b",epsp_b
             !print*, "epsp ",epsp
             !print*, "epse",epse
             !print*,"Balpha",Balpha
             !print*,"alpha",alpha
             
             

            xE        = inputmat(1)
            xnu       = inputmat(2)  !Poison    
            !xsigy0    = inputmat(3)  !Yield point     
            q1        = inputmat(6)
            q2        = inputmat(7)
            q3        = inputmat(8)
            f_0       = inputmat(9)
            f_n       = inputmat(10)
            s_n       = inputmat(11)
            f_f       = inputmat(12)
            E_n       = inputmat(13)     

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c Constants 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
            xmu  = xE/(2.0*(1.0+xnu))    !(G)Shear modulus 
            xk   = xE/(3.0*(1-2.0*xnu))  !(K) Compression Modulus
            call Ident1(xid,3)
            pi=4.D0*DATAN(1.D0)
            tol_conv_cor=1e-15


          !str_case=0

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c Initial values to start iteration
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            PHI=(q/sig_0)**2+2.0*f*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f**2)
         if (str_case==0) then
c-----------------------------------------------------------------------
c          THREE DIMENSIONAL GEOMETRIES
c-----------------------------------------------------------------------
            print*, "initial phi", PHI,"p_tr", p,  "q_tr", q, 
     &            "fn", f,"sig_0 print",sig_0
c  return direction (only depends on the trial step            )
            dsig_e=sig-Trace(sig)*xid
            print*, "dsig_entry",dsig
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
                   print*, "initial phi", PHI,"p_tr", p,  "q_tr", q, 
     &            "fn", f,"sig_0 print",sig_0
c  return direction (only depends on the trial step            )
            dsig_e=sig-Trace(sig)*xid
            print*, "dsig_entry",dsig
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
            print*,"Befo it", "D_eps_p", D_eps_p, "D_eps_q",D_eps_q
                        
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

c   10     if ((conv_D_eps_p .gt. tol_conv_cor) .or. 
c     $           (conv_D_eps_q.gt. tol_conv_cor )) then       
          do i=1,40
            iterartion=iterartion+1
            print*, "iterartion", iterartion
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
            !sig_0= contrac_2nd_2nd(sig,d_epsp)/((1.0-f)*D_epsp_b)
            !<--------Update subsequent yield stress !I review the code a never evaluates H I thing that is evaluate phi  or ealuate sig_0 solving the equation for this parameter
            p=p_tr+xk*D_eps_p
            q=q_tr-3.0*xmu*D_eps_q

            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b-E_n)/
     &              s_n)**2) !Chu and Needleman parameter A
            d_f=(1.0-f)*D_eps_p+A_c*D_epsp_b
            f=f+d_f  
            print*, "p",p,"q",q, "Ac", A_c, "epsp b", epsp_b, "f", f     
            print*, "Deltas", " D_epsp_b",D_epsp_b, "Df", d_f 
c     &        "D_epsp", D_epsp,       
          else
c-----------------------------------------------------------------------
c                       PLAIN STRAIN CASE
c-----------------------------------------------------------------------
      !D_epsp_b=0.01
c Update Hydrostatic and equivalent stress

            D_epsp_b=(-p*D_eps_p+q*D_eps_q)/(1.0-f)/sig_0
            epsp_b=epsp_b+D_epsp_b
            !sig_0= contrac_2nd_2nd(sig,d_epsp)/((1.0-f)*D_epsp_b)
            !<--------Update subsequent yield stress !I review the code a never evaluates H I thing that is evaluate phi  or ealuate sig_0 solving the equation for this parameter
            
            p=p_tr-xk*(D_eps3-D_eps_p)
            q=-3.0*xmu*D_eps_q+(q_tr**2+6.0*xmu*se(3,3)*D_eps3+  
     &         4.0*xmu**2*D_eps3**2)**0.5   !here q^e is q_tr 
            !p=p_tr+xk*D_eps_p
            !q=q_tr-3.0*xmu*D_eps_q
            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b-E_n)/
     &              s_n)**2) !Chu and Needleman parameter A
            d_f=(1.0-f)*D_eps_p+A_c*D_epsp_b
            f=f+d_f  
            print*, "p",p,"q",q, "Ac", A_c, "epsp b", epsp_b, "f", f     
            print*, "Deltas", " D_epsp_b",D_epsp_b, "Df", d_f 
c     &        "D_epsp", D_epsp,       

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
            print*,"dg_",dg_dq, dg_dp, d2g_dpdq, d2g_d2p,d2g_d2q
            !d2g_dqd_alpha, alpha1 = f, alpha2 = epsp (epsilon plastic)
            !alpha1 = f
            d_sig0_d_f=(-p*D_eps_p+q*D_eps_q)/((1-f)**2*D_epsp_b)
            d_sig0_d_epsp_b=-(-p*D_eps_p+q*D_eps_q)/((1-f)*D_epsp_b**2)  !Option 1 using the equation where is sig_0, D_epsp_b
            NU=0.1
c            d_sig0_d_epsp_b= 3.0*xmu*NU*(sig_0/sig+3.0*xmu*epsp_b/sig)
c     &       **(NU-1)/(1.0-NU*(sig_0/sig+3.0*xmu*epsp_b/sig)**(NU-1)) !option 2
            !d2g_dpdf=-3.0*q1*q2*p*sinh(-3.0*q2*p/(2.0*sig_0))/sig_0  !replaced in the next line
            d2g_dpdf=-3.0*q1*q2*(p/sig_0*sinh(-3.0*q2*p/(2.0*sig_0))-
     &        f*p/sig_0**2*d_sig0_d_f*sinh(-q2*3.0*p/(2.0*sig_0))+
     &        f*p/sig_0*cosh(-q2*3.0*p/(2.0*sig_0))*q2*3.0*p/
     &        (2.0*sig_0**2)*d_sig0_d_f)

            !d2g_dp_d_epsp=0.0

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
            !d2g_dq_d_epsp=0.0   <---------------------Therically this is not required as the other epsp derivatives 
            !d2g_dq_d_epsp_b= d2g_dqdf*df_d_epsp_b  !<---------Previous
            d2g_dq_d_epsp_b= d2g_dqdsig_0*d_sig0_d_epsp_b
c            print*,"dg_ pat2", d2g_dpdf, df_d_epsp_b, d2g_dp_d_epsp_b,
c     &       d2g_dqdf,d2g_dq_d_epsp_b 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
c    !dH_alpha_dDeps_p = df_dDeps_p + d_epsp_dDeps_p + d_epsp_b_dDeps_p
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  !Compute df_dDeps_p
            !dDf_dDeps_p=(1.0-f)/(1.0+D_eps_p) !<---Verify at end Delta_f or only f  
            dDf_dDeps_p=(1.0-f)            
            !dDf_dp=A_c*(-D_eps_p/((1.0-f)*sig_0))  
            dDf_d_epsp_b=dA_d_epsp_b*D_epsp_b+ A_c*dDepsp_b_d_epsp_b
            d_epsp_b_dp=-D_eps_p/((1.0-f)*sig_0)
            dDf_dp=dDf_d_epsp_b*d_epsp_b_dp
            !dDepsp_dDeps_p=(1.0/3.0)*xid <-----Useless
            !dDepsp_dp=contrac_2nd_2nd((D_eps_q*xid),-3.0/2.0*dsig/q**2) <-----Useless
            dDepsp_b_dDeps_p=-p/((1.0-f)*sig_0)
            dDepsp_b_dp=-D_eps_p/((1.0-f)*sig_0)
c            print*,"dD_dh ", dDf_dDeps_p,dDf_dp,dDepsp_b_dDeps_p,  
c     &          dDepsp_b_dp
            dDf_df= -D_eps_p  
            c_f_f = 1.0/(Kron_d(1,1)-dDf_df)
            !c_f_epsp=0.0 !inv_T_2D(xid)
            c_f_epsp_b=1.0/(Kron_d(1,3)-dDf_d_epsp_b)
            df_dDeps_p=c_f_f*(dDf_dDeps_p+xk*dDf_dp) 
     &            +c_f_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp) !<---I'm doing the second term 0, since epsp this is not a H_beta  
c     &            +c_f_epsp*(dDepsp_dDeps_p+xk*dDepsp_dp)*0
c  Compute  d_epsp_dDeps_p
            !c_epsp_f=inv_T_2D(Kron_d(2,1)-1.0/3.0*xid*(-1.0/(f-1.0)+
c     &               (A_c*D_epsp_b+d_f)/(f-1.0)**2)) !(xid) <-----verify in my expression I have explictly f_(n+1) but here I have fn, and most of the part is a 2nd order tensor which does not agrre with Kronecker evaluation
            !c_epsp_epsp=1/Kron_d(2,2) ! The result should be a 4th grade tensor which does not agree with the 0th grade of the Kronecker evaluation 
            !c_epsp_epsp_b=inv_T_2D(Kron_d(2,3)-(1.0/3.0*xid*A_c/ 
c     &       (f-1.0)))!<-----verify in my expression I have explictly f_(n+1) but here I have fn. The result should be a 2th grade tensor which does not agree with the 0th grade of the Kronecker evaluation 
            !d_epsp_dDeps_p=c_epsp_f*(dDf_dDeps_p+xk*dDf_dp)*0
     &      !      +c_epsp_epsp*(dDepsp_dDeps_p+xk*dDepsp_dp)*0
     &      !      +c_epsp_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp)*0 !<---I'm doing this all zero , since epsp this is not a H_beta
c  Compute  d_epsp_b_dDeps_p
            c_epsp_b_f =1.0/(Kron_d(3,1)-dDepsp_b_df) 
            !c_epsp_b_epsp =1/(Kron_d(3,2)-(f-1.0)/A_c)
            c_epsp_b_epsp_b=1/(Kron_d(3,3)-dDepsp_b_d_epsp_b)
c            print*, "cf_ ",c_f_f ,c_f_epsp_b ,"c_epsp_b ", c_epsp_b_f,
c     &          c_epsp_b_epsp_b
            !d_epsp_b_dDeps_p=-p/(1-f)/sig_0
            !The problems are reflected here  <--------------------------------------
            d_epsp_b_dDeps_p=c_epsp_b_f*(dDf_dDeps_p+xk*dDf_dp)
     &            +c_epsp_b_epsp_b*(dDepsp_b_dDeps_p+xk*dDepsp_b_dp) !<---I'm doing the second term 0, since epsp this is not a H_beta
c     &            +c_epsp_b_epsp*(dDepsp_dDeps_p+xk*dDepsp_dp)*0             
c   Finally compute dH_alpha_dDeps_p
            dH_alpha_dDeps_p= df_dDeps_p + d_epsp_b_dDeps_p 
c            print*,"df_dDeps_p", df_dDeps_p,"d_epsp_b_dDeps_p",
c     &            d_epsp_b_dDeps_p, "dH_alpha_dDeps_p",dH_alpha_dDeps_p
c     &                + d_epsp_dDeps_p                     !Still miss to compute second term
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     dH_alpha_dDeps_q = df_dDeps_q + d_epsp_dDeps_q+ d_epsp_b_dDeps_q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Compute df_dDeps_q
            !dDf_dDeps_q=-3.0*Trace(n_a) !Trace of n must be 0
            dDf_dDeps_q=0
            d_epsp_b_dq=D_eps_q/((1-f)*sig_0)
            !dDf_dq=A_c * D_eps_q/(1.0-f)/sig_0
            dDf_dq=dDf_d_epsp_b*d_epsp_b_dq
            !dDepsp_dDeps_q=nxi_a   
            !dDepsp_dq= contrac_2nd_2nd (D_eps_q*xid, -3.0*dsig/
c     &       (2.0*q**2))         <--------------Useless
            dDepsp_b_dDeps_q= q/((1.0-f)*sig_0)
            dDepsp_b_dq     = D_eps_q/((1.0-f)*sig_0)
c           print*,"df_dD", dDf_dDeps_q,dDf_dq,dDepsp_dq,dDepsp_b_dDeps_q 
            df_dDeps_q=c_f_f*(dDf_dDeps_q-3.0*xmu*dDf_dq)
     &         +c_f_epsp_b*(dDepsp_b_dDeps_q-3.0*xmu*d_epsp_b_dq)
c     &        +c_f_epsp*(dDepsp_dDeps_q-3*xmu*dDepsp_dq)
            !!Compute d_epsp_dDeps_q
c            d_epsp_dDeps_q=c_epsp_f*(dDf_dDeps_q-3*xmu*dDf_dq)
c     &        +c_epsp_epsp*(dDepsp_dDeps_q-3*xmu*dDepsp_dq)
c     &        +c_epsp_epsp_b*(dDepsp_b_dDeps_q-3*xmu*d_epsp_b_dq)
c Compute d_epsp_b_dDeps_q
            d_epsp_b_dDeps_q=c_epsp_b_f*(dDf_dDeps_q-3.0*xmu*dDf_dq)
     &        +c_epsp_b_epsp_b*(dDepsp_b_dDeps_q-3.0*xmu*d_epsp_b_dq)
c     &        +c_epsp_b_epsp*(dDepsp_dDeps_q-3*xmu*dDepsp_dq)
c            print*, "df_dDeps_q", df_dDeps_q, "d_epsp_b_dDeps_q",
c     &            d_epsp_b_dDeps_q
c   Finally compute dH_alpha_dDeps_q 
            dH_alpha_dDeps_q = df_dDeps_q + d_epsp_b_dDeps_q
c     &               +d_epsp_dDeps_q
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    dPHI_d
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            dPHI_dp =dg_dp !q1*f*sinh(-3.0*q2*p/(2.0*sig_0))
            dPHI_dq=dg_dq  !2*(q/sig_0**2)
c  Compute all  dPHI_dH_alpha = dPHI_df, dPHI_d_epsp 
            !dPHI_df = 2.0*q1*cosh(-3.0*q2*p/(2.0*sig_0))-2.0*q3*f !Incomplete expression
            dPHI_df = -2*q**2/sig_0**3*d_sig0_d_f + 2.0*q1*(cosh(-3.0*
     &           q2*p/(2.0*sig_0))+f*sinh(-3.0*q2*p/(2.0*sig_0))*
     &           (3.0*q2*p/(2.0*sig_0**2))*d_sig0_d_f)-2.0*q3*f
            !sign_l=(3.0/2.0)**0.5*dsig/(contrac_2nd_2nd(dsig,dsig)) <---useles
c            dPHI_d_epsp =dPHI_dp*(1.0/3.0*contrac_4th_2nd(C_e,xid))
c     &            +dPHI_dq*sign_l           !<-------Verify the right operation for the required contractions between Identity matrix
c            dPHI_d_epsp_b= (-2.0*q**2/sig_0**3+3.0*f*q1*q2*p*sinh
c     &       (-3.0/2.0*q2*p/sig_0)/sig_0**2)*(-(p*D_eps_p+q*D_eps_q)/
c     &       ((1.0-f)*D_epsp_b**2))                      !Incomplete expression

c            dPHI_d_epsp_b= (-2.0*q**2/sig_0**3)*d_sig0_d_epsp_b+2*q1*
c     &        (df_d_epsp_b *cosh(-3.0/2.0*q2*p/sig_0) + f*sinh
c     &        (-3.0/2.0*q2*p/sig_0)*(3.0/2.0*q2*p/sig_0**2)*
c     &        d_sig0_d_epsp_b)-2.0*f*q3*df_d_epsp_b                         ! Option 2 considering depending on epsp_b

            dPHI_d_epsp_b= (-2.0*q**2/sig_0**3)*d_sig0_d_epsp_b+2*q1*
     &        (f*sinh(-3.0/2.0*q2*p/sig_0)*(3.0/2.0*q2*p/sig_0**2)*
     &        d_sig0_d_epsp_b)                         ! Option 1 considering f = cte
     
            PHI=(q/sig_0)**2+2.0*f*q1*cosh(-3.0*q2*p/(2.0*sig_0))
     &                -(1.0+q3*f**2) 
c            print*, "dPHI_d",dPHI_dp,dPHI_dq,dPHI_df,dPHI_d_epsp_b   
            print*, "phi with Deltas", PHI,  "q_tr", q,"p_tr", p, 
     &            "fn", f,"sig_0 print",sig_0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
c Consolidation of Elastoplastic equation constants A11, A12, A21, A22, 
c b1, b2 and solution of LSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   
            sum1=d2g_dqdf*df_dDeps_p+d2g_dq_d_epsp_b*d_epsp_b_dDeps_p
c     $            +d2g_dq_d_epsp*d_epsp_dDeps_p
            sum2=d2g_dpdf*df_dDeps_p+d2g_dp_d_epsp_b*d_epsp_b_dDeps_p
c     $             +d2g_dp_d_epsp*d_epsp_dDeps_p
            A11=dg_dq+D_eps_p*(xk*d2g_dpdq+sum1)+D_eps_q*(xk*d2g_d2p
     &          +sum2)
            sum3=d2g_dqdf*df_dDeps_q+d2g_dq_d_epsp_b*d_epsp_b_dDeps_q
            sum4=d2g_dpdf*df_dDeps_q+d2g_dp_d_epsp_b*d_epsp_b_dDeps_q
c            A12= dg_dp+D_eps_p*(-3.0*xmu*d2g_d2q+sum1)
c     &             + D_eps_p*(-3*xmu*d2g_dpdq+sum2)  !The article has this error, the derivative is not with respect to dDeps_q
            A12= dg_dp+D_eps_p*(-3.0*xmu*d2g_d2q+sum3)
     &             + D_eps_q*(-3*xmu*d2g_dpdq+sum4)  
            sum32=dPHI_df*df_dDeps_p+dPHI_d_epsp_b*d_epsp_b_dDeps_p
c     $              +dPHI_d_epsp*d_epsp_dDeps_p
            A21=xk*dPHI_dp+sum32
            sum42=dPHI_df*df_dDeps_q+dPHI_d_epsp_b*d_epsp_b_dDeps_q
c     $            +dPHI_d_epsp*d_epsp_dDeps_q
            A22=-3.0*xmu*dPHI_dq+sum42
            !b1=-D_eps_p*dg_dq-D_eps_q*dg_dq  !The article has this error, in thes econd term the derivative is not with respect to d_q
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
            print*, "system to solve"
            call print_matrix(M,2,2)
            call print_matrix(b_v,2,1)
            M_inv = inv_T_2D(M)
            c_corr = matmul(inv_T_2D(M), b_v)
            print*, "correction", c_corr

c Update volumetric and deviatoric plastic strain increments
            print*,"Befo it ", "D_eps_p", D_eps_p, "D_eps_q",D_eps_q
            D_eps_p=D_eps_p+c_corr(1) 
            D_eps_q=D_eps_q+c_corr(2)
            print*,"Post it ", "D_eps_p", D_eps_p, "D_eps_q",D_eps_q
c Evaluate step convergence  
            conv_D_eps_p=abs(D_eps_p-D_eps_pn)
            conv_D_eps_q=abs(D_eps_p-D_eps_pn)
            print*, "convergence af"," p",conv_D_eps_p,"q", conv_D_eps_q
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
            print*,"Befo it ", "D_eps_p", D_eps_p, "D_eps_q",D_eps_q
            D_eps_p=D_eps_p+c_corr_ps(1) 
            D_eps_q=D_eps_q+c_corr_ps(2)
            D_eps3=D_eps_q+c_corr_ps(3)

            print*,"Post it ", "D_eps_p", D_eps_p, "D_eps_q",D_eps_q
c Evaluate step convergence  
            conv_D_eps_p=abs(D_eps_p-D_eps_pn)
            conv_D_eps_q=abs(D_eps_p-D_eps_pn)
            conv_D_eps3=abs(D_eps3-D_eps3n)
            print*, "convergence af"," p",conv_D_eps_p,"q", conv_D_eps_q
            D_eps_pn=D_eps_p
            D_eps_qn=D_eps_q
            D_eps3n=D_eps3

         end if
         end do
c        goto 10
c        endif

c  Update hydrostatic and equivalent stresses 
            !I don't agree to make this update 
c Update stress and strain tensors
            
            print*, "sig bef", sig  
            sig=-p*xid+(2.0/3.0)*q*n_a
c            sig=sig-xk*D_eps_p-2*xmu*D_eps_q*n_a   !<-------------------Try this option for sigma update
c            print*, "n_a with step", n_a
c            print*, "sig aft", sig       
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
         !if (str_case==1) then

         !end if


c=======================================================================
c                 Compute Material tangent stiffness
c=======================================================================
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Additional derivatives 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            df_dp=1.0  !<--------Finish
            d_epsp_dp=1.0
            !d_epsp_b_dp=1.0 !it was moved above
            df_dq=1.0  !<--------Finish
            d_epsp_dq=1.0
            !d_epsp_b_dq=1.0 !it was moved above
            dou_sum1=(D_eps_p*d2g_dqdf+D_eps_q*d2g_dpdf)*(c_f_f*
     &         df_dDeps_p+c_f_epsp_b*d_epsp_b_dDeps_p)+
     &         (D_eps_p*d2g_dq_d_epsp_b+D_eps_q*d2g_dp_d_epsp_b)* 
     &         (c_epsp_b_f*df_dDeps_p+
     &         c_epsp_b_epsp_b*d_epsp_b_dDeps_p)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Consolidation of constants to compute LSE to solve D: A11, A12, 
c  A21, A22, B11, B12, B21, B22 and solution of LSE 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c            dou_sum1_P1=(D_eps_p*d2g_dqdf+D_eps_q*d2g_dpdf)
c            dou_sum1_P2=(/c_f_f*df_dDeps_p+c_f_epsp*d_epsp_dDeps_p+
c     &      c_f_epsp_b*d_epsp_b_dDeps_p,  c_epsp_f*df_dDeps_p +
c     &       c_epsp_epsp*d_epsp_dDeps_p+ c_epsp_epsp_b*d_epsp_b_dDeps_p,
c     &        c_epsp_b_f*df_dDeps_p+c_epsp_b_epsp*d_epsp_dDeps_p+
c     &        c_epsp_b_epsp_b*d_epsp_b_dDeps_p/)   
c     &       A11_D=dg_dq+dou_sum1_P1*dou_sum1_P2
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

       M = reshape((/A11_D, A12_D,
     $               A21_D, A22_D/),shape(M), order=(/2,1/))
        b1_D=contrac_2nd_2nd(B11_D*xid+B12_D*n_a,D_sig)
        b2_D=contrac_2nd_2nd(B21_D*xid+B22_D*n_a,D_sig)
        b_v=(/b1_D,b2_D/)
c       print*, "system to solve"
c       call print_matrix(M,2,2)
c       call print_matrix(b_v,2,1)
        M_inv = inv_T_2D(M)
        dD_eps_pq = matmul(inv_T_2D(M), b_v)
       print*, "partial delta eps_p, partial delta eps_p", c_corr

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Find ATS using patial(D_eps_p) and partial(D_eps_q)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       dD_eps_p=dD_eps_pq(1)
       dD_eps_q=dD_eps_pq(2)
       call I4sym(I4sym_r)
       !
       dn_d_sig=(1.0/q)*(3.0/2.0*I4sym_r)-1.0/2.0*diadic_prod_T2_T2(xid,
     &       xid)- diadic_prod_T2_T2(n_a,n_a)
       !E_p=
       M_4T=1.0/3.0*diadic_prod_T2_T2(dD_eps_p*inv_T_2D(D_sig),xid)
     &       +diadic_prod_T2_T2(dD_eps_q*inv_T_2D(D_sig),n_a)
     &       +D_eps_q*dn_d_sig
       call I4dikdjl(I4dikdjl_r)
       call I4dildjk(I4dildjk_r)  !<---------------------------Note that in the artiicle there is a little difference with the indices
       C_e=2.0*xmu*I4dikdjl_r-(xk-2.0/3.0*xmu)*I4dildjk_r
       call invT4(C_e,C_e_inv)
       call invT4(M_4T+C_e_inv,ATS)
           
      end subroutine plas_corre_GTN

      