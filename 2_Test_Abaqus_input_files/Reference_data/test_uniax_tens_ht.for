      program test_uniax_tens_ht
        use tensor_operations

        implicit none
        double precision::  sig, eps, H, C, ome,alp, gam, sig_0, epsp_b,
     &      epsp_bn, D_epsp_b,A_c, f, Teps
        !Parameters GTN mode
        double precision:: xE, xnu, sigy0, q1, q2, q3, f_0, f_n, s_n, 
     &   f_f, E_n, NU, mu, q_el, kap
        double precision :: d_sig_0_d_epsp_b, dsig_deps, df_deps, 
     &        depsp_b_deps
        double precision :: dp_deps_v, df_deps_v, depsp_b_deps_v, epsv_y
        double precision, dimension(13):: mat_param
        double precision, dimension(10):: epsv_l, p_l, f_l
        double precision:: l0, uf, Teps_f, D_u ,u, Teps1, D_eps, eps_fi,
     &   l, lsi,lf, l_y
        double precision:: pi, steps
        double precision:: ph, eps_vf, epsv, epsv1, D_epsv, eps_y
        integer :: i, stat
        

        mat_param(1)   = 300.0*200.0              ! xE
        mat_param(2)   = 0.33                 ! xnu 
        mat_param(3)   = 200.0 !200           ! xsigy0
        mat_param(4)   = 50.0*mat_param(3);   ! 50   ! xH
        mat_param(5)   = 100.0*mat_param(3);  ! 10  ! xh
        mat_param(6)   = 1.5                  ! q1 
        mat_param(7)   = 1.0                  ! q2
        mat_param(8)   = 1.5                  ! q1=q3 Aricle G.Vadillo
        mat_param(9)   = 0.0                ! f_0
        mat_param(10)  = 0.04                  ! f_n
        mat_param(11)  = 0.1                  ! s_n
        mat_param(12)  = 0.2025               ! f_f
        mat_param(13)  = 0.3                  ! E_n

        xE       = mat_param(1)
        xnu      = mat_param(2) 
        sigy0    = mat_param(3)
        q1       = mat_param(6) 
        q2       = mat_param(7)
        q3       = q1**2
        f        = mat_param(9)
        f_n      = mat_param(10)
        s_n      = mat_param(11)
        f_f      = mat_param(12)
        E_n      = mat_param(13)
        mu       = xE/(2.0*(1.0+xnu))
        kap      = xE/(3.0*(1.0-2.0*xnu))
        q_el     = -0.5*(kap-2.0/3.0*mu)/(kap+1.0/3.0*mu)
        NU       = 0.1

        pi       = 4.D0*DATAN(1.D0)

c=======================================================================
c        Hydrostatic tension
c=======================================================================
c test on an 8-node brick element
c Equal displacement increments are used so that a logarithmic volumetric strain of 40 per cent is reached in 10 increments

        f=0.04
        D_u=0.2/10000.0

        l0=0.004  !0.25 Works !Dimenision initial test piece
        eps_vf=0.4  !0.4 Final volumetric strain
c Without include the strain until yield point
        epsv_l=(/0.042, 0.084, 0.126, 0.166, 0.207, 0.246, 0.285, 0.324, 
     &          0.362, 0.400/) 
        p_l=(/201.4, 172.8, 148.7, 128.8, 112.1, 98.1, 86.3, 76.3, 
     &        67.8, 60.4/) 
        f_l=(/0.0728, 0.1130, 0.1524, 0.1911, 0.2287, 0.2645, 0.2980, 
     &       0.3290, 0.3575, 0.3838/)
        uf=(DEXP(eps_vf/3)-1)*l0
        lf=(DEXP(eps_vf/3)-1)*l0+uf !Final length 
        steps=uf/D_u
        print*, "steps" , steps
        eps_fi=(lf-l0)/l0 !Final engineering strain
        u=0.0
        i=0
        eps_y=sigy0/xE !Strain at yield point
        l_y=eps_y*l0+l0  !Length at yied point
        epsv_y=3.0*dlog(1+eps_y)  !Volumetric strain
        !epsv=0.00  !Initial Volumetric strain
        epsv=epsv_y !<-----------------------------Alternative
        ph=sigy0 !0.1 !Initial hidrostatic stress. Here we supouse that start a yield point
        sig_0=sigy0 !Equivalent tensile flow stress 
        epsp_b=0.0 
        print*, "eps_y", eps_y,"eps_vol_y",epsv ,"Length at yied point",  
     &         l, "Final lengt", lf, "eps max"
     &       ,eps_fi

        open(unit=1, iostat=stat, file='data_test_ht_vol.csv', 
     $         status='old')
             if (stat == 0) close(1, status='delete')
        open(1, file = 'data_test_ht_vol.csv', status = 'new')
   20   if (u<uf) then 
            i =i+ 1

            u=D_u*i
            !epsv1=3.0*dlog(1+u/l0)
            epsv1=3.0*dlog(1+u/l0+eps_y) !<-----------------------------Alternativ
            D_epsv=epsv1-epsv 

            !ome=-3.0*ph/sig_0 Aravas
          ome=3.0*ph/sig_0 !Z.L. Zhang
          alp=0.5*q1*q2*f*sinh(0.5*q2*ome)
          !gam=q1*cosh(0.5*q2*ome)-q3*f !Aravas
          gam=q1*cosh(-1.5*q2*sig/sig_0)-q3*f !Z.L. Zhang
          print*, "ome", ome, "alp", alp,"gam" , gam
          d_sig_0_d_epsp_b=3.0*mu*NU*(sig_0/sigy0+3.0*mu*epsp_b/sigy0)
     &              **(NU-1.0)/(1.0-NU*(sig_0/sigy0+
     &              3.0*mu*epsp_b/sigy0)**(NU-1))
          A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b-E_n)/
     &           s_n)**2)
c            H=d_sig_0_d_epsp_b*alp*ome/(1.0-f)*(alp*ome-sig_0*A_c*gam)-
c     &         3.0*sig_0*(1.0-f)*alp*gam  !Aravas
          H=alp*ome/(1.0-f)*(d_sig_0_d_epsp_b*alp*ome-sig_0*A_c*gam)-
     &         3.0*sig_0*(1.0-f)*alp*gam   !Z.L. Zhang
          !print*, d_sig_0_d_epsp_b, "A_c", A_c, "H", H, "A_c", A_c


          !dp_deps_v=-kap*H/(H+9*alp**2*kap) Aravas
          dp_deps_v=kap*H/(H+9*alp**2*kap) !Z.L. Zhang
          depsp_b_deps_v=3*alp**2*kap/(H+9*alp**2*kap)*ome/(1-f)
          df_deps_v=9*alp**2*kap/(H+9*alp**2*kap)*(1-f+1.0/3.0*A_c*
     &          ome/(1-f))
          ph =ph+ dp_deps_v*D_epsv
          f =f+ df_deps_v*D_epsv
          epsp_bn =epsp_b
          epsp_b =epsp_b+ depsp_b_deps_v*D_epsv
          D_epsp_b=epsp_b-epsp_bn
          sig_0=sig_0+d_sig_0_d_epsp_b*D_epsp_b
          epsv=epsv1
            
c Figure 7, hydrostatic stress as a function of the logarithmic volumetric strain
c Figure 8, void volume fraction as a function of the logarithmic volumetric strain
            !output data into a file as follow: 
              !epsv       = volumetric strain
              !ph = realation hidrostatic stress to yield stress 
              !f         = void volume fraction
              !epsp_b    = equivalent plastic strain
            write(1,*) epsv, ph, f, epsp_b
            goto 20
        end if
        close(1)
        
c=======================================================================
c        Hydrostatic tension
c=======================================================================

      end program test_uniax_tens_ht
        
