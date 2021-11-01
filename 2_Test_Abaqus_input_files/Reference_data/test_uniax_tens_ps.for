      program test_uniax_tens_ps
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
c        initial values
        
c=======================================================================
c                     Unixial tension Plane stress
c=======================================================================
c test on a plane stress element
c the direction of stressing is at 45” with the X-axis
c Test test is displacement controlled, with equal increments of 2.5 per cent the original length (I understand that this is the Abaqus job)

cThe set of non-linear equations (32)-(34) is integrated using a forward Euler scheme with equal
cstrain increments of 1/1000 of the yield strain        

        D_u=0.2/10000.0
        
        l0=0.004
        Teps_f=0.38  !1.0 maximum true strain 
        lf=(DEXP(Teps_f)-1)*l0+l0
        eps_fi=(lf-l0)/l0
        steps=(lf-l0)/D_u
        print*, "steps" , steps
        
        i=0
        eps_y=sigy0/xE !Strain at yield point 
        l=0
        l=l+eps_y*l0+l0  !Length at yied point
        lsi=l !longitud od start of iteration
        Teps=dlog(1+eps_y)  !This must be the logarithmic strain defined as eps_t=ln (1+eps)
        print*, "eps_y", eps_y,"Teps_y",Teps , "Length at yied point",  
     &       l, "Final lengt", lf, "eps max"
     &     ,eps_fi
        sig=sigy0 !This must be the macroscopic axial true stress defined as sig_t=sig*(1+eps)
        sig_0=sigy0
        epsp_b=0.0
        open(unit=1, iostat=stat, file='data_test_ut_ps.csv', 
     $         status='old')
             if (stat == 0) close(1, status='delete')
        open(1, file = 'data_test_ut_ps.csv', status = 'new')
        !u=0.0
   10   if (l<lf) then 
            i =i+ 1

            u=D_u*i
            l=lsi+u
            print*, "l",l
            Teps1=dlog(1+u/l0+eps_y)
            D_eps=Teps1-Teps
            ome=sig/sig_0
            alp=0.5*q1*q2*f*sinh(0.5*q2*ome)
            gam=q1*cosh(0.5*q2*ome)-q3*f
            !print*, "ome", ome, "alp", alp,"gam" , gam
            d_sig_0_d_epsp_b=3.0*mu*NU*(sig_0/sigy0+3.0*mu*epsp_b/sigy0)
     &              **(NU-1.0)/(1.0-NU*(sig_0/sigy0+
     &              3.0*mu*epsp_b/sigy0)**(NU-1))
            A_c=f_n/(s_n*(2.0*pi)**0.5)*DEXP(-1.0/2.0*((epsp_b-E_n)/
     &           s_n)**2)
            H=d_sig_0_d_epsp_b*(ome**2+alp*ome)*(ome**2+alp*ome-sig_0*
     &             A_c*gam)-3.0*sig_0*(1.0-f)*alp*gam
            C=1.0+alp/ome+H/(xE*ome*(ome+alp))
            !print*, d_sig_0_d_epsp_b, "A_c", A_c, "H", H, "C", A_c

            dsig_deps=H/(C*ome*(ome+alp))
            df_deps=3.0*alp*(1.0-f)/(C*ome)+A_c*(ome+alp)/(C*(1.0-f))
            depsp_b_deps=(ome+alp)/(C*(1.0-f))

            sig =sig+ dsig_deps*D_eps
            f =f+ df_deps*D_eps
            epsp_bn =epsp_b
            epsp_b =epsp_b+ depsp_b_deps*D_eps
            D_epsp_b=epsp_b-epsp_bn
            sig_0=sig_0+d_sig_0_d_epsp_b*D_epsp_b
            Teps=Teps1
c In Figure 4, the stress is plotted as a function of the logarithmic strain
c In Figure 5, The void volume fraction,f, is plotted as a function of strain 
            !output data into a file as follow: 
              !eps       = strain
              !sig/sigy0 = realtion stress 
              !f         = void volume fraction
              !epsp_b    = equivalent plastic strain
            write(1,*) Teps, sig/sigy0, f, epsp_b
            goto 10
        end if
        close(1)

      end program test_uniax_tens_ps
        
