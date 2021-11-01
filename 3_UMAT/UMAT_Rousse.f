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
      DIMENSION EPLAS(6), SE(NTENS), SDE(NTENS), S(NTENS), DEP(NTENS), 
     2 EET(NTENS), DDE(NTENS,NTENS), DDP(NTENS,NTENS),
     3 SA(NTENS),HEF(2,4), CC(2,2), FQ(4), DW(6), DJ(6), FP(4)
c======================================================================
c              Additional variables used in this UMAT file
c=======================================================================
      PARAMETER (ZERO = 0.D0,OP5 = 1.5D0,ONE = 1.0D0,TWO = 2.0D0,
     1     THREE = 3.0D0,FOUR = 4.0D0,SIX = 6.0D0,
     2     TOL = 1.D-14,MNI= 100)
       integer:: options(3)    
       DATA DW/1.0,1.0,1.0,0.5,0.5,0.5/
       DATA DJ/1.0,1.0,1.0,0.0,0.0,0.0/

C -----------------------------------------------------------

        call UEXTERNALDB(1,0,TIME,DTIME,KSTEP,KINC)
        WRITE(16,*) stran(1:3),stress, statev,KINC,KINC,KINC
c------------------MATERIAL PARAMETERS----------------------------------
        xE = props(1)
        POISSON = props(2)
        YS = props(3)
        S1 = props(4)
        D = props(5)
        F0 = props(6)


        options(1) = 0
        options(3) = 0

c ---Assign 3D or plane stress program
  
        if (NTENS==6) then
           options(2) = 0
        elseif (NTENS==3) then
           options(2) = 1
        else  
          WRITE(6,2)
 2       FORMAT(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',
     &          'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
       end if

c Compute the values of elastic properties for UMAT
        c3k = xE/(one -two*poisson)
        c2g=xE/(1.0+poisson)
        g=c2g/2.0
        c3g=three*g
        clame=(c3k-c2g)/three
        c1k=c3k/three
C Calculate a number of ordered pair of hardening data
        nvalue=(nprops/2)-3
c recover the value of elastic(eet)and plastic(eplas) strain
c tensor at the begining of increment
        call rotsig(statev(1), drot,eet,2,ndi,nshr)
        call rotsig(statev(1+ntens),drot,eplas,2,ndi,nshr)

c recover the value of equivalent plastic strain (eeqpt),
c void volume fraction (ft) and relative density (rho)
c at the begining of increment
        eeqpt = statev(1 + 2*ntens)
        ft = statev(2 + 2*ntens)
        rhot = statev(3 + 2*ntens)
 
        if(eeqpt.eq.0.0) then
        ft = f0
        rhot = one
        end if      

c define stiffness matrix,[k] or jacobian matrix
c for elastic case(dde)
        do 200 k1 = 1, ntens
          do 100 k2 = 1, ntens
             dde(k2,k1) = 0.0
             ddp(k2,k1) = 0.0
             ddsdde(k2,k1) = 0.0
  100     continue
  200  continue
          do 400 k1 = 1, ndi
            do 300 k2 = 1, ndi
                dde(k2,k1) = clame
  300       continue
            dde(k1,k1) = c2g + clame
  400     continue
        do 500 k1 = ndi + 1, ntens
          dde(k1,k1) = g
  500 continue
      
c compute trial stress tensor from elastic strains
c increment(dstran)
       

c compute trial stress tensor from elastic strains
c increment(dstran)

       do 700 k1 = 1, ntens
         do 600 k2 = 1, ntens
           stress(k2)= stress(k2)+ dde(k2,k1)*dstran(k1)
  600 continue
  700 continue       
c calculate equivalent values(seq) and mean value(sm) of trial
c stress tensor(stress)
       call sinv(stress,sm,seq,ndi,nshr)
       
c determine subsequent yield stress(seqo)
        call khard(h,dh,eeqpt,props(7),nvalue)

c check if the stress value is beyond the yield surface
       sigdif=(seq/rhot)-h+ft*s1*d*exp(sm/(s1*rhot)) 
        if ( sigdif .le. zero ) then
c  in case of elastic (not yielding)
c store elastic(eet), plastic(eplas) strain tensor
c in state variable array
        do 800 k1= 1, ntens
        statev(k1)=eet(k1)
        statev(k1+ntens) =eplas(k1)
  800      continue
  
c store value of equivalent plastic strain (eeqpt),
c void volume fraction (ft) and relative density (rho)
c in state variable array
        statev(1 + 2*ntens) = eeqpt
        statev(2 + 2*ntens) = ft
        statev(3 + 2*ntens) = rhot

c set jacobain matrix(ddsdde) of current increment
c = jacobian matrix of elastic part(dde)
        ddsdde = dde
        else
          
c in case of plastic (yielding)

        if(time(1).eq.0.0.and.(seq.ne.0.0.and.sigdif.gt.10000.0))then
        pnewdt = 0.5
        print*,"The subroutines finish due to wrong setting of parameters"
        go to 5000
        endif

c define new parameters to store the trial stress tensor(se)
        se = stress

c Determine mean value(SME)and equivalent value(SEQE) of
c trial stress tensor(SE)
        call sinv(se,sme,seqe,ndi,nshr)
c Determine hydrostatic pressure(PE)
        pe=-sme

C Calculate deviatoric component(SDE) of SE
        sde=se
        do 1000 i=1,ndi
           sde(i) =sde(i)-sme
 1000   continue
c Set the initial value of equivalent(DEEQP) and
c mean value(DEMP) of plastic strain increment tensor
        demp = 0.D0
        deeqp = 0.D0
c Set the initial value of void volume fraction(DF)
        DF=0.D0
c Set the initial value of iteration number(ITER)
        ITER =0
c---------------------------------------------------------------
c starting iterative process for computing the increment values
c of mean value(cm) and equivalent value(ceq) of plastic strain
c increment tensor
c---------------------------------------------------------------

c update iteration number
 1500   iter=iter+1
c
c update equivalent value of plastic strain(eeqp) for
c current iteration
        eeqp=eeqpt+deeqp
c calculate subsequent yield stress(h) and slop on the hardening
c curve(dh) for current iteration
        call khard(h,dh,eeqp,props(7),nvalue)
c
c compute mean(sm) and equivalent values(seq) of stress tensor
c for current iteration using eq.(3.21)
        p=pe+c1k*demp
        seq=seqe-c3g*deeqp
        if(seq.le.0.0) seq=0.d0

c determine the increment of void volume fraction(df)
c for current iteration using newton-raphson method
         call kfi(ft, demp, df)
c update the values of void volume fraction(f)
c for current iteration
        f=ft + df
        rho=(one-f)/(one-f0)
        fzeb=(one-f0)/((one-f)*(one-f))
        epr=exp(-p/(rho*s1))
        depr =d*epr

        ff0=(seq/rho)-h+f*s1*depr
        fq0old=one/rho
        fp0old=-(f/rho)*depr
        fem=-dh
        ff=seq*fzeb+depr*(s1-f*p*fzeb)

        hef(1,1)=zero
        hef(1,2)=one
        hef(1,3)=zero
        hef(1,4)=zero
      
        hef(2,1)=one-f
        hef(2,2)=zero
        hef(2,3)=zero
        hef(2,4)=zero
        
        emh1=zero
        emh2=zero
        fh1=zero
        fh2=-demp

        divid=(one-emh1)*(one-fh2)-(emh2*fh1)
        cc(1,1)=(one-fh2)/divid
        cc(2,2)=(one-emh1)/divid
        cc(1,2)=emh2/divid
        cc(2,1)=fh1/divid
        
        h1ep=0.d0
        h2ep=0.d0
        h1eq=0.d0
        h2eq=0.d0
        
        do 1600 i = 1,2
          h1ep = h1ep + cc(1,i)*(hef(i,1) + c1k*hef(i,3))
          h2ep = h2ep + cc(2,i)*(hef(i,1) + c1k*hef(i,3))
          h1eq = h1eq + cc(1,i)*(hef(i,2) - c3g*hef(i,4))
 1600     h2eq = h2eq+cc(2,i)*(hef(i,2)-c3g*hef(i,4))
c        
        fp0= -(f/rho)*depr
        fp(1) = (f/(s1*rho*rho))*depr
        fp(2) = zero
        fp(3) = zero

        fp(4) = -(one/(s1*rho*rho))*depr
     1         *(s1*(rho +(f/(one-f0)))-rho*f*p*fzeb)
        
        fq0 =one/rho
        fq(1) = zero
        fq(2) = zero
        fq(3) = zero
        fq(4) = fzeb
        
        A11 = fq0+demp*(c1k*fq(1)+fq(3)*h1ep+fq(4)*h2ep)
     1   + deeqp*(c1k*fp(1)+fp(3)*h1ep+fp(4)*h2ep)
        
        A12 = fp0+demp*(-c3g*fq(2)+fq(3)*h1eq+fq(4)*h2eq)
     1    + deeqp*(-c3g*fp(2) + fp(3)*h1eq + fp(4)*h2eq)
        
        B11 = - demp*fq0 - deeqp*fp0
        A21 = c1k*fp0old + fem*h1ep + ff*h2ep
        A22= - c3g*fq0old + fem*h1eq + ff*h2eq
        B2=- ff0
        
c Determine the correction values for volumetric(CM) and
c deviatoric(CEQ) component of plastic strain increment tensor
c for current iteration
        DET = A11*A22 - A12*A21
        cm = (A22*B1 - A12*B2)/DET
        ceq = (A11*B2 - A21*B1)/DET
c update mean value(demp) and equivalent value(deeqp) of
c plastic strain increment tensor for current iteration
        demp = demp+ cm
        deeqp = deeqp + ceq

c Checking whether the number of current iteration exceed the
c maximum allowed number
        if (iter.lt.mni) goto 1700
        WRITE(102,*)'THE ITERATION CANNOT CONVERGE'
        STOP

c Start new iteration, if the increment values of volumetric or
c deviatoric part of plastic strain tensor is not less than
c the tolerance
 1700   if((abs(cm).gt.tol).or.(abs(ceq).gt.tol)) goto 1500
c Update final values of hydrostatic pressure(p)and equivalent
c value(seq) of stress tensor(p) for current increment
        p = pe + c1k*demp
        seq = seqe - c3g*deeqp

c Update the final values of stress tensor(s) and
c plastic strain increment(dep) for current increment
        s = 0.d0
        dep = 0.d0
        s = seq*(sde/seqe)
        dep = 1.5*(sde/seqe)*deeqp
        do 1900 i= 1, ndi
        s(i) = s(i) - p
        dep(i) = dep(i) + demp/three
 1900   continue
c store the stress tensor of current increment to variable
c ‘stress’
        stress = s

c update plastic strain tensor(eplas) at the end of increment
        do 2000 k1 = 1, ntens
          if (k1.le.ndi) then
          eplas(k1) = eplas(k1) + dep(k1)
          eplas(k1) = eplas(k1) + two*dep(k1)
          else
          end if
 2000   continue
  
c update elastic strain tensor(eet) at the end of time step
c using eq.(3.35)
        eet = eet - dep
  
c update the final values of equivalent plastic strain(eeqpt)
c and void volume fraction (ft) for current increment
        eeqpt =eeqpt+ deeqp
        ft = ft + df
  
  
c consistent tangent modulus for plastic part
c---------------------------------------------------------------
  
c determine the flow direction(sa) using eq.(3.20)
        do 2200 i = 1, ntens
        sa(i) =op5*(sde(i)/seqe)
 2200 continue

        h1ep = zero
        h1eq = zero
        h1p = zero
        h1q = zero
        h2ep = zero
        h2eq = zero

        h2p = zero
        h2q = zero
        
        do 2300 i = 1,2
        h1ep = h1ep + cc(1,i)*hef(i,1)
        h1eq = h1eq + cc(1,i)*hef(i,2)
        h1p = h1p + cc(1,i)*hef(i,3)
        h1q = h1q + cc(1,i)*hef(i,4)
        
        h2ep = h2ep + cc(2,i)*hef(i,1)
        h2eq = h2eq + cc(2,i)*hef(i,2)
        h2p = h2p + cc(2,i)*hef(i,3)
        h2q = h2q + cc(2,i)*hef(i,4)
        
 2300   continue
        
        AA11 = fq0 + (demp*fq(3) + deeqp*fp(3))*h1ep
     &   + (demp*fq(4) + deeqp*fp(4))*h2ep
        AA12 = fp0 + (demp*fq(3) + deeqp*fp(3))*h1eq
     & + (demp*fq(4) + deeqp*fp(4))*h2eq
        
        AA21 = fem*h1ep + ff*h2ep
        AA22 = fem*h1eq + ff*h2eq
        
        BB11=(demp/three)*(fq(1)+ fq(3)*h1p + fq(4)*h2p) +
     &      (deeqp/three)*(fp(1)+ fp(3)*h1p + fp(4)*h2p)
        BB12 = -(demp)* (fq(2) + fq(3)*h1q + fq(4)*h2q) -
     &   (deeqp)*(fp(2) + fp(3)*h1q + fp(4)*h2q)
        BB21 = (fp0 + fem*h1p + ff*h2p)/three
        BB22 = -(fq0 + fem*h1q + ff*h2q)
C Determine coefficients in Eq.(3.66) using formula given
C in Appendix III

        DIV=(AA11 + C3K*BB11)*(AA22 + C3G*BB22) -
     &      (AA12 + C3G*BB12)*(AA21 + C3K*BB21)
        
        C11 = ((AA22 + C3G*BB22)*BB11 - (AA12 + C3G*BB12)*BB21)/DIV
        C21 = ((AA11 + C3K*BB11)*BB21 - (AA21 + C3K*BB21)*BB11)/DIV
        C12 = ((AA22 + C3G*BB22)*BB12 - (AA12 + C3G*BB12)*BB22)/DIV
        C22 = ((AA11 + C3K*BB11)*BB22 - (AA21 + C3K*BB21)*BB12)/DIV

C Evaluate the four coefficients given by Zhange using Eq.(3.70)
        d0 = c2g*(seq/seqe)
        d1 = c1k-(c2g/three)*(seq/seqe)-three*c1k*c1k*c11
        d2 = ((four*g*g*deeqp)/seqe)-four*g*g*c22
        d3 = -(c2g)*c1k*c12
        d4 = -(six*g)*c1k*c21

        do 2600 i = 1, ntens
           do 2500 j = 1, ntens
           dij = 0.0
           if(i.eq.j) dij = 1.0
              ddp(i,j) = dij*dw(i)*d0 + d1*dj(i)*dj(j) +
     &        d2*sa(i)*sa(j) + d3*dj(i)*sa(j) +
     &        d4*sa(i)*dj(j)

 2500   continue
 2600   continue
c set jacobain matrix(ddsdde) of current increment =
c jacobian matrix of plastic part(ddp)
        ddsdde = ddp
c store elastic(eet), plastic(eplas) and equivalent value of
c plastic strain(eeqpt)in state variable array
        do 2800 k1 = 1, ntens
          statev(k1) = eet(k1)
          statev(k1 + ntens) = eplas(k1)
 2800   continue
        statev(1+ 2*ntens)= eeqpt
        statev(2 + 2*ntens) = ft
        statev(3 + 2*ntens) = rhot
      endif
      print*,"End of the routine, look the msg and dat files"
 5000    return
      end subroutine UMAT


      subroutine khard(syield,hard,eqplas,table,nvalue)
        include 'aba_param.inc'
        dimension table(2,nvalue)

c set yield stress to last value of table, hardening to zero
        syield=table(1,nvalue)
        hard=0.0

c if more than one entry, search table
        if(nvalue.gt.1) then
          do 10 k1=1,nvalue-1
            eqpl1=table(2,k1+1)
            if(eqplas.lt.eqpl1) then
              eqpl0=table(2,k1)
              if(eqpl1.le.eqpl0) then
                 write(6,1)
 1             format(//,30x,'***error - plastic strain must be ',
     1                   'entered in ascending order')
                 call xit
              endif
C Current Yield Stress And Hardening
              Deqpl=Eqpl1-Eqpl0
              Syiel0=Table(1,K1)
              Syiel1=Table(1,K1+1)
              Dsyiel=Syiel1-Syiel0
              Hard=Dsyiel/Deqpl
              
              Syield=Syiel0+(Eqplas-Eqpl0)*Hard 
            Goto 20
            Endif
  10    Continue
  20    Continue  
        Endif
        Return
      End Subroutine Khard


      subroutine kfi(ft, demp, df)
        include 'aba_param.inc'
        parameter (one = 1.d0, tol = 1.d-14, mni = 500)
        
c set the initial value of iteration number (i) and
c the increment of void volume fraction(df)
        i = 0
        df = 0.d0
c---------------------------------------------------------------
c        starting iterative process for finding df
c        by newton-raphson method
c---------------------------------------------------------------
c update void volume fraction(f) for current iteration
 1000   f= ft + df
        
c determine the correction values(cf) of void volume fraction
c increment using newton-raphson method
        f1=df-(one-f)*demp
        f2=one+demp
        cf=-f1/f2
c update the value of void volume fraction increment(df) for
c current iteration
        df=df+cf
c update the iteration number
        i = i + 1
c checking whether the number of current iteration exceed
c the maximum allowed number or not
        if (i .lt. mni) goto 1100
           write(102,*)'kbi cannot converge'
        stop
c checking whether the correction value(cf) is greater than the tolerance

          tcf = cf*1000.
 1100     if (abs(tcf) .gt. tol) goto 1000
                  

       end subroutine kfi

       SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
        !
          INCLUDE 'ABA_PARAM.INC'
          DIMENSION TIME(2),VNLEPSEQ_OUT(8)
          CHARACTER*256 JOBNAME,OUTDIR
         
          CALL GETJOBNAME(JOBNAME,LENJOBNAME)
          CALL GETOUTDIR(OUTDIR,LENOUTDIR)
        !
          IF(LOP.EQ.0) THEN
             OPEN(16,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.feh',
     &              STATUS='REPLACE')
             OPEN(17,FILE=TRIM(OUTDIR)//'/'//TRIM(JOBNAME)//'.fep',
     &              STATUS='REPLACE')
          ELSE IF(LOP.EQ.3) THEN
             CLOSE(16)
          END IF
          RETURN
        END SUBROUTINE UEXTERNALDB

       
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
        print*, "Enter load general"
        
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
          print*, "End load 2"
      
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
          print*, "End load 4"
      
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
