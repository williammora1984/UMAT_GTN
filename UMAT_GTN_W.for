c======================================================================
c    UMAT_for the Gurson, Tvergaard and Needleman (GTN) damage model
c             coded by  Wiliam Fernando Mora Pulido   
c=======================================================================
c  This is a user defined subroutine UMAT for Abaqus, is and inter-
c face that calls the material subroutine kGTN based on the on the 
c algorithm proposed by N.Aravas and Z.L.Zhang
c     
c-----------------------------------------------------------------------
c Version 0.9.1
c coded by: W. Mora Oct 2021
c-----------------------------------------------------------------------

c      NDI         : number of direct components of DDSDDE, DDSDDT, and DRPLDE
c      NSHR        : number of engineering shear components of DDSDDE, DDSDDT, and DRPLDE
c      NTENS       : NDI + NSHR: Size of the stress or strain component array
c                      For solids   NTENS = 6
c                      For plane stress/strain: NTENS = 5
c      DDSDDE      : Algorithmic tangent stiffness (Jacobian) 
c      EET         : Elastic strain tensor
c      EPLAS       : Plastic strain tensor
c      STRAN       : strain from previous increment
c      stress      : stress  
c      DSTRAN      : increment of strain for current iteration
c      DDSDDE      : stiffness matrix
c
c      PROPS(NPROPS): Array with material property data
c      NPROPS       : number of material properties. In CAE > Property > User Material
c                      NPROPS = 13          

c======================================================================
c                    Standard declaration UMAT
c=======================================================================

      SUBROUTINE umat_GTN_W(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

       INCLUDE 'ABA_PARAM.INC'

c======================================================================
c                    Standard variables UMAT
c=======================================================================
       implicit none
       character*80 cmname
       integer:: ndi,nshr,ntens,nstatv,nprops,noel,npt,
     &  layer, kspt, kstep, kinc
       double precision:: dtime,temp,dtemp,pnewdt,celent,sse,spd,scd
     &  ,rpl,drpldt
     &  double precision:: stress(ntens),statev(nstatv),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     &  props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
       !
c======================================================================
c              Additional variables used in this UMAT file
c=======================================================================       

      double precision ::sdv(21)
      integer :: ttype 
      double precision ::sdvup(21)
       
c      sdv   : internal state variables input to material routine
c      ttype : type of tangent stiffnes 
c              0: analytical 
c              1: numerical
c      sdvup : updated internal state varialbes in the material routine
       ntens  = 6 
       nprops = 13
       nstatv = 3
c MATERIAL PARAMETERS

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
c       props(12)  ! :  0.2025              !f_f
c       props(13)  ! : 0.1                  !E_n

C INTERNAL STATE VARIABLES
c       statev(1:6)    : plastic strain
c       statev(7:12)   : Back stress
c       statev(13)     : Drag stress, alpha,
c       statev(14:19)  : Elastic strain
c       statev(20)     : void fraction 
c       statev(21)     : Microstrain eps_b  

c Get the updated rotated strain
       CALL ROTSIG(STATEV(1+NTENS),DROT,EPLAS,2,NDI,NSHR)
       CALL ROTSIG(STATEV(1), DROT,EET, 2,NDI,NSHR)
       
       sdv        = statev(1:21)
       sdv(1:6)   = EPLAS
       sdv(14:19) = EET
       ttype=0

        call kGTN (STRAN+DSTRAN,DSTRAN,sdv,ttype, props, stress, ddsdde, 
     &         sdvup)

       sdvl=sdvup
       END SUBROUTINE umat_GTN_W
       
   