C -----------------------------------------------------------------
C User Subroutine Umat for ABAQUS - Hypoelastic material
C Jaumann rate
C
C By Lizt - School of Aerospace Engeerning, Tsinghua
C 2022-05-07 
C -----------------------------------------------------------------

    Subroutine Umat(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      Implicit None
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

C
      parameter(
     1      zero = 0.0d0
     1       one = 1.0d0
     2       two = 2.0d0
     3      half = 0.5d0  )

      real*8 E, nu, mu, lamda

C Inputs
      E = props(1)      ! Young's modulus
      nu = props(2)     ! Poisson's ratio

C Lame's parameters
      mu = E/two/(one + nu) 
      lamda = mu*E/(one + mu)/(one - two*mu)

C Calculate Jacobi Matrix of Material
      do i = 1,NDI
        do j = 1,NDI
            DDSDDE(i,j) = lamda
        end do
        DDSDDE(i,i) = lamda + two*mu
      end do

      do i = 1,NSHR
        DDSDDE(NDI+i,NDI+i) = mu
      end do

C Update Stress
      STRESS = STRESS + matmul(DDSDDE,DSTRAN)

      RETURN
      END