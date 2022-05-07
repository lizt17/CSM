C -----------------------------------------------------------------
C User Subroutine Vumat for ABAQUS - Hypoelastic material
C Green - Naghdi rate
C
C By Lizt - School of Aerospace Engeerning, Tsinghua
C 2022-05-07 
C -----------------------------------------------------------------

      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
      implicit none
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
      real*8 E, nu, mu, lamda
      dimension C(ndir+nshr,ndir+nshr)
C
      parameter(
     1      zero = 0.0d0
     1       one = 1.0d0
     2       two = 2.0d0
     3      half = 0.5d0  )

C Inputs
      E = props(1)      ! Young's modulus
      nu = props(2)     ! Poisson's ratio

C Lame's parameters
      mu = E/two/(one + nu) 
      lamda = mu*E/(one + mu)/(one - two*mu)

C Calculate Jacobi Matrix of Material
      C = zero
      do i = 1,ndir
            do j = 1,ndir
                  C(i,j) = lamda
            end do
            C(i,i) = lamda + two*mu
      end do
      do i = 1,nshr
            C(ndir+i,ndir+i) = two*mu
      end do

C Update Stress
      do 100 km = 1,nblock
      
            stressNew(km,:) = stressOld(km,:) 
     *             + matmul(C,strainInc(km,:))
      
  100 continue

      return
      end