! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MaterialMatrixVpModule

  USE parallel_mod
  USE DIMENSIONS_MOD
  USE READ_INPUT_MOD
  use StressSolveVpModule
  USE KinematicsModule

  IMPLICIT  NONE

  PRIVATE
  PUBLIC :: material_matrix_vp, defrate

CONTAINS

      SUBROUTINE material_matrix_vp(&
  &   type, stif, dndx, dndy, dndz, gvel, scale,&
  &   qr5x5, wts, epseff, dtime, incr)

!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER   ::  type, incr
      REAL(RK)  ::  dtime
      REAL(RK)  ::  stif(TVEC, TVEC, el_sub1:el_sup1)
      REAL(RK)  ::  dndx(0:nnpe, el_sub1:el_sup1), dndy(0:nnpe, el_sub1:el_sup1)
      REAL(RK)  ::  dndz(0:nnpe, el_sub1:el_sup1), gvel(0:kdim1, el_sub1:el_sup1)
      REAL(RK)  ::  scale(el_sub1:el_sup1)
      REAL(RK)  ::  qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  ::  wts(0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  ::  epseff(el_sub1:el_sup1)
      ! tm268_M13:
!      REAL(RK)  ::  eqplas(el_sub1:el_sup1)      
!
!     Locals:
!    
      INTEGER i, j, m_el
      REAL(RK)&
  &   d(0:DIMS1, 0:DIMS1, el_sub1:el_sup1),&
  &   d_vec(0:TVEC1, el_sub1:el_sup1)
      REAL(RK) cmu(el_sub1:el_sup1), temp_k(el_sub1:el_sup1)
      REAL(RK) dsdt_iso(el_sub1:el_sup1), state_iso(el_sub1:el_sup1)
      REAL(RK) cconst(TVEC)
!
!     Data:
!
      DATA cconst /1.0d0, 3.0d0, 4.0d0, 4.0d0, 4.0d0/
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
      call defrate(d, dndx, dndy, dndz, gvel)

      ! tm268_M13:
!      call eff_def(epseff, eqplas, d, dtime, m_el)
      call eff_def(epseff, d, dtime, m_el)

! hr-tm
! we aren't going through this portion to update for two-phase compatibility
      if (type .eq. ANISOTROPIC_VP) then
         call par_quit('Error  :     > ANISOTROPIC_VP is no longer implemented.')
      else if (type .eq. ISOTROPIC_VP) then
!
!dbg    Note the hardwired temperature and state variable.
!
        temp_k    = 273.0 + 400.0
        state_iso = 20.0e06

        call isotropic(epseff, temp_k, state_iso, cmu, dsdt_iso)

        stif = 0.0
        stif(1, 1, :) = cconst(1) * cmu
        stif(2, 2, :) = cconst(2) * cmu
        stif(3, 3, :) = cconst(3) * cmu
        stif(4, 4, :) = cconst(4) * cmu
        stif(5, 5, :) = cconst(5) * cmu

      endif

      scale = 0.0
      do i = 1, TVEC
        scale = scale + stif(i, i, :) / cconst(i)
      enddo
      
      scale = scale / 5.0

      END SUBROUTINE material_matrix_vp

SUBROUTINE isotropic(dii, t, s_pa, visc, dsdt)
!
!     Return viscosity for isotropic problem.
!
!     Apparently, this routine was disabled and was replaced by a
!     linear viscous problem with unit viscosity.
!
!----------------------------------------------------------------------
!
USE READ_INPUT_MOD
USE IntrinsicTypesModule, RK=>REAL_KIND
!
IMPLICIT NONE
!
!
!     compute the viscosity for an element based on the effective
!     strain rate.  uses the fit of dawson to hart's model for pure
!     aluminum.
!
!     input arguments -
!
!     dii    - the effective strain rate for the element
!     t      - element temperature (deg f)
!     s_pa   - the value of the state variable (in Pascals)
!
REAL(RK) dii(el_sub1:el_sup1), t(el_sub1:el_sup1)
REAL(RK) s(el_sub1:el_sup1), s_pa(el_sub1:el_sup1)
!
!     output arguments -
!
!     visc   - the resulting viscosity
!     dsdt   - the derivative of the state variable
!
REAL(RK) visc(el_sub1:el_sup1), dsdt(el_sub1:el_sup1)
!
!     modifications -
!
!     the effective viscosity computation has been changed to
!     visc = (2nd invariant of sigma) / 3*(2nd invariant of strain rate)
!     so as to be compatible with definition of strain rate invariant
!     by ajb (9/25/87).
!
!     note - units are mpa, kjoule/mole-k, etc.
!
!     constants -
!
!     REAL(RK) a0
!     PARAMETER ( a0     = 9.64d52     )
!
REAL(RK) log_a0
PARAMETER ( log_a0 = 122.0       )

REAL(RK) c0, qpr, ge, m, f0, qr, smallm, lambda, mp, n

PARAMETER ( c0     = 6.19d-9     )
PARAMETER ( qpr    = 1.45d4      )
PARAMETER ( ge     = 24.20d3     )
PARAMETER ( m      = 7.80        )
PARAMETER ( f0     = 2.12d19     )
PARAMETER ( qr     = qpr         )
PARAMETER ( smallm = 5.0         )
PARAMETER ( lambda = 0.14        )
PARAMETER ( mp     = 3.5         )
PARAMETER ( n      = 6.0         )
!
REAL(RK) a(el_sub1:el_sup1)
REAL(RK) dstar(el_sub1:el_sup1)
REAL(RK) sigma(el_sub1:el_sup1)
REAL(RK) sigmap(el_sub1:el_sup1)
REAL(RK) sigmav(el_sub1:el_sup1)
!
!     Locals:
!
INTEGER i
!
!----------------------------------------------------------------------
!
!
!     Place the state variable in Pascals.
!
visc = 1.0
!
RETURN

s = s_pa / 1.0d6

!     compute flow stress in viscous (friction) element.

!     a = a0 * exp(-qpr/t)

a = log_a0 + (-qpr/t)
a = exp(a)

sigmav = ge * (dii/a)**(1.0/m)
!
!     compute flow stress in the plastic element.
!
dstar  = f0 * ((s/ge)**smallm) * exp(-qr/t)
sigmap = s * exp( - (dstar/dii)**lambda )

!     total flow stress

sigma = sigmav + sigmap

!     compute the derivative of the state vaiable

dsdt = c0*s*dii*((ge/s)**mp)*((sigmap/s)**n)

!     Convert to Pascals

!     sigma = s * dii**0.10
sigma = sigma * 1.0e6

!     compute viscosity

visc = sigma / (3.d0*dii)

print *,'sigma', minval(sigma), maxval(sigma)
print *,'visc', minval(visc), maxval(visc)

RETURN
END
      
END MODULE MaterialMatrixVpModule
