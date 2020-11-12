! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE KINEMATICS_MOD
!
! Kinematic calculations
!
! Contains subroutines:
! VEL_GRADIENT: Compute velocity gradient
! EFF_DEF: Compute effective deformation rate
! DEFRATE: Compute deformation rate tensor
!
! From libf95:
!
USE LIBF95, RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE READ_INPUT_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: VEL_GRADIENT
PUBLIC :: EFF_DEF
PUBLIC :: DEFRATE
!
CONTAINS
    !
    SUBROUTINE VEL_GRADIENT(VGRAD, DNDX, DNDY, DNDZ, GVEL)
    !
    ! Compute velocity gradient.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VGRAD:
    ! DNDX:
    ! DNDY:
    ! DNDZ:
    ! GVEL:
    !
    REAL(RK), INTENT(OUT) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DNDX(0:NNPE,  EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DNDY(0:NNPE,  EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DNDZ(0:NNPE,  EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !    
    INTEGER :: I
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    !
    !---------------------------------------------------------------------------
    !
    VGRAD = 0.0D0
    !
    DO I = 0, NNPE
        !
        I1 = 3 * I
        I2 = I1 + 1
        I3 = I2 + 1
        !
        VGRAD(0, 0, :) = VGRAD(0, 0, :) + DNDX(I, :) * GVEL(I1, :)
        VGRAD(0, 1, :) = VGRAD(0, 1, :) + DNDY(I, :) * GVEL(I1, :)
        VGRAD(0, 2, :) = VGRAD(0, 2, :) + DNDZ(I, :) * GVEL(I1, :)
        VGRAD(1, 0, :) = VGRAD(1, 0, :) + DNDX(I, :) * GVEL(I2, :)
        VGRAD(1, 1, :) = VGRAD(1, 1, :) + DNDY(I, :) * GVEL(I2, :)
        VGRAD(1, 2, :) = VGRAD(1, 2, :) + DNDZ(I, :) * GVEL(I2, :)
        VGRAD(2, 0, :) = VGRAD(2, 0, :) + DNDX(I, :) * GVEL(I3, :)
        VGRAD(2, 1, :) = VGRAD(2, 1, :) + DNDY(I, :) * GVEL(I3, :)
        VGRAD(2, 2, :) = VGRAD(2, 2, :) + DNDZ(I, :) * GVEL(I3, :)
        !
    END DO
    !
    END SUBROUTINE VEL_GRADIENT
    !
    !===========================================================================
    !
    SUBROUTINE EFF_DEF(EPSEFF, D, DTIME, M)
    !
    ! Compute the effective deformation rate and accumulated deformation
    !   (strain) from the deformation rate tensor.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! EPSEFF:
    ! D:
    ! DTIME:
    ! M:
    REAL(RK), INTENT(OUT) :: EPSEFF(0:(M - 1))
    REAL(RK), INTENT(IN) :: D(0:DIMS1, 0:DIMS1, 0:(M - 1))
    REAL(RK), INTENT(IN) :: DTIME
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    REAL(RK), PARAMETER :: TWOTHIRDS = 2.0D0 / 3.0D0
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    ! Effective deformation rate
    DO I = 0, M - 1
        !
        EPSEFF(I) = DSQRT(TWOTHIRDS * (D(0, 0, I) * D(0, 0, I) + D(1, 1, I) * &
            & D(1, 1, I) + D(2, 2, I) * D(2, 2, I) + 2.0 * ( D(0, 1, I) * &
            & D(0, 1, I) + D(0, 2, I) * D(0, 2, I) + D(1, 2, I) * D(1, 2, I))))
        !
    END DO
    !
    END SUBROUTINE EFF_DEF
    !
    !===========================================================================
    !
    SUBROUTINE DEFRATE(D, DNDX, DNDY, DNDZ, GVEL)
    !
    ! Calculate deformation rate tensor
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D:
    ! DNDX:
    ! DNDY:
    ! DNDZ:
    ! GVEL:
    !
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    REAL(RK) :: DIVV(EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    D = 0.0D0
    !
    DO I = 0, NNPE
        !
        I1 = 3 * I
        I2 = I1 + 1
        I3 = I2 + 1
        !
        D(0, 0, :) = D(0, 0, :) + DNDX(I, :) * GVEL(I1, :)
        D(1, 1, :) = D(1, 1, :) + DNDY(I, :) * GVEL(I2, :)
        D(2, 2, :) = D(2, 2, :) + DNDZ(I, :) * GVEL(I3, :)
        D(1, 0, :) = D(1, 0, :) + DNDX(I, :) * GVEL(I2, :) &
            & + DNDY(I, :) * GVEL(I1, :)
        D(2, 0, :) = D(2, 0, :) + DNDX(I, :) * GVEL(I3, :)&
            & + DNDZ(I, :) * GVEL(I1, :)
        D(2, 1, :) = D(2, 1, :) + DNDY(I, :) * GVEL(I3, :)&
            & + DNDZ(I, :) * GVEL(I2, :)
        !
    END DO
    !
    D(1, 0, :) = 0.5D0 * D(1, 0, :)
    D(2, 0, :) = 0.5D0 * D(2, 0, :)
    D(2, 1, :) = 0.5D0 * D(2, 1, :)
    !
    D(0, 1, :) = D(1, 0, :)
    D(0, 2, :) = D(2, 0, :)
    D(1, 2, :) = D(2, 1, :)
    !
    DIVV = D(0, 0, :) + D(1, 1, :) + D(2, 2, :)
    DIVV = DIVV / 3.0D0
    !
    D(0, 0, :) = D(0, 0, :) - DIVV
    D(1, 1, :) = D(1, 1, :) - DIVV
    D(2, 2, :) = D(2, 2, :) - DIVV
    !
    END SUBROUTINE DEFRATE
    !
END MODULE KINEMATICS_MOD
