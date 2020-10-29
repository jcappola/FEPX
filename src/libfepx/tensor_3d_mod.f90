! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE TENSOR_3D_MOD
!
! This module handles operations on 3D tensors.
!
! Contains subroutines:
! TENSOR3DCOMPOSE: Form matrix from deviatoric, skew, or spherical parts. If no
!   parts are passed, the matix is zeroed.
! TESNRO3DDECOMPOSE:
!
! In both subroutines, we use the basis advocated by Tome' and Kocks, 1985.
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND, IK=>INTEGER_KIND, &
    & LK=>LOGICAL_KIND
USE CONSTANTS_MOD, ONLY: ZERO=>RK_ZERO, ONE=>RK_ONE, TWO=>RK_TWO, &
    & HALF=>RK_ONE_HALF, THIRD=>RK_ONE_THIRD, ROOT_2=>RK_ROOT_2, &
    & ROOT_6=>RK_ROOT_6
!
IMPLICIT NONE
!
! Private
!
PRIVATE 
!
! Decomposition conventions.
!
INTEGER, PARAMETER :: DECOMP_MPSIM=0
INTEGER, PARAMETER :: DECOMP_FEMEVPS=1
INTEGER :: DECOMP_DFLT = DECOMP_MPSIM
!INTEGER :: DECOMP_DFLT = DECOMP_FEMEVPS
!
! Constants.
!
REAL(RK), PARAMETER :: SQ2_I = ONE / ROOT_2
REAL(RK), PARAMETER :: SQ6_I = ONE / ROOT_6
REAL(RK), PARAMETER :: TWOSQ6_I = TWO * ROOT_6
!
! Public
!
PUBLIC :: TENSOR3DDECOMPOSE
PUBLIC :: TENSOR3DCOMPOSE
!
PUBLIC :: DECOMP_MPSIM
PUBLIC :: DECOMP_FEMEVPS
!
CONTAINS 
    !
    SUBROUTINE TENSOR3DCOMPOSE(MAT, DEV, SKW, SPH)
    !
    ! Form  matrix from deviatoric, skew, or spherical parts. If no parts are
    !   passed, the matix is zeroed.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Resulting array of 3x3 matrices
    ! DEV: Arry of 5-vec representing symmetric, traceless portion
    ! SKW: Array of 3-vec representing axial vectors of the skew pportion
    ! SPH: Array of scalars representing one third of the trace
    !
    REAL(RK), INTENT(OUT) :: MAT(:, :, :)
    REAL(RK), INTENT(IN), OPTIONAL :: DEV(:, :)
    REAL(RK), INTENT(IN), OPTIONAL :: SKW(:, :)
    REAL(RK), INTENT(IN), OPTIONAL :: SPH(:)
    !
    !---------------------------------------------------------------------------
    !
    MAT = ZERO
    !
    IF (PRESENT(DEV)) THEN
        !
        MAT(1, 1, :) = -DEV(1, :) * SQ2_I - DEV(2, :) * SQ6_I
        MAT(2, 2, :) = DEV(1, :) * SQ2_I - DEV(2, :) * SQ6_I
        MAT(3, 3, :) = DEV(2, :) * TWOSQ6_I
        !
        MAT(2, 3, :) = DEV(3, :) * SQ2_I
        MAT(3, 2, :) = MAT(2, 3, :)
        !
        MAT(1, 3, :) = DEV(4, :) * SQ2_I
        MAT(3, 1, :) = MAT(1, 3, :)
        !
        MAT(1, 2, :) = DEV(5, :) * SQ2_I
        MAT(2, 1, :) = MAT(1, 2, :)
        !
    END IF
    !
    IF (PRESENT(SKW)) THEN
        !
        MAT(2, 3, :) = MAT(2, 3, :) - SKW(1, :)
        MAT(3, 2, :) = MAT(3, 2, :) + SKW(1, :)
        !
        MAT(1, 3, :) = MAT(1, 3, :) + SKW(2, :)
        MAT(3, 1, :) = MAT(3, 1, :) - SKW(2, :)
        !
        MAT(1, 2, :) = MAT(1, 2, :) - SKW(3, :)
        MAT(2, 1, :) = MAT(2, 1, :) + SKW(3, :)
        !
    END IF
    !
    IF (PRESENT(SPH)) THEN
        !
        MAT(1, 1, :) = MAT(1, 1, :) + SPH
        MAT(2, 2, :) = MAT(2, 2, :) + SPH
        MAT(3, 3, :) = MAT(3, 3, :) + SPH
        !
    END IF
    !
    END SUBROUTINE TENSOR3DCOMPOSE
    !
    !===========================================================================
    !
    SUBROUTINE TENSOR3DDECOMPOSE(MAT, DEV, SKW, SPH, DECOMP)
    !
    ! Decompose matrix into deviatoric, skew, or spherical parts.
    !
    ! Note that the three components of the decomposition are orthogonal, but
    !   the underlying basis is not orthonormal due to scaling.
    !
    ! Note: dev(2) = sqrt(3/2)*mat(3,3), when mat is deviatoric.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MAT: Array of 3x3 matrices
    ! DEV: Array of 5-vec representing symmetric portion (MPSIM convention)
    ! SKW: Array of 3-vec representing skew portion (MPSIM convention)
    ! SPH: Array of scalars representing one third of the trace
    ! DECOMP: Flag indicating decomposition convention
    !
    REAL(RK), INTENT(IN) :: MAT(:, :, :)
    REAL(RK), INTENT(OUT), OPTIONAL :: DEV(:, :)
    REAL(RK), INTENT(OUT), OPTIONAL :: SKW(:, :)
    REAL(RK), INTENT(OUT), OPTIONAL :: SPH(:)
    INTEGER, INTENT(IN), OPTIONAL :: DECOMP
    !
    ! Locals:
    !
    INTEGER ::  DECOMP_CONV
    !
    !---------------------------------------------------------------------------
    !
    DECOMP_CONV = DECOMP_DFLT
    !
    IF (PRESENT(DECOMP)) THEN
        !
        DECOMP_CONV = DECOMP
        !
    END IF
    !
    IF (PRESENT(DEV)) THEN
        !
        SELECT CASE(DECOMP_CONV)
        !
        CASE (DECOMP_MPSIM)
            !
            DEV(1, :) = (MAT(2, 2, :) - MAT(1, 1, :)) * SQ2_I
            DEV(2, :) = (MAT(3, 3, :) + MAT(3, 3, :) - MAT(1, 1, :) - &
                & MAT(2, 2, :)) * SQ6_I
            !
            DEV(3, :) = SQ2_I * (MAT(2, 3, :) + MAT(3, 2, :))
            DEV(4, :) = SQ2_I * (MAT(1, 3, :) + MAT(3, 1, :))
            DEV(5, :) = SQ2_I * (MAT(1, 2, :) + MAT(2, 1, :))
            !
        CASE (DECOMP_FEMEVPS)
            !
            DEV(1,:) = ( MAT(1,1,:) - MAT(2,2,:) )*SQ2_I
            DEV(2,:) = ( MAT(3,3,:) + MAT(3,3,:) - MAT(1,1,:) - MAT(2,2,:) ) * &
                & SQ6_I
            !
            DEV(3, :) = SQ2_I * (MAT(1, 2, :) + MAT(2, 1, :))
            DEV(4, :) = SQ2_I * (MAT(1, 3, :) + MAT(3, 1, :))
            DEV(5, :) = SQ2_I * (MAT(2, 3, :) + MAT(3, 2, :))
            !
        CASE DEFAULT
        !
        ! Return error status
        !
        END SELECT
        !
    END IF
    !
    IF (PRESENT(SKW)) THEN
        !
        SELECT CASE(DECOMP_CONV)
        !
        CASE (DECOMP_MPSIM)
            !
            SKW(1, :) = HALF * (MAT(3, 2, :) - MAT(2, 3, :))
            SKW(2, :) = HALF * (MAT(1, 3, :) - MAT(3, 1, :))
            SKW(3, :) = HALF * (MAT(2, 1, :) - MAT(1, 2, :))
            !
        CASE (DECOMP_FEMEVPS)
            !
            SKW(1, :) = HALF * (MAT(2, 1, :) - MAT(1, 2, :))
            SKW(2, :) = -HALF * (MAT(1, 3, :) - MAT(3, 1, :))
            SKW(3, :) = HALF * (MAT(3, 2, :) - MAT(2, 3, :))
            !
        CASE DEFAULT
            !
            ! Return error status
            !
        END SELECT
        !
    END IF
    !
    IF (PRESENT(SPH)) THEN
        !
        SPH = (MAT(1, 1, :) + MAT(2, 2, :) + MAT(3, 3, :)) * THIRD
        !
    END IF
    !
    END SUBROUTINE TENSOR3DDECOMPOSE
    !
END MODULE TENSOR_3D_MOD
