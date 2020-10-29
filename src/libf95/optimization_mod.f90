! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE OPTIMIZATION_MOD
!
!  Optimization routines
!
! Contains subroutines:
! MINSEARCHGOLDEN:
! SOLVENEWTON1D:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
USE CONSTANTS_MOD, GCONJ=>RK_GOLD_CONJ
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: MINSEARCHGOLDEN
PUBLIC :: SOLVENEWTON1D
PUBLIC :: NEWTON_RES
PUBLIC :: NEWTON_DER
!
INTEGER, PARAMETER :: NEWTON_RES = 0
INTEGER, PARAMETER :: NEWTON_DER = 1
!
CONTAINS
    !
    SUBROUTINE MINSEARCHGOLDEN(F, D, INTERVAL, FMIN, XMIN, TOLX)
    !
    ! Use golden section search to find a minimum.
    !
    ! In the future, more sophisticated convergence criteria could be specified.
    !   The current implementation only accepts an absolute tolerance on x, but
    !   it would be nice to have relative tolerances both on x and on the value.
    !   For that reason, the tolerance arguemnt is made optional so that other
    !   tolerance settings can be easily incorporated.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! F: The function to be minimized
    ! D: is a vector, data (fixed parameters) to be passed to f
    ! INTERVAL: A 2-vector giving endpoints of interval containing min; this is
    !   updated in this routine
    ! FMIN: A scalar, the minimum value found
    ! XMIN: A scalar, the location of the minimizing point
    ! TOLX: A scalar, the absolute tolerance on x
    !
    INTERFACE
        !
        FUNCTION F(X, D) RESULT(RES)
        !
        USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
        !
        REAL(RK), INTENT(IN) :: X
        REAL(RK), INTENT(IN) :: D(:)
        REAL(RK) :: RES
        !
        END FUNCTION F
        !
    END INTERFACE
    !
    REAL(RK), INTENT(IN) :: D(:)
    REAL(RK), INTENT(IN OUT) :: INTERVAL(2)
    REAL(RK), INTENT(OUT) :: FMIN
    REAL(RK), INTENT(OUT) :: XMIN
    REAL(RK), INTENT(IN) :: TOLX
    OPTIONAL :: TOLX
    !
    ! Locals:
    !
    REAL(RK) :: L_K
    REAL(RK) :: R_K
    REAL(RK) :: AK
    REAL(RK) :: BK
    REAL(RK) :: DK
    REAL(RK) :: FA
    REAL(RK) :: FB
    !
    !---------------------------------------------------------------------------
    !
    L_K = INTERVAL(1)
    R_K = INTERVAL(2)
    DK = R_K - L_K
    AK = R_K - GCONJ * DK
    BK = L_K + GCONJ * DK
    FA = F(AK, D)
    FB = F(BK, D)
    FMIN = MIN(FA, FB)
    !
    ! Loop until tolerance is reached.
    ! Interval is reduced by a constant factor each time.
    !
    DO WHILE (DK > TOLX)
        !
        IF (FA <= FB) THEN
            !
            ! L_K = L_K
            R_K = BK
            DK = R_K - L_K
            AK = R_K - GCONJ * DK
            BK = AK
            FB = FA
            FA = F(AK, D)
            !
        ELSE
            !
            L_K = AK
            ! R_K = R_K
            DK = R_K - L_K
            AK = BK
            BK = L_K + GCONJ * DK
            FA = FB
            FB = F(BK, D)
            !
        END IF
        !
    END DO
    !
    ! Update output variables on return.
    !
    FMIN = MIN(FA, FB)
    IF (FA <= FB) THEN
        !
        INTERVAL(1) = L_K
        XMIN = AK
        INTERVAL(2) = BK
        !
    ELSE
        !
        INTERVAL(1) = AK
        XMIN = BK
        INTERVAL(2) = R_K
        !
    END IF
    !
    END SUBROUTINE MINSEARCHGOLDEN
    !
    !===========================================================================
    !
    SUBROUTINE SOLVENEWTON1D(F, D, X, TOLX, STATUS)
    !
    ! Newton solver for nonlinear equations
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! F: is a function the RESidual to make zero; it also returns the derivative
    !   of the RESidual on request
    ! D: An array, fixed parameters to pass to f and FPRIME
    ! X: A scalar on input, the initial guess; on ouput, the solution, if
    !   converged
    ! TOLX: A scalar, the solution tolerance on x
    ! STATUS: Flag indicating success(==0)|failure(/=0)
    !
    INTERFACE
        !
        FUNCTION F(X, D, DER) RESULT(RES)
        !
        ! X: The point to evaluate
        ! D: Array of fixed data
        ! DER: Flag indicating return value of RESidual or derivative
        !
        USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
        !
        REAL(RK), INTENT(IN) :: X
        REAL(RK), INTENT(IN) :: D(:)
        INTEGER :: DER
        REAL(RK) :: RES
        !
        END FUNCTION F
        !
    END INTERFACE
    !
    REAL(RK), INTENT(IN) :: D(:)
    REAL(RK), INTENT(IN OUT) :: X
    REAL(RK), INTENT(IN) :: TOLX
    INTEGER, INTENT(OUT) :: STATUS
    !
    ! Locals:
    ! (Hardwire max iterations to prevent infinite loops)
    !
    INTEGER :: ITERMAX = 1000
    INTEGER :: N_INCREASED = 3
    INTEGER :: I
    REAL(RK) :: RES
    REAL(RK) :: FPRIME
    REAL(RK) :: DX
    REAL(RK) :: DXOLD
    REAL(RK) :: N_INC
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 1  ! so that failure is set iterations reaching itermax
    !
    N_INC = 0
    DXOLD = HUGE(DX)
    !
    DO I = 1, ITERMAX
        !
        RES = F(X, D, NEWTON_RES)
        FPRIME = F(X, D, NEWTON_DER)
        !
        IF (FPRIME == RK_ZERO) THEN
            !
            ! Zero derivative
            !
            STATUS = 2
            !
            EXIT
            !
        END IF
        !
        DX = RES/FPRIME
        X = X - DX
        DX = ABS(DX)
        !
        IF (DX < TOLX) THEN
            !
            ! Converged
            !
            STATUS = 0
            !
            EXIT
            !
        END IF
        !
        ! Check for divergence
        !
        IF (DX < DXOLD) THEN
            !
            N_INC = 0
            !
        ELSE
            !
            N_INC = N_INC + 1
            !
            IF (N_INC >= N_INCREASED) THEN
                !
                ! Diverged
                !
                STATUS = 1
                !
                EXIT
                !
            END IF
            !
        END IF
        !
        DXOLD = DX
        !
    END DO
    !
    END SUBROUTINE SOLVENEWTON1D
    !
END MODULE OPTIMIZATION_MOD
