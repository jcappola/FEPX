! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE PARALLEL_MATRIX_MOD
!
! Parallel matrix operations
!
! Contains subroutines:
! SPARSE_MATVEC_EBE: Matrix times vector, element by element
! GEN_MATRIX_VECTOR_MULT: Matrix times vector
! GEN_MATRIX_MULT: Matrix times matrix
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
USE GATHER_SCATTER
!
IMPLICIT NONE
!
CONTAINS
    !
    SUBROUTINE SPARSE_MATVEC_EBE(RES, SOL, TEMP1, TEMP2, GSTIF, BCS, NNPE, &
        & NSUB, NSUP, ESUB, ESUP, DTRACE, NP)
    !
    !  Matrix vector multiply, element by element.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! RES:
    ! SOL:
    ! TEMP1:
    ! TEMP2:
    ! GSTIF:
    ! BCS:
    ! NNPE:
    ! NSUB:
    ! NSUP:
    ! ESUB:
    ! ESUP:
    ! DTRACE:
    ! NP:
    !
    REAL(RK) :: RES(NSUB:NSUP)
    REAL(RK) :: SOL(NSUB:NSUP)
    REAL(RK) :: TEMP1(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: TEMP2(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: GSTIF(0:(NNPE - 1), 0:(NNPE - 1), ESUB:ESUP)
    LOGICAL :: BCS(NSUB:NSUP)
    INTEGER :: NNPE
    INTEGER :: NSUB
    INTEGER :: NSUP
    INTEGER :: ESUB
    INTEGER :: ESUP
    TYPE(TRACE) :: DTRACE
    INTEGER :: NP(0:(NNPE-1),ESUB:ESUP)
    !
    !  Locals:
    !
    INTEGER :: IER
    INTEGER :: I
    INTEGER :: J
    INTEGER :: IDUMMY
    !
    !---------------------------------------------------------------------------
    !
    CALL PART_GATHER(TEMP1, SOL, NP, DTRACE)
    !
    TEMP2 = 0.0_RK
    !
    CALL GEN_MATRIX_VECTOR_MULT(TEMP2, GSTIF, TEMP1, IDUMMY, IDUMMY, IDUMMY, &
        & IDUMMY, IER)
    !
    RES = 0.0_RK
    !
    CALL PART_SCATTER(RES, TEMP2, NP, .FALSE., DTRACE)
    !
    WHERE (BCS)
        !
        RES = 0.0_RK
        !
    END WHERE
    !
    RETURN
    !
    END SUBROUTINE SPARSE_MATVEC_EBE
    !
    !===========================================================================
    !
    SUBROUTINE GEN_MATRIX_VECTOR_MULT(Y, A, X, I1, I2, I3, I4, IER)
    !
    !  Array of matrices times array of vectors: y(i) = A(i)*x(i)
    !
    ! Note that the arrays x,y,a are all assumed shape arrays and that they are
    !   therefore dimensioned 1:p, 1:q, etc., even though they may be
    !   dimensioned differently in the calling routine.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! Y:
    ! A:
    ! X:
    ! I1:
    ! I2:
    ! I3:
    ! I4:
    ! IER:
    !
    REAL(RK) :: Y(:,:)
    REAL(RK) :: A(:,:,:)
    REAL(RK) :: X(:,:)
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    INTEGER :: I4
    INTEGER :: IER
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: AD1
    INTEGER :: AD2
    INTEGER :: AD3
    !
    !---------------------------------------------------------------------------
    !
    AD1 = UBOUND(A, 1)
    AD2 = UBOUND(A, 2)
    AD3 = UBOUND(A, 3)
    !
    DO K = 1, AD3
        !
        DO J = 1, AD2
            !
!           DO I = 1, AD1
            Y(:, K) = Y(:, K) + A(:, J, K) * X(J, K)
!           END DO
            !
        END DO
        !
    END DO
    !
    IER = 0
    !
    RETURN
    !
    END SUBROUTINE GEN_MATRIX_VECTOR_MULT
    !
    !===========================================================================
    !
    SUBROUTINE GEN_MATRIX_MULT(C, A, B, I1, I2, IER)
    !
    ! Accumulate Array of matrices times array of matrices: c(i) = a(i) * b(i)
    !
    ! Note: c is not zeroed in this subroutine, so the values are accumulated;
    !   in all the calling routines, however, the c array is zeroed before
    !   calling this routine.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! C: Array of result matrices (c=a*b)
    ! A, B: Array of input matrices
    ! I1, I2: Not used (MPK: So why are they here?)
    ! IER: Return status (0 = okay)
    !
    REAL(RK) :: C(:,:,:)
    REAL(RK) :: A(:,:,:)
    REAL(RK) :: B(:,:,:)
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: IER
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: M
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDC
    INTEGER :: LTA
    INTEGER :: LTB
    INTEGER :: N
    !
    !----------------------------------------------------------------------
    !
    LDA = UBOUND(A, 1)
    LDB = UBOUND(B, 1)
    LDC = UBOUND(C, 1)
    LTA = UBOUND(A, 2)
    LTB = UBOUND(B, 2)
    N = UBOUND(A, 3)
    !
    DO M = 1, N
        !
        DO J = 1, LTB
            !
            DO K = 1, LTA
                !
!               DO I = 1,LDC
                C(:, J, M) = C(:, J, M) + A(:, K, M) * B(K, J, M)
!               END DO
                !
            END DO
            !
        END DO
        !
    END DO
    !
    IER = 0
    !
    RETURN
    !
    END SUBROUTINE GEN_MATRIX_MULT
    !
END MODULE PARALLEL_MATRIX_MOD
