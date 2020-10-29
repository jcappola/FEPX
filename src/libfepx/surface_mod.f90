! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE SURFACE_MOD
!
! Module for using surface information.
!
! Contains function:
! ALLOCATE_SURFACE_SECTION: Allocate space to enough elements
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY:RK=>REAL_KIND
!
! From libparallel
!
USE GATHER_SCATTER
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: NDIM3 = 3
!
! The type `surface_section' is for use in parallel codes, and contains a
!   section of the surface mesh.
!
TYPE SURFACE_SECTION
    !
    ! TYPE: Currently a number indicating nodes per element
    ! NSEL: Number of surface elements in this section
    ! SEMIN to SEMAX: Range of surface element numbers
    ! CONN: Connectivity
    ! CONN3D: DOF connectivity
    ! ECONN: Elemental connectivity to be used for gather/scatter stress
    ! CRDS: Coordinate array
    ! ELEM:
    ! TR: The trace data structure for gather/scatter operations
    ! ETR: Elemental trace for stress
    !
    INTEGER :: TYPE, SEMIN, SEMAX
    INTEGER, POINTER, DIMENSION(:,:) :: CONN
    INTEGER, POINTER, DIMENSION(:,:) :: CONN3D
    INTEGER, POINTER, DIMENSION(:,:) :: ECONN
    REAL(RK), POINTER, DIMENSION(:,:) :: CRDS
    INTEGER, POINTER, DIMENSION(:) :: ELEM
    TYPE(TRACE) :: TR
    TYPE(TRACE) :: ETR
    !
END TYPE SURFACE_SECTION
!
CONTAINS
    !
    FUNCTION ALLOCATE_SURFACE_SECTION(TYPE, SEMIN, SEMAX, SURF) RESULT(STATUS)
    !
    ! Allocate space to enough elements.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TYPE: Surface element type, now only q6 is available
    ! SEMIN, SEMAX: Range of surface element numbers
    ! SURF: The surface to be allocated
    !
    INTEGER :: TYPE
    INTEGER :: SEMIN
    INTEGER :: SEMAX
    TYPE(SURFACE_SECTION) :: SURF
    INTEGER :: STATUS
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 0
    !
    !if (type /= 4) then
    ! tsh for 3-node triangle
    !if (type /= 3) then
    ! tsh for 6-node triangle
    IF (TYPE /= 6) THEN
        !
        STATUS = 1
        !
        RETURN
    END IF
    SURF%TYPE = TYPE ! TYPE=6
    !
    IF (SEMAX >= SEMIN) THEN
        !
        SURF%SEMIN = SEMIN
        SURF%SEMAX = SEMAX
        !
    ELSE
        !
        STATUS = 2
        !
        RETURN
        !
    END IF
    !
    ALLOCATE(SURF%CONN(0:(TYPE-1), SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%CONN3D(0:(3*TYPE-1), SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%ECONN(0:5, SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%CRDS(0:(NDIM3*TYPE-1), SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%ELEM(SEMIN:SEMAX), STAT = STATUS)
    !
    RETURN
    !
    END FUNCTION ALLOCATE_SURFACE_SECTION
    !
END MODULE SURFACE_MOD
