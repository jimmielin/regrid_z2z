!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: regrid_z2z_mod.F90
!
! !DESCRIPTION: Module REGRID\_Z2Z\_MOD is a quick dependency-free module to
!  vertically interpolate data from and to arbitrary pressure grids.
!  Regridding is supported for concentration / flux quantities.
!
! !INTERFACE:
!
module Regrid_Z2Z_Mod
!
! !USES:
!
    implicit none
    private
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: Map_Z2Z            ! 1-D Column variant
    !public :: Map_Z2Z_Chunk     ! 2-D Chunk variant

    interface Map_Z2Z
        module procedure Map_Z2Z_R8R8
    end interface Map_Z2Z
!
! !MODULE VARIABLES:
!
    ! Built in precision definitions
    integer, parameter :: r8 = kind(real(0.0, 8))
    integer, parameter :: r4 = kind(real(0.0, 4))
!
! !REVISION HISTORY:
!  25 May 2020 - H.P. Lin   - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_Z2Z_R8R8
!
! !DESCRIPTION: Subroutine Map\_Z2Z\_R8R8 is a vertical arbitrary pressure grid to
!  arbitrary pressure grid conservative mapping. Input and output data have double
!  precision. This is a 1-D column operation.
!
!\\
!\\
! !INTERFACE:
!
    subroutine Map_Z2Z_R8R8(IZ, OZ, IPEDGE, OPEDGE, IQ, OQ, missval)
!
! !INPUT PARAMETERS:
!
        ! Input and output vertical dimensions
        integer,  intent(in) :: IZ, OZ

        ! Input and output vertical pressure edges. The bottommost later is 1.
        ! If you want to flip the atmosphere, flip the input data first.
        real(r8), intent(in) :: IPEDGE(IZ+1), OPEDGE(IZ+1)

        ! Data (quantity) on input pressure grid
        real(r8), intent(in) :: IQ(IZ)
!
! !OUTPUT PARAMETERS:
!
        ! Regridded quantity on output pressure grid
        real(r8), intent(out) :: OQ(OZ)
!
! !OPTIONAL ARGUMENTS:
!
        ! Missing fill value
        real(r8), intent(in), optional :: missval
!
! !REMARKS:
!  Performs linear weighted interpolation on the pressure grid.
!
! !REVISION HISTORY:
!  25 May 2020 - H.P. Lin   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer               :: Li, Lo
        integer               :: Lo_below, Lo_above
        real(r8), allocatable :: ZLevRegridWgt(:,:)     ! Regrid weight (OZ, IZ) diml.

        ! Initialize
        if(present(missval)) then
            OQ = missval
        else
            OQ = 0.0_r8
        endif

        ! Allocate regridding weights. The weights can be reused if regridding
        ! in a chunk, which is much more efficient as the handles do not need
        ! to be re-computed at every iteration.
        !              IZ
        !             ----
        !   OQ(Lo) =  \    IQ(Li) * ZLevRegridWgt(Li, Lo)
        !             /
        !             ----
        !            Li = 1
        allocate(ZLevRegridWgt(IZ, OZ)) ! Column major. The read is done more often so optimized for that sum
        ZLevRegridWgt = 0.0_r8

        ! Sift through original input levels for weight generation (computing the
        ! pressure overlap between PEDGEs in input and output)
        !
        ! The weight generation loop is INTENTIONALLY separate from the matrix multiply.
        ! This is to allow reusing of the weights.
        Lo_above = 1
        do Li = 1, IZ
            ! If we have already traversed to the topmost layer, then leave
            if(Lo_above .eq. OZ) then
                exit
            endif

            ! Look if my pressure (IPEDGE(Li)) is above which Lo?
            do Lo = Lo_above, OZ
                if( (IPEDGE(Li) .ge. OPEDGE(Lo) .and. IPEDGE(Li+1) .le. OPEDGE(Lo  ) ) .or. &
                    (IPEDGE(Li) .le. OPEDGE(Lo) .and. IPEDGE(Li)   .ge. OPEDGE(Lo+1) )) then
                    Lo_above = Lo
                    ! There is an overlap. Compute the overlap which is the weight that this input lev
                    ! contributes to the output
                    !                      -- OPEDGE(Lo+1) --     ^ low pressure  (up there)
                    !   -- IPEDGE(Li+1) -- ...
                    !   //////////////////
                    !                      --  OPEDGE(Lo)  --     v high pressure (ground)
                    !   --  IPEDGE(Li)  --
                    ! Full overlap edge case: output layer completely covers input
                    !                      -- OPEDGE(Lo+1) --
                    !   -- IPEDGE(Li+1) --
                    !   --  IPEDGE(Li)  --
                    !                      --  OPEDGE(Lo)  --
                    ! No overlap edge case: opedge is way below or way up
                    !   -- IPEDGE(Li+1) -- -- OPEDGE(Lo+2) --   |                      -- OPEDGE(Lo+1) --
                    !   \\\\\\\\\\\\\\\\\\                      |                      --  OPEDGE(Lo)  --
                    !   --  IPEDGE(Li)  --                      |
                    !                      -- OPEDGE(Lo+1) --   |   -- IPEDGE(Li+1) --
                    !                      --  OPEDGE(Lo)  --   |   --  IPEDGE(Li)  --
                    ZLevRegridWgt(Li, Lo) = (min(OPEDGE(Lo), IPEDGE(Li)) - max(IPEDGE(Li+1), OPEDGE(Lo+1))) / &
                                            (IPEDGE(Li) - IPEDGE(Li+1))

                    write(6, *) "Regird_Z2Z debug: Weighting of Input L = ", Li, " to Output Lo = ", Lo, ZLevRegridWgt(Li, Lo)
                endif
            enddo
        enddo

        ! Apply regridding weights. The matrix multiply could have been more efficient
        do Lo = 1, OZ
            OQ(Lo) = sum(IQ(:) * ZLevRegridWgt(:, Lo))
        enddo

        ! Return
    end subroutine Map_Z2Z_R8R8
end module Regrid_Z2Z_Mod

