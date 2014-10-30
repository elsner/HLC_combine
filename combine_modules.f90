!
!    Copyright (C) 2013, 2014 Franz Elsner <f.elsner@ucl.ac.uk>
!
!    This file is part of HLC_combine.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program. If not, see <http://www.gnu.org/licenses/>.
!
MODULE COMBINE_MODULES

  USE MPI
  USE HEALPIX_TYPES

  IMPLICIT NONE


CONTAINS

!-----------------------------------------------------------------------
! Allocate memory
!-----------------------------------------------------------------------

  SUBROUTINE ALLOCATE_ARRAYS(noise_lfi__pix_IQU_channel,&
       & noise_hfi__pix_IQU_channel, hlc_map__pix_IQU,&
       & hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
       & w8rings_hfi__ring_IQU, filenames__channel_nr, zbounds__lat,&
       & nr_channels, nside, npix, lmax, mmax, nside_lfi, nside_hfi)

    USE PIX_TOOLS,     ONLY: NSIDE2NPIX

    INTEGER(I4B),                           INTENT(IN)    :: nside_lfi
    INTEGER(I4B),                           INTENT(IN)    :: nside_hfi
    CHARACTER(LEN=*),  DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(OUT)   :: filenames__channel_nr
    INTEGER(I4B),      DIMENSION(1:),       INTENT(OUT)   :: nr_channels
    INTEGER(I4B),      DIMENSION(1:),       INTENT(OUT)   :: nside
    INTEGER(I4B),      DIMENSION(1:),       INTENT(OUT)   :: npix
    INTEGER(I4B),      DIMENSION(1:),       INTENT(OUT)   :: lmax
    INTEGER(I4B),      DIMENSION(1:),       INTENT(OUT)   :: mmax
    REAL(DP),          DIMENSION(1:),       INTENT(OUT)   :: zbounds__lat
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(OUT)   :: w8rings_lfi__ring_IQU
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(OUT)   :: w8rings_hfi__ring_IQU
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(OUT)   :: hlc_map__pix_IQU
    REAL(DP),          DIMENSION(:,:,:),    ALLOCATABLE,&
         &                                  INTENT(OUT)   :: noise_lfi__pix_IQU_channel
    REAL(DP),          DIMENSION(:,:,:),    ALLOCATABLE,&
         &                                  INTENT(OUT)   :: noise_hfi__pix_IQU_channel
    COMPLEX(DPC),      DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(OUT)   :: hlc_weights__l_channel

    INTEGER(I4B)                                          :: max_lmax
    INTEGER(I4B)                                          :: npix_lfi
    INTEGER(I4B)                                          :: npix_hfi

    max_lmax = 3*nside_hfi - 1
    npix_lfi = NSIDE2NPIX(nside_lfi)
    npix_hfi = NSIDE2NPIX(nside_hfi)

    nr_channels(1) = 3
    nr_channels(2) = 6
    nr_channels(3) = 9
    nside(1:3)     = nside_lfi
    nside(4:9)     = nside_hfi
    npix(1:3)      = npix_lfi
    npix(4:9)      = npix_hfi
    lmax           = max_lmax
    mmax           = max_lmax

    ALLOCATE(filenames__channel_nr(nr_channels(3)+1,10000))
    ALLOCATE(w8rings_lfi__ring_IQU(1:2*nside_lfi,1:1))
    ALLOCATE(w8rings_hfi__ring_IQU(1:2*nside_hfi,1:1))
    ALLOCATE(hlc_map__pix_IQU(0:npix_hfi-1,1:1))
    ALLOCATE(noise_lfi__pix_IQU_channel(0:npix_lfi-1,1:1,1:nr_channels(1)))
    ALLOCATE(noise_hfi__pix_IQU_channel(0:npix_hfi-1,1:1,1:nr_channels(2)))
    ALLOCATE(hlc_weights__l_channel(0:max_lmax,1:nr_channels(3)))
    filenames__channel_nr      = ''
    zbounds__lat               = (/ -1.0_DP, 1.0_DP /)
    w8rings_lfi__ring_IQU      = 0.0_DP
    w8rings_hfi__ring_IQU      = 0.0_DP
    hlc_map__pix_IQU           = 0.0_DP
    noise_lfi__pix_IQU_channel = 0.0_DP
    noise_hfi__pix_IQU_channel = 0.0_DP
    hlc_weights__l_channel     = 0.0_DPC


  END SUBROUTINE ALLOCATE_ARRAYS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read auxiliary files
!-----------------------------------------------------------------------

  SUBROUTINE READ_AUX_FILES(hlc_weights__l_channel,&
       & w8rings_lfi__ring_IQU, w8rings_hfi__ring_IQU,&
       & filenames__channel_nr, nr_maps, nside,&
       & lmax, mmax, nr_channels, thread_id, dir_w8ring,&
       & file_hlc_coeff, file_filenamelist)

    USE FITSTOOLS,     ONLY: READ_DBINTAB

    CHARACTER(LEN=*),                       INTENT(IN)    :: dir_w8ring
    CHARACTER(LEN=*),                       INTENT(IN)    :: file_hlc_coeff
    CHARACTER(LEN=*),                       INTENT(IN)    :: file_filenamelist
    INTEGER(I4B),                           INTENT(IN)    :: thread_id
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nr_channels
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nside
    CHARACTER(LEN=*),  DIMENSION(1:,1:),    INTENT(INOUT) :: filenames__channel_nr
    INTEGER(I4B),      DIMENSION(1:),       INTENT(INOUT) :: lmax
    INTEGER(I4B),      DIMENSION(1:),       INTENT(INOUT) :: mmax
    INTEGER(I4B),                           INTENT(OUT)   :: nr_maps
    REAL(DP),          DIMENSION(1:,1:),    INTENT(OUT)   :: w8rings_lfi__ring_IQU
    REAL(DP),          DIMENSION(1:,1:),    INTENT(OUT)   :: w8rings_hfi__ring_IQU
    COMPLEX(DPC),      DIMENSION(0:,1:),    INTENT(OUT)   :: hlc_weights__l_channel

    CHARACTER(LEN=150)                                    :: filename
    CHARACTER(LEN=5)                                      :: nside_string
    LOGICAL(LGT)                                          :: anynull
    INTEGER(I4B)                                          :: ierror
    REAL(DP)                                              :: nullval

    IF (thread_id .EQ. 0) THEN

       WRITE(nside_string,'(I5.5)') nside(1)
       filename = TRIM(dir_w8ring) // '/weight_ring_n'&
            & // nside_string // '.fits'
       CALL READ_DBINTAB(filename, w8rings_lfi__ring_IQU, 2*nside(1),&
            & 1, nullval, anynull)
       w8rings_lfi__ring_IQU = 1.0_DP + w8rings_lfi__ring_IQU

       WRITE(nside_string,'(I5.5)') nside(4)
       filename = TRIM(dir_w8ring) // '/weight_ring_n'&
            & // nside_string // '.fits'
       CALL READ_DBINTAB(filename, w8rings_hfi__ring_IQU, 2*nside(4),&
            & 1, nullval, anynull)
       w8rings_hfi__ring_IQU = 1.0_DP + w8rings_hfi__ring_IQU

       CALL INPUT_HLC_WEIGHTS(hlc_weights__l_channel,&
            & lmax, mmax, nr_channels, file_hlc_coeff)

    ENDIF

    CALL INPUT_FILENAMES(filenames__channel_nr, nr_maps,&
         & file_filenamelist, nr_channels(3))

    CALL MPI_BCAST(lmax, SIZE(lmax), MPI_INTEGER, 0, MPI_COMM_WORLD,&
         & ierror)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", ierror)

    CALL MPI_BCAST(mmax, SIZE(mmax), MPI_INTEGER, 0, MPI_COMM_WORLD,&
         & ierror)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", ierror)

    CALL MPI_BCAST(w8rings_lfi__ring_IQU, SIZE(w8rings_lfi__ring_IQU),&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", ierror)

    CALL MPI_BCAST(w8rings_hfi__ring_IQU, SIZE(w8rings_hfi__ring_IQU),&
         & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", ierror)

    CALL MPI_BCAST(hlc_weights__l_channel, SIZE(hlc_weights__l_channel),&
         & MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
    CALL CHECK_ERROR_FLAG("MPI_BCAST", ierror)


  END SUBROUTINE READ_AUX_FILES

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read input/output filenames
!-----------------------------------------------------------------------

  SUBROUTINE INPUT_FILENAMES(filenames__channel_nr, nr_maps,&
       & file_filenamelist, nr_channels)

    CHARACTER(LEN=*),                       INTENT(IN)    :: file_filenamelist
    INTEGER(I4B),                           INTENT(IN)    :: nr_channels
    CHARACTER(LEN=*),  DIMENSION(1:,1:),    INTENT(OUT)   :: filenames__channel_nr
    INTEGER(I4B),                           INTENT(OUT)   :: nr_maps

    INTEGER(I4B)                                          :: i, j
    INTEGER(I4B)                                          :: counter
    CHARACTER(LEN=150)                                    :: name

    OPEN(UNIT=10, FILE=file_filenamelist, STATUS='OLD',&
         & ACTION='READ', FORM='FORMATTED')
    REWIND 10

    counter = 0
    DO i = 1,UBOUND(filenames__channel_nr, 2)
       DO j = 1,nr_channels+1
          READ(10,'(A)',END=1000) name
          filenames__channel_nr(j,i) = name
          counter = counter + 1
       ENDDO
    ENDDO

    1000 CLOSE(10)

    CALL CHECK_ERROR_FLAG("INPUT_FILENAMES, file format mismatch",&
         & MOD(counter, nr_channels + 1))

    nr_maps = counter/(nr_channels + 1)

    WRITE(*,'(3X, I5, " simulations will be processed")') nr_maps


  END SUBROUTINE INPUT_FILENAMES

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read harmonic weights from file
!-----------------------------------------------------------------------

  SUBROUTINE INPUT_HLC_WEIGHTS(hlc_weights__l_channel, lmax, mmax,&
       & nr_channels, file_weights)

    CHARACTER(LEN=*),                       INTENT(IN)    :: file_weights
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nr_channels
    INTEGER(I4B),      DIMENSION(1:),       INTENT(INOUT) :: lmax
    INTEGER(I4B),      DIMENSION(1:),       INTENT(INOUT) :: mmax
    COMPLEX(DPC),      DIMENSION(0:,1:),    INTENT(OUT)   :: hlc_weights__l_channel

    INTEGER(I4B)                                          :: i, j
    REAL(DP)                                              :: c1_read
    REAL(DP)                                              :: c2_read
    REAL(DP)                                              :: c3_read
    REAL(DP)                                              :: c4_read
    REAL(DP)                                              :: c5_read
    REAL(DP)                                              :: c6_read
    REAL(DP)                                              :: c7_read
    REAL(DP)                                              :: c8_read
    REAL(DP)                                              :: c9_read

    OPEN(UNIT=10, FILE=file_weights, STATUS='OLD',&
         & ACTION='READ', FORM='FORMATTED')
    REWIND 10

    DO i = 0,MAXVAL(lmax)
       READ(10,*,END=1000) c1_read, c2_read, c3_read, c4_read, c5_read,&
            & c6_read, c7_read, c8_read, c9_read
       hlc_weights__l_channel(i,:) = DCMPLX((/ c1_read, c2_read,&
            & c3_read, c4_read, c5_read, c6_read, c7_read, c8_read,&
            & c9_read /), 0.0_DP)
    ENDDO

    1000 CLOSE(10)

    DO i = 1,nr_channels(3)
       lmax(i) = 0
       mmax(i) = 0
       DO j = UBOUND(hlc_weights__l_channel(:,i), 1)-1,&
            & LBOUND(hlc_weights__l_channel(:,i), 1),-1
          IF (hlc_weights__l_channel(j,i) .NE. 0.0_DPC) THEN
             lmax(i) = j
             mmax(i) = j
             EXIT
          ENDIF
       ENDDO
    ENDDO

    lmax(1:3) = MAXVAL(lmax(1:3))
    mmax(1:3) = MAXVAL(mmax(1:3))
    lmax(4:9) = MAXVAL(lmax(4:9))
    mmax(4:9) = MAXVAL(mmax(4:9))

    hlc_weights__l_channel(0:1,:) = 0.0_DPC


  END SUBROUTINE INPUT_HLC_WEIGHTS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Load noise realizations from file
!-----------------------------------------------------------------------

  SUBROUTINE LOAD_NOISE_MAPS(noise_lfi__pix_IQU_channel,&
       & noise_hfi__pix_IQU_channel, filenames__channel_nr,&
       & nr_channels, npix, map_nr)

    USE PIX_TOOLS,     ONLY: CONVERT_INPLACE, NEST2RING

    CHARACTER(LEN=*),  DIMENSION(1:,1:),    INTENT(IN)    :: filenames__channel_nr
    INTEGER(I4B),                           INTENT(IN)    :: map_nr
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: npix
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nr_channels
    REAL(DP),          DIMENSION(0:,1:,1:), TARGET,&
         &                                  INTENT(IN)    :: noise_lfi__pix_IQU_channel
    REAL(DP),          DIMENSION(0:,1:,1:), TARGET,&
         &                                  INTENT(IN)    :: noise_hfi__pix_IQU_channel

    INTEGER(I4B)                                          :: i
    REAL(DP),          DIMENSION(:,:),      POINTER       :: noise__pix_IQU

    DO i = 1,nr_channels(3)

       NULLIFY(noise__pix_IQU)

       IF (i .LE. nr_channels(1)) THEN
          noise__pix_IQU&
               & => noise_lfi__pix_IQU_channel(:,:,i)
       ELSE
          noise__pix_IQU&
               & => noise_hfi__pix_IQU_channel(:,:,i-nr_channels(1))
       ENDIF

       CALL READ_MAP(noise__pix_IQU, filenames__channel_nr(i,map_nr),&
            & npix(i), 1)

    ENDDO

    NULLIFY(noise__pix_IQU)


  END SUBROUTINE LOAD_NOISE_MAPS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Combine to harmonic ILC map
!-----------------------------------------------------------------------

  SUBROUTINE COMBINE_MAPS(hlc_map__pix_IQU,&
       & noise_lfi__pix_IQU_channel, noise_hfi__pix_IQU_channel,&
       & hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
       & w8rings_hfi__ring_IQU, zbounds__lat, nr_channels, nside,&
       & lmax, mmax)

    USE ALM_TOOLS,     ONLY: MAP2ALM, ALM2MAP

    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nr_channels
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nside
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: lmax
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: mmax
    REAL(DP),          DIMENSION(1:)   ,    INTENT(IN)    :: zbounds__lat
    REAL(DP),          DIMENSION(1:,1:),    TARGET,&
         &                                  INTENT(IN)    :: w8rings_lfi__ring_IQU
    REAL(DP),          DIMENSION(1:,1:),    TARGET,&
         &                                  INTENT(IN)    :: w8rings_hfi__ring_IQU
    REAL(DP),          DIMENSION(0:,1:,1:), TARGET,&
         &                                  INTENT(IN)    :: noise_lfi__pix_IQU_channel
    REAL(DP),          DIMENSION(0:,1:,1:), TARGET,&
         &                                  INTENT(IN)    :: noise_hfi__pix_IQU_channel
    COMPLEX(DPC),      DIMENSION(0:,1:),    INTENT(IN)    :: hlc_weights__l_channel
    REAL(DP),          DIMENSION(0:,1:),    INTENT(OUT)   :: hlc_map__pix_IQU

    INTEGER(I4B)                                          :: i, j
    REAL(DP),          DIMENSION(:,:),      POINTER       :: w8rings__ring_IQU
    REAL(DP),          DIMENSION(:,:),      POINTER       :: noise__pix_IQU
    COMPLEX(DPC),      DIMENSION(:,:,:),    POINTER       :: alm__TEB_l_m
    COMPLEX(DPC),      DIMENSION(:,:,:),    TARGET,&
         &                                  ALLOCATABLE   :: alm_lfi__TEB_l_m
    COMPLEX(DPC),      DIMENSION(:,:,:),    TARGET,&
         &                                  ALLOCATABLE   :: alm_hfi__TEB_l_m
    COMPLEX(DPC),      DIMENSION(:,:,:),    ALLOCATABLE   :: hlc_alm__TEB_l_m

    ALLOCATE(alm_lfi__TEB_l_m(1:1,0:lmax(1),0:mmax(1)))
    ALLOCATE(alm_hfi__TEB_l_m(1:1,0:lmax(4),0:mmax(4)))
    ALLOCATE(hlc_alm__TEB_l_m(1:1,0:lmax(4),0:mmax(4)))
    hlc_map__pix_IQU = 0.0_DP
    alm_lfi__TEB_l_m = 0.0_DPC
    alm_hfi__TEB_l_m = 0.0_DPC
    hlc_alm__TEB_l_m = 0.0_DPC

    DO i = 1,nr_channels(3)

       NULLIFY(w8rings__ring_IQU)
       NULLIFY(noise__pix_IQU)
       NULLIFY(alm__TEB_l_m)

       IF (i .LE. nr_channels(1)) THEN
          w8rings__ring_IQU&
               & => w8rings_lfi__ring_IQU
          noise__pix_IQU&
               & => noise_lfi__pix_IQU_channel(:,:,i)
          alm__TEB_l_m&
               & => alm_lfi__TEB_l_m
       ELSE
          w8rings__ring_IQU&
               & => w8rings_hfi__ring_IQU
          noise__pix_IQU&
               & => noise_hfi__pix_IQU_channel(:,:,i-nr_channels(1))
          alm__TEB_l_m&
               & => alm_hfi__TEB_l_m
       ENDIF

       CALL MAP2ALM(nside(i), lmax(i), mmax(i), noise__pix_IQU(:,1),&
            & alm__TEB_l_m, zbounds__lat, w8rings__ring_IQU)

       DO j = 0,lmax(i)
          hlc_alm__TEB_l_m(1,j,0:mmax(i))&
               & = hlc_alm__TEB_l_m(1,j,0:mmax(i))&
               & + hlc_weights__l_channel(j,i) * alm__TEB_l_m(1,j,:)
       ENDDO

    ENDDO

    CALL ALM2MAP(nside(4), MAXVAL(lmax), MAXVAL(mmax),&
         & hlc_alm__TEB_l_m, hlc_map__pix_IQU(:,1))

    NULLIFY(w8rings__ring_IQU)
    NULLIFY(noise__pix_IQU)
    NULLIFY(alm__TEB_l_m)
    DEALLOCATE(alm_lfi__TEB_l_m)
    DEALLOCATE(alm_hfi__TEB_l_m)
    DEALLOCATE(hlc_alm__TEB_l_m)


  END SUBROUTINE COMBINE_MAPS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read map from file
!-----------------------------------------------------------------------

  SUBROUTINE READ_MAP(map__pix_IQU, input_file, npix, nIQU)

    USE FITSTOOLS,     ONLY: INPUT_MAP, GETSIZE_FITS
    USE PIX_TOOLS,     ONLY: CONVERT_INPLACE, NEST2RING

    CHARACTER(LEN=*),                       INTENT(IN)    :: input_file
    INTEGER(I4B),                           INTENT(IN)    :: npix
    INTEGER(I4B),                           INTENT(IN)    :: nIQU
    REAL(DP),          DIMENSION(0:,1:),    INTENT(INOUT) :: map__pix_IQU

    CHARACTER(LEN=80), DIMENSION(1)                       :: header
    INTEGER(I4B)                                          :: npix_file
    INTEGER(I4B)                                          :: order_file

    npix_file = GETSIZE_FITS(input_file, ordering=order_file)

    CALL CHECK_ERROR_FLAG("READ_MAP, size mismatch", npix_file - npix)

    CALL INPUT_MAP(input_file, map__pix_IQU, npix, nIQU,&
         & 0.0_DP, header)

    IF (order_file .EQ. 2) THEN
       CALL CONVERT_INPLACE(NEST2RING, map__pix_IQU)
    ENDIF


  END SUBROUTINE READ_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Write map to file
!-----------------------------------------------------------------------

  SUBROUTINE OUTPUT_MAP(filenames__channel_nr, map__pix_IQU,&
          & nside, npix, lmax, mmax, nr_channels, map_nr)

    USE MISC_UTILS,    ONLY: ASSERT_NOT_PRESENT
    USE HEAD_FITS,     ONLY: WRITE_MINIMAL_HEADER
    USE FITSTOOLS,     ONLY: WRITE_BINTAB

    CHARACTER(LEN=*),  DIMENSION(1:,1:),    INTENT(IN)    :: filenames__channel_nr
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nside
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: npix
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: lmax
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: mmax
    INTEGER(I4B),                           INTENT(IN)    :: map_nr
    INTEGER(I4B),      DIMENSION(1:),       INTENT(IN)    :: nr_channels
    REAL(DP),          DIMENSION(0:,1:),    INTENT(IN)    :: map__pix_IQU

    CHARACTER(LEN=150)                                    :: filename
    CHARACTER(LEN=80), DIMENSION(1:80)                    :: header

    CALL WRITE_MINIMAL_HEADER(header, 'MAP', append=.FALSE.,&
         & nside=MAXVAL(nside), ordering='RING', coordsys='G',&
         & nlmax=MAXVAL(lmax), polar=.FALSE., nmmax=MAXVAL(mmax))

    filename = TRIM(filenames__channel_nr(nr_channels(3)+1,map_nr))

    CALL ASSERT_NOT_PRESENT(filename)
    CALL WRITE_BINTAB(REAL(map__pix_IQU, KIND=SP), MAXVAL(npix), 1,&
         & header, 80, filename)


  END SUBROUTINE OUTPUT_MAP

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compute loop boundaries for each thread
!-----------------------------------------------------------------------

  SUBROUTINE GET_LOOP_BOUNDS(loop_min, loop_max, loop_start,&
       & loop_stop, thread_id, nr_chunks_total)

    INTEGER(I4B),                           INTENT(IN)    :: loop_start
    INTEGER(I4B),                           INTENT(IN)    :: loop_stop
    INTEGER(I4B),                           INTENT(IN)    :: thread_id
    INTEGER(I4B),                           INTENT(IN)    :: nr_chunks_total
    INTEGER(I4B),                           INTENT(OUT)   :: loop_min
    INTEGER(I4B),                           INTENT(OUT)   :: loop_max

    INTEGER(I4B)                                          :: nr_iterations
    REAL(DP)                                              :: index_fraction
    REAL(DP)                                              :: index_min
    REAL(DP)                                              :: index_max

    nr_iterations = loop_stop - loop_start + 1

    IF ((thread_id .LT. 0) .OR. (nr_chunks_total .LT. thread_id - 1)&
         & .OR. (nr_iterations .LT. nr_chunks_total)) THEN
       WRITE(*,'(3X, "Error: Loop indices")')
       STOP
    ENDIF

    index_fraction = REAL(nr_iterations, KIND=DP)&
         & / REAL(nr_chunks_total, KIND=DP)

    index_min = index_fraction * REAL(thread_id, KIND=DP)
    index_max = index_fraction * REAL(thread_id + 1, KIND=DP)

    loop_min = loop_start + INT(index_min, KIND=I4B)
    loop_max = loop_start + INT(index_max, KIND=I4B) - 1

    IF (thread_id .EQ. nr_chunks_total - 1) THEN
       loop_max = loop_stop
    ENDIF

    IF (loop_max .LT. loop_min) THEN
       WRITE(*,'(3X, "Error: Loop indices")')
       STOP
    ENDIF


  END SUBROUTINE GET_LOOP_BOUNDS

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Check error flag
!-----------------------------------------------------------------------

  SUBROUTINE CHECK_ERROR_FLAG(code, error)

    CHARACTER(LEN=*),                       INTENT(IN)    :: code
    INTEGER(I4B),                           INTENT(IN)    :: error

    IF (error .NE. 0) THEN
       WRITE(*,'(/3X, "Error detected by ", A"; errorcode ", I)')&
            & code, error
       STOP
    ENDIF


  END SUBROUTINE CHECK_ERROR_FLAG

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Clock program
!-----------------------------------------------------------------------

  SUBROUTINE CLOCK(thread_id)

    USE MISC_UTILS,    ONLY: WALL_CLOCK_TIME

    INTEGER(I4B),                           INTENT(IN)    :: thread_id

    LOGICAL(LGT),                           SAVE          :: firstcall = .TRUE.
    REAL(SP),                               SAVE          :: wct_start
    REAL(SP),                               SAVE          :: wct_stop

    IF (thread_id .NE. 0) RETURN

    IF (firstcall) THEN
       CALL WALL_CLOCK_TIME(wct_start)
       firstcall = .FALSE.
    ELSE
       CALL WALL_CLOCK_TIME(wct_stop)
       WRITE(*,'(3X, "Time: ", F8.1, " min")')&
            & (wct_stop - wct_start) / 60.0_SP
    ENDIF


  END SUBROUTINE CLOCK

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Syncronize and then abort program
!-----------------------------------------------------------------------

  SUBROUTINE BURN

    INTEGER(I4B)                                          :: ierror

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    CALL MPI_ABORT(MPI_COMM_WORLD, ierror, ierror)


  END SUBROUTINE BURN

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Cleanup
!-----------------------------------------------------------------------

  SUBROUTINE DEALLOCATE_ARRAYS(noise_lfi__pix_IQU_channel,&
       & noise_hfi__pix_IQU_channel, hlc_map__pix_IQU,&
       & hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
       & w8rings_hfi__ring_IQU, filenames__channel_nr)

    CHARACTER(LEN=*),  DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(INOUT) :: filenames__channel_nr
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(INOUT) :: w8rings_lfi__ring_IQU
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(INOUT) :: w8rings_hfi__ring_IQU
    REAL(DP),          DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(INOUT) :: hlc_map__pix_IQU
    REAL(DP),          DIMENSION(:,:,:),    ALLOCATABLE,&
         &                                  INTENT(INOUT) :: noise_lfi__pix_IQU_channel
    REAL(DP),          DIMENSION(:,:,:),    ALLOCATABLE,&
         &                                  INTENT(INOUT) :: noise_hfi__pix_IQU_channel
    COMPLEX(DPC),      DIMENSION(:,:),      ALLOCATABLE,&
         &                                  INTENT(INOUT) :: hlc_weights__l_channel

    DEALLOCATE(filenames__channel_nr)
    DEALLOCATE(w8rings_lfi__ring_IQU)
    DEALLOCATE(w8rings_hfi__ring_IQU)
    DEALLOCATE(hlc_map__pix_IQU)
    DEALLOCATE(noise_lfi__pix_IQU_channel)
    DEALLOCATE(noise_hfi__pix_IQU_channel)
    DEALLOCATE(hlc_weights__l_channel)


  END SUBROUTINE DEALLOCATE_ARRAYS

!-----------------------------------------------------------------------

END MODULE COMBINE_MODULES
