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
PROGRAM COMBINE_MAIN

  USE MPI
  USE HEALPIX_TYPES
  USE COMBINE_MODULES

  IMPLICIT NONE

  CHARACTER(LEN=150)                                      :: dir_w8ring
  CHARACTER(LEN=150)                                      :: file_hlc_coeff
  CHARACTER(LEN=150)                                      :: file_filenamelist
  CHARACTER(LEN=150),  DIMENSION(:,:),      ALLOCATABLE   :: filenames__channel_nr
  INTEGER(I4B)                                            :: nr_maps
  INTEGER(I4B)                                            :: map_nr
  INTEGER(I4B)                                            :: nside_lfi
  INTEGER(I4B)                                            :: nside_hfi
  INTEGER(I4B)                                            :: loop_min, loop_max
  INTEGER(I4B)                                            :: ierror
  INTEGER(I4B)                                            :: thread_id, nr_cpus
  INTEGER(I4B),        DIMENSION(1:3)                     :: nr_channels
  INTEGER(I4B),        DIMENSION(1:9)                     :: nside
  INTEGER(I4B),        DIMENSION(1:9)                     :: npix
  INTEGER(I4B),        DIMENSION(1:9)                     :: lmax
  INTEGER(I4B),        DIMENSION(1:9)                     :: mmax
  REAL(DP),            DIMENSION(1:2)                     :: zbounds__lat
  REAL(DP),            DIMENSION(:,:),      ALLOCATABLE   :: w8rings_lfi__ring_IQU
  REAL(DP),            DIMENSION(:,:),      ALLOCATABLE   :: w8rings_hfi__ring_IQU
  REAL(DP),            DIMENSION(:,:),      ALLOCATABLE   :: hlc_map__pix_IQU
  REAL(DP),            DIMENSION(:,:,:),    ALLOCATABLE   :: noise_lfi__pix_IQU_channel
  REAL(DP),            DIMENSION(:,:,:),    ALLOCATABLE   :: noise_hfi__pix_IQU_channel
  COMPLEX(DPC),        DIMENSION(:,:),      ALLOCATABLE   :: hlc_weights__l_channel

!-----------------------------------------------------------------------
! Parse input parameters
!-----------------------------------------------------------------------

  NAMELIST / input_parameters / dir_w8ring, file_hlc_coeff,&
       & file_filenamelist

  OPEN(UNIT=7, FILE="input.cfg", STATUS='OLD', ACTION='READ')
  READ(7, NML=input_parameters)
  CLOSE(7)

  nside_lfi = 1024
  nside_hfi = 2048

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Initialize MPI
!-----------------------------------------------------------------------

  CALL MPI_INIT(ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nr_cpus, ierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, thread_id, ierror)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Allocate memory
!-----------------------------------------------------------------------

  CALL ALLOCATE_ARRAYS(noise_lfi__pix_IQU_channel,&
       & noise_hfi__pix_IQU_channel, hlc_map__pix_IQU,&
       & hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
       & w8rings_hfi__ring_IQU, filenames__channel_nr, zbounds__lat,&
       & nr_channels, nside, npix, lmax, mmax, nside_lfi, nside_hfi)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Read auxiliary files
!-----------------------------------------------------------------------

  CALL READ_AUX_FILES(hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
       & w8rings_hfi__ring_IQU, filenames__channel_nr, nr_maps, nside,&
       & lmax, mmax, nr_channels, thread_id, dir_w8ring,&
       & file_hlc_coeff, file_filenamelist)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! Compute loop boundaries for each MPI thread
!-----------------------------------------------------------------------

  CALL GET_LOOP_BOUNDS(loop_min, loop_max, 1, nr_maps,&
       & thread_id, nr_cpus)

!-----------------------------------------------------------------------


  CALL CLOCK(thread_id)

!-----------------------------------------------------------------------
! Sample maps
!-----------------------------------------------------------------------

  DO map_nr = loop_min,loop_max
     WRITE(*,'(/20X, "--- Map ", 1X, I4, " ---"/)') map_nr


   !--------------------------------------------------------------------
   ! Load noise realizations
   !--------------------------------------------------------------------

     CALL LOAD_NOISE_MAPS(noise_lfi__pix_IQU_channel,&
          & noise_hfi__pix_IQU_channel, filenames__channel_nr,&
          & nr_channels, npix, map_nr)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Combine to harmonic ILC map
   !--------------------------------------------------------------------

     CALL COMBINE_MAPS(hlc_map__pix_IQU,&
          & noise_lfi__pix_IQU_channel, noise_hfi__pix_IQU_channel,&
          & hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
          & w8rings_hfi__ring_IQU, zbounds__lat, nr_channels, nside,&
          & lmax, mmax)

   !--------------------------------------------------------------------



   !--------------------------------------------------------------------
   ! Output residual noise maps
   !--------------------------------------------------------------------

     CALL OUTPUT_MAP(filenames__channel_nr, hlc_map__pix_IQU,&
          & nside, npix, lmax, mmax, nr_channels, map_nr)

   !--------------------------------------------------------------------


  ENDDO

!-----------------------------------------------------------------------

  CALL CLOCK(thread_id)


!-----------------------------------------------------------------------
! Cleanup
!-----------------------------------------------------------------------

  CALL DEALLOCATE_ARRAYS(noise_lfi__pix_IQU_channel,&
       & noise_hfi__pix_IQU_channel, hlc_map__pix_IQU,&
       & hlc_weights__l_channel, w8rings_lfi__ring_IQU,&
       & w8rings_hfi__ring_IQU, filenames__channel_nr)

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! MPI Cleanup
!-----------------------------------------------------------------------

  CALL MPI_FINALIZE(ierror)

!-----------------------------------------------------------------------

END PROGRAM COMBINE_MAIN
