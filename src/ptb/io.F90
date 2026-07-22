! This file is part of xtb.
!
! Copyright (C) 2024 xtb developers
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Export of PTB density matrix in the atomic-orbital basis.
module xtb_ptb_io
#if WITH_TBLITE
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_symbols, only: to_symbol
   use tblite_basis_type, only: basis_type, cgto_type
   implicit none
   private

   !> In-memory .npy payload for one member of a .npz (ZIP) archive.
   type :: zip_member
      character(len=:), allocatable :: name
      integer(int8), allocatable :: bytes(:)
   end type zip_member

   public :: write_ptb_matrix_npy
   public :: write_ptb_ao_order
   public :: write_ptb_matrix_npz_csr

contains

   !> Write a dense matrix to a NumPy .npy file (format version 1.0).
   !>
   !> The matrix is stored in its native column-major (Fortran) order, so it is
   !> read back without transposition by numpy.load. The row/column AO ordering
   !> is that of write_ptb_ao_order (tblite's m = -l..+l per shell).
   !>
   !> Args:
   !>   filename: Path of the .npy file to create.
   !>   mat: Matrix in the spherical AO basis.
   subroutine write_ptb_matrix_npy(filename, mat)
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: mat(:, :)

      integer :: unit
      character(len=40) :: shapebuffer

      open (newunit=unit, file=filename, access='stream', form='unformatted', &
         & status='replace')
      write (shapebuffer, '(a,i0,a,i0,a)') "(", size(mat, 1), ", ", size(mat, 2), ")"
      call write_npy_header(unit, trim(shapebuffer))
      write (unit) mat
      close (unit)
   end subroutine write_ptb_matrix_npy

   !> Write the atomic-orbital ordering of the exported matrices to a NumPy .npy
   !> file as an integer array of shape (nao, 3).
   !>
   !> Row iao holds (atom, l, m) for the iao-th row/column of the density matrix.
   !> This makes the AO layout explicit so a consumer can build the permutation
   !> to any target program's convention without assuming ours.
   !>
   !> Args:
   !>   filename: Path of the .npy file to create.
   !>   bas: Persistent PTB basis set.
   subroutine write_ptb_ao_order(filename, bas)
      character(len=*), intent(in) :: filename
      type(basis_type), intent(in) :: bas

      integer :: unit, iao, iat, ishg, angmom
      integer(int64), allocatable :: ao_order(:, :)
      character(len=40) :: shapebuffer

      allocate (ao_order(bas%nao, 3))
      do iao = 1, bas%nao
         iat = bas%ao2at(iao)
         ishg = bas%ao2sh(iao)
         angmom = bas%cgto(ishg - bas%ish_at(iat), iat)%ang
         ao_order(iao, 1) = int(iat, int64)
         ao_order(iao, 2) = int(angmom, int64)
         ao_order(iao, 3) = int(iao - bas%iao_sh(ishg) - 1 - angmom, int64)
      end do

      open (newunit=unit, file=filename, access='stream', form='unformatted', &
         & status='replace')
      write (shapebuffer, '(a,i0,a)') "(", bas%nao, ", 3)"
      call write_npy_header(unit, trim(shapebuffer), descr="<i8")
      write (unit) ao_order
      close (unit)
   end subroutine write_ptb_ao_order

   !> Write a matrix as a SciPy CSR sparse matrix in a .npz file, loadable via
   !> scipy.sparse.load_npz. Only elements with abs(value) above the threshold
   !> are stored. The row/column AO ordering is that of write_ptb_ao_order.
   !>
   !> Args:
   !>   filename: Path of the .npz file to create.
   !>   mat: Matrix in the spherical AO basis.
   !>   threshold: Elements with abs(value) <= threshold are dropped.
   subroutine write_ptb_matrix_npz_csr(filename, mat, threshold)
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: mat(:, :)
      real(wp), intent(in) :: threshold

      integer :: unit, nao, irow, icol, pos
      integer(int64) :: nnz
      integer(int64), allocatable :: indptr(:)
      integer(int64), allocatable :: indices(:)
      real(wp), allocatable :: values(:)
      type(zip_member), allocatable :: members(:)

      nao = size(mat, 1)

      nnz = 0
      do irow = 1, nao
         do icol = 1, nao
            if (abs(mat(irow, icol)) > threshold) nnz = nnz + 1
         end do
      end do

      allocate (indptr(nao + 1), indices(nnz), values(nnz))
      indptr(1) = 0
      pos = 0
      do irow = 1, nao
         do icol = 1, nao
            if (abs(mat(irow, icol)) > threshold) then
               pos = pos + 1
               indices(pos) = int(icol - 1, int64)
               values(pos) = mat(irow, icol)
            end if
         end do
         indptr(irow + 1) = int(pos, int64)
      end do

      allocate (members(5))
      members(1) = npy_member_char("format.npy", "csr")
      members(2) = npy_member_int64("shape.npy", [int(nao, int64), int(nao, int64)])
      members(3) = npy_member_real("data.npy", values)
      members(4) = npy_member_int64("indices.npy", indices)
      members(5) = npy_member_int64("indptr.npy", indptr)

      open (newunit=unit, file=filename, access='stream', form='unformatted', &
         & status='replace')
      call write_zip_archive(unit, members)
      close (unit)
   end subroutine write_ptb_matrix_npz_csr

   !> Write a NumPy .npy version 1.0 header for a given dtype descriptor and
   !> shape tuple to a stream-access unit positioned at the start of the file.
   subroutine write_npy_header(unit, shapestr, descr, fortran_order)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: shapestr
      character(len=*), intent(in), optional :: descr
      logical, intent(in), optional :: fortran_order

      character(len=:), allocatable :: header, dtype, order
      integer :: total_len, pad, header_len

      dtype = "<f8"
      if (present(descr)) dtype = descr
      order = "True"
      if (present(fortran_order)) then
         if (.not. fortran_order) order = "False"
      end if

      header = "{'descr': '"//dtype//"', 'fortran_order': "//trim(order)// &
         & ", 'shape': "//shapestr//", }"

      total_len = 10 + len(header) + 1
      pad = mod(64 - mod(total_len, 64), 64)
      header = header//repeat(" ", pad)//char(10)
      header_len = len(header)

      write (unit) int(-109, int8)
      write (unit) "NUMPY"
      write (unit) int(1, int8), int(0, int8)
      write (unit) int(header_len, int16)
      write (unit) header
   end subroutine write_npy_header

   !> Serialize a real(wp) array into an in-memory .npy payload.
   function npy_member_real(name, array) result(member)
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: array(:)
      type(zip_member) :: member

      integer :: unit
      character(len=40) :: shapebuffer

      open (newunit=unit, status='scratch', access='stream', form='unformatted')
      write (shapebuffer, '(a,i0,a)') "(", size(array), ",)"
      call write_npy_header(unit, trim(shapebuffer))
      write (unit) array
      call slurp_scratch(unit, member%bytes)
      close (unit)
      member%name = name
   end function npy_member_real

   !> Serialize an integer(int64) array into an in-memory .npy payload.
   function npy_member_int64(name, array) result(member)
      character(len=*), intent(in) :: name
      integer(int64), intent(in) :: array(:)
      type(zip_member) :: member

      integer :: unit
      character(len=40) :: shapebuffer

      open (newunit=unit, status='scratch', access='stream', form='unformatted')
      write (shapebuffer, '(a,i0,a)') "(", size(array), ",)"
      call write_npy_header(unit, trim(shapebuffer), descr="<i8")
      write (unit) array
      call slurp_scratch(unit, member%bytes)
      close (unit)
      member%name = name
   end function npy_member_int64

   !> Serialize a short ASCII string into an in-memory .npy payload (a
   !> zero-dimensional byte-string array).
   function npy_member_char(name, string) result(member)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: string
      type(zip_member) :: member

      integer :: unit
      character(len=40) :: descrbuffer

      open (newunit=unit, status='scratch', access='stream', form='unformatted')
      write (descrbuffer, '(a,i0)') "|S", len(string)
      call write_npy_header(unit, "()", descr=trim(descrbuffer), fortran_order=.false.)
      write (unit) string
      call slurp_scratch(unit, member%bytes)
      close (unit)
      member%name = name
   end function npy_member_char

   !> Read the entire contents of a scratch stream unit into a byte array.
   subroutine slurp_scratch(unit, bytes)
      integer, intent(in) :: unit
      integer(int8), allocatable, intent(out) :: bytes(:)

      integer(int64) :: nbytes

      inquire (unit=unit, size=nbytes)
      allocate (bytes(nbytes))
      read (unit, pos=1) bytes
   end subroutine slurp_scratch

   !> Write an uncompressed (stored) ZIP archive containing the given members.
   subroutine write_zip_archive(unit, members)
      integer, intent(in) :: unit
      type(zip_member), intent(in) :: members(:)

      integer :: imember, nmember
      integer(int64), allocatable :: local_offset(:)
      integer(int32), allocatable :: crc(:)
      integer(int64) :: central_start, central_size
      integer(int64) :: offset

      nmember = size(members)
      allocate (local_offset(nmember), crc(nmember))

      offset = 0
      do imember = 1, nmember
         local_offset(imember) = offset
         crc(imember) = crc32(members(imember)%bytes)
         call write_local_header(unit, members(imember), crc(imember))
         write (unit) members(imember)%bytes
         offset = offset + 30_int64 + len(members(imember)%name) &
            & + size(members(imember)%bytes)
      end do

      central_start = offset
      do imember = 1, nmember
         call write_central_header(unit, members(imember), crc(imember), &
            & local_offset(imember))
      end do

      central_size = 0
      do imember = 1, nmember
         central_size = central_size + 46_int64 + len(members(imember)%name)
      end do

      call write_end_of_central_directory(unit, nmember, central_size, central_start)
   end subroutine write_zip_archive

   !> Write a ZIP local file header followed by nothing (data written by caller).
   subroutine write_local_header(unit, member, crc)
      integer, intent(in) :: unit
      type(zip_member), intent(in) :: member
      integer(int32), intent(in) :: crc

      write (unit) int(z'04034b50', int32)
      write (unit) int(20, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) crc
      write (unit) int(size(member%bytes), int32)
      write (unit) int(size(member%bytes), int32)
      write (unit) int(len(member%name), int16)
      write (unit) int(0, int16)
      write (unit) member%name
   end subroutine write_local_header

   !> Write a ZIP central-directory file header for one member.
   subroutine write_central_header(unit, member, crc, local_offset)
      integer, intent(in) :: unit
      type(zip_member), intent(in) :: member
      integer(int32), intent(in) :: crc
      integer(int64), intent(in) :: local_offset

      write (unit) int(z'02014b50', int32)
      write (unit) int(20, int16)
      write (unit) int(20, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) crc
      write (unit) int(size(member%bytes), int32)
      write (unit) int(size(member%bytes), int32)
      write (unit) int(len(member%name), int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(0, int32)
      write (unit) int(local_offset, int32)
      write (unit) member%name
   end subroutine write_central_header

   !> Write the ZIP end-of-central-directory record closing the archive.
   subroutine write_end_of_central_directory(unit, nmember, central_size, central_start)
      integer, intent(in) :: unit
      integer, intent(in) :: nmember
      integer(int64), intent(in) :: central_size
      integer(int64), intent(in) :: central_start

      write (unit) int(z'06054b50', int32)
      write (unit) int(0, int16)
      write (unit) int(0, int16)
      write (unit) int(nmember, int16)
      write (unit) int(nmember, int16)
      write (unit) int(central_size, int32)
      write (unit) int(central_start, int32)
      write (unit) int(0, int16)
   end subroutine write_end_of_central_directory

   !> Compute the CRC-32 (IEEE 802.3, reflected) checksum of a byte array, as
   !> required by the ZIP format.
   pure function crc32(bytes) result(checksum)
      integer(int8), intent(in) :: bytes(:)
      integer(int32) :: checksum

      integer :: ibyte, ibit
      integer(int32) :: crc, lowbit
      !> Reflected CRC-32 polynomial 0xEDB88320 as a signed 32-bit integer.
      integer(int32), parameter :: polynomial = -306674912_int32

      crc = not(0_int32)
      do ibyte = 1, size(bytes)
         crc = ieor(crc, iand(int(bytes(ibyte), int32), 255_int32))
         do ibit = 1, 8
            lowbit = iand(crc, 1_int32)
            crc = ishft(crc, -1)
            if (lowbit /= 0) crc = ieor(crc, polynomial)
         end do
      end do
      checksum = not(crc)
   end function crc32

#endif
end module xtb_ptb_io
