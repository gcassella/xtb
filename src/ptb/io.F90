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

!> Export of PTB wavefunction data in the atomic-orbital basis: the density
!> matrix, the overlap matrix and the vDZP basis set (in NWChem format).
!> Together these are sufficient to reconstruct the real-space single-particle
!> density rho(r) = sum_{mu,nu} P_{mu,nu} phi_mu(r) phi_nu(r).
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

   !> Double factorial (2*l-1)!! for l = 0..7, see OEIS A001147.
   real(wp), parameter :: double_factorial(0:7) = &
      & [1.0_wp, 1.0_wp, 3.0_wp, 15.0_wp, 105.0_wp, 945.0_wp, 10395.0_wp, 135135.0_wp]

   public :: write_ptb_matrix_npy
   public :: write_ptb_matrix_npz_csr
   public :: write_ptb_basis_nwchem
   public :: nwchem_primitive_coeff
   public :: primitive_normalizer

   !> Angular momentum labels used by the NWChem basis format, indexed by l.
   character(len=1), parameter :: nwchem_angmom_label(0:6) = &
      & ["S", "P", "D", "F", "G", "H", "I"]

contains

   !> Write a dense matrix to a NumPy .npy file (format version 1.0).
   !>
   !> The matrix is stored in its native column-major (Fortran) order, so it is
   !> read back without transposition by numpy.load. The AO ordering follows
   !> tblite's spherical-harmonic convention (m = -l..+l per shell), which
   !> differs from NWChem's; a consumer combining these matrices with the
   !> exported NWChem basis must reorder rows/columns accordingly.
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

   !> Write a matrix as a SciPy CSR sparse matrix in a .npz file, loadable via
   !> scipy.sparse.load_npz. Only elements with abs(value) above the threshold
   !> are stored. The AO ordering caveat of write_ptb_matrix_npy applies.
   !>
   !> The archive uses 32-bit ZIP size fields, so an individual stored member
   !> must stay below 2 GB (not ZIP64); this is not a concern for PTB systems.
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

   !> Write the vDZP basis set in NWChem format, faithful to the normalized AO
   !> basis in which the exported density and overlap matrices are expressed.
   !>
   !> Each atom is emitted under a unique tag (element symbol + one-based index,
   !> e.g. "O1", "H2"); a consumer must apply the same tags to the geometry. The
   !> contraction coefficients are scaled to unit self-overlap so that NWChem's
   !> default primitive renormalization reproduces the functions unchanged.
   !>
   !> The exported matrices use tblite's spherical-harmonic component ordering,
   !> which differs from NWChem's; a consumer reading this basis through NWChem
   !> must reorder the matrix rows/columns accordingly.
   !>
   !> Args:
   !>   unit: Open, writable file unit.
   !>   mol: Molecular structure data.
   !>   bas: Persistent PTB basis set (must use base vDZP exponents, expscal = 1).
   !>   aonorm: Per-atomic-orbital normalization factors from PTB.
   subroutine write_ptb_basis_nwchem(unit, mol, bas, aonorm)
      integer, intent(in) :: unit
      type(structure_type), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      real(wp), intent(in) :: aonorm(:)

      integer :: iat, ish, ishg, iprim, angmom
      real(wp) :: shellnorm
      character(len=:), allocatable :: tag
      type(cgto_type) :: cgto

      write (unit, '(a)') "# PTB vDZP basis set in NWChem format"
      write (unit, '(a)') "# per-atom tags (element symbol + atom index) must be"
      write (unit, '(a)') "# matched by the geometry block of the consumer"
      write (unit, '(a)') 'basis "ao basis" spherical'
      do iat = 1, mol%nat
         tag = atom_tag(mol%num(mol%id(iat)), iat)
         do ish = 1, bas%nsh_at(iat)
            ishg = bas%ish_at(iat) + ish
            cgto = bas%cgto(ish, iat)
            angmom = cgto%ang
            shellnorm = aonorm(bas%iao_sh(ishg) + 1)
            write (unit, '(a,4x,a)') tag, nwchem_angmom_label(angmom)
            do iprim = 1, cgto%nprim
               write (unit, '(4x,es24.16,4x,es24.16)') cgto%alpha(iprim), &
                  & nwchem_primitive_coeff(cgto%alpha(iprim), cgto%ang, &
                  & cgto%coeff(iprim), shellnorm)
            end do
         end do
      end do
      write (unit, '(a)') "end"
   end subroutine write_ptb_basis_nwchem

   !> Contraction coefficient written to the NWChem basis for one primitive.
   !>
   !> NWChem stores bare primitive coefficients and re-applies the primitive
   !> normalization on read, so this divides out the normalization tblite folded
   !> into the coefficient and includes the per-shell factor relating the raw
   !> contracted function to the normalized AO in which P and S are expressed.
   !>
   !> Args:
   !>   alpha: Primitive Gaussian exponent.
   !>   angmom: Angular momentum of the shell.
   !>   coeff: Contraction coefficient as stored by tblite (normalized primitive).
   !>   shellnorm: Per-shell AO normalization factor.
   pure function nwchem_primitive_coeff(alpha, angmom, coeff, shellnorm) result(bare)
      real(wp), intent(in) :: alpha
      integer, intent(in) :: angmom
      real(wp), intent(in) :: coeff
      real(wp), intent(in) :: shellnorm
      real(wp) :: bare

      bare = coeff / primitive_normalizer(alpha, angmom) * shellnorm
   end function nwchem_primitive_coeff

   !> Normalization factor tblite folds into a primitive Gaussian's contraction
   !> coefficient, used here to recover the bare primitive.
   pure function primitive_normalizer(alpha, angmom) result(normfac)
      real(wp), intent(in) :: alpha
      integer, intent(in) :: angmom
      real(wp) :: normfac

      normfac = (2.0_wp * alpha / pi)**0.75_wp * sqrt(4.0_wp * alpha)**angmom &
         & / sqrt(double_factorial(angmom))
   end function primitive_normalizer

   !> Build a unique per-atom tag from the element symbol and atom index, e.g. "O1".
   pure function atom_tag(znum, iat) result(tag)
      integer, intent(in) :: znum
      integer, intent(in) :: iat
      character(len=:), allocatable :: tag

      character(len=16) :: idxbuffer

      write (idxbuffer, '(i0)') iat
      tag = trim(to_symbol(znum))//trim(idxbuffer)
   end function atom_tag

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
   !> zero-dimensional byte-string array, as SciPy expects for the CSR tag).
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
   !> This produces a valid .npz file readable by numpy/scipy.
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
