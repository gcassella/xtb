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
   use, intrinsic :: iso_fortran_env, only: int8, int16
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_symbols, only: to_symbol
   use tblite_basis_type, only: basis_type, cgto_type
   implicit none
   private

   !> Double factorial (2*l-1)!! for l = 0..7, see OEIS A001147.
   real(wp), parameter :: double_factorial(0:7) = &
      & [1.0_wp, 1.0_wp, 3.0_wp, 15.0_wp, 105.0_wp, 945.0_wp, 10395.0_wp, 135135.0_wp]

   public :: write_ptb_matrix_npy
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

   !> Write the vDZP basis set in NWChem format, faithful to the normalized AO
   !> basis in which the exported density and overlap matrices are expressed.
   !>
   !> One block is emitted per element under the plain element symbol. The PTB
   !> basis carries no charge dependence (base vDZP exponents, expscal = 1), so
   !> all atoms of the same element share an identical basis. The contraction
   !> coefficients are scaled to unit self-overlap so that NWChem's default
   !> primitive renormalization reproduces the functions unchanged.
   !>
   !> The exported matrices use tblite's spherical-harmonic component ordering,
   !> which differs from NWChem's for p shells (tblite m = -1, 0, +1 versus
   !> NWChem/PySCF px, py, pz); a consumer must reorder the p rows/columns of
   !> the matrices accordingly.
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

      integer :: isp, iat, ish, ishg, iprim, angmom
      real(wp) :: shellnorm
      character(len=:), allocatable :: symbol
      type(cgto_type) :: cgto

      do isp = 1, mol%nid
         iat = first_atom_of_species(mol, isp)
         symbol = trim(to_symbol(mol%num(isp)))
         write (unit, '(a)') "#BASIS SET: PTB vDZP"
         do ish = 1, bas%nsh_at(iat)
            ishg = bas%ish_at(iat) + ish
            cgto = bas%cgto(ish, iat)
            angmom = cgto%ang
            shellnorm = aonorm(bas%iao_sh(ishg) + 1)
            write (unit, '(a,4x,a)') symbol, nwchem_angmom_label(angmom)
            do iprim = 1, cgto%nprim
               write (unit, '(4x,es24.16,4x,es24.16)') cgto%alpha(iprim), &
                  & nwchem_primitive_coeff(cgto%alpha(iprim), cgto%ang, &
                  & cgto%coeff(iprim), shellnorm)
            end do
         end do
      end do
      write (unit, '(a)') "END"
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

   !> Index of the first atom belonging to a given species.
   pure function first_atom_of_species(mol, isp) result(iat)
      type(structure_type), intent(in) :: mol
      integer, intent(in) :: isp
      integer :: iat

      do iat = 1, mol%nat
         if (mol%id(iat) == isp) return
      end do
      iat = 0
   end function first_atom_of_species

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

#endif
end module xtb_ptb_io
