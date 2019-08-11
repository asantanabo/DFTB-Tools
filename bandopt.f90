module global

! Number of atoms in a unit cell
  integer :: Natom
! Number of cells on three direction of supercell
  integer :: Ncell(3)
! number of total atoms
  integer :: Ntotal
! number of k-space scanning path
  integer :: npath
! number of total scanning
!  integer :: nktotal
! grid for each path
  integer, allocatable :: nkgrid(:)
! symmetric points (start and end)
  double precision, allocatable :: kpath(:, :, :)

! parameters of crystal
  double precision :: length_a, length_b, length_c
  double precision :: angle_alpha, angle_beta, angle_gamma

! molar mass of atoms
!  double precision, allocatable :: molarmass(:)
! name of atoms in supercell
!  character*3, allocatable :: chatom(:)

! molar mass list for elements
  integer :: atom_number
  parameter (atom_number = 110)
  double precision :: atom_mass(atom_number)
  parameter (atom_mass = (/1.00790d0, 4.00260d0, 6.9410d0, 9.01220d0, 10.8110d0, &
                           12.01070d0, 14.00670d0, 15.99940d0, 18.99840d0, 20.17970d0, &
                           22.98970d0, 24.3050d0, 26.98150d0, 28.08550d0, 30.97380d0, &
                           32.0650d0, 35.4530d0, 39.09830d0, 39.9480d0, 40.0780d0, &
                           44.95590d0, 47.8670d0, 50.94150d0, 51.99610d0, 54.9380d0, &
                           55.8450d0, 58.69340d0, 58.93320d0, 63.5460d0, 65.390d0, &
                           69.7230d0, 72.640d0, 74.92160d0, 78.960d0, 79.9040d0, 83.80d0, &
                           85.46780d0, 87.620d0, 88.90590d0, 91.2240d0, 92.90640d0, 95.940d0, &
                           98.0d0, 101.070d0, 102.90550d0, 106.420d0, 107.86820d0, 112.4110d0, &
                           114.8180d0, 118.710d0, 121.760d0, 126.90450d0, 127.60d0, 131.2930d0, &
                           132.90550d0, 137.3270d0, 138.90550d0, 140.1160d0, 140.90770d0, &
                           144.240d0, 1450d0, 150.360d0, 151.9640d0, 157.250d0, 158.92530d0, &
                           162.50d0, 164.93030d0, 167.2590d0, 168.93420d0, 173.040d0, &
                           174.9670d0, 178.490d0, 180.94790d0, 183.840d0, 186.2070d0, 190.230d0, &
                           192.2170d0, 195.0780d0, 196.96650d0, 200.590d0, 204.38330d0, 207.20d0, &
                           208.98040d0, 209.0d0, 210.0d0, 222.0d0, 223.0d0, 226.0d0, 227.0d0, &
                           231.03590d0, 232.03810d0, 237.0d0, 238.02890d0, 243.0d0, 244.0d0, &
                           247.0d0, 247.0d0, 251.0d0, 252.0d0, 257.0d0, 258.0d0, 259.0d0, 261.0d0, &
                           262.0d0, 262.0d0, 264.0d0, 266.0d0, 268.0d0, 272.0d0, 277.0d0/)) 
  character*3 :: atom_name(atom_number)
  parameter (atom_name = (/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
                           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'K ', 'Ar', 'Ca', &
                           'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Ni', 'Co', 'Cu', 'Zn', &
                           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
                           'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
                           'Sb', 'I ', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
                           'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
                           'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
                           'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Pa', &
                           'Th', 'Np', 'U ', 'Am', 'Pu', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
                           'Md', 'No', 'Rf', 'Lr', 'Db', 'Bh', 'Sg', 'Mt', 'Rg', 'Hs'/))

end module

program band

  call initial()
  call calculatepath()

end program

! initial information of supercell (Natom, Ncell, etc.)
subroutine initial()

  use global

  integer :: ierror
  character*5 :: ch

  nktotal = 0
  open(11, file = 'input')
  do
    do
      read(11, '(a5)', iostat = ierror) ch
      if (ch == '#CELL' .or. ch == '#PATH' .or. ierror /= 0) exit
    end do
    if (ch == '#CELL') then
      read(11, *) Natom, Ncell(1), Ncell(2), Ncell(3)
    else if (ch == '#PATH') then
      read(11, *) npath
      allocate(nkgrid(npath))
      allocate(kpath(npath, 2, 3))
      do i = 1, npath
        read(11, *) kpath(i, 1, :), kpath(i, 2, :), nkgrid(i)
      end do
    else
      exit
    end if
  end do
  close(11)

!write(*, *) Ncell
!write(*, *) nkgrid
!write(*, *) nktotal
!write(*, *) npath
!do i = 1, npath
!write(*, '(6f7.3)') kpath(i, 1, :), kpath(i, 2, :)
!end do
!stop

  Ntotal = Natom * Ncell(1) * Ncell(2) * Ncell(3)

!  allocate(molarmass(Natom))
!  allocate(chatom(Ntotal))

end subroutine

subroutine calculatepath()

  use global

! index of each atom in supercell for reorganization of Hessian matrix:
! 1, 2, 3 for cell position, 4 for index in unit cell
  integer :: position(Ntotal, 4)
! number of Hessian elements in hessian output file
  integer :: Nperline
  parameter (Nperline = 4)
! position of reference cell
  integer :: refer(3)!, t_position(3)
  integer :: index1, index2, index3, index4, ios
  integer :: ncount, nkcount
  integer :: index_order(Natom)
  double precision :: getmass, r
  double precision :: fxyz(Ntotal, 3)
  double precision :: minr(Natom, Ntotal)
  double precision :: fxyz1(3), fxyz2(3)
  double precision :: t_hessian(3)
  double precision :: V(3 * Natom, 3 * Ntotal)
  double precision :: molarmass(Natom)
  double precision :: pi
  parameter (pi = 3.1415926535898)
  double precision :: kpoint(2, 3)
  double complex :: phase
  double complex :: freq(3 * Natom)
!  double complex :: leigenvector(3 * Natom, 3 * Natom), reigenvector(3 * Natom, 3 * Natom)
  double complex :: B(3 * Natom, 3 * Natom)
  character*3 :: chatom(Ntotal)
  character*20 :: ch

  do i = 1, 3
    refer(i) = int((Ncell(i)+1)/2)
  end do
!write(*, *) refer
!stop

!  position = 0
  open(11, file = 'mol.cif')
  read(11, *)
  read(11, *) ch(1:14), length_a
  read(11, *) ch(1:14), length_b
  read(11, *) ch(1:14), length_c
  read(11, *) ch(1:17), angle_alpha
  read(11, *) ch(1:16), angle_beta
  read(11, *) ch(1:17), angle_gamma
  angle_alpha = angle_alpha * pi / 180.0d0
  angle_beta = angle_beta * pi / 180.0d0
  angle_gamma = angle_gamma * pi / 180.0d0
!write(*, *) length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma
!stop
  do i = 1, 6
    read(11, *)
  end do
  ncount = 1
  do i1 = 1, Ncell(1); do i2 = 1, Ncell(2); do i3 = 1, Ncell(3)
    do i4 = 1, Natom
      read(11, *) chatom(ncount), fxyz(ncount, :)
      position(ncount, 1) = i1
      position(ncount, 2) = i2
      position(ncount, 3) = i3
      position(ncount, 4) = i4
      if(i1 * i2 * i3 == 1) molarmass(ncount) = getmass(chatom(i4))
      ncount = ncount + 1
    end do
  end do; end do; end do
  close(11)

!  do i = 1, Natom
!    write(*, *) i, molarmass(i)
!  end do
!stop
!  do i = 1, Ntotal
!    write(*, '(i4, 3f10.5)') i, fxyz(i, :)
!  end do
!stop

  index1 = (((refer(1) - 1) * Ncell(2) + refer(2) - 1) * Ncell(3) + refer(3) - 1) * Natom
  do i = 1 + index1, Natom + index1
    do j = 1, Ntotal
      call neighbour(fxyz(i, :), fxyz(j, :), minr(i - index1, j))
    end do
  end do

!do i = 1, Natom
!do j = 1, Ntotal
!write(*, '(f4.1)', advance = 'no') minr(i, j)
!end do
!write(*, *)
!end do
!stop

  open(11, file = 'hessian.out')
  do i1 = 1, Ntotal; do j1 = 1, 3
    if (position(i1, 1) == refer(1) .and. position(i1, 2) == refer(2) .and. position(i1, 3) == refer(3)) then
!write(*, *) i1, 'pass1'
      nread = 0
      do i2 = 1, Ntotal; do j2 = 1, 3
        nread = nread + 1
        index1 = j1 + (position(i1, 4) - 1) * 3
        index2 = j2 + (i2 - 1) * 3
!write(*, *) index1, index2
        read(11, '(f16.10)', advance = 'no') V(index1, index2)
        V(index1, index2) = V(index1, index2) * 26424605.20760d0 / &
                            sqrt(molarmass(position(i1, 4)) * molarmass(position(i2, 4)))
        if(nread == Nperline) then
          read(11, *)
          nread = 0
        end if
      end do; end do
      if (mod(3 * Ntotal, Nperline) /= 0) read(11, *)
    else
!write(*, *) i1, 'pass2'
      do i2 = 1, (Ntotal * 3 - 1) / Nperline + 1
        read(11, *)
      end do
    end if
  end do; end do
  close(11)

!stop
!  do i = 1, 6
!    write(*, '(6f20.5)') V(i, 1:6)
!  end do
!stop

  open(33, file = 'kpoints.dat')
  open(22, file = 'bandstructure.dat')
  kpoint = 0.0d0
  nkcount = 1
  write(33, '(i10, f20.5)') nkcount, 0.0d0
  write(33, '(i10, f20.5)') nkcount, 3500.0d0
  write(33, *)
  do i = 1, npath
    do j = 1, nkgrid(i)
!write(*, *) nkcount
      do k = 1, 3
        kpoint(2, k) = (j - 1.0d0) / (nkgrid(i) - 1.0d0) * (kpath(i, 2, k) - kpath(i, 1, k)) + kpath(i, 1, k)
      end do
      if(i > 1 .and. j == 1 .and. &
         kpoint(1, 1) == kpoint(2, 1) .and. &
         kpoint(1, 2) == kpoint(2, 2) .and. &
         kpoint(1, 3) == kpoint(2, 3)) cycle
      B = cmplx(0.0, 0.0)
      index1 = 0
      index4 = (((refer(1) - 1) * Ncell(2) + refer(2) - 1) * Ncell(3) + refer(3) - 1) * Natom
      do i1 = 1, Natom; do j1 = 1, 3
        index1 = index1 + 1
        index2 = 0
        fxyz1(:) = fxyz(i1 + index4, :)
        do i2 = 1, Ntotal; do j2 = 1, 3
!write(*, '(3f10.5)') fxyz1
          index2 = index2 + 1
          index3 = (position(i2, 4) - 1) * 3 + j2
!write(*, *) index1, index2, index3
!write(*, *) minr(i1, i2)
          ncount_atom = 0
          do k1 = 1, 3; do k2 = 1, 3; do k3 = 1, 3
            fxyz2(1) = fxyz(i2, 1) + 2.0d0 - k1
            fxyz2(2) = fxyz(i2, 2) + 2.0d0 - k2
            fxyz2(3) = fxyz(i2, 3) + 2.0d0 - k3
            if(dabs(r(fxyz1, fxyz2) - minr(i1, i2)) < 1.0e-6) then
              t_hessian(1) = position(i2, 1) - refer(1) + (2.0d0 - k1) * Ncell(1)
              t_hessian(2) = position(i2, 2) - refer(2) + (2.0d0 - k2) * Ncell(2)
              t_hessian(3) = position(i2, 3) - refer(3) + (2.0d0 - k3) * Ncell(3)
              ncount_atom = ncount_atom + 1
!write(*, '(3f10.5, f10.6)') fxyz2, r(fxyz1, fxyz2)
            end if
          end do; end do; end do 
!write(*, '(3i4)') 2 - k1, 2 - k2, 2 - k3
!write(*, '(2i4, 3f7.3)') i1 + index4, i2, t_hessian
!          if(ncount_atom == 1) then
          phase = cdexp(dcmplx(0.0d0, - 2.0d0 * pi * &
                           (t_hessian(1) * kpoint(2, 1) + &
                            t_hessian(2) * kpoint(2, 2) + &
                            t_hessian(3) * kpoint(2, 3))))
          B(index1, index3) = B(index1, index3) + phase * V(index1, index2)
!write(*, '(i3, f10.5)') ncount_atom, minr(i1, i2)
!if(ncount_atom > 2) then 
!write(*, '(3f20.5)') V((i1 - 1) * 3 + 1, ((i2 - 1) * 3 +1):((i2 - 1) * 3 +3))
!write(*, '(3f20.5)') V((i1 - 1) * 3 + 2, ((i2 - 1) * 3 +1):((i2 - 1) * 3 +3))
!write(*, '(3f20.5)') V((i1 - 1) * 3 + 3, ((i2 - 1) * 3 +1):((i2 - 1) * 3 +3))
!end if
!          end if
        end do; end do
      end do; end do
      write(22, '(i5)', advance = 'no') nkcount
      call diag(3 * Natom, B, freq)
      call selection_sort(Natom * 3, freq)
      do k = 1, 3 * Natom
        write(22, '(f14.6)', advance = 'no') real(freq(k))
      end do
!write(*, '(6f14.6)') (aimag(freq(k)), k = 1, 6)
      write(22, *)
      nkcount = nkcount + 1
      kpoint(1, :) = kpoint(2, :)
    end do
    write(33, '(i10, f20.5)') nkcount, 0.0d0
    write(33, '(i10, f20.5)') nkcount, 3500.0d0
    write(33, *)
  end do
  close(22)
  close(33)

end subroutine

subroutine neighbour(fxyz1, fxyz2, minr)

  double precision, intent(in) :: fxyz1(3), fxyz2(3)
  double precision, intent(out) :: minr
  double precision :: fxyz2_re(3)
  double precision :: r

!write(*, '(6f10.5)') fxyz1, fxyz2
!write(*, '(6f10.5)') fxyz1, fxyz2
!stop
  minr = r(fxyz1, fxyz2)
!write(*, *) minr
!stop
  do i1 = 1, 3; do i2 = 1, 3; do i3 = 1, 3
    fxyz2_re(1) = fxyz2(1) + 2.0d0 - i1
    fxyz2_re(2) = fxyz2(2) + 2.0d0 - i2
    fxyz2_re(3) = fxyz2(3) + 2.0d0 - i3
!write(*, *) 2 - i1, 2 - i2, 2 - i3, r(fxyz1, fxyz2_re)
    if(r(fxyz1, fxyz2_re) < minr) minr = r(fxyz1, fxyz2_re)
!write(*, *) minr
  end do; end do; end do
!write(*, *)
!stop

end subroutine

subroutine diag(n_M, M, lambda)!, l_eigenvector, r_eigenvector)

  integer, intent(in) :: n_M
  double complex, intent(in) :: M(n_M, n_M)
  double complex, intent(out) :: lambda(n_M)
! left eigenvector
  double complex :: l_eigenvector(n_M, n_M)
! right eigenvector
  double complex :: r_eigenvector(n_M, n_M)
  integer :: info, lwmax
  parameter (lwmax = 1000000)
  integer :: lwork
  double precision :: rwork(2 * n_M)
  double complex :: work(lwmax)


  lwork = -1
  call zgeev('V', 'V', n_M, M, n_M, lambda, l_eigenvector, n_M, r_eigenvector, n_M, work, lwork, rwork, info) 

  lwork = min(lwmax, int(work(1)))

  call zgeev('V', 'V', n_M, M, n_M, lambda, l_eigenvector, n_M, r_eigenvector, n_M, work, lwork, rwork, info)

  do i = 1, n_M
    lambda(i) = zsqrt(lambda(i))
  end do

end subroutine

subroutine selection_sort(n_a, a)

  integer, intent(in) :: n_a
  double complex :: a(n_a)
  integer :: i, minindex
  double precision :: tmpdp
  double precision :: tmpreal(n_a)
  double complex :: tempcmplx

  do i = 1, n_a
    tmpreal(i) = real(a(i))
  end do

  do i = 1, n_a - 1
    minindex = minloc(tmpreal(i:n_a), 1) + i - 1
    if (tmpreal(i) > tmpreal(minindex)) then
      tempdp = tmpreal(i)
      tmpreal(i) = tmpreal(minindex)
      tmpreal(minindex) = tempdp
      tempcmplx = a(i)
      a(i) = a(minindex)
      a(minindex) = tempcmplx
    end if
  end do

end subroutine

double precision function getmass(char_atom)

  use global

  character*3, intent(in) :: char_atom

  do i = 1, atom_number
    if(char_atom == atom_name(i)) exit
  end do

  getmass = atom_mass(i)

end function

double precision function r(fxyz1, fxyz2)

  use global

  double precision, intent(in) :: fxyz1(3), fxyz2(3)
  double precision :: cxyz(3)

  cxyz(1) = (fxyz2(1) - fxyz1(1)) * length_a
  cxyz(2) = (fxyz2(2) - fxyz1(2)) * length_b
  cxyz(3) = (fxyz2(3) - fxyz1(3)) * length_c

  r = dsqrt(cxyz(1) ** 2 + cxyz(2) ** 2 + cxyz(3) ** 2 + &
            2 * cxyz(1) * cxyz(2) * dcos(angle_gamma) + &
            2 * cxyz(1) * cxyz(3) * dcos(angle_beta) + &
            2 * cxyz(2) * cxyz(3) * dcos(angle_alpha))

end function
