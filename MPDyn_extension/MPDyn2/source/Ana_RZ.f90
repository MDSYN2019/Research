! ############################
! ## SUBROUTINE LIST 
! ## -- RZ 
! ## -- RZG 
! ## -- RD 
! ## -- RDdisk 
! ## -- RDdiskG 
! ## -- ElecDensProfile 
! ## -- chargedensity 
! ## -- RZCyl 
! ## -- RZRes 
! ############################


!######################################################################
!######################################################################


subroutine RZ

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : AtomName, ResidName

implicit none

integer :: ndmin, ndmax
integer :: i , j, k, NumF
real(8), dimension(:,:), allocatable :: Pz
integer, dimension(NumRZ) :: NumPart
integer, dimension(N) :: Ibelong
integer :: iz
real(8) :: zz, dd, InvDR
integer :: NumChar
character(len=4) :: AIName, RIName
character(len=4) :: AName
real(8), dimension(3) :: a, b
real(8) :: Area, AvArea

   InvDR = 1.d0 / DRgrid

   ndmax = nint( Xmax * InvDR )
   ndmin = nint( Xmin * InvDR )

   allocate( Pz(ndmin:ndmax,NumRZ) )

   Ibelong = 0
   NumPart = 0

   do j = 1 , NumRZ

     AIName = RZAtom(j)
     RIName = RZResi(j)

     NumChar = len( trim(AIName) )

     if( NumChar == 1 ) then

       do i = 1 , N
         if( ResidName(i) == RIName ) then
           Aname = AtomName(i)
           if( Aname(1:1) == trim(AIName) ) then
             Ibelong(i) = j
             NumPart(j) = NumPart(j) + 1
           end if
         end if
       end do

     else

       do i = 1 , N
         if( ResidName(i) == RIName ) then
           if( AtomName(i) == AIName ) then
             if(Ibelong(i) /= 0) then
               NumPart(Ibelong(i)) = NumPart(Ibelong(i)) - 1
             end if
             Ibelong(i) = j
             NumPart(j) = NumPart(j) + 1
           end if
         end if
       end do

     end if

   end do

   do i = 1 , NumRZ
     if(NumPart(i)<=0) then
       write(*,*) 'ERROR : check the defined name of atom',i
       call Finalize
     end if
   end do

   Pz = 0.d0
   NumF = 0
   AvArea = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

       a  = H(:,1)
       b  = H(:,2)

       call CArea(a,b,Area)

       AvArea = AvArea + Area

       do k = 1, N

         if(Ibelong(k)==0) cycle
         iz = nint(R(RZdir,k)*InvDR)
         if(iz<ndmin.or.iz>ndmax) cycle
         Pz(iz,Ibelong(k)) = Pz(iz,Ibelong(k)) + 1.d0

       end do

     end do

   end do

!   AvArea = AvArea / dble(NumF)

   do i = 1 , NumRZ

     open(18,file=RZfile(i),status='unknown')

     write(18,'(a)') '# Probability density profile '
     if(RZdir==1) then
       write(18,'(a)') '# x[A]   Rho[A^-3] '
     else if(RZdir==2) then
       write(18,'(a)') '# y[A]   Rho[A^-3] '
     else if(RZdir==3) then
       write(18,'(a)') '# z[A]   Rho[A^-3] '
     end if

     do iz = ndmin, ndmax
       zz = iz*DRgrid
       dd = Pz(iz,i)*InvDR/AvArea
       write(18,'(f9.2,e15.6)') zz, dd
     end do

     close(18)

   end do

end subroutine RZ


!######################################################################
!######################################################################


subroutine RZG

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : ResidName, Mass

implicit none

integer :: ndmin, ndmax
integer :: i , j, k, NumF
real(8), dimension(:,:), allocatable :: Pz
integer, dimension(NumRZ) :: NumPart, NumA
integer :: iz, ii, i0, kk, ll
real(8) :: zz, dd, InvDR
character(len=4) :: RIName
real(8), dimension(3) :: a, b
real(8) :: Area, AvArea, Mg, Rg

   InvDR = 1.d0 / DRgrid

   ndmax = nint( Xmax * InvDR )
   ndmin = nint( Xmin * InvDR )

   allocate( Pz(ndmin:ndmax,NumRZ) )

   NumPart = 0

   do j = 1 , NumRZ

     RIName = RZResi(j)

     ii = 1
p1:  do i = 1 , NumSpec
       if(ResidName(ii) == RIName) then
         NumPart(j) = i
         NumA(j) = ii-1
         exit p1
       end if
       if(i==NumSpec) then
         write(*,*) "error : selected molecule is not found!"
         write(*,*) RIName
         call Finalize
       end if
       ii = ii + NumMol(i)*NumAtm(i)
     end do p1

   end do

   Pz = 0.d0
   NumF = 0
   AvArea = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

       a  = H(:,1)
       b  = H(:,2)

       call CArea(a,b,Area)

       AvArea = AvArea + Area

       do k = 1, NumRZ

         ii = NumA(k)
         i0 = NumPart(k)

         do kk = 1, NumMol(i0)
           Rg = 0.d0
           Mg = 0.d0
           do ll = 1, NumAtm(i0)
             ii = ii + 1
             Rg = Rg + R(RZdir,ii) * Mass(ii)
             Mg = Mg + Mass(ii)
           end do
           Rg = Rg / Mg
           iz = nint(Rg*InvDR)
           if(iz<ndmin.or.iz>ndmax) cycle
           Pz(iz,k) = Pz(iz,k) + 1.d0
         end do

       end do

     end do

   end do

   do i = 1 , NumRZ

     open(18,file=RZfile(i),status='unknown')

     write(18,'(a)') '# Probability density profile '
     if(RZdir==1) then
       write(18,'(a)') '# x[A]   Rho[A^-3] '
     else if(RZdir==2) then
       write(18,'(a)') '# y[A]   Rho[A^-3] '
     else if(RZdir==3) then
       write(18,'(a)') '# z[A]   Rho[A^-3] '
     end if

     do iz = ndmin, ndmax
       zz = iz*DRgrid
       dd = Pz(iz,i)*InvDR/AvArea
       write(18,'(f9.2,e15.6)') zz, dd
     end do

     close(18)

   end do

end subroutine RZG


!######################################################################
!######################################################################


subroutine RD

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : AtomName, ResidName, Mass
use UnitExParam, only : pi

implicit none

integer :: ndmin, ndmax
integer :: i, j, k, ii, NumF, Numa, Numb
real(8), dimension(:,:), allocatable :: Pz
integer, dimension(NumRZ) :: NumPart
integer, dimension(N) :: Ibelong
integer :: iz, ngridpoint
real(8) :: dd, InvDR, Xmax2, Xmin2, R2, R1
integer :: NumChar
character(len=4) :: AIName, RIName
character(len=4) :: AName
real(8) :: rl, rs, vl, InvMasG
real(8), dimension(3) :: Rg

   InvDR = 1.d0 / DRgrid

   if(Xmin < 0.) then
     write(*,*) "RANGE should be positive values"
     write(*,*) "Minimum is set to 0"
     Xmin = 0.d0
   end if
   if(Xmax <= 0.) then
     write(*,*) "Error: RANGE should be positive values"
     stop
   end if

   ndmax = nint( Xmax * InvDR )
   ndmin = nint( Xmin * InvDR )

   Xmax2 = Xmax ** 2
   Xmin2 = Xmin ** 2

   ngridpoint = int( (Xmax - Xmin)*InvDR )

   allocate( Pz(ngridpoint,NumRZ) )

   Ibelong = 0
   NumPart = 0

   do j = 1 , NumRZ

     AIName = RZAtom(j)
     RIName = RZResi(j)

     NumChar = len( trim(AIName) )

     if( NumChar == 1 ) then

       do i = 1 , N
         if( ResidName(i) == RIName ) then
           Aname = AtomName(i)
           if( Aname(1:1) == trim(AIName) ) then
             Ibelong(i) = j
             NumPart(j) = NumPart(j) + 1
           end if
         end if
       end do

     else

       do i = 1 , N
         if( ResidName(i) == RIName ) then
           if( AtomName(i) == AIName ) then
             if(Ibelong(i) /= 0) then
               NumPart(Ibelong(i)) = NumPart(Ibelong(i)) - 1
             end if
             Ibelong(i) = j
             NumPart(j) = NumPart(j) + 1
           end if
         end if
       end do

     end if

   end do

   do i = 1 , NumRZ
     if(NumPart(i)==0) then
       write(*,*) 'ERROR : check the defined name of atom'
       call Finalize
     end if
   end do

   Numb = 0
   ii = 0
   do i = 1, NumSpec
     if(i==NComp) then
       Numa = ii+1
       Numb = ii+NumMol(i)*NumAtm(i)
       exit
     end if
     ii = ii + NumMol(i)*NumAtm(i)
   end do
   InvMasG = 0.d0
   if( Numb /= 0 ) then
     do i = Numa, Numb
       InvMasG = InvMasG + Mass(i)
     end do
     InvMasG = 1.d0 / InvMasG
   end if

   Pz(:,:) = 0.d0
   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

! ## correction of COM when needed
       if( Numb /= 0 ) then
         Rg(:) = 0.d0
         do k = Numa, Numb
           Rg(:) = Rg(:) + Mass(k)*R(:,k)
         end do
         Rg(:) = Rg(:) * InvMasG
         do k = 1, N
           R(:,k) = R(:,k) - Rg(:)
         end do
         call PBC
         Rg(:) = 0.d0
         do k = Numa, Numb
           Rg(:) = Rg(:) + Mass(k)*R(:,k)
         end do
         Rg(:) = Rg(:) * InvMasG
         do k = 1, N
           R(:,k) = R(:,k) - Rg(:)
         end do
         call PBC
       end if

       do k = 1, N

         if(Ibelong(k)==0) cycle
         R2 = dot_product(R(:,k),R(:,k))
         if((R2 > Xmax2).or.(R2 < Xmin2)) cycle
         R1 = sqrt(R2) - Xmin
         iz = int(R1*InvDR) + 1
         Pz(iz,Ibelong(k)) = Pz(iz,Ibelong(k)) + 1.d0

       end do

     end do

   end do

   do i = 1 , NumRZ

     open(18,file=RZfile(i),status='unknown')

     write(18,'(a)') '# Probability density profile '
     write(18,'(a)') '# r[A]   Rho[A^-3] '

     do iz = 1, ngridpoint
       rl = Xmin + iz*DRgrid
       rs = rl - DRgrid
       vl = 4.d0 * pi / 3.d0 * (rl ** 3 - rs ** 3)
       dd = Pz(iz,i)/(vl*NumF)
       write(18,'(f9.2,e15.6)') rs+0.5*DRgrid, dd
     end do

     close(18)

   end do

end subroutine RD


!######################################################################
!######################################################################


subroutine RDdisk

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : AtomName, ResidName
use UnitExParam, only : pi

implicit none

integer :: ndmin, ndmax
integer :: i , j, k, NumF
real(8), dimension(:,:), allocatable :: Pz
integer, dimension(NumRZ) :: NumPart
integer, dimension(N) :: Ibelong
integer :: iz, ix, iy, ii
real(8) :: dd, InvDR, Xmax2, Xmin2, R2, R1
integer :: NumChar, ngridpoint
character(len=4) :: AIName, RIName
character(len=4) :: AName
real(8) :: rl, rs, vl, aveH

   InvDR = 1.d0 / DRgrid

   if(Xmin < 0.) then
     write(*,*) "RANGE should be positive values"
     write(*,*) "Minimum is set to 0"
     Xmin = 0.d0
   end if
   if(Xmax <= 0.) then
     write(*,*) "Error: RANGE should be positive values"
     stop
   end if

   ndmax = nint( Xmax * InvDR )
   ndmin = nint( Xmin * InvDR )

   Xmax2 = Xmax ** 2
   Xmin2 = Xmin ** 2

   ngridpoint = int( (Xmax - Xmin)*InvDR )

   allocate( Pz(ngridpoint,NumRZ) )

   Ibelong = 0
   NumPart = 0

   do j = 1 , NumRZ

     AIName = RZAtom(j)
     RIName = RZResi(j)

     NumChar = len( trim(AIName) )

     if( NumChar == 1 ) then

       do i = 1 , N
         if( ResidName(i) == RIName ) then
           Aname = AtomName(i)
           if( Aname(1:1) == trim(AIName) ) then
             Ibelong(i) = j
             NumPart(j) = NumPart(j) + 1
           end if
         end if
       end do

     else

       do i = 1 , N
         if( ResidName(i) == RIName ) then
           if( AtomName(i) == AIName ) then
             if(Ibelong(i) /= 0) then
               NumPart(Ibelong(i)) = NumPart(Ibelong(i)) - 1
             end if
             Ibelong(i) = j
             NumPart(j) = NumPart(j) + 1
           end if
         end if
       end do

     end if

   end do

   do i = 1 , NumRZ
     if(NumPart(i)==0) then
       write(*,*) 'ERROR : check the defined name of atom'
       call Finalize
     end if
   end do

   ii = 0
   do i = 1, 3
     if(i==RZdir) cycle
     ii = ii + 1
     if(ii==1) then
       ix = i
     else if(ii==2) then
       iy = i
     end if
   end do


   Pz(:,:) = 0.d0
   NumF = 0
   aveH = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------
       call CellTransform

       aveH = aveH + H(3,3)

       do k = 1, N

         if(Ibelong(k)==0) cycle
         R2 = R(ix,k)*R(ix,k) + R(iy,k)*R(iy,k)
         if((R2 > Xmax2).or.(R2 < Xmin2)) cycle
         R1 = sqrt(R2) - Xmin
         iz = int(R1*InvDR) + 1
         Pz(iz,Ibelong(k)) = Pz(iz,Ibelong(k)) + 1.d0

       end do

     end do

   end do

   do i = 1 , NumRZ

     open(18,file=RZfile(i),status='unknown')

     write(18,'(a)') '# Probability density profile in 2D'
     if(RZdir==1) then
       write(18,'(a)') '# Ryz[A]   Rho[A^-3] '
     else if(RZdir==2) then
       write(18,'(a)') '# Rxz[A]   Rho[A^-3] '
     else if(RZdir==3) then
       write(18,'(a)') '# Rxy[A]   Rho[A^-3] '
     end if

     do iz = 1, ngridpoint
       rl = Xmin + iz*DRgrid
       rs = rl - DRgrid
       vl = 2.d0 * pi * (rl ** 2 - rs ** 2)
       dd = Pz(iz,i)/(vl*aveH)
       write(18,'(f9.2,e15.6)') rs+0.5*DRgrid, dd
     end do

     close(18)

   end do

end subroutine RDdisk


!######################################################################
!######################################################################


subroutine RDdiskG

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H
use AtomParam, only : ResidName, Mass
use UnitExParam, only : pi

implicit none

integer :: ndmin, ndmax
integer :: i , j, k
real(8), dimension(:,:), allocatable :: Pz
integer, dimension(NumRZ) :: NumPart, NumA
integer :: iz, ii, i0, kk, ll, ix, iy, ngridpoint
real(8) :: dd, InvDR, R1, R2
character(len=4) :: RIName
real(8) :: aveH, Mg
real(8) :: Rgx, Rgy, Xmax2, Xmin2, rs, rl, vl

   InvDR = 1.d0 / DRgrid

   if(Xmin < 0.) then
     write(*,*) "RANGE should be positive values"
     write(*,*) "Minimum is set to 0"
     Xmin = 0.d0
   end if
   if(Xmax <= 0.) then
     write(*,*) "Error: RANGE should be positive values"
     stop
   end if

   ndmax = nint( Xmax * InvDR )
   ndmin = nint( Xmin * InvDR )

   Xmax2 = Xmax ** 2
   Xmin2 = Xmin ** 2

   ngridpoint = int( (Xmax - Xmin)*InvDR )

   allocate( Pz(ngridpoint,NumRZ) )

   NumPart = 0

   do j = 1 , NumRZ

     RIName = RZResi(j)

     ii = 1
p1:  do i = 1 , NumSpec
       if(ResidName(ii) == RIName) then
         NumPart(j) = i
         NumA(j) = ii-1
         exit p1
       end if
       if(i==NumSpec) then
         write(*,*) "error : selected molecule is not found!"
         write(*,*) RIName
         call Finalize
       end if
       ii = ii + NumMol(i)*NumAtm(i)
     end do p1

   end do

   ii = 0
   do i = 1, 3
     if(i==RZdir) cycle
     ii = ii + 1
     if(ii==1) then
       ix = i
     else if(ii==2) then
       iy = i
     end if
   end do

   Pz = 0.d0
   aveH = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------
       call CellTransform

       aveH = aveH + H(3,3)

       do k = 1, NumRZ

         ii = NumA(k)
         i0 = NumPart(k)

         do kk = 1, NumMol(i0)
           Rgx = 0.d0
           Rgy = 0.d0
           Mg = 0.d0
           do ll = 1, NumAtm(i0)
             ii = ii + 1
             Rgx = Rgx + R(ix,ii) * Mass(ii)
             Rgy = Rgy + R(iy,ii) * Mass(ii)
             Mg = Mg + Mass(ii)
           end do
           Rgx = Rgx / Mg
           Rgy = Rgy / Mg
           R2 = Rgx*Rgx + Rgy*Rgy
           if((R2 > Xmax2).or.(R2 < Xmin2)) cycle
           R1 = sqrt(R2) - Xmin
           iz = int(R1*InvDR) + 1
           Pz(iz,k) = Pz(iz,k) + 1.d0
         end do

       end do

     end do

   end do

   do i = 1 , NumRZ

     open(18,file=RZfile(i),status='unknown')

     write(18,'(a)') '# Probability density profile in 2D'
     if(RZdir==1) then
       write(18,'(a)') '# Ryz[A]   Rho[A^-3] '
     else if(RZdir==2) then
       write(18,'(a)') '# Rxz[A]   Rho[A^-3] '
     else if(RZdir==3) then
       write(18,'(a)') '# Rxy[A]   Rho[A^-3] '
     end if

     do iz = 1, ngridpoint
       rl = Xmin + iz*DRgrid
       rs = rl - DRgrid
       vl = 2.d0 * pi * (rl ** 2 - rs ** 2)
       dd = Pz(iz,i)/(vl*aveH)
       write(18,'(f9.2,e15.6)') rs+0.5*DRgrid, dd
     end do

     close(18)

   end do

end subroutine RDdiskG


!######################################################################
!######################################################################


subroutine ElecDensProfile

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H

implicit none

integer :: ndmin, ndmax
integer :: i, j, k, NumF
real(8), dimension(:), allocatable :: Ed, Es
integer :: iz
real(8) :: zz, dd, InvDR
real(8), dimension(3) :: a, b
real(8) :: Area

   InvDR = 1.d0 / DRgrid

   ndmax = nint( Xmax * InvDR )
   ndmin = nint( Xmin * InvDR )

   allocate( Ed(ndmin:ndmax) )
   allocate( Es(ndmin:ndmax) )

   call chargedensity

   Ed = 0.d0
   NumF = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       Es(:) = 0.d0

       NumF = NumF + 1
!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

       a  = H(:,1)
       b  = H(:,2)

       call CArea(a,b,Area)

       do k = 1, N

         iz = nint(R(RZdir,k)*InvDR)
         if(iz<ndmin.or.iz>ndmax) cycle
         Es(iz) = Es(iz) + elecnum(k)

       end do

       Area = InvDR / Area

       Ed(:) = Ed(:) + Es(:) * Area

     end do

   end do

   open(18,file='./Analy/ElectronDensity.data',status='unknown')

     write(18,'(a)') '# Electron density profile '

     if(RZdir==1) then
       write(18,'(a)') '# x[A]  Rho[e/A^3] '
     else if(RZdir==2) then
       write(18,'(a)') '# y[A]  Rho[e/A^3] '
     else if(RZdir==3) then
       write(18,'(a)') '# z[A]  Rho[e/A^3] '
     end if

     do iz = ndmin, ndmax

       zz = iz*DRgrid
       dd = dble(Ed(iz))/dble(NumF)

       write(18,'(f6.2,f10.5)') zz, dd

     end do

   close(18)

end subroutine ElecDensProfile


!######################################################################
!######################################################################


subroutine chargedensity

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : ec
use NonbondParam, only : Charge
use AtomParam, only : AtomName

implicit none

character(len=4) :: Aname
integer :: i

   allocate( elecnum(N) )

   do i = 1 , N

     Aname = AtomName(i)

     if(Aname(1:1) == 'C') then
       elecnum(i) = 6.d0
     else if(Aname(1:1) == 'H') then
       elecnum(i) = 1.d0
     else if(Aname(1:1) == 'O') then
       elecnum(i) = 8.d0
     else if(Aname(1:1) == 'N') then
       elecnum(i) = 7.d0
     else if(Aname(1:1) == 'P') then
       elecnum(i) = 15.d0
     else if(Aname == 'SOD') then
       elecnum(i) = 11.d0
     else if(Aname == 'CLA') then
       elecnum(i) = 17.d0
     else
       write(*,*) 'ERROR : missing parameter ',Aname
       write(*,*) 'add your parameter for atom,', Aname
       write(*,*) 'in subroutine chargedensity in Ana_RZ.f90'
       call Finalize
     end if

   end do

   do i = 1 , N

     elecnum(i) = elecnum(i) + Charge(i) / sqrt(ec)

   end do

end subroutine chargedensity


!######################################################################
!######################################################################


subroutine RZCyl

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : AtomName, ResidName

implicit none

integer, parameter :: ndz = 400
integer :: i, j, k, TotalSteps
integer, dimension(-ndz:ndz,NumRZ) :: RzIdentIC
integer, dimension(-ndz:ndz,NumRZ) :: RgIdentIC
integer, dimension(-ndz:ndz,NumRZ) :: RzIdentOC
integer, dimension(-ndz:ndz,NumRZ) :: RgIdentOC
real(8), dimension(3) :: Rg, Rij
integer :: iz, ig
real(8) :: zz, R2
integer, dimension(NumRZ) :: NinC
integer, dimension(NumRZ) :: NumPart
integer, dimension(N) :: Ibelong
real(8) :: rz1, rz2

   Ibelong = 0
   NumPart = 0

   do i = 1 , N

     do j = 1 , NumRZ

       if((AtomName(i) == RZAtom(j)).and.(ResidName(i) == RZResi(j))) then

         Ibelong(i) = j
         NumPart(j) = NumPart(j) + 1

       end if

     end do

   end do

   do i = 1 , NumRZ

     if(NumPart(i)==0) then

       write(*,*) 'ERROR : check the defined name of atom'
       call Finalize

     end if

   end do

   RzIdentIC = 0
   RzIdentOC = 0
   RgIdentIC = 0
   RgIdentOC = 0

   NinC = 0

   TotalSteps = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     TotalSteps = TotalSteps + NTrjStep(i)

     do j = 1 , NTrjStep(i)

!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

       call COM_SpecComponent(Rg,Kcomp)

       do k = 1, N

         if(Ibelong(k)==0) cycle

         Rij = R(:,k) - Rg
         R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2)

         iz = nint( R(3,k) * 10.d0)
         ig = nint((R(3,k)-Rg(3)) * 10.d0)

         if(R2 < RadCyl2) then
           if(abs(iz)<400) RzIdentIC(iz,Ibelong(k)) = RzIdentIC(iz,Ibelong(k)) + 1
           if(abs(ig)<400) RgIdentIC(ig,Ibelong(k)) = RgIdentIC(ig,Ibelong(k)) + 1
           NinC(Ibelong(k)) = NinC(Ibelong(k)) + 1
         else
           if(abs(iz)<400) RzIdentOC(iz,Ibelong(k)) = RzIdentOC(iz,Ibelong(k)) + 1
           if(abs(ig)<400) RgIdentOC(ig,Ibelong(k)) = RgIdentOC(ig,Ibelong(k)) + 1
         end if

       end do

     end do

   end do

   do i = 1 , NumRZ

     open(41,file=RZfileI(i),status='unknown')
     open(42,file=RZfileO(i),status='unknown')

     do j = -ndz,ndz

       zz = dble(j)*0.1d0

       rz1 = dble(RzIdentIC(j,i)) / dble(NinC(i)) * 1.d+02
       rz2 = dble(RgIdentIC(j,i)) / dble(NinC(i)) * 1.d+02
       write(41,'(f6.1,2f10.6)') zz, rz1, rz2

       rz1 = dble(RzIdentOC(j,i)) / dble(TotalSteps*NumPart(i)-NinC(i))*1.d+2
       rz2 = dble(RgIdentOC(j,i)) / dble(TotalSteps*NumPart(i)-NinC(i))*1.d+2
       write(42,'(f6.1,2f10.6)') zz, rz1, rz2

     end do

     close(41)
     close(42)

   end do

end subroutine RZCyl


!######################################################################
!######################################################################


subroutine RZRes

use Numbers, only : N
use Configuration, only : R
use ParamAnalyze
use AtomParam , only : ResidNum, Mass

implicit none

integer, parameter :: ndz = 400
integer :: i , j, k, l
character(len=72) :: Ch
real(8), dimension(-ndz:ndz,300) :: RzIdent
real(8), dimension(-ndz:ndz,300) :: RgIdent
real(8), dimension(3) :: Rg
real(8) :: SumM
integer :: iz
real(8) :: Sekz, Sekg, zz

   RzIdent = 0.
   RgIdent = 0.

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call CellTransform

! ## CHANGE the origin of axes

       if(cWhat=='RZResG') then

         call COM_SpecComponent(Rg,Kcomp)

         do k = 1 , N

           R(:,k) = R(:,k) - Rg(:)

         end do

       end if
! ####

       do l = 1 , NumResRZ

         Rg   = 0.
         SumM = 0.

         do k = 1 , N

           if(ResidNum(k) == ResRZ(l)) then

             Rg   = Rg   + Mass(k) * R(:,k)
             SumM = SumM + Mass(k)
             iz   = nint( R(3,k) * 10.d0 )
             if(abs(iz)<400) RzIdent(iz,l) = RzIdent(iz,l) + Mass(k)

           end if

         end do

         Rg = Rg / SumM
         iz = nint(Rg(3)*10.d0)
         if(abs(iz)<400) RgIdent(iz,l) = RgIdent(iz,l) + 1.d0

       end do

     end do

   end do

   do i = 1 , NumResRZ

     if(cWhat=='RZRes') then

       if(ResRZ(i)<10) then

         write(Ch,'(a,i1,a)') './Analy/RZ000',ResRZ(i),'.dat'

       else if(ResRZ(i)<100) then

         write(Ch,'(a,i2,a)') './Analy/RZ00' ,ResRZ(i),'.dat'

       else if(ResRZ(i)<1000) then

         write(Ch,'(a,i3,a)') './Analy/RZ0'  ,ResRZ(i),'.dat'

       else if(ResRZ(i)<10000) then

         write(Ch,'(a,i4,a)') './Analy/RZ'   ,ResRZ(i),'.dat'

       else

         write(*,*) 'too large Residue number'
         call Finalize

       end if

     else if(cWhat=='RZResG') then

        if(ResRZ(i)<10) then

         write(Ch,'(a,i1,a)') './Analy/RZG000',ResRZ(i),'.dat'

       else if(ResRZ(i)<100) then

         write(Ch,'(a,i2,a)') './Analy/RZG00' ,ResRZ(i),'.dat'

       else if(ResRZ(i)<1000) then

         write(Ch,'(a,i3,a)') './Analy/RZG0'  ,ResRZ(i),'.dat'

       else if(ResRZ(i)<10000) then

         write(Ch,'(a,i4,a)') './Analy/RZG'   ,ResRZ(i),'.dat'

       else

         write(*,*) 'too large Residue number'
         call Finalize

       end if

     end if

     open(10,file=Ch,status='unknown')

     Sekz = 0.d0
     Sekg = 0.d0

     do j = -ndz,ndz

       Sekz = Sekz + RzIdent(j,i)
       Sekg = Sekg + RgIdent(j,i)

     end do

     do j = -ndz,ndz

       zz = dble(j)*0.1d0
       write(10,'(f6.1,2f10.6)') zz, RzIdent(j,i)/Sekz*1.d2, RgIdent(j,i)/Sekg*1.d2

     end do

     close(10)

   end do

end subroutine RZRes
