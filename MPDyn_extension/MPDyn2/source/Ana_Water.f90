! ############################
! ## SUBROUTINE LIST 
! ## -- PoreWaterTrajectroy 
! ## -- PoreWaterDensity 
! ## -- OrientationOfWater 
! ## -- CylinderWater 
! ## -- CountWaterInPore 
! ## -- WaterCOM 
! ## -- Water_Diffusion_Zaxis 
! ## -- Water_Exchange_Correlation 
! ############################


!######################################################################
!######################################################################


subroutine PoreWaterTrajectroy

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use AtomParam, only : ResidName
use TimeParam, only : Timeps

implicit none

integer, parameter :: Steps = 10
integer :: i, j, k, l, Count, ii
real(8), dimension(3) :: Rg, ScRg
real(8), dimension(:,:), allocatable :: RgL

integer, parameter :: ndz = 400
real(8) :: RzProt, zz
integer, dimension(-ndz:ndz) :: IzProt
integer :: TotalStep, ISteps, iz
integer :: Numa, Numb, Wcomp, Numwa, Numwb
integer, dimension(:), allocatable :: ListN
real(8), dimension(:,:), allocatable :: Rw
integer :: NumSearch
real(8), dimension(:,:,:), allocatable :: CurrRwg

   open(1,file='./Analy/RzProteinG.dat',status='unknown')
   open(2,file='./Analy/NumWatPore.dat',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   ii = 1
   do i = 1 , NumSpec

     if( ResidName(ii) == 'TIP3' ) then
       Wcomp = i
       Numwa = ii
       Numwb = ii + NumMol(i)*NumAtm(i) - 1
     end if

     ii = ii + NumMol(i)*NumAtm(i)

   end do

   allocate( PoreWater(NumMol(Wcomp)) )
   allocate( Rw(3,NumMol(Wcomp)) )

   PoreWater = 0

   ISteps = 0
   IzProt = 0

   TotalStep = 0

   do i = 1 , NJobs

     TotalStep = TotalStep + NTrjStep(i)

   end do

   allocate(RgL(3,TotalStep/Steps))

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

       ISteps = ISteps + 1

!     -------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     -------------------

       call COM_SpecComponent(Rg,Kcomp)

       call InversMatrix(H,InvH)

       ScRg = matmul( InvH, Rg )

       call CellTransform

       Rg = matmul( H, ScRg )

       iz = nint(Rg(3)*10.)
       IzProt(iz) = IzProt(iz) + 1

       if(mod(ISteps,Steps)==0) RgL(:,ISteps/Steps) = Rg

       call CountWaterInPore( ScRg,Count,Wcomp,Numwa )

       write(2,'(f9.4,i8)') Timeps, Count

     end do

   end do

   close(2)

   NumSearch = 0

   do i = 1 , NumMol(Numwa)

     if(PoreWater(i)/=0) then

       NumSearch = NumSearch + 1
       PoreWater(i) = NumSearch

     end if

   end do

   allocate( ResultFile(NumSearch) )
   allocate( ListN(NumSearch) )
   allocate( CurrRwg(3,NumSearch,TotalStep/Steps) )

   NumSearch = 0

   do i = 1 , NumMol(Numwa)

     if(PoreWater(i)/=0) then

       NumSearch = NumSearch + 1
       ListN(NumSearch) = i

     end if

   end do

   do i = 1 , NumSearch

     j = ListN(i)

     if(j<10) then

       write(ResultFile(i),'(a,i1,a)') 'WatTraj000',j,'.dat'

     else if(j<100) then

       write(ResultFile(i),'(a,i2,a)') 'WatTraj00' ,j,'.dat'

     else if(j<1000) then

       write(ResultFile(i),'(a,i3,a)') 'WatTraj0'  ,j,'.dat'

     else if(j<10000) then

       write(ResultFile(i),'(a,i4,a)') 'WatTraj'   ,j,'.dat'

     else

       write(*,*) 'too large particles'
       call Finalize

     end if

   end do


   ISteps = 0

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

       if(mod(j,Steps)==0) then

         ISteps = ISteps + 1

         call CellTransform

         do k = Numwa , Numwb
           R(:,k) = R(:,k) - RgL(:,ISteps)
         end do

         call WaterCOM(Rw,Wcomp,Numwa)

         do k = 1 , NumSearch

           l = ListN(k)
           CurrRwg(:,k,ISteps) = Rw(:,l)

         end do

       end if

     end do

   end do

   do i = 1 , NumSearch

     open(18,file=ResultFile(i),status='unknown')

     do j = 1, ISteps

!       if(j/=1) then
!         Rij = CurrRwg(:,i,j) - CurrRwg(:,i,j-1)
!         if((Rij(1) > 10.).or.(Rij(2) > 10.).or.(Rij(3) > 10.)) then
!           write(18,'(a)') ' = = = '
!           write(18,'(3f8.3)') CurrRwg(:,i,j)
           write(18,'(a1,3f8.3)') 'O',CurrRwg(:,i,j)
!       else
!           write(18,'(3f8.3)') CurrRwg(:,i,j)
!       end if

     end do

     close(18)

   end do

   do i = -400, 400

     zz = dble(i) * 0.1
     RzProt = dble(IzProt(i)) / dble(TotalStep)
     write(1,'(f5.1,f10.5)') zz,RzProt*1.d+2

   end do

close(1)

end subroutine PoreWaterTrajectroy


!######################################################################
!######################################################################


subroutine PoreWaterDensity

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use AtomParam, only : ResidName

implicit none

integer :: i, j, k, l, ii
real(8), dimension(3) :: Rg

real(8) :: zo, zh
integer :: TotalStep, ix, iy, iz
integer :: Numa, Numb, Numwa, Numwb
integer, dimension(:,:,:), allocatable :: NumGridO
integer, dimension(:,:,:), allocatable :: NumGridH
integer, dimension(3) :: MaxGrid

   open(1,file='./Analy/GridDensityWaterO.dat',status='unknown')
   open(2,file='./Analy/GridDensityWaterH.dat',status='unknown')

   Numa = 0
   Numb = NumMol(1) * NumAtm(1)

   if(Kcomp /= 1) then

     do i = 2, Kcomp

       Numa = Numb
       Numb = Numb + NumMol(i)*NumAtm(i)

     end do

   end if

   Numa = Numa + 1

   ii = 1
   do i = 1 , NumSpec

     if( ResidName(ii) == 'TIP3' ) then
       Numwa = ii
       Numwb = ii + NumMol(i)*NumAtm(i) - 1
     end if

     ii = ii + NumMol(i)*NumAtm(i)

   end do

   MaxGrid(1) = nint(RadCyl)
   MaxGrid(2) = nint(RadCyl)
   MaxGrid(3) = nint(Zsh)

   allocate( NumGridO( -MaxGrid(1) : MaxGrid(1), &
   &                   -MaxGrid(2) : MaxGrid(2), &
   &                   -MaxGrid(3) : MaxGrid(3) ) )
   allocate( NumGridH( -MaxGrid(1) : MaxGrid(1), &
   &                   -MaxGrid(2) : MaxGrid(2), &
   &                   -MaxGrid(3) : MaxGrid(3) ) )

   NumGridO = 0
   NumGridH = 0

   TotalStep = 0

   do i = 1 , NJobs

     TotalStep = TotalStep + NTrjStep(i)

   end do

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

       call COM_SpecComponent(Rg,Kcomp)

       do k = 1 , N

         R(:,k) = R(:,k) - Rg

       end do

       call PBC

       do k = Numwa, Numwb

         l = mod((k-Numa+1),3)

         ix = nint( R(1,k) )
         iy = nint( R(2,k) )
         iz = nint( R(3,k) )

         if((abs(ix)>MaxGrid(1)).or.(abs(iy)>MaxGrid(2))&
         &                      .or.(abs(iz)>MaxGrid(3))) cycle

         if(l==1) then
           NumGridO(ix,iy,iz) = NumGridO(ix,iy,iz) + 1
         else
           NumGridH(ix,iy,iz) = NumGridH(ix,iy,iz) + 1
         end if

       end do

     end do

   end do

   write(1,'(a,/,a/)') &
   & '! ## the grid density of the oxygen atom of water molecules',&
   & '! ## x, y, z, density(%) , form=(3i4,f12.5) '

   write(2,'(a,/,a/)') &
   & '! ## the grid density of the hydrogen atoms of water molecules',&
   & '! ## x, y, z, density(%) , form=(3i4,f12.5) '

   do i = -MaxGrid(1) , MaxGrid(1)

     do j = -MaxGrid(2) , MaxGrid(2)

       do k = -MaxGrid(3) , MaxGrid(3)

         zo = dble(NumGridO(i,j,k)) / dble(TotalStep) * 1.d2
         zh = dble(NumGridH(i,j,k)) / dble(TotalStep) * 1.d2

         write(1,'(3i4,f12.5)') i,j,k,zo
         write(2,'(3i4,f12.5)') i,j,k,zh

       end do

     end do

   end do

close(1)
close(2)

end subroutine PoreWaterDensity


!######################################################################
!######################################################################


subroutine OrientationOfWater

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use UnitExParam, only : pi
use CellParam, only : H, InvH
use AtomParam, only : ResidName

implicit none

integer :: i, j, ii
real(8), dimension(3) :: Rg, ScRg

real(8), parameter :: drz = 1.d0  ! [A]
real(8), dimension(-40:40) :: RhoInCyl, RhoOutCyl
integer :: TotalStep
real(8) :: SumArea, AveArea, CylArea
integer :: Numa, Numb, Wcomp
real(8) :: Area, zz, cs
real(8), dimension(3) :: a, b

   if(cWhat=='WatOrient') then

     open(42,file='./Analy/DensWater.dat',status='unknown')
     open(44,file='./Analy/OriWater.dat',status='unknown')

   else if(cWhat=='WatOrientG') then

     open(42,file='./Analy/DensWaterG.dat',status='unknown')
     open(44,file='./Analy/OriWaterG.dat',status='unknown')
 
   else if(cWhat=='WatOriCyl') then

     open(41,file='./Analy/DensWater_inCyl.dat',status='unknown')
     open(42,file='./Analy/DensWater_outCyl.dat',status='unknown')
     open(43,file='./Analy/OriWater_inCyl.dat',status='unknown')
     open(44,file='./Analy/OriWater_outCyl.dat',status='unknown')

   else if(cWhat=='WatOriCylG') then

     open(41,file='./Analy/DensWater_inCylG.dat',status='unknown')
     open(42,file='./Analy/DensWater_outCylG.dat',status='unknown')
     open(43,file='./Analy/OriWater_inCylG.dat',status='unknown')
     open(44,file='./Analy/OriWater_outCylG.dat',status='unknown')

   end if

   if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

     Numa = 0
     Numb = NumMol(1) * NumAtm(1)

     if(Kcomp /= 1) then

       do i = 2, Kcomp

         Numa = Numb
         Numb = Numb + NumMol(i)*NumAtm(i)

       end do

     end if

     Numa = Numa + 1

   end if

   ii = 1

   do i = 1 , NumSpec

     if( ResidName(ii) == 'TIP3' ) then
       Wcomp = i
     end if

     ii = ii + NumMol(i)*NumAtm(i)

   end do

   TotalStep = 0

   NInCyl  = 0
   NOutCyl = 0

   NTotalInCyl  = 0
   NTotalOutCyl = 0

   SumArea = 0.d0

   WcsInCyl  = 0.d0
   WcsOutCyl = 0.d0

   do i = 1 , NJobs

     call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       if((cWhat=='WatOrient').or.(cWhat=='WatOrientG')) then

         ScRg = 0.d0

       else

         call COM_SpecComponent(Rg,Kcomp)

         call InversMatrix(H,InvH)

         ScRg = matmul( InvH, Rg )

       end if

       call CellTransform

       a = H(:,1)
       b = H(:,2)
       call CArea(a,b,Area)

       SumArea = SumArea + Area

       call CylinderWater( ScRg,drz,Wcomp )

     end do

   end do

   if((cWhat=='WatOrient').or.(cWhat=='WatOrientG')) then

     CylArea = RadCyl2 * pi

   else

     CylArea = 0.d0

   end if

   AveArea = SumArea / dble(TotalStep) - CylArea

   RhoInCyl  = 0.d0
   RhoOutCyl = 0.d0

   if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

     do i = -40, 40

       if(NInCyl(i)==0) cycle

       RhoInCyl(i)  = dble(NInCyl(i))  / ( CylArea * drz ) / TotalStep

       do j = -20, 20

         WcsInCyl(i,j)  = WcsInCyl(i,j)  / dble(NInCyl(i))

       end do

     end do

     write(6,'(a,3f8.2/)') 'Rg= ', Rg

   end if

   do i = -40, 40

     if(NOutCyl(i)==0) cycle

     RhoOutCyl(i) = dble(NOutCyl(i)) / ( AveArea * drz ) / TotalStep

     do j = -20, 20

       WcsOutCyl(i,j) = WcsOutCyl(i,j) / dble(NOutCyl(i))

     end do

   end do

   if((cWhat=='WatOrient').or.(cWhat=='WatOrientG')) then

     write(6,'(a,f10.0)') 'Number of water ',&
     &                     real(NTotalInCyl)/TotalStep

   else if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

     write(6,'(a,f10.0)') 'Number of water (in cylinder) ',&
     &                     real(NTotalInCyl)/TotalStep
     write(6,'(a,f10.0)') '           (out of  cylinder) ',&
     &                     real(NTotalOutCyl)/TotalStep

   end if

   do i = -40, 40

     zz = dble(i) * drz

     write(42,'(f5.1,f10.7)') zz,RhoOutCyl(i)

     if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

       write(41,'(f5.1,f10.7)') zz,RhoInCyl(i)

     end if

     do j = -20, 20

       cs = dble(j)/20.d0

       if(j==-20.or.j==20) then

          WcsOutCyl(i,j) = WcsOutCyl(i,j)*2.

         if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

           WcsInCyl(i,j)  = WcsInCyl(i,j) *2.

         end if

       end if

       write(44,'(f5.1,f6.2,f10.5)') zz,cs,WcsOutCyl(i,j)*1.d+2

       if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

         write(43,'(f5.1,f6.2,f10.5)') zz,cs,WcsInCyl(i,j)*1.d+2

       end if

     end do

   end do

close(41)
close(42)
close(43)
close(44)

end subroutine OrientationOfWater


!######################################################################
!######################################################################


subroutine CylinderWater(ScRg,drz,Wcomp)

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use AtomParam, only : Mass

implicit none

integer :: i, j, ii, jj, j1, j2, j3, iz, ics
integer :: Wcomp
integer :: Numa
real(8), dimension(3) :: ScRg, Sij, Rij
real(8), dimension(3,NumMol(Wcomp)) :: Rw,ScRw,Ori
real(8) :: drz, SumM, dvc, R2

   Numa = 0

   do i = 1 , Wcomp-1

     Numa = Numa + NumMol(i)*NumAtm(i)

   end do

   do i = 1 , NumMol(Wcomp)

     Rw(:,i) = 0.d0
     SumM    = 0.d0

     ii = (i-1)*NumAtm(Wcomp) + Numa

     do j = 1 , NumAtm(Wcomp)

       jj = ii + j

       Rw(:,i) = Rw(:,i) + Mass(jj)*R(:,jj)
       SumM    = SumM    + Mass(jj)

     end do

     Rw(:,i) = Rw(:,i) / SumM
     ScRw(:,i) = matmul( InvH, Rw(:,i) )

     j1 = ii + 1
     j2 = ii + 2
     j3 = ii + 3

     Ori(:,i) = R(:,j2) + R(:,j3) - 2.d0 * R(:,j1)
     dvc = dot_product( Ori(:,i), Ori(:,i) )
     dvc = sqrt( dvc )
     Ori(:,i) = Ori(:,i) / dvc

   end do

   if((cWhat=='WatOrient').or.(cWhat=='WatOrientG')) then

   do i = 1, NumMol(Wcomp)

     Sij = ScRw(:,i) - ScRg
     Sij = Sij - nint(Sij)
     Rij = matmul(H,Sij)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2)

     if(cWhat=='WatOrient' ) iz = nint( Rw(3,i) / drz )
     if(cWhat=='WatOrientG') iz = nint( Rij(3)  / drz )

     ics = nint(Ori(3,i)*20)

     NTotalOutCyl = NTotalOutCyl + 1
     NOutCyl(iz)  = NOutCyl(iz)  + 1
     WcsOutCyl(iz,ics) = WcsOutCyl(iz,ics) + 1.

   end do

   else if((cWhat=='WatOriCyl').or.(cWhat=='WatOriCylG')) then

   do i = 1, NumMol(Wcomp)

     Sij = ScRw(:,i) - ScRg
     Sij = Sij - nint(Sij)
     Rij = matmul(H,Sij)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2)

     if(cWhat=='WatOriCyl' ) iz = nint( Rw(3,i) / drz )
     if(cWhat=='WatOriCylG') iz = nint( Rij(3)  / drz )

     ics = nint(Ori(3,i)*20)

     if(R2 <= RadCyl2) then

       NTotalInCyl = NTotalInCyl + 1
       NInCyl(iz)  = NInCyl(iz)  + 1
       WcsInCyl(iz,ics) = WcsInCyl(iz,ics) + 1.

     else

       NTotalOutCyl = NTotalOutCyl + 1
       NOutCyl(iz)  = NOutCyl(iz)  + 1
       WcsOutCyl(iz,ics) = WcsOutCyl(iz,ics) + 1.

     end if

   end do

   end if

end subroutine CylinderWater


!######################################################################
!######################################################################


subroutine CountWaterInPore(ScRg,Count,Wcomp,Numa)

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use AtomParam, only : Mass

implicit none

integer :: i, j, ii, jj
integer :: Numa, Count, Wcomp
real(8), dimension(3) :: ScRg, Sij, Rij
real(8), dimension(3,NumMol(Wcomp)) :: Rw,ScRw
real(8) :: SumM, R2

   Numa = Numa - 1

   do i = 1 , NumMol(Wcomp)

     Rw(:,i) = 0.d0
     SumM    = 0.d0

     ii = (i-1)*NumAtm(Wcomp) + Numa

     do j = 1 , NumAtm(Wcomp)

       jj = ii + j

       Rw(:,i) = Rw(:,i) + Mass(jj)*R(:,jj)
       SumM    = SumM    + Mass(jj)

     end do

     Rw(:,i) = Rw(:,i) / SumM
     ScRw(:,i) = matmul( InvH, Rw(:,i) )

   end do

   Count = 0

   do i = 1, NumMol(Wcomp)

     Sij = ScRw(:,i) - ScRg
     Sij = Sij - nint(Sij)
     Rij = matmul(H,Sij)
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2)

     if((R2 <= RadCyl2).and.(abs(Rij(3)) < Zsh)) then
       PoreWater(i)=1
       Count = Count + 1
     end if

   end do


end subroutine CountWaterInPore


!######################################################################
!######################################################################


subroutine WaterCOM(Rw,Wcomp,Numa)

use Numbers, only : NumMol, NumAtm
use Configuration, only : R
use AtomParam, only : Mass

implicit none

integer :: i, j, ii, jj
integer :: Numa, Wcomp
real(8), dimension(3,NumMol(Wcomp)) :: Rw
real(8) :: SumM

   Numa = Numa - 1

   do i = 1 , NumMol(Wcomp)

     Rw(:,i) = 0.d0
     SumM    = 0.d0

     ii = (i-1)*NumAtm(Wcomp) + Numa

     do j = 1 , NumAtm(Wcomp)

       jj = ii + j
       Rw(:,i) = Rw(:,i) + Mass(jj)*R(:,jj)
       SumM    = SumM    + Mass(jj)

     end do

     Rw(:,i) = Rw(:,i) / SumM

   end do


end subroutine WaterCOM


!######################################################################
!######################################################################


subroutine Water_Diffusion_Zaxis

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use AtomParam, only : Mass
use TimeParam, only : Timeps

implicit none

character(len = 72) :: FileName
integer :: i, j, k, l, lz, Ncw, Nmw, Numa, ll
integer :: lbb, Step
real(8), dimension(:,:,:), allocatable :: ScomW, Ht, InvHt
real(8), dimension(:,:), allocatable :: RcomW
real(8), dimension(:,:), allocatable :: disZ, disXY
integer, dimension(:), allocatable :: zcount
real(8) :: MassW

   Ncw = NumAtm(Kcomp)
   Nmw = NumMol(Kcomp)

   Numa = 0
   do i = 1, NumSpec
     if(i == Kcomp) exit
     Numa = Numa + NumMol(i) * NumAtm(i)
   end do

   MassW = 0.d0
   do i = Numa+1, Numa+Ncw
     MassW = MassW + Mass(i)
   end do

   if(mod(Nsnap,BlockStep) /= 0) then
     write(*,*) 'ERROR : BlockStep'
     write(*,*) Nsnap, BlockStep
     call Finalize
   end if

   lbb  = (BlockStep - ilength) + 1
   lz = maxdz - mindz

   allocate( ScomW(3,Nmw,BlockStep) )
   allocate( RcomW(3,Nmw) )
   allocate( Ht(3,3,BlockStep) )
   allocate( InvHt(3,3,BlockStep) )

   allocate( disZ(lz,ilength) )
   allocate( disXY(lz,ilength) )
   allocate( zcount(lz) )

   do i = 1, lz

     j = nint( ( i + mindz ) * deltaZ )

     if(j<10) then
       write(FileName,'(a,i1,a)') './Analy/MSD_W_Z=',j,'.dat'
     else if(j<100) then
       write(FileName,'(a,i2,a)') './Analy/MSD_W_Z=',j,'.dat'
     else
       write(*,*) 'ERROR : too large Z value'
       write(*,*) j
       call Finalize
     end if

     open(i,file=trim(FileName),status='unknown')

   end do

   Step = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       call Transform

       Step = Step + 1

       Ht(:,:,Step)    = H(:,:)
       InvHt(:,:,Step) = InvH(:,:)

       ll = Numa
       do k = 1, Nmw
         RcomW(:,k) = 0.d0
         do l = 1, Ncw
           ll = ll + 1
           RcomW(:,k) = RcomW(:,k) + Mass(ll) * R(:,ll)
         end do
         RcomW(:,k) = RcomW(:,k) / MassW
       end do

       do k = 1, Nmw
         ScomW(:,k,Step) = matmul( InvH, RcomW(:,k) )
       end do

       if(Step == BlockStep) then
         call MSD_WZ
         Step = 0
       end if

     end do

   end do

Contains

   subroutine MSD_WZ

   implicit none

   integer :: ii, jj, ini, iz
   real(8), dimension(3,ilength) :: Rw
   real(8), dimension(3) :: Si, Ri, Rf
   real(8) :: pref

     disZ  = 0.d0
     disXY = 0.d0
     zcount = 0

     do ini = 1, lbb, interv

       do k = 1, Nmw

         H(:,:) = Ht(:,:,ini)
         RcomW(:,k) = matmul( H, ScomW(:,k,ini) )
         iz = nint(abs(RcomW(3,k)/3.))
         if((iz > maxdz).or.(iz <= mindz)) cycle
         iz = iz - mindz
         zcount(iz) = zcount(iz) + 1
         Rw(:,1) = RcomW(:,k)

         do ii = ini+1, ini+ilength-1

           Rf(:) = Rw(:,ii-ini)

           Si = ScomW(:,k,ii) - ScomW(:,k,ii-1)
           Si = Si - nint( Si )
           Ri(1) = Ht(1,1,ii)*Si(1) + Ht(1,2,ii)*Si(2) + Ht(1,3,ii)*Si(3)
           Ri(2) = Ht(2,1,ii)*Si(1) + Ht(2,2,ii)*Si(2) + Ht(2,3,ii)*Si(3)
           Ri(3) = Ht(3,1,ii)*Si(1) + Ht(3,2,ii)*Si(2) + Ht(3,3,ii)*Si(3)

           Rw(:,ii-ini+1) = Rf(:) + Ri(:)

         end do

         do ii = 2, ilength

           Ri(:) = Rw(:,ii) - Rw(:,1)
           disZ(iz,ii-1)  = disZ(iz,ii-1)  + Ri(3) * Ri(3)
           disXY(iz,ii-1) = disXY(iz,ii-1) + Ri(1) * Ri(1) + Ri(2) * Ri(2)

         end do

       end do

     end do

     do ii = 1, lz

       Timeps = 0.d0

       pref = 1.d0 / dble(zcount(ii))

       do jj = 1, ilength - 1

         timeps = dtime * jj

         write(ii,'(f9.3,2f10.4)') timeps, disZ(ii,jj)*pref, disXY(ii,jj)*pref

       end do

       write(ii,'(a)') ' = = = '

     end do

   end subroutine MSD_WZ

end subroutine Water_Diffusion_Zaxis


!######################################################################
!######################################################################


subroutine Water_Exchange_Correlation

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : R
use ParamAnalyze
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2
use AtomParam, only : AtomName, ResidName
use TimeParam, only : Timeps

implicit none

integer :: i, j, Ncw, Nmw, Numa, Numb
integer :: lbb, Step
integer, dimension(:,:), allocatable :: HdW
integer, dimension(:), allocatable :: ID

open(51,file='./Analy/LifeTimeW.dat')

   Ncw = NumAtm(Kcomp)
   Nmw = NumMol(Kcomp)

   Numa = 0
   do i = 1, NumSpec
     if(i == Kcomp) exit
     Numa = Numa + NumMol(i) * NumAtm(i)
   end do

   Numb = Numa + Nmw * Ncw

   if(mod(Nsnap,BlockStep) /= 0) then
     write(*,*) 'ERROR : BlockStep'
     write(*,*) Nsnap, BlockStep
     call Finalize
   end if

   lbb  = (BlockStep - ilength) + 1

   allocate( Hdw(Nmw,BlockStep) )
   allocate( ID(Nmw) )

   allocate( NumGR_I(NumGR) )
   allocate( ListGR_I(N,NumGR) )

   NumGR_I = 0
   ListGR_I = 0

   do i = 1, NumGR
     do j = 1, N
       if((AtomName(j)==GRAtomI(i)).and.(ResidName(j)==GRResiI(i))) then
         NumGR_I(i) = NumGR_I(i) + 1
         ListGR_I(NumGR_I(i),i) = j
       end if
     end do
   end do

   print *, 'NumGR_I=',NumGR_I

   Step = 0
   Hdw = 0

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ------------------

       call InversMatrix(H,InvH)

       Step = Step + 1

       call NeibourCount

       Hdw(:,Step) = ID(:)

       if(Step == BlockStep) then
         call LifeTimeCo
         Step = 0
         Hdw = 0
       end if

     end do

   end do

Contains

   subroutine NeibourCount

   implicit none

   integer :: k, l, ii, jj, m
   real(8), dimension(3,N) :: ScR
   real(8), dimension(3) :: Sij, Rij
   real(8) :: R2

      do k = 1 , N

        ScR(:,k) = matmul( InvH , R(:,k) )

      end do

      ID = 0

      do k = 1, NumGR

        do l = 1, NumGR_I(k)

          ii = ListGR_I(l,k)

          do jj = Numa+1, Numb, Ncw

            Sij = ScR(:,ii) - ScR(:,jj)
            Sij = Sij - nint( Sij )
            Rij = matmul( H, Sij )

            R2 = dot_product( Rij, Rij )

            if(R2 < Rcutoff2) then

              m = (jj - Numa)/Ncw + 1
              ID(m) = 1

            end if

          end do

        end do

      end do

   end subroutine NeibourCount


   subroutine LifeTimeCo

   implicit none

   integer :: ini, k, ii, jj, Nsample
   integer, dimension(ilength) :: Corr
   real(8) :: xx

     Corr = 0
     Nsample = 0

     do ini = 1, lbb, interv

       do k = 1, Nmw

         if(Hdw(k,ini)==1) then

           Nsample = Nsample + 1

inner:     do ii = ini+1, ini+ilength-1

             if(Hdw(k,ii)==1) then
               Corr(ii-ini) = Corr(ii-ini) + 1
             else
               exit inner
             end if

           end do inner

         end if

       end do

     end do

     write(51,'(a)') '    0.000   1.00000'

     do jj = 1, ilength - 1

       Timeps = dtime * jj
       xx = dble(Corr(jj))/dble(Nsample)

       write(51,'(f9.3,f10.5)') timeps, xx

     end do

   end subroutine LifeTimeCo

end subroutine Water_Exchange_Correlation
