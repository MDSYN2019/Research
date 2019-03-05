! ############################
! ## SUBROUTINE LIST 
! ## -- PairList 
! ## -- PairListTICR
! ## -- PairListCG
! ## -- PairListMacro
! ## -- PairListCGMacro
! ## -- PairListMacroMacro
! ############################


!#####################################################################
!#####################################################################


! ********************************************
! **  making a pair list between particles  **
! ********************************************

subroutine PairList

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QMaster, Qstdout, QRigidBody, QPathInt, ForceField, QMacro
use CommonMPI
use CommonPI
use RBparam
use BookParam, only : MaxPair, Npair, ListIJ
use CellParam, only : H, CellShft
use CutoffParam, only : Rbook2

implicit NONE

integer :: i, j, k, Nas
integer :: IRB

real(8) :: Six, Siy, Siz
real(8) :: Rix, Riy, Riz
real(8) :: Sx, Sy, Sz
real(8) :: Rx, Ry, Rz
integer :: Nx, Ny, Nz

real(8), dimension(3,N) :: ScR

real(8) :: R2

integer :: NProcsTemp, MyRankTemp

   if(ForceField(1:2)=='CG') Return

   if(QMacro) then
     call PairListMacro
     Return
   end if

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
     call PBC
   end if

   call ScaledCoordinate(ScR)

   Npair = 0

   Nas = NProcsTemp - MyRankTemp

   if(QRigidBody) then

     do i = Nas, N, NProcsTemp

       IRB = AtomUnitNum(i)
       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(IRB == AtomUnitNum(j)) cycle

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

       do j = i+1 , N, 2

         if(IRB == AtomUnitNum(j)) cycle

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   else

     do i = Nas, N, NProcsTemp

       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

       do j = i+1 , N, 2

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   end if

   if(Qstdout.and.QMaster) write(*,*) 'Npair=',Npair

   if( Npair > MaxPair ) then
     if(QMaster) then
       write(*,*)  'ERROR : Npair exceeds MaxPair! '
#ifndef BMONI
       write(11,*) 'ERROR : Npair exceeds MaxPair! '
#endif
     end if
     call Finalize
   end if

end subroutine PairList


!#####################################################################
!#####################################################################


! ********************************************
! **  making a pair list between particles  **
! ********************************************

subroutine PairListTICR

use Numbers, only : N
use CommonBlocks, only : QMaster, Qstdout, QPathInt
use Configuration, only : R
use CommonMPI
use CommonPI
use RBparam
use FEparam
use BookParam, only : MaxPair, Npair, ListIJ
use CellParam, only : H, CellShft
use CutoffParam, only : Rbook2

implicit NONE

integer :: i, j, k, Nas
integer :: IRB, Nx, Ny, Nz

real(8) :: Six, Siy, Siz
real(8) :: Rix, Riy, Riz
real(8) :: Sx, Sy, Sz
real(8) :: Rx, Ry, Rz

real(8), dimension(3,N) :: ScR

real(8) :: R2

integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
     call PBC
   end if

   call ScaledCoordinate(ScR)

   Npair = 0

   Nas = NProcsTemp - MyRankTemp

   do i = Nas, N, NProcsTemp

     if(.not.QTICRcalc(i)) cycle

     IRB = AtomUnitNum(i)

     Six = ScR(1,i)
     Siy = ScR(2,i)
     Siy = ScR(3,i)
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     do j = i-2 , 1, -2

       if(.not.QTICRcalc(j)) cycle

       if(IRB == AtomUnitNum(j)) cycle

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rbook2 ) then

         Npair           = Npair + 1
         ListIJ(1,Npair) = i
         ListIJ(2,Npair) = j
         ListIJ(3,Npair) = k

       end if

     end do

     do j = i+1 , N, 2

       if(.not.QTICRcalc(j)) cycle

       if(IRB == AtomUnitNum(j)) cycle

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rbook2 ) then

         Npair           = Npair + 1
         ListIJ(1,Npair) = i
         ListIJ(2,Npair) = j
         ListIJ(3,Npair) = k

       end if

     end do

   end do

   if(Qstdout.and.QMaster) write(*,*) 'Npair=',Npair

   if( Npair > MaxPair ) then

     if(QMaster) then
       write(*,*)  'ERROR : Npair exceeds MaxPair! '
#ifndef BMONI
       write(11,*) 'ERROR : Npair exceeds MaxPair! '
#endif
     end if

     call Finalize

   end if

end subroutine PairListTICR


!#####################################################################
!#####################################################################


! ********************************************
! **  making a pair list between particles  **
! ********************************************

subroutine PairListCG

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QMaster, QRigidBody, Qstdout, QMacro
use CommonMPI
use CommonPI
use RBparam
use NoLJparam
use CGdata
use BookParam, only : MaxPair, Npair, ListIJ
use CellParam, only : H, CellShft

implicit NONE

integer :: i, j, k, Nas, non, check_bonded
integer :: IRB, itype, jtype, Nx, Ny, Nz

#ifdef GEN
real(8) :: Six, Siy, Siz
real(8) :: Sx, Sy, Sz
real(8), dimension(3,N) :: ScR
#else
real(8) :: clhx, clhy, clhz
#endif
real(8) :: Rix, Riy, Riz
real(8) :: Rx, Ry, Rz

real(8) :: R2

external check_bonded

   if(QMacro) then
     call PairListCGMacro
     Return
   end if

#ifdef GEN
   call ScaledCoordinate(ScR)
#else
   clhx = H(1,1)*0.5d0
   clhy = H(2,2)*0.5d0
   clhz = H(3,3)*0.5d0
#endif

   Npair = 0

   Nas = NProcs - MyRank

   if(QRigidBody) then

     do i = Nas, N, NProcs

       IRB = AtomUnitNum(i)
       itype = NBAtomType(i)
       non   = NumNoLJ(i)
#ifdef GEN
       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
#endif
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(IRB == AtomUnitNum(j)) cycle

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)

#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         endif

       end do

       do j = i+1 , N, 2

         if(IRB == AtomUnitNum(j)) cycle

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   else

     do i = Nas, N, NProcs

       non   = NumNoLJ(i)
       itype = NBAtomType(i)
#ifdef GEN
       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
#endif
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

       do j = i+1 , N, 2

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   end if

   if(QMaster.and.Qstdout) write(*,*) 'Npair=',Npair

   if( Npair > MaxPair ) then

     if(QMaster) then
       write(*,*)  'ERROR : Npair exceeds MaxPair! '
#ifndef BMONI
       write(11,*) 'ERROR : Npair exceeds MaxPair! '
#endif
     end if

     call Finalize

   end if

end subroutine PairListCG


!#####################################################################
!#####################################################################


subroutine PairListMacro

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QMaster, Qstdout, QRigidBody, QPathInt, ForceField
use CommonMPI
use CommonPI
use RBparam
use BookParam, only : MaxPair, Npair, ListIJ
use CellParam, only : H, CellShft
use CutoffParam, only : Rbook2
use CGball, only : Rlstmacro2, NmacList, MaxPairMac, IDsphere, Listmac

implicit NONE

integer :: i, j, k, Nas
integer :: IRB

real(8) :: Six, Siy, Siz
real(8) :: Rix, Riy, Riz
real(8) :: Sx, Sy, Sz
real(8) :: Rx, Ry, Rz
integer :: Nx, Ny, Nz

real(8), dimension(3,N) :: ScR

real(8) :: R2

integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
     call PBC
   end if

   call ScaledCoordinate(ScR)

   Npair = 0

   Nas = NProcsTemp - MyRankTemp

   if(QRigidBody) then

     do i = Nas, N, NProcsTemp

       if(IDsphere(i)/=0) cycle

       IRB = AtomUnitNum(i)
       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(IRB == AtomUnitNum(j)) cycle
         if(IDsphere(j)/=0) cycle

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

       do j = i+1 , N, 2

         if(IRB == AtomUnitNum(j)) cycle
         if(IDsphere(j)/=0) cycle

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   else

     do i = Nas, N, NProcsTemp

       if(IDsphere(i)/=0) cycle

       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(IDsphere(j)/=0) cycle

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

       do j = i+1 , N, 2

         if(IDsphere(j)/=0) cycle

         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else 
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbook2 ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   end if

   call PairListMacrovsAtom(ScR)

   if(Qstdout.and.QMaster) then
     write(*,*) 'Npair   =',Npair
     write(*,*) 'NmacList=',NmacList
   end if

   if( Npair > MaxPair ) then
     if(QMaster) then
       write(*,*)  'ERROR : Npair exceeds MaxPair! '
#ifndef BMONI
       write(11,*) 'ERROR : Npair exceeds MaxPair! '
#endif
     end if
     call Finalize
   end if

end subroutine PairListMacro


!#####################################################################
!#####################################################################


! ********************************************
! **  making a pair list between particles  **
! ********************************************

subroutine PairListCGMacro

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QMaster, QRigidBody, Qstdout
use CommonMPI
use RBparam
use NoLJparam
use CGdata
use BookParam, only : MaxPair, Npair, ListIJ
use CellParam, only : H, CellShft
use CGball, only : NumSphere, Rlstmacro2, NmacList, MaxPairMac, IDsphere, Listmac

implicit NONE

integer :: i, j, k, Nas, non, check_bonded
integer :: IRB, itype, jtype, Nx, Ny, Nz

#ifdef GEN
real(8) :: Six, Siy, Siz
real(8) :: Sx, Sy, Sz
real(8), dimension(3,N) :: ScR
#else
real(8) :: clhx, clhy, clhz
#endif
real(8) :: Rix, Riy, Riz
real(8) :: Rx, Ry, Rz

real(8) :: R2

external check_bonded

#ifdef GEN
   call ScaledCoordinate(ScR)
#else
   clhx = H(1,1)*0.5d0
   clhy = H(2,2)*0.5d0
   clhz = H(3,3)*0.5d0
#endif

   Npair = 0

   Nas = NProcs - MyRank

   if(QRigidBody) then

     do i = Nas, N, NProcs

       if(IDsphere(i)/=0) cycle

       IRB = AtomUnitNum(i)
       itype = NBAtomType(i)
       non   = NumNoLJ(i)
#ifdef GEN
       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
#endif
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(IRB == AtomUnitNum(j)) cycle

         if(IDsphere(j)/=0) cycle

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)

#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         endif

       end do

       do j = i+1 , N, 2

         if(IRB == AtomUnitNum(j)) cycle

         if(IDsphere(j)/=0) cycle

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   else

     do i = Nas, N, NProcs

       if(IDsphere(i)/=0) cycle

       non   = NumNoLJ(i)
       itype = NBAtomType(i)
#ifdef GEN
       Six = ScR(1,i)
       Siy = ScR(2,i)
       Siz = ScR(3,i)
#endif
       Rix = R(1,i)
       Riy = R(2,i)
       Riz = R(3,i)

       do j = i-2 , 1, -2

         if(check_bonded(i,j,non) == 1) cycle

         if(IDsphere(j)/=0) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

       do j = i+1 , N, 2

         if(check_bonded(i,j,non) == 1) cycle

         if(IDsphere(j)/=0) cycle

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
#ifdef GEN
         Sx = Six - ScR(1,j)
         Sy = Siy - ScR(2,j)
         Sz = Siz - ScR(3,j)
         if(Sx>0.5) then
           Nx = -9
         else if(Sx<-0.5) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Sy>0.5) then
           Ny = -3
         else if(Sy<-0.5) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Sz>0.5) then
           Nz = -1
         else if(Sz<-0.5) then
           Nz =  1
         else
           Nz =  0
         end if
#else
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
#endif
         k = Nx + Ny + Nz
         Rx = Rx + CellShft(1,k)
         Ry = Ry + CellShft(2,k)
         Rz = Rz + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if ( R2 <= Rbk2(itype,jtype) ) then
           Npair           = Npair + 1
           ListIJ(1,Npair) = i
           ListIJ(2,Npair) = j
           ListIJ(3,Npair) = k
         end if

       end do

     end do

   end if

#ifdef GEN
   call PairListMacrovsCGParticle(ScR,Nas)
   if(NumSphere>1) call PairListMacroMacro(Nas)
#else
   call PairListMacrovsCGParticle(clhx,clhy,clhz,Nas)
   if(NumSphere>1) call PairListMacroMacro(clhx,clhy,clhz,Nas)
#endif
   if(Qstdout.and.QMaster) then
     write(*,*) 'Npair   =',Npair
     write(*,*) 'NmacList=',NmacList
   end if

   if( Npair > MaxPair ) then
     if(QMaster) then
       write(*,*)  'ERROR : Npair exceeds MaxPair! '
#ifndef BMONI
       write(11,*) 'ERROR : Npair exceeds MaxPair! '
#endif
     end if
     call Finalize
   end if

end subroutine PairListCGMacro


!#####################################################################
!#####################################################################


subroutine PairListMacrovsAtom(ScR)

use Numbers, only : N
use Configuration, only : R
use CommonMPI, only : NProcs
use CellParam, only : H, CellShft
use CommonBlocks, only : QMaster
use CGball, only : Rlstmacro2, NmacList, MaxPairMac, Listmac, Nox, IdOx, &
&   NumSphere, IdSph

implicit NONE

integer :: i, j, k, ii, jj, Nl
real(8) :: Six, Siy, Siz
real(8) :: Rix, Riy, Riz
real(8) :: Sx, Sy, Sz
real(8) :: Rx, Ry, Rz
integer :: Nx, Ny, Nz
real(8), dimension(3,N) :: ScR
real(8) :: R2, RC2

   RC2 = Rlstmacro2(1,1)

   do ii = 1, NumSphere

     i = IdSph(ii)

     Six = ScR(1,i)
     Siy = ScR(2,i)
     Siz = ScR(3,i)
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     Nl = 0

     do jj = 1, Nox

       j = IdOx(jj)

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else 
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= RC2 ) then
         Nl               = Nl + 1
         Listmac(1,Nl,ii) = j
         Listmac(2,Nl,ii) = k
       end if

     end do

     NmacList(ii) = Nl

     if( NmacList(ii) > MaxPairMac ) then
       if(QMaster) then
         write(*,*)  'ERROR : NmacList exceeds MaxPairMac! '
#ifndef BMONI
         write(11,*) 'ERROR : NmacList exceeds MaxPairMac! '
#endif
       end if
       call Finalize
     end if

   end do

end subroutine PairListMacrovsAtom


!#####################################################################
!#####################################################################


#ifdef GEN
subroutine PairListMacrovsCGParticle(ScR,Nas)
#else
subroutine PairListMacrovsCGParticle(clhx,clhy,clhz,Nas)
#endif

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use CommonBlocks, only : QMaster
use CGdata
use CellParam, only : H, CellShft
use CGball, only : NumSphere, IdSph, Rlstmacro2, NmacList, MaxPairMac, IDsphere,&
&  Listmac, ItypeSph

implicit NONE

integer :: i, j, k, Nas, jtype, ii, Nl, itype
integer :: Nx, Ny, Nz

#ifdef GEN
real(8) :: Six, Siy, Siz
real(8) :: Sx, Sy, Sz
real(8), dimension(3,N) :: ScR
#else
real(8) :: clhx, clhy, clhz
#endif
real(8) :: Rix, Riy, Riz
real(8) :: Rx, Ry, Rz
real(8) :: R2

   NmacList(:) = 0

   do ii = 1, NumSphere

     i = IdSph(ii)
     itype = ItypeSph(ii)

#ifdef GEN
     Six = ScR(1,i)
     Siy = ScR(2,i)
     Siz = ScR(3,i)
#endif
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     Nl = 0

     do j = Nas, N, NProcs

       if(IDsphere(j)/=0) cycle
       jtype = NBAtomType(j)

       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
#ifdef GEN
       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else 
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
#else
       if(Rx>clhx) then
         Nx = -9
       else if(Rx<-clhx) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Ry>clhy) then
         Ny = -3
       else if(Ry<-clhy) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Rz>clhz) then
         Nz = -1
       else if(Rz<-clhz) then
         Nz =  1
       else
         Nz =  0
       end if
#endif
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rlstmacro2(jtype,itype) ) then
         Nl               = Nl + 1
         Listmac(1,Nl,ii) = j
         Listmac(2,Nl,ii) = k
       end if

     end do

     NmacList(ii) = Nl

     if( NmacList(ii) > MaxPairMac ) then
       if(QMaster) then
         write(*,*)  'ERROR : NmacList exceeds MaxPairMac! '
#ifndef BMONI
         write(11,*) 'ERROR : NmacList exceeds MaxPairMac! '
#endif
       end if
       call Finalize
     end if

   end do

end subroutine PairListMacrovsCGParticle


!#####################################################################
!#####################################################################


#ifdef GEN
subroutine PairListMacroMacro(Nas)
#else
subroutine PairListMacroMacro(clhx,clhy,clhz,Nas)
#endif

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use CommonBlocks, only : QMaster
use CGdata
#ifdef GEN
use CellParam, only : H, InvH, CellShft
#else
use CellParam, only : H, CellShft
#endif
use CGball, only : NumSphere, Rbkmm2, MaxPairMM, NumPairMM, ListMM, IdSph, ItypeSph

implicit NONE

integer :: i, j, k, Nas, ii, itype, jtype
integer :: Nx, Ny, Nz

#ifdef GEN
real(8) :: Six, Siy, Siz
real(8) :: Sx, Sy, Sz
real(8), dimension(3,NumSphere) :: ScRmc
#else
real(8) :: clhx, clhy, clhz
#endif
real(8), dimension(3,NumSphere) :: Rmc
real(8) :: Rix, Riy, Riz
real(8) :: Rx, Ry, Rz
real(8) :: R2

   NumPairMM = 0

   do ii = 1, NumSphere
     i = IdSph(ii)
     Rmc(:,ii) = R(:,i)
   end do

#ifdef GEN
   do i = 1, NumSphere
     ScRmc(1,i) = InvH(1,1)*Rmc(1,i) + InvH(1,2)*Rmc(2,i) + InvH(1,3)*Rmc(3,i)
     ScRmc(2,i) = InvH(2,1)*Rmc(1,i) + InvH(2,2)*Rmc(2,i) + InvH(2,3)*Rmc(3,i)
     ScRmc(3,i) = InvH(3,1)*Rmc(1,i) + InvH(3,2)*Rmc(2,i) + InvH(3,3)*Rmc(3,i)
   end do
#endif

   do i = Nas, NumSphere, NProcs

     itype = ItypeSph(i)
#ifdef GEN
     Six = ScRmc(1,i)
     Siy = ScRmc(2,i)
     Siz = ScRmc(3,i)
#endif
     Rix = Rmc(1,i)
     Riy = Rmc(2,i)
     Riz = Rmc(3,i)

     do j = i-2 , 1, -2

       jtype = ItypeSph(j)
       Rx = Rix - Rmc(1,j)
       Ry = Riy - Rmc(2,j)
       Rz = Riz - Rmc(3,j)
#ifdef GEN
       Sx = Six - ScRmc(1,j)
       Sy = Siy - ScRmc(2,j)
       Sz = Siz - ScRmc(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else 
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
#else
       if(Rx>clhx) then
         Nx = -9
       else if(Rx<-clhx) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Ry>clhy) then
         Ny = -3
       else if(Ry<-clhy) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Rz>clhz) then
         Nz = -1
       else if(Rz<-clhz) then
         Nz =  1
       else
         Nz =  0
       end if
#endif
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rbkmm2(itype,jtype) ) then
         NumPairMM        = NumPairMM + 1
         ListMM(1,NumPairMM) = i
         ListMM(2,NumPairMM) = j
         ListMM(3,NumPairMM) = k
       end if

     end do

     do j = i+1 , NumSphere, 2

       jtype = ItypeSph(j)
       Rx = Rix - Rmc(1,j)
       Ry = Riy - Rmc(2,j)
       Rz = Riz - Rmc(3,j)
#ifdef GEN
       Sx = Six - ScRmc(1,j)
       Sy = Siy - ScRmc(2,j)
       Sz = Siz - ScRmc(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else 
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
#else
       if(Rx>clhx) then
         Nx = -9
       else if(Rx<-clhx) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Ry>clhy) then
         Ny = -3
       else if(Ry<-clhy) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Rz>clhz) then
         Nz = -1
       else if(Rz<-clhz) then
         Nz =  1
       else
         Nz =  0
       end if
#endif
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rbkmm2(itype,jtype) ) then
         NumPairMM        = NumPairMM + 1
         ListMM(1,NumPairMM) = i
         ListMM(2,NumPairMM) = j
         ListMM(3,NumPairMM) = k
       end if

     end do

   end do

   if( NumPairMM > MaxPairMM ) then
     if(QMaster) then
       write(*,*)  'ERROR : NumPairMM exceeds MaxPairMM! '
#ifndef BMONI
       write(11,*) 'ERROR : NumPairMM exceeds MaxPairMM! '
#endif
     end if
     call Finalize
   end if

end subroutine PairListMacroMacro
