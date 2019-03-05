! ############################
! ## SUBROUTINE LIST 
! ## -- SaveParam 
! ## -- StoreRestart 
! ## -- Write_Config 
! ## -- Write_Velocity 
! ## -- Write_Bath 
! ## -- Write_Time 
! ## -- Write_PDB 
! ## -- Write_CONECT 
! ## -- Write_ConfigDP 
! ############################


!######################################################################
!######################################################################


! ******************************
! ** Output for sequential MD **
! ******************************

subroutine SaveParam

use CommonBlocks, only : QMaster, QPathInt, QResForm
use IOparam, only : DirectoryName

implicit none

   if(QPathInt) then

     call RnmPI

     call VnmPI

     call RVBathPI

   end if

   if(QMaster) then

     if(QResForm) then
       open(1,file=trim(DirectoryName)//'restart.dat',form='formatted',status='unknown')
     else
       open(1,file=trim(DirectoryName)//'restart.dat',form='unformatted',status='unknown')
     end if

     call Write_Config(0)

     call Write_Velocity

     call Write_Bath

     call Write_Time

     close(1)

   end if

end subroutine SaveParam


!######################################################################
!######################################################################


subroutine StoreRestart

use CommonBlocks, only : Job_name, QMaster, QPathInt, QResForm
use IOparam, only : NtrjF, DirectoryName

implicit none

character(len=100) :: RestartFile

   if(QPathInt) then

     call RnmPI

     call VnmPI

     call RVBathPI

   end if

   if(QMaster) then

     if(NtrjF<10) then

       write(RestartFile,'(3a,i1)') &
       & trim(DirectoryName),trim(adjustl(Job_name)),'_restart.00',NtrjF

     else if(NtrjF<100) then

       write(RestartFile,'(3a,i2)') &
       & trim(DirectoryName),trim(adjustl(Job_name)),'_restart.0',NtrjF

     else

       write(RestartFile,'(3a,i3)') &
       & trim(DirectoryName),trim(adjustl(Job_name)),'_restart.',NtrjF

     end if

     if(QResForm) then
       open(1,file=trim(RestartFile),form='formatted',status='unknown')
     else
       open(1,file=trim(RestartFile),form='unformatted',status='unknown')
     end if

     call Write_Config(1)

     call Write_Velocity

     call Write_Bath

     call Write_Time

     close(1)

   end if

end subroutine StoreRestart


!######################################################################
!######################################################################


! **************************
! **  Final Configration  **
! **************************

subroutine Write_Config(ii)

use Numbers, only : N
use CommonBlocks, only : Job_name, QRigidBody, QPathInt, QResForm, QPDB
use Configuration, only : R
use CommonPI
use RBparam, only : NumRB, R_RB, Quaternion
use IOparam, only : NtrjF, DirectoryName
use AtomParam, only : MolName, ResidNum, DefaultAtomName, DefaultResidName

implicit none

integer :: i, ii, j
character(len=100) :: RestartFile
character(len=10) :: ModelName

   if(QPathInt) then
     do i = 1 , N
       R(:,i) = Rnm(:,i,1)
     end do
   end if

   if(QResForm) then

     ModelName = MolName(1)

     write(1,'(2a)') '* ',ModelName
     write(1,'(a)') '* '
     write(1,*) N


     if(N>100000) then
       do i = 1 , N
         write(1,'(2i9,2(x,a4),/3d24.16)')  &
         &     i , ResidNum(i) , DefaultResidName(i) ,  &
         &     DefaultAtomName(i) , R(:,i)
       end do
     else
       do i = 1 , N
         write(1,'(2i5,2(x,a4),/3d24.16)')  &
         &     i , ResidNum(i) , DefaultResidName(i) ,  &
         &     DefaultAtomName(i) , R(:,i)
       end do
     end if

   else

     write(1) N

     write(1) ResidNum
     write(1) DefaultResidName
     write(1) DefaultAtomName
     write(1) R

   end if

   if(QPathInt) then
     call Write_Conf_Bead
   end if

   if(QPDB) call Write_PDB

   if(QRigidBody) then

     close(1)

     if(ii == 0) then

       if(QResForm) then
         open(1,file=trim(DirectoryName)//'restartRB.dat',form='formatted')
       else
         open(1,file=trim(DirectoryName)//'restartRB.dat',form='unformatted')
       end if

     else if(ii == 1) then

       if(NtrjF<10) then
         write(RestartFile,'(3a,i1)') &
         & trim(DirectoryName),trim(adjustl(Job_name)),'_restartRB.00',NtrjF
       else if(NtrjF<100) then
         write(RestartFile,'(3a,i2)') &
         & trim(DirectoryName),trim(adjustl(Job_name)),'_restartRB.0',NtrjF
       else
         write(RestartFile,'(3a,i3)') &
         & trim(DirectoryName),trim(adjustl(Job_name)),'_restartRB.',NtrjF
       end if

       if(QResForm) then
         open(1,file=trim(RestartFile),form='formatted',status='unknown')
       else
         open(1,file=trim(RestartFile),form='unformatted',status='unknown')
       end if

     else if(ii == 2) then

       if(QResForm) then
         open(1,file=trim(DirectoryName)//'restartRB.converted',form='formatted')
       else
         open(1,file=trim(DirectoryName)//'restartRB.converted',form='unformatted')
       end if

     end if

     if(QResForm) then

       do i = 1 , NumRB
         write(1,'(3d24.16)') ( R_RB(j,i) , j = 1 , 3 )
       end do

       do i = 1 , NumRB
         write(1,'(4d24.16)') ( Quaternion(j,i) , j = 1 , 4 )
       end do

     else

       write(1) R_RB
       write(1) Quaternion

     end if

   end if

end subroutine Write_Config


!#####################################################################
!#####################################################################


subroutine Write_Conf_Bead

use Numbers, only : N
use CommonBlocks, only : QResForm
use CommonPI, only : Nbead, Rnm

implicit none

integer :: i, j

  if(QResForm) then

     do j = 2, Nbead
       do i = 1, N
         write(1,'(3d24.16)') Rnm(:,i,j)
       end do
     end do

   else

     do j = 2, Nbead
       do i = 1, N
         write(1) Rnm(:,i,j)
       end do
     end do

   end if

end subroutine Write_Conf_Bead


!#####################################################################
!#####################################################################


! *********************
! **  Final Velocity **
! *********************

subroutine Write_Velocity

use Numbers, only : N
use CommonBlocks, only : QPathInt, QResForm, QRigidBody
use Configuration, only : Vel
use CommonPI, only : Nbead, Vnm
use RBparam, only : NumRB, V_RB, Lmoment, QLinear, QSingle

implicit NONE

integer :: i, j

   if(QPathInt) then

     if(QResForm) then

       do j = 1, Nbead
         do i = 1, N
           write(1,'(3d24.16)') Vnm(:,i,j)
         end do
       end do

     else

       write(1) Vnm

     end if

   else if(QRigidBody) then

     if(QResForm) then

       do i = 1 , NumRB
         write(1,'(3d24.16)') ( V_RB(j,i) , j = 1 , 3 )
       end do

       do i = 1 , NumRB
         if(QLinear(i)) Lmoment(1,i) = 0.d0
         if(QSingle(i)) Lmoment(:,i) = 0.d0
         write(1,'(3d24.16)') ( Lmoment(j,i) , j = 1 , 3 )
       end do

     else

       write(1) V_RB

       do i = 1 , NumRB
         if(QLinear(i)) Lmoment(1,i) = 0.d0
         if(QSingle(i)) Lmoment(:,i) = 0.d0
       end do

       write(1) Lmoment

     end if

   else

     if(QResForm) then

       do i = 1 , N
         write(1,'(3d24.16)') ( Vel(j,i) , j = 1 , 3 )
       end do

     else

       write(1) Vel

     end if

   end if

end subroutine Write_Velocity


!#####################################################################
!#####################################################################


! **************************
! ** Final Bath Parameter **
! **************************

subroutine Write_Bath

use Numbers, only : N
use CommonBlocks, only : QPBC, QThermostat, cThermostatMethod, &
&   QResForm, QPathInt
use CommonPI
use BathParam
use CellParam, only : H

implicit NONE

integer :: i, j, k, l
real(8), parameter :: zero = 0.d0

! ## Thermostat

   if(QThermostat) then

     if((cThermostatMethod == 'NHC').or. &
     &  (cThermostatMethod == 'NH')) then

       if(QResForm) then

         do i = 1 , NHchain
           write(1,'(2d23.16)') Rss(i), Vss(i)
         end do

       else

         do i = 1 , NHchain
           write(1) Rss(i), Vss(i)
         end do

       end if

     else if(cThermostatMethod == 'MNHC') then

       if(QResForm) then

         do i = 1 , NumMNHC
           do j = 1 , NHchain
             write(1,'(2d23.16)') RMNHC(j,i), VMNHC(j,i)
           end do
         end do

       else

         do i = 1 , NumMNHC
           do j = 1 , NHchain
             write(1) RMNHC(j,i), VMNHC(j,i)
           end do
         end do

       end if

     end if

   end if

   if(QPathInt) then

     if(QResForm) then

       do k = 2, Nbead
         do j = 1, NHchain
           do i = 1, N
             do l = 1, 3
               write(1,'(2d23.16)') Rbath(l,i,j,k), Vbath(l,i,j,k)
             end do
           end do
         end do
       end do

     else

       do k = 2, Nbead
         do j = 1, NHchain
           do i = 1, N
             do l = 1, 3
               write(1) Rbath(l,i,j,k), Vbath(l,i,j,k)
             end do
           end do
         end do
       end do

     end if

   end if

! ## Cell matrix 

   if(QPBC) then

     if(QResForm) then

       write(1,'(3d23.16)') (H(1,i), i = 1, 3)
       write(1,'(3d23.16)') (H(2,i), i = 1, 3)
       write(1,'(3d23.16)') (H(3,i), i = 1, 3)
       write(1,'(3d23.16)') (Vg(1,i), i = 1, 3)
       write(1,'(3d23.16)') (Vg(2,i), i = 1, 3)
       write(1,'(3d23.16)') (Vg(3,i), i = 1, 3)

     else

       write(1) H
       write(1) Vg

     end if

   end if

end subroutine Write_Bath


!#####################################################################
!#####################################################################


! **************************
! **  Write Time **
! **************************

subroutine Write_Time

use CommonBlocks, only : SimMethod, QResForm
use CommonHMC, only : TimeMC
use TimeParam, only : Timeps

implicit NONE

   if(SimMethod == 'HMC') then

     if(QResForm) then
       write(1,'(i12)') TimeMC
     else
       write(1) TimeMC
     end if

   else

     if(QResForm) then
       write(1,'(f12.4)') Timeps
     else
       write(1) Timeps
     end if

   end if

end subroutine Write_Time


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_PDB

use Numbers, only : N, NumSpec, NumMol
use CommonBlocks, only : QPBC
use Configuration, only : R
use UnitExParam, only : pi
use IOparam, only : DirectoryName
use CellParam, only : H
use AtomParam, only : MolName, ResidNum, DefaultAtomName, DefaultResidName
use TimeParam, only : Timeps

implicit NONE

integer :: i, ir
character(len=4) :: NameA, RName
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0

real(8), dimension(3) :: va, vb, vc
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca

open(2,file=trim(DirectoryName)//'final.pdb')

   write(2,'(a)') &
   & 'TITLE     The Last Structure by a Molecular Dynamics Run'

   write(2,'(a)') &
   & 'REMARK   1 The Simulated System is composed with'

   do i = 1 , NumSpec

     write(2,'(a21,i1,a3,a10,i5)') &
     &     'REMARK   1 Component ',i,' : ',MolName(i),NumMol(i)

   end do

   write(2,'(a)')         'REMARK   2 '
   write(2,'(a,f11.3,a)') 'REMARK   2  Time = ',Timeps,' ps'
   write(2,'(a)')         'REMARK   2 '

   if(QPBC) then

     va = H(:,1)
     vb = H(:,2)
     vc = H(:,3)

     LLa = sqrt( dot_product(va,va) )
     LLb = sqrt( dot_product(vb,vb) )
     LLc = sqrt( dot_product(vc,vc) )
     Aab = acos( dot_product(va,vb) / LLa / LLb ) * 180. / pi
     Abc = acos( dot_product(vb,vc) / LLb / LLc ) * 180. / pi
     Aca = acos( dot_product(vc,va) / LLc / LLa ) * 180. / pi

     write(2,'(a,3f9.3,3f7.2)') 'CRYST1',LLa,LLb,LLc,Abc,Aca,Aab

   end if

   do i = 1, N

     NameA = DefaultAtomName(i)
     RName = DefaultResidName(i)
     if(ResidNum(i)<100000) then
       ir = ResidNum(i)
     else
       ir = mod(ResidNum(i),100000)
     end if

     if(DefaultResidName(i)=='HSD') RName='HIS'
     if(NameA(4:4)/=' ') then

       write(2,'(a4,i7,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,DefaultAtomName(i),RName,ir,  &
       &  R(:,i),zero,zero,NameA(1:1)

     else

       write(2,'(a4,i7,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,DefaultAtomName(i),RName,ir, &
       & R(:,i),zero,zero,NameA(1:1)

     end if

   end do

   write(2,'(a)') 'END'

close(2)

end subroutine Write_PDB


!######################################################################
!######################################################################


subroutine Write_CONECT

use BondedParam, only : NumBond, BondI, BondJ

implicit none

integer :: i,j,k

   if(NumBond/=0) then
     open(2,file='bond.pdb')
     do k = 1, NumBond
       i = BondI(k)
       j = BondJ(k)
       if(i>100000.or.j>100000) cycle
       write(2,'(a6,2i5)') 'CONECT',i,j
     end do
     write(2,'(a)') 'END'
     close(2)
   end if

end subroutine Write_CONECT


!######################################################################
!######################################################################


subroutine Write_ConfigDP

use Numbers, only : NumSpec, NumMol, NumAtm
use TimeParam, only : Timeps
use Configuration
use CommonDPD
use RBparam, only : R_RB, V_RB, Quaternion, Lmoment

implicit none

integer :: i, j, k, l
!
! #####################################
! #        write config to disk       #
! #####################################
!
open(1,file='restartDP.dat',form='formatted',status='unknown')

   l = 0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do j = 1 , NumMol(i)

         write(1,'(6f10.6)') R_RB(:,j),V_RB(:,j)
         write(1,'(7f10.6)') Quaternion(:,j),Lmoment(:,j)

       end do

       l = l + NumMol(i) * NumAtm(i)

     else

       do j = 1, NumMol(i)

         do k = 1 , NumAtm(i)

           l = l + 1
           write(1,'(6f10.6)') R(:,l),Vel(:,l)

         end do

       end do

     end if

   end do

   write(1,'(f10.6)') SlideGap

   write(1,'(f15.3)') Timeps

close(1)

end subroutine Write_ConfigDP
