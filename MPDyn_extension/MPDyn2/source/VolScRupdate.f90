!! ############################
! ## SUBROUTINE LIST 
! ## -- VolScRupdate 
! ## -- VolScRupdateRB 
! ## -- VolScVcorrect 
! ## -- VolScVcorrectRB 
! ############################


!######################################################################
!######################################################################


subroutine VolScRupdate

use Numbers, only : N
use CellParam, only : H, InvH, Vsc_Rate
use Configuration, only : R

implicit none

integer :: i
real(8), dimension(3,N) :: ScR

   call ScaledCoordinate(ScR)
   do i = 1, 3
     H(i,i) = H(i,i) + Vsc_Rate(i)
   end do
   call InversMatrix(H,InvH)
   do i = 1, N
     R(:,i) = matmul(H,ScR(:,i))
   end do

end subroutine VolScRupdate


!######################################################################
!######################################################################


subroutine VolScRupdateRB

use CellParam, only : H, InvH, Vsc_Rate
use RBparam, only : R_RB, NumRB

implicit none

integer :: i
real(8), dimension(3,NumRB) :: ScR

   do i = 1, NumRB
     ScR(:,i) = matmul(InvH,R_RB(:,i))
   end do
   do i = 1, 3
     H(i,i) = H(i,i) + Vsc_Rate(i)
   end do
   call InversMatrix(H,InvH)
   do i = 1, NumRB
     R_RB(:,i) = matmul(H,ScR(:,i))
   end do

end subroutine VolScRupdateRB


!######################################################################
!######################################################################


subroutine VolScVcorrect

use CommonBlocks, only : QCorrectCutoff
use CellParam, only : H, InvH, Volume
use TailCorrect

implicit none

real(8) :: det
external det

   call BcastRH
   call InversMatrix(H,InvH)
   call TransCellList
   Volume = det(H)
   if(QCorrectCutoff) then
     Virial_co = 0.d0
     Virial_co(1,1) = CorrectV / (3.d0*Volume)
     Virial_co(2,2) = Virial_co(1,1)
     Virial_co(3,3) = Virial_co(1,1)
     Ene_LJ_co = CorrectE / Volume
   end if

end subroutine VolScVcorrect


!######################################################################
!######################################################################


subroutine VolScVcorrectRB

use CommonBlocks, only : QCorrectCutoff
use CellParam, only : H, InvH, Volume
use TailCorrect

implicit none

real(8) :: det
external det

   call BcastRgQuatH
   call InversMatrix(H,InvH)
   call TransCellList
   Volume = det(H)
   if(QCorrectCutoff) then
     Virial_co = 0.d0
     Virial_co(1,1) = CorrectV / (3.d0*Volume)
     Virial_co(2,2) = Virial_co(1,1)
     Virial_co(3,3) = Virial_co(1,1)
     Ene_LJ_co = CorrectE / Volume
   end if

end subroutine VolScVcorrectRB
