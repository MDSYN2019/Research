subroutine GetAcc(Acc,iflag)

use Numbers, only : N
use CommonBlocks, only : QPBC, QRigidBody
use EAM_param, only : Frc_EAM
use EwaldParam, only : Frc_Eksp
use OptConstraintParam, only : Frc_OptC
use NonbondParam, only : Frc_Ersp, Frc_NBlong, Frc_NBshrt, Frc_Elec
use BondedParam, only : Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro
use AtomParam, only : InvMass

implicit none

integer :: iflag, i
real(8), dimension(3,N) :: Acc

   if(QPBC) then

     if(iflag==1) then
       do i = 1, N
         Acc(:,i) = ( Frc_Bond (:,i) &
         &          + Frc_Angle(:,i) &
         &          + Frc_UB   (:,i) &
         &          + Frc_Dihed(:,i) &
         &          + Frc_Impro(:,i) &
         &          + Frc_OptC (:,i) )
       end do
     else if(iflag==2) then
       do i = 1, N
         Acc(:,i) = ( Frc_EAM   (:,i) &
         &          + Frc_Ersp  (:,i) &
         &          + Frc_NBshrt(:,i) )
       end do
     else if(iflag==3) then
       do i = 1, N
         Acc(:,i) = ( Frc_Eksp  (:,i) &
         &          + Frc_NBlong(:,i) )
       end do
     else
       write(*,*) 'ERROR: GetAcc'
       call Finalize
     end if

   else ! isolated system

     if(iflag==1) then
       do i = 1 , N
         Acc(:,i) = ( Frc_Bond (:,i) &
         &          + Frc_Angle(:,i) &
         &          + Frc_UB   (:,i) &
         &          + Frc_OptC (:,i) )
       end do
     else if(iflag==2) then
       do i = 1 , N
         Acc(:,i) = ( Frc_Dihed (:,i) &
         &          + Frc_Impro (:,i) &
         &          + Frc_NBshrt(:,i) )
       end do
     else if(iflag==3) then
       do i = 1 , N
         Acc(:,i) = ( Frc_Elec  (:,i) &
         &          + Frc_NBlong(:,i) )
       end do
     else
       write(*,*) 'ERROR: GetAcc'
       call Finalize
     end if

   end if

   if(.not.QRigidBody) then
     do i = 1, N
       Acc(:,i) = Acc(:,i) * InvMass(i)
     end do
   end if

end subroutine GetAcc
