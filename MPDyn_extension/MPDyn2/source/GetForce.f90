! ############################
! ## SUBROUTINE LIST 
! ## -- GetForce
! ## -- GetForceIso
! ############################


!######################################################################
!######################################################################


subroutine GetForce(i,istep)

use CommonBlocks, only : QTICR, ForceField, QCoulomb, QMacro
use Numbers, only : N
use CellParam, only : H
use CommonMPI, only : NProcs, MyRank
use CGball, only : NumSphere

implicit none

integer :: i, istep, Nas
#ifdef GEN
real(8), dimension(3,N) :: ScR
#endif
real(8) :: clhx, clhy, clhz

if(ForceField(1:2) == 'CG') then

   if(i == 0) then

     if(QMacro.and.(istep==0)) then
       Nas = NProcs - MyRank
#ifdef GEN
       call ScaledCoordinate(ScR)
       call PairListMacrovsCGParticle(ScR,Nas)
       if(NumSphere>1) call PairListMacroMacro(ScR,Nas)
#else
       clhx = H(1,1)*0.5d0
       clhy = H(2,2)*0.5d0
       clhz = H(3,3)*0.5d0
       call PairListMacrovsCGParticle(clhx,clhy,clhz,Nas)
       if(NumSphere>1) call PairListMacroMacro(clhx,clhy,clhz,Nas)
#endif
     end if

     call Force_OptC
     call Force_HamC
     call Force_Bond_CG
     call Force_Angle_CG
     call Force_Dihed_CG
     call Force_Improper_CG

     call Force_nonbond_CG_long(istep)
     call Force_nonbond_CG_short
     if(QCoulomb) call Force_kspace

   else if(i == 1) then

     call Force_OptC
     call Force_HamC
     call Force_Bond_CG
     call Force_Angle_CG
     call Force_Dihed_CG
     call Force_Improper_CG

   else if(i == 2) then

     call Force_nonbond_CG_short

   else if(i == 3) then

     call Force_nonbond_CG_long(istep)
     if(QCoulomb) call Force_kspace

   end if

else

   if(i == 0) then

     call Force_OptC
     call Force_HamC
     call Force_Bond
     call Force_Angle
     call Force_UB
     call Force_Dihedral
     call Force_Improper
     if(QTICR) then
       call PairListTICR
     else
       call PairList                !  list-up pairs of atoms
     end if
     call Force_rspace            !  calculation of interaction in the real space
     call Force_EAM
     if(QCoulomb) call Force_kspace  !  calculation of interaction in the reciprocal space

   else if(i == 1) then

     call Force_OptC
     call Force_HamC
     call Force_Bond
     call Force_Angle
     call Force_UB
     call Force_Dihedral
     call Force_Improper

   else if(i == 2) then

     call Force_rspace            !  calculation of interaction in real space
     call Force_EAM

   else if(i == 3) then

     if(QCoulomb) call Force_kspace  !  calculation of interaction in reciprocal space

   end if

end if

end subroutine GetForce


!######################################################################
!######################################################################


subroutine GetForceIso(i)

use CommonBlocks, only : ForceField

implicit none

integer :: i

if(ForceField(1:2) == 'CG') then

   if(i == 0) then

     call Force_OptC
     call Force_HamC
     call Force_Bond_CG
     call Force_Angle_CG
     call Force_Dihed_CG
     call Force_Improper_CG

     call Force_nonbond_CG_long(0)
     call Force_nonbond_CG_short

   else if(i == 1) then

     call Force_OptC
     call Force_HamC
     call Force_Bond_CG
     call Force_Angle_CG
     call Force_Dihed_CG
     call Force_Improper_CG

   else if(i == 2) then

     call Force_nonbond_CG_short

   else if(i == 3) then

     call Force_nonbond_CG_long(0)

   end if

else

   if(i == 0) then

     call Force_OptC
     call Force_HamC
     call Force_Bond
     call Force_Angle
     call Force_UB
     call Force_Dihedral
     call Force_Improper
     call Force_iso

   else if(i == 1) then

     call Force_OptC
     call Force_HamC
     call Force_Bond
     call Force_Angle
     call Force_UB

   else if(i == 2) then

     call Force_Dihedral
     call Force_Improper

   else if(i == 3) then

     call Force_iso

   end if

end if

end subroutine GetForceIso


!######################################################################
!######################################################################


subroutine GetForceVV

use CommonBlocks, only : ForceField

implicit none

   if(ForceField=='F11') then

     call ForceDPD_Original

   else if(ForceField=='F22') then

     call ForceDPD_Smooth22

   else if(ForceField=='F23') then

     call ForceDPD_Smooth23

   else if(ForceField=='Morse') then

     call ForceDPD_Attractive

   end if

end subroutine GetForceVV


!######################################################################
!######################################################################


subroutine GetForceVVSC

use CommonBlocks, only : ForceField

implicit none

   if(ForceField=='F11') then

     call ForceDPD_Original_SC

   else if(ForceField=='F22') then

     call ForceDPD_Smooth22_SC

   else if(ForceField=='F23') then

     call ForceDPD_Smooth23_SC

   else if(ForceField=='Morse') then

     call ForceDPD_Attractive_SC

   end if


end subroutine GetForceVVSC


!######################################################################
!######################################################################


subroutine GetForceLowePeters

use CommonBlocks, only : ForceField

implicit none

   if(ForceField=='F11') then

     call ForceDPD_Original_Lowe

   else if(ForceField=='F12') then

     call ForceDPD_Smooth12_Lowe

   else if(ForceField=='F21') then

     call ForceDPD_Smooth21_Lowe

   else if(ForceField=='F22') then

     call ForceDPD_Smooth22_Lowe

   else if(ForceField=='F23') then

     call ForceDPD_Smooth23_Lowe

   else if(ForceField=='Morse') then

     call ForceDPD_Attractive_Lowe

   end if

end subroutine GetForceLowePeters
