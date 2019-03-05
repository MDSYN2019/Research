! ############################
! ## SUBROUTINE LIST 
! ## -- EnergyMinimize 
! ############################


!######################################################################
!######################################################################


!*******************************************
!*     Minimize the energy                 *
!*     Steepest Decent Algorithm           *
!*******************************************

subroutine EnergyMinimize

use Numbers, only : N, Nf, NfT, NfR
use CommonBlocks, only : QMaster, QPBC, QRigidBody, QSHAKE, &
&   ForceField, QCoulomb
use Configuration, only : R
use RBparam, only : NumRB, R_RB, F_RB, Torque, Quaternion
use Eminparam
use CommonMPI
use SHAKEparam, only : R_o
use UnitExParam, only : pi, cvol
use EwaldParam, only : Frc_Eksp, Ene_Eksp, Ene_Eslf
use OptConstraintParam, only : NHam, Frc_OptC, Ene_OptC
use NonbondParam, only : Frc_Elec, Frc_NBlong, Frc_NBshrt, Ene_Elec, &
&   Ene_Ersp, Ene_LJ, Ene_NBlong, Ene_NBshrt, Ene_ELlong, Ene_ELshrt
use BondedParam, only : Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro, &
& Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro


implicit NONE

integer :: Itry, i
real(8) :: sumF2, RmsF, RmsF_o
real(8) :: sumT2, RmsT, RmsT_o
real(8) :: deltaR, deltaQ, dQmax
real(8) :: deltaE, Energy, Energy_Pre

real(8), dimension(3,N) :: Force
real(8), dimension(3,N) :: F_o
real(8), dimension(3,NumRB) :: F_RBo
real(8), dimension(3,NumRB) :: Torqueo
real(8), dimension(3,NumRB) :: R_RBo
real(8), dimension(4) :: qq, qqt


   if(QMaster) then

     write( 6,'(a)') 'Energy minimization has just started'
#ifndef BMONI
     write(11,'(5x,a)') '===== Energy Minimizaion ====='
#endif

   end if

   deltaR = dRmax

   dQmax  = pi / 56.d0
   deltaQ = dQmax

   if(NHam/=0) call Rot_FixedPoint

! ## FORCE -----------------------------------

   call Force_OptC
   call Force_HamC

   if(ForceField(1:2) == 'CG') then

     call Force_Bond_CG
     call Force_Angle_CG
     call Force_Dihed_CG
     call Force_Improper_CG

     call Force_nonbond_CG_long(0)
     call Force_nonbond_CG_short
     if(QCoulomb) call Force_kspace

     Ene_Elec = Ene_ELlong + Ene_ELshrt + Ene_Eksp + Ene_Eslf/NProcs
     Frc_Elec = Frc_Eksp + Frc_NBlong + Frc_NBshrt
     Ene_LJ = Ene_NBlong + Ene_NBshrt

   else

     call Force_Bond
     call Force_Angle
     call Force_UB
     call Force_Dihedral
     call Force_Improper

     if(QPBC) then
       call PairList
       call Force_rcut
     else
       call Force_iso
     end if

   end if

!----------------------------------------------------------------------
   call SumEnergy(Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro, &
   &              Ene_LJ, Ene_Elec, Ene_Ersp, Ene_Eksp, Ene_OptC)
!----------------------------------------------------------------------

   Energy = Ene_Bond  + Ene_Angle + Ene_UB  &
   &      + Ene_Dihed + Ene_Impro           &
   &      + Ene_LJ    + Ene_Elec  + Ene_OptC

   Energy_Pre = Energy

   Force = Frc_Bond  + Frc_Angle + Frc_UB  &
   &     + Frc_Dihed + Frc_Impro           &
   &     + Frc_Elec  + Frc_OptC

   call SumFrc( Force )

! ----------------------------------------------

   if(QMaster) then

     write( 6,'(10x,a,e18.8,a)')                               &
     &  'Initial Energy : ', Energy_Pre * cvol , ' [kcal / mol]'
#ifndef BMONI
     write(11,'(10x,a,e18.8,a)')                               &
     &  'Initial Energy : ', Energy_Pre * cvol , ' [kcal / mol]'
#endif
     sumF2 = 0.d0

     if(QRigidBody) then

       sumT2 = 0.d0

       call Force_Div_Component(Force,F_RB,Torque)

       do i = 1 , NumRB

         sumF2 = sumF2 + dot_product( F_RB(:,i), F_RB(:,i) )
         sumT2 = sumT2 + dot_product( Torque(:,i), Torque(:,i) )

       end do

       RmsF = sqrt( sumF2 / NfT )
       RmsT = sqrt( sumT2 / NfR )

     else

       do i = 1 , N

         sumF2 = sumF2 + dot_product( Force(:,i), Force(:,i) )

       end do

       RmsF = sqrt( sumF2 / Nf )

     end if

   end if

! ------------------------
! ## Minimization Routine
! ------------------------
   do Itry = 1 , MinTry

     if(QMaster) then

       R_o = R

       RmsF_o = RmsF

       if(QRigidBody) then

         R_RBo = R_RB
         F_RBo = F_RB
         Torqueo = Torque

         RmsT_o = RmsT

         do i = 1 , NumRB

           R_RB(:,i) = R_RB(:,i) + deltaR * F_RB(:,i) / RmsF

           qq = Quaternion(:,i)

           qqt(1) = - qq(2) * Torque(1,i) - qq(3) * Torque(2,i) - qq(4) * Torque(3,i)
           qqt(2) =   qq(1) * Torque(1,i) - qq(4) * Torque(2,i) + qq(3) * Torque(3,i)
           qqt(3) =   qq(4) * Torque(1,i) + qq(1) * Torque(2,i) - qq(2) * Torque(3,i)
           qqt(4) = - qq(3) * Torque(1,i) + qq(2) * Torque(2,i) + qq(1) * Torque(3,i)

           Quaternion(:,i) = Quaternion(:,i) + 0.5d0 * qqt(:) / RmsT * deltaQ

         end do

         call IntraMolVec

       else

         F_o = Force

         do i = 1 , N

           R(:,i) = R(:,i) + deltaR * Force(:,i) / RmsF

         end do

         if(QSHAKE) call SHAKE

       end if

     end if

     call BcastR

! ## FORCE -----------------------------------
     call Force_OptC
     call Force_HamC

     if(ForceField(1:2) == 'CG') then

       call Force_Bond_CG
       call Force_Angle_CG
       if(mod(Itry,10)==0) then
         call Force_nonbond_CG_long(0)
         if(QCoulomb) call Force_kspace
       end if
       call Force_nonbond_CG_short

       Ene_Elec = Ene_ELlong + Ene_ELshrt + Ene_Eksp + Ene_Eslf/NProcs
       Frc_Elec = Frc_Eksp + Frc_NBlong + Frc_NBshrt
       Ene_LJ = Ene_NBlong + Ene_NBshrt

     else

       call Force_Bond
       call Force_Angle
       call Force_UB
       call Force_Dihedral
       call Force_Improper
       if(QPBC) then
         if(mod(Itry,10)==0) call PairList
         call Force_rcut
       else
         call Force_iso
       end if

      end if

!-------------------------------------------------------------------------
      call SumEnergy(Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro, &
      &              Ene_LJ, Ene_Elec, Ene_Ersp, Ene_Eksp, Ene_OptC)
!-------------------------------------------------------------------------

     if(QMaster) then

       Energy = Ene_Bond  + Ene_Angle + Ene_UB  &
       &      + Ene_Dihed + Ene_Impro           &
       &      + Ene_LJ    + Ene_Elec  + Ene_OptC

     end if

     Force = Frc_Bond  + Frc_Angle + Frc_UB  &
     &     + Frc_Dihed + Frc_Impro           &
     &     + Frc_Elec  + Frc_OptC

     call SumFrc( Force )

! ----------------------------------------------

     if(QMaster) then

       deltaE = Energy - Energy_Pre

       sumF2 = 0.d0

       if(QRigidBody) then

         sumT2 = 0.d0

         call Force_Div_Component(Force,F_RB,Torque)

         do i = 1 , NumRB

           sumF2 = sumF2 + dot_product( F_RB(:,i), F_RB(:,i) )
           sumT2 = sumT2 + dot_product( Torque(:,i), Torque(:,i) )

         end do

         RmsF = sqrt( sumF2 / NfT )
         RmsT = sqrt( sumT2 / NfR )

       else

         do i = 1 , N
           sumF2 = sumF2 + dot_product( Force(:,i), Force(:,i) )
         end do

         RmsF = sqrt( sumF2 / Nf )

       end if

       write(6,'(a,i3,a,3e12.4/8x,a,3e12.4)') 'Step ',Itry,&
       &                          'Energy= ',Energy*cvol,Ene_LJ*cvol,Ene_Elec*cvol,&
       &                          'st,bn,ub',Ene_Bond*cvol,Ene_Angle*cvol,Ene_UB*cvol

       if(QRigidBody) then
         write(6,'(8x,a,e12.4)')    'Rms F, Torque =',RmsF * cvol,RmsT * cvol
       else
         write(6,'(8x,a,e12.4)')    'Rms F =',RmsF * cvol
       end if

     end if

     call BcastE(deltaE,Energy)

     if( abs(deltaE/Energy) < dev_relative ) exit

     if(QMaster) then

       if( Energy < Energy_Pre ) then

         Energy_Pre = Energy
         deltaR     = deltaR * 1.2d0
         deltaR     = min(deltaR , dRmax)

         if(QRigidBody) then

           deltaQ = deltaQ * 1.2d0
           deltaQ = min(deltaQ , dQmax)

         end if

       else

         if(QRigidBody) then
           R_RB   = R_RBo
           F_RB   = F_RBo
           Torque = Torqueo
           RmsT   = RmsT_o
           deltaQ = deltaQ * 0.5d0
         end if

         Force  = F_o
         R      = R_o
         RmsF   = RmsF_o
         deltaR = deltaR * 0.5d0

       end if

     end if

   end do

   if(QMaster) then

     write( 6,'(10x,a,e18.8,a)')                              &
     &     'Final   Energy : ', Energy * cvol , ' [kcal / mol]'
#ifndef BMONI
     write(11,'(10x,a,e18.8,a)')                              &
     &     'Final   Energy : ', Energy * cvol , ' [kcal / mol]'
#endif
     write( 6,*) 'Enegy minimization has just finished'
#ifndef BMONI
     write(11,'(/)')
#endif
   end if

end subroutine EnergyMinimize
