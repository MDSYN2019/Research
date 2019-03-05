! ############################
! ## SUBROUTINE LIST 
! ## -- CalcTemp 
! ## -- CalcTempDP 
! ## -- CalcTempPI 
! ## -- BathTemp 
! ############################


!#####################################################################
!#####################################################################


! *****************************************************
! **  subroutine for calculating temparature         **
! **  this subroutine is coded for a system          **
! **  including constraints.                         **
! *****************************************************

subroutine CalcTemp

use Numbers, only : N, Nf, NfT, NfR
use CommonBlocks, only : QRigidBody
use Configuration, only : Vel
use RBparam, only : NumRB, V_RB, QSingle, RBType, NumRBAtom, MassRB, &
&   Lmoment, InvInertiaRB, QLinear
use UnitExParam, only : kb
use AtomParam, only : Mass
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Ene_kinT, Ene_kinR, Pkinp

implicit NONE

integer :: i, k, MyType
real(8) :: SumM

   Pkinp = 0.d0

   if(QRigidBody) then

     SumM = 0.d0

     k = 0

     do i = 1 , NumRB

       if(QSingle(i)) then

         k = k + 1

         Pkinp(1,1) = Pkinp(1,1) + Mass(k) * V_RB(1,i) * V_RB(1,i)
         Pkinp(1,2) = Pkinp(1,2) + Mass(k) * V_RB(1,i) * V_RB(2,i)
         Pkinp(1,3) = Pkinp(1,3) + Mass(k) * V_RB(1,i) * V_RB(3,i)
         Pkinp(2,2) = Pkinp(2,2) + Mass(k) * V_RB(2,i) * V_RB(2,i)
         Pkinp(2,3) = Pkinp(2,3) + Mass(k) * V_RB(2,i) * V_RB(3,i)
         Pkinp(3,3) = Pkinp(3,3) + Mass(k) * V_RB(3,i) * V_RB(3,i)

       else

         MyType = RBType(i)
         k = k + NumRBAtom(MyType)

         Pkinp(1,1) = Pkinp(1,1) + MassRB(MyType) * V_RB(1,i) * V_RB(1,i)
         Pkinp(1,2) = Pkinp(1,2) + MassRB(MyType) * V_RB(1,i) * V_RB(2,i)
         Pkinp(1,3) = Pkinp(1,3) + MassRB(MyType) * V_RB(1,i) * V_RB(3,i)
         Pkinp(2,2) = Pkinp(2,2) + MassRB(MyType) * V_RB(2,i) * V_RB(2,i)
         Pkinp(2,3) = Pkinp(2,3) + MassRB(MyType) * V_RB(2,i) * V_RB(3,i)
         Pkinp(3,3) = Pkinp(3,3) + MassRB(MyType) * V_RB(3,i) * V_RB(3,i)

         if(QLinear(i)) then
           SumM = SumM + Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) &
           &           + Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType)
         else
           SumM = SumM + Lmoment(1,i) * Lmoment(1,i) * InvInertiaRB(1,MyType) &
           &           + Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) &
           &           + Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType)
         end if

       end if

     end do

   else

     do i = 1 , N

       Pkinp(1,1) = Pkinp(1,1) + Mass(i) * Vel(1,i) * Vel(1,i)
       Pkinp(1,2) = Pkinp(1,2) + Mass(i) * Vel(1,i) * Vel(2,i)
       Pkinp(1,3) = Pkinp(1,3) + Mass(i) * Vel(1,i) * Vel(3,i)
       Pkinp(2,2) = Pkinp(2,2) + Mass(i) * Vel(2,i) * Vel(2,i)
       Pkinp(2,3) = Pkinp(2,3) + Mass(i) * Vel(2,i) * Vel(3,i)
       Pkinp(3,3) = Pkinp(3,3) + Mass(i) * Vel(3,i) * Vel(3,i)

     end do

   end if

   Pkinp(2,1) = Pkinp(1,2)
   Pkinp(3,1) = Pkinp(1,3)
   Pkinp(3,2) = Pkinp(2,3)

   if(QRigidBody) then

     Ene_kinT = Pkinp(1,1) + Pkinp(2,2) + Pkinp(3,3)
     Ene_kinR = SumM

     Temp_Translation = Ene_kinT / ( dble(NfT) * kb )
     Temp_Rotation    = Ene_kinR / ( dble(NfR) * kb )

     Ene_kin = Ene_kinT + Ene_kinR

   else

     Ene_kin = Pkinp(1,1) + Pkinp(2,2) + Pkinp(3,3)

   end if

   Temp = Ene_kin / ( dble(Nf) * kb )

end subroutine CalcTemp


!#####################################################################
!#####################################################################


! *****************************************************
! **  subroutine for calculating temparature         **
! **  this subroutine is coded for a system          **
! **  including constraints.                         **
! *****************************************************

! ##  INPUT  : Vel, V_RB, Lmoment
! ##  OUTPUT : Temp, Ene_kin

subroutine CalcTempDP

use Numbers, only : Nf, NumSpec, NumMol, NumAtm
use CommonDPD
use RBparam, only : MassRB, R_RB, V_RB, Lmoment, InvInertiaRB
use Configuration
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation, &
&   Ene_kin, Ene_kinT, Ene_kinR, Pkinp

implicit NONE

integer :: i, j, k, l
real(8) :: SumM, VxSheared

   Pkinp = 0.d0

   SumM = 0.d0

   l = 0

! #####################
   if(QSheared) then
! #####################

     do i = 1 , NumSpec

       if(ColloidFlag(i)) then

         do j = 1 , NumMol(i)

           VxSheared  = V_RB(1,j) - ShearRate * R_RB(2,j)
           Pkinp(1,1) = Pkinp(1,1) + MassRB(1) * VxSheared * VxSheared
           Pkinp(1,2) = Pkinp(1,2) + MassRB(1) * VxSheared * V_RB(2,j)
           Pkinp(1,3) = Pkinp(1,3) + MassRB(1) * VxSheared * V_RB(3,j)
           Pkinp(2,2) = Pkinp(2,2) + MassRB(1) * V_RB(2,j) * V_RB(2,j)
           Pkinp(2,3) = Pkinp(2,3) + MassRB(1) * V_RB(2,j) * V_RB(3,j)
           Pkinp(3,3) = Pkinp(3,3) + MassRB(1) * V_RB(3,j) * V_RB(3,j)

           SumM = SumM + Lmoment(1,j) * Lmoment(1,j) * InvInertiaRB(1,1) &
           &           + Lmoment(2,j) * Lmoment(2,j) * InvInertiaRB(2,1) &
           &           + Lmoment(3,j) * Lmoment(3,j) * InvInertiaRB(3,1)

         end do

         l = l + NumMol(i)*NumAtm(i)

       else

         do j = 1 , NumMol(i)

           do k = 1 , NumAtm(i)

             l = l + 1

             VxSheared  = Vel(1,l) - ShearRate * R(2,l)
             Pkinp(1,1) = Pkinp(1,1) + VxSheared * VxSheared
             Pkinp(1,2) = Pkinp(1,2) + VxSheared * Vel(2,l)
             Pkinp(1,3) = Pkinp(1,3) + VxSheared * Vel(3,l)
             Pkinp(2,2) = Pkinp(2,2) + Vel(2,l)  * Vel(2,l)
             Pkinp(2,3) = Pkinp(2,3) + Vel(2,l)  * Vel(3,l)
             Pkinp(3,3) = Pkinp(3,3) + Vel(3,l)  * Vel(3,l)

           end do

         end do

       end if

     end do

! #####################
   else
! #####################

     do i = 1 , NumSpec

       if(ColloidFlag(i)) then

         do j = 1 , NumMol(i)

           Pkinp(1,1) = Pkinp(1,1) + MassRB(1) * V_RB(1,j) * V_RB(1,j)
           Pkinp(1,2) = Pkinp(1,2) + MassRB(1) * V_RB(1,j) * V_RB(2,j)
           Pkinp(1,3) = Pkinp(1,3) + MassRB(1) * V_RB(1,j) * V_RB(3,j)
           Pkinp(2,2) = Pkinp(2,2) + MassRB(1) * V_RB(2,j) * V_RB(2,j)
           Pkinp(2,3) = Pkinp(2,3) + MassRB(1) * V_RB(2,j) * V_RB(3,j)
           Pkinp(3,3) = Pkinp(3,3) + MassRB(1) * V_RB(3,j) * V_RB(3,j)

           SumM = SumM + Lmoment(1,j) * Lmoment(1,j) * InvInertiaRB(1,1) &
           &           + Lmoment(2,j) * Lmoment(2,j) * InvInertiaRB(2,1) &
           &           + Lmoment(3,j) * Lmoment(3,j) * InvInertiaRB(3,1)

         end do

         l = l + NumMol(i)*NumAtm(i)

       else

         do j = 1 , NumMol(i)

           do k = 1 , NumAtm(i)

             l = l + 1

             Pkinp(1,1) = Pkinp(1,1) + Vel(1,l) * Vel(1,l)
             Pkinp(1,2) = Pkinp(1,2) + Vel(1,l) * Vel(2,l)
             Pkinp(1,3) = Pkinp(1,3) + Vel(1,l) * Vel(3,l)
             Pkinp(2,2) = Pkinp(2,2) + Vel(2,l) * Vel(2,l)
             Pkinp(2,3) = Pkinp(2,3) + Vel(2,l) * Vel(3,l)
             Pkinp(3,3) = Pkinp(3,3) + Vel(3,l) * Vel(3,l)

           end do

         end do

       end if

     end do

! #####################
   end if
! #####################

   Pkinp(2,1) = Pkinp(1,2)
   Pkinp(3,1) = Pkinp(1,3)
   Pkinp(3,2) = Pkinp(2,3)

   Ene_kinT = Pkinp(1,1) + Pkinp(2,2) + Pkinp(3,3)
   Ene_kinR = SumM

   Ene_kin = Ene_kinT + Ene_kinR

   Temp = Ene_kin / ( dble(Nf) )

   if(QColloid) then

     Temp_Translation = Ene_kinT / ( dble(Nf - 3*NumColloid) )
     Temp_Rotation    = Ene_kinR / ( dble(3*NumColloid) )

   end if

end subroutine CalcTempDP


!#####################################################################
!#####################################################################


! *****************************************************
! **  subroutine for calculating temparature         **
! **  this subroutine is coded for a system          **
! **  including constraints.                         **
! *****************************************************

subroutine CalcTempPI

use Numbers, only : N, Nf
use CommonPI
use UnitExParam, only : kb
use ThermoData, only : Temp, Ene_kin, Pkinp

implicit NONE

integer :: i, j

   Pkinp = 0.d0

   if(QMasterPI) then

     do j = IniBead, FinBead

       do i = 1 , N

         Pkinp(1,1) = Pkinp(1,1) + FictMass(i,j) * Vnm(1,i,j) * Vnm(1,i,j)
         Pkinp(1,2) = Pkinp(1,2) + FictMass(i,j) * Vnm(1,i,j) * Vnm(2,i,j)
         Pkinp(1,3) = Pkinp(1,3) + FictMass(i,j) * Vnm(1,i,j) * Vnm(3,i,j)
         Pkinp(2,2) = Pkinp(2,2) + FictMass(i,j) * Vnm(2,i,j) * Vnm(2,i,j)
         Pkinp(2,3) = Pkinp(2,3) + FictMass(i,j) * Vnm(2,i,j) * Vnm(3,i,j)
         Pkinp(3,3) = Pkinp(3,3) + FictMass(i,j) * Vnm(3,i,j) * Vnm(3,i,j)

       end do

     end do

   end if

! -----------------
   call SumTempPI
! -----------------

   Pkinp(2,1) = Pkinp(1,2)
   Pkinp(3,1) = Pkinp(1,3)
   Pkinp(3,2) = Pkinp(2,3)

   Ene_kin = Pkinp(1,1) + Pkinp(2,2) + Pkinp(3,3)

   Temp = Ene_kin / ( dble( Nf * Nbead ) * kb )

end subroutine CalcTempPI


!######################################################################
!######################################################################


subroutine BathTemp

use CommonBlocks, only : cBarostatMethod
use BathParam, only : Baro_kin, Mp, Vg

implicit none

integer :: i, j

   if( cBarostatMethod == 'PR' .or.  cBarostatMethod == 'ST' ) then

     Baro_kin = 0.d0

     do j = 1 , 3
       do i = 1 , 3
         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)
       end do
     end do

   else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then

     Baro_kin = 0.d0

     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)  ! Wg*Tr(Vg^t*Vg)
     end do

   else if( cBarostatMethod == 'AN' ) then

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

   end if

end subroutine BathTemp
