! ############################
! ## SUBROUTINE LIST 
! ## -- IntegrEOM_NVT_PI
! ## -- IntegrEOM_NPT_PI
! ## -- IntegrEOM_iso_PI
! ############################


!######################################################################
!######################################################################


! ***********************************************************
! ** ####### Microcanonical or Canonical Ensemble  ######  **
! ** #######    <<<< Constraint Dynamics >>>>      ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine IntegrEOM_NVT_PI

use Numbers, only : N
use CommonBlocks, only : QMaster, QCorrectCutoff, &
&   QHeatCap, QDelCellMove
use Configuration, only : R
use EAM_param, only : Frc_EAM, Vir_EAM, Ene_EAM
use CommonPI
use EwaldParam, only : Frc_Eksp, Vir_Eksp, Ene_Eksp
use OptConstraintParam, only : Frc_OptC, Vir_OptC, Ene_OptC
use NonbondParam, only : Frc_Ersp, Ene_LJ, Ene_Ersp, Vir_Ersp
use BondedParam, only : &
&   Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro, &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : Volume
use TailCorrect
use TimeParam, only : Nstep, ixc, ixv, lp, lk, isampleHC, BookFreq, Timeps, &
&   deltat, dt2, irs
use ThermoData, only : Virial

implicit none

integer :: i, j, istep, iref_step
real(8), dimension(3,N,Nbead) :: F_fast, Fnm_fast
real(8), dimension(3,N,Nbead) :: F_mode, Fnm_mode
real(8), dimension(3,N,Nbead) :: F_slow, Fnm_slow
real(8), dimension(3,3) :: V_fast, V_mode, V_slow
real(8) :: E_fast, E_mode, E_slow
real(8) :: dt2ref, dtp2, dtk2

real(8), dimension(3,N,Nbead) :: FrcPI
real(8) :: CVtrace, Dummy
real(8) :: HessCVTrace
real(8), dimension(3) :: Ftmp, Rtmp
integer :: ioutput

   dt2ref = 0.5d0 * deltat_ref

   dtp2 = dt2 * lp
   dtk2 = dt2 * lk

   if(QHeatCap) then
     ioutput = 1000
     ioutput = ioutput / isampleHC
   end if

!   if(QMaster.and.QOpFix) then
!     call ConstPrepare
!   end if

   call GetForce_Ref

   call NormalModeTransform

   do i = 1, N
     R(:,i) = Rnm(:,i,1)
   end do
   call PairList

   Rcentroid = R

   F_fast = 0.d0
   F_mode = 0.d0
   F_slow = 0.d0
   Virial = 0.d0

   do j = IniBead, FinBead

     do i = 1, N
       R(:,i) = Rpi(:,i,j)
     end do

     call GetForce(1,0)
     call GetForce(2,0)
     call GetForce(3,0)

     do i = 1, N

       F_fast(:,i,j) = ( Frc_Bond (:,i) &
     &                 + Frc_Angle(:,i) &
     &                 + Frc_UB   (:,i) &
     &                 + Frc_Dihed(:,i) &
     &                 + Frc_Impro(:,i) &
     &                 + Frc_OptC (:,i) ) * InvP

       F_mode(:,i,j) = ( Frc_EAM  (:,i) &
     &                 + Frc_Ersp (:,i) ) * InvP

       F_slow(:,i,j) = Frc_Eksp(:,i) * InvP

     end do

     Virial = Virial + Vir_Bond  + Vir_Angle + Vir_UB   &
     &               + Vir_Dihed + Vir_Impro + Vir_OptC &
     &               + Vir_Ersp  + Vir_EAM   + Vir_Eksp

   end do

   call GetNormalModeForce( F_fast, Fnm_fast )
   call GetNormalModeForce( F_mode, Fnm_mode )
   call GetNormalModeForce( F_slow, Fnm_slow )

   call SumFrcPI( Fnm_fast )
   call SumFrcPI( Fnm_mode )
   call SumFrcPI( Fnm_slow )

   call SumVir( Virial )

   if(QCorrectCutoff) then

     Virial_co = 0.d0

     Virial_co(1,1) = CorrectV / (3.d0*Volume)
     Virial_co(2,2) = Virial_co(1,1)
     Virial_co(3,3) = Virial_co(1,1)

     Ene_LJ_co = CorrectE / Volume

   end if

   Virial = Virial - Virial_co

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

! ## Multiple Time Step
! ---------------------------------------------------
     if((QMaster).and.(mod(istep-1,lk) == 0)) then
!
! ## thermostat velocities and thermostat positions
!    ------------------------------
       call ThermostatPIcent(lk)
!    ------------------------------
     end if
!
! ----------------------------------
! ## update the particle velocities
! ----------------------------------
!
     if(QMasterPI) then

       if(mod(istep-1,lk) == 0) then

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ( dtk2 * Fnm_slow(:,i,j) &
             &                         + dtp2 * Fnm_mode(:,i,j) &
             &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

           end do

         end do

       else if(mod(istep-1,lp) == 0) then

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ( dtp2 * Fnm_mode(:,i,j) &
             &                         + dt2  * Fnm_fast(:,i,j) )* InvFictMass(i,j)

           end do

         end do

       else

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2 * Fnm_fast(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if

! ## Multiple time step

       do iref_step = 1, Ncref
!       ------------------------
         call ThermostatPInonc
!       ------------------------
! -------------------------------
! update the particle positions
! -------------------------------
         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2ref * Fnm_ref(:,i,j) * InvFictMass(i,j)

           end do

         end do

         do j = IniBead, FinBead

           do i = 1, N

             Rnm(:,i,j) = Rnm(:,i,j) + deltat_ref * Vnm(:,i,j)

           end do

         end do

!         if(QOpFix) call AddConstR

         call GetForce_Ref

         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2ref * Fnm_ref(:,i,j) * InvFictMass(i,j)

           end do

         end do

!         if(QOpFix) call AddConstV(istep)

         call ThermostatPInonc

       end do

     end if

     if(QMaster.and.mod(istep,BookFreq)==0) then
! ## periodic boundary condition
       call PBC_PI
     end if

     call NormalModeTransform

     if(QMaster) then
       do i = 1, N
         R(:,i) = Rnm(:,i,1)
       end do
     end if
     call BcastR
     Rcentroid = R

! ## get the new force
!   ---------------------------------------------
! ## Book keeping 
     if( mod(istep,BookFreq) == 0 ) then
       R = Rcentroid
       call PairList
     end if

! ## dt_short
     E_fast = 0.d0
     V_fast = 0.d0

     do j = IniBead, FinBead

       do i = 1, N
         R(:,i) = Rpi(:,i,j)
       end do

       call GetForce(1,0)

       do i = 1, N

         F_fast(:,i,j) = ( Frc_Bond (:,i) &
         &               + Frc_Angle(:,i) &
         &               + Frc_UB   (:,i) &
         &               + Frc_Dihed(:,i) &
         &               + Frc_Impro(:,i) &
         &               + Frc_OptC (:,i) ) * InvP

       end do

       E_fast = E_fast + Ene_Bond  + Ene_Angle + Ene_UB &
       &      + Ene_Dihed + Ene_Impro + Ene_OptC

       V_fast = V_fast + Vir_Bond  + Vir_Angle + Vir_UB &
       &      + Vir_Dihed + Vir_Impro + Vir_OptC

     end do

     call GetNormalModeForce( F_fast, Fnm_fast )

     call SumFrcPI( Fnm_fast )

! ---------------------------------------------
! Multiple time step

     if(mod(istep,lp) == 0) then

       E_mode = 0.d0
       V_mode = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call GetForce(2,0)

         do i = 1, N

           F_mode(:,i,j) = ( Frc_EAM  (:,i) &
           &               + Frc_Ersp (:,i) ) * InvP

         end do

         E_mode = E_mode + Ene_LJ + Ene_EAM + Ene_Ersp

         V_mode = V_mode + Vir_Ersp + Vir_EAM

       end do

       call GetNormalModeForce( F_mode, Fnm_mode )

       call SumFrcPI( Fnm_mode )

     end if

! ------------------------------------------------------

     if(mod(istep,lk) == 0) then

       E_slow = 0.d0
       V_slow = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call GetForce(3,0)

         do i = 1, N

           F_slow(:,i,j) = Frc_Eksp(:,i) * InvP

         end do

         E_slow = E_slow + Ene_Eksp

         V_slow = V_slow + Vir_Eksp

       end do

       call GetNormalModeForce( F_slow, Fnm_slow )

       call SumFrcPI( Fnm_slow )

     end if

     if(QMasterPI) then

       if(mod(istep,lk) == 0) then

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ( dtk2 * Fnm_slow(:,i,j) &
             &                         + dtp2 * Fnm_mode(:,i,j) &
             &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

           end do

         end do

       else if(mod(istep,lp) == 0) then

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ( dtp2 * Fnm_mode(:,i,j) &
             &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

           end do

         end do

       else

         do j = IniBead, FinBead

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2 * Fnm_fast(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if

     end if

!       if(QOpFix) call AddConstV(istep)

     if(QMaster.and.(mod(istep,lk) == 0)) then
! ## thermostat velocities and thermostat positions
       call ThermostatPIcent(lk)
! ## remove the cell-momentum
       if(QDelCellMove) call Elim_CellMove

     end if
! ------------------------------------------------------

!   - save parameters ---------------------------------------

     if(QHeatCap.and.(mod(istep,isampleHC) == 0)) then

       Virial = V_fast + V_mode + V_slow
       Virial = Virial * InvP

       call SumVir( Virial )

       FrcPI = F_fast + F_mode + F_slow

       call SumFrcPI( FrcPI )

       EnePI = E_fast + E_mode + E_slow
       EnePI = EnePI * InvP

       if(QMaster) then
         do i = 1, N
           R(:,i) = Rnm(:,i,1)
         end do
       end if
       call BcastR
       if(.not.QMaster) then
         do i = 1, N
           Rnm(:,i,1) = R(:,i)
         end do
       end if

       CVtrace = 0.d0

       if(QMasterPI) then

         do j = IniBead, FinBead

           do i = 1, N

             Ftmp(:) = - FrcPI(:,i,j)
             Rtmp(:) = Rpi(:,i,j) - Rnm(:,i,1)

             CVtrace = CVtrace + dot_product(Rtmp, Ftmp)

           end do

         end do

       end if

       HessCVTrace = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call Calc_Hessian( HessCVTrace, Dummy )

       end do

       HessCVTrace = HessCVTrace * InvP * 0.5d0

       call SumHC(EnePI, CVtrace, HessCVTrace, Dummy, Dummy)

       if(QMaster) then
         call HeatCapacity(EnePI,CVtrace,Dummy,HessCVTrace,Dummy,istep/isampleHC,ioutput)
       end if

     end if

     if( mod(istep,lk) == 0 ) then
       EnePI  = E_fast + E_mode + E_slow
       Virial = V_fast + V_mode + V_slow
       call Print_Energy_NV_PI(istep)
     end if

     if( mod(istep,ixc) == 0 ) call Print_Config_PI
     if( mod(istep,ixv) == 0 ) call Print_Velocity_PI
     if( mod(istep,irs) == 0 ) call SaveParam

!   ---------------------------------------------------------

   end do

end subroutine IntegrEOM_NVT_PI


!######################################################################
!######################################################################


! ***********************************************************
! ** #######     Isothermal-Isobaric Ensemble      ######  **
! ** #######    <<<< Flexible chain model >>>>     ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine IntegrEOM_NPT_PI

use Numbers, only : N
use CommonBlocks, only : QMaster, QHeatCap, QCorrectCutoff, &
&   cBarostatMethod, ForceField
use Configuration, only : R
use EAM_param, only : Frc_EAM, Vir_EAM, Ene_EAM
use CommonPI
use BathParam, only : Vg
use EwaldParam, only : Frc_Eksp, Vir_Eksp, Ene_Eksp
use OptConstraintParam, only : NHam, Frc_OptC, Vir_OptC, Ene_OptC
use NonbondParam, only : Frc_Ersp, Ene_LJ, Ene_Ersp, Vir_Ersp
use BondedParam, only : &
&   Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro, &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : H, InvH, Volume
use TailCorrect
use TimeParam, only : Nstep, ixc, ixv, lp, lk, isampleHC, BookFreq, Timeps, &
&   deltat, dt2, irs
use ThermoData, only : Virial

implicit none

real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8), dimension(3,3) :: tempoH
real(8), dimension(3)   :: eigenValue
real(8), dimension(3)   :: fc1, fc2
real(8), dimension(3)   :: tempor, tempov
integer :: i, j, istep, iref_step
real(8) :: cf, arg2
real(8) :: poly
real(8) :: bb
real(8) :: Anaa, Anaa2, ClSc
real(8), dimension(3) :: aa, aa2, aa3
real(8), dimension(3) :: arg3, poly3
real(8), dimension(3) :: bb3

real(8), dimension(3,N,Nbead) :: F_fast, Fnm_fast
real(8), dimension(3,N,Nbead) :: F_mode, Fnm_mode
real(8), dimension(3,N,Nbead) :: F_slow, Fnm_slow

real(8) :: E_fast, E_mode, E_slow
real(8), dimension(3,3) :: V_fast, V_mode, V_slow

real(8) :: dt2ref, dtp2, dtk2

integer :: ioutput

real(8) :: det
External det

if(ForceField(1:3) == 'EAM') then
open(13,file='CellMatrix.dat',form='unformatted')
end if

   dt2ref = 0.5d0 * deltat_ref

   dtp2 = dt2 * lp
   dtk2 = dt2 * lk

   if(QHeatCap) then
     ioutput = 1000
     ioutput = ioutput / isampleHC
   end if

!   if(QMaster.and.QOpFix) then
!     call ConstPrepare
!   end if

   call GetForce_Ref

   call NormalModeTransform

   do i = 1, N
     R(:,i) = Rnm(:,i,1)
   end do
   call PairList

   Rcentroid = R

   F_fast = 0.d0
   F_mode = 0.d0
   F_slow = 0.d0
   Virial = 0.d0

   do j = IniBead, FinBead

     do i = 1, N
       R(:,i) = Rpi(:,i,j)
     end do

     call GetForce(1,0)
     call GetForce(2,0)
     call GetForce(3,0)

     do i = 1, N

       F_fast(:,i,j) = ( Frc_Bond (:,i) &
       &               + Frc_Angle(:,i) &
       &               + Frc_UB   (:,i) &
       &               + Frc_Dihed(:,i) &
       &               + Frc_Impro(:,i) &
       &               + Frc_OptC (:,i) ) * InvP

       F_mode(:,i,j) = ( Frc_EAM  (:,i) &
       &               + Frc_Ersp (:,i) ) * InvP

       F_slow(:,i,j) =   Frc_Eksp(:,i) * InvP

     end do

     Virial = Virial + Vir_Bond  + Vir_Angle + Vir_UB   &
     &               + Vir_Dihed + Vir_Impro + Vir_OptC &
     &               + Vir_Ersp  + Vir_EAM   + Vir_Eksp

   end do

   call GetNormalModeForce( F_fast, Fnm_fast )
   call GetNormalModeForce( F_mode, Fnm_mode )
   call GetNormalModeForce( F_slow, Fnm_slow )

   call SumFrcPI( Fnm_fast )
   call SumFrcPI( Fnm_mode )
   call SumFrcPI( Fnm_slow )

   call SumVir( Virial )

   if(QCorrectCutoff) then

     Virial_co = 0.d0

     Virial_co(1,1) = CorrectV / (3.d0*Volume)
     Virial_co(2,2) = Virial_co(1,1)
     Virial_co(3,3) = Virial_co(1,1)

     Ene_LJ_co = CorrectE / Volume

   else

     Virial_co = 0.d0
     Ene_LJ_co = 0.d0

   end if

   Virial = Virial*InvP - Virial_co

! ### start MD time evolution ###

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

! ## Master - Slave
     if((QMaster).and.(mod(istep-1,lk) == 0)) then
! ## baro & thermostat
!     -----------------------------------------
       call BaroThermostatPIcent(lk)
!     -----------------------------------------
     end if

! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)
     if(QMasterPI) then

       if(mod(istep-1,lk) == 0) then

         do j = IniBead, FinBead

         do i = 1 , N

           Vnm(:,i,j) = Vnm(:,i,j) + ( dtk2 * Fnm_slow(:,i,j) &
           &                         + dtp2 * Fnm_mode(:,i,j) &
           &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

         end do

         end do

       else if(mod(istep-1,lp) == 0) then

! ## update the particle velocities

         do j = IniBead, FinBead

         do i = 1 , N

           Vnm(:,i,j) = Vnm(:,i,j) + ( dtp2 * Fnm_mode(:,i,j) &
           &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

         end do

         end do

       else

! ## update the particle velocities
         do j = IniBead, FinBead

         do i = 1 , N

           Vnm(:,i,j) = Vnm(:,i,j) + dt2 * Fnm_fast(:,i,j) * InvFictMass(i,j)

         end do

         end do

       end if

! ## Multiple time step

       do iref_step = 1, Ncref
!       ------------------------
         call ThermostatPInonc
!       ------------------------
! -------------------------------
! update the particle positions
! -------------------------------
         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2ref * Fnm_ref(:,i,j) * InvFictMass(i,j)

           end do

         end do

         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Rnm(:,i,j) = Rnm(:,i,j) + deltat_ref * Vnm(:,i,j)

           end do

         end do

!         if(QOpFix) call AddConstR

         call GetForce_Ref

         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2ref * Fnm_ref(:,i,j) * InvFictMass(i,j)

           end do

         end do

!         if(QOpFix) call AddConstV(istep)

         call ThermostatPInonc

       end do

     end if

     if(QMaster) then

       if( ( cBarostatMethod == 'PR' ) .or. &
       &   ( cBarostatMethod == 'ST' ) ) then

! ## update the particle positions
!       -------------------------------------
         call Jacobi(Vg,eigenVec,eigenValue) ! diagonalize Vg matrix
!       -------------------------------------

         do i = 1 , 3

           cf     = exp( dt2 * eigenValue(i) )
           fc1(i) = cf * cf

           arg2   = ( eigenValue(i) * dt2 ) * ( eigenValue(i) * dt2 )
           poly   = ((( e9 * arg2 + e7 ) * arg2 + e5 ) * arg2 + e3) * arg2 + 1.d0

           fc2(i) = cf * poly * deltat

         end do

         TeigenVec = Transpose(eigenVec)

         do i = 1 , N

           tempor = matmul( TeigenVec, Rnm(:,i,1) )      ! cg^t*r
           tempov = matmul( TeigenVec, Vnm(:,i,1) )    ! cg^t*v
           tempor = tempor * fc1 + tempov * fc2      ! Ie*cg^t*r + Is*cg^t*v*dt

           Rnm(:,i,1) = matmul( eigenVec, tempor )       ! cg*(Ie*cg^t*r + Is*cg^t*v*dt)

         end do

! ## update H
         tempoH = matmul( TeigenVec, H )      ! cg^t*H

         do i = 1 , 3

           tempoH(:,i) = tempoH(:,i) * fc1    ! Ie*cg^t*H

         end do

         H = matmul( eigenVec, tempoH )       ! cg*Ie*cg^t*H

       else if( ( cBarostatMethod == 'A3'  ).or.& ! Anisotropic Andersen
       &        ( cBarostatMethod == 'A2' ) ) then

         do i = 1 , 3

           aa(i)  = exp( dt2 * Vg(i,i) )

         end do

         aa2  = aa * aa

         do i = 1 , 3

           arg3(i) = ( Vg(i,i) * dt2 ) * ( Vg(i,i) * dt2 )

         end do

         poly3 = ((( e9 * arg3 + e7 ) * arg3 + e5 ) * arg3 + e3) * arg3 + 1.d0
         aa3   = aa * poly3
         bb3   = aa3 * deltat

         do i = 1, N
           Rnm(:,i,1) = Rnm(:,i,1) * aa2 + Vnm(:,i,1) * bb3
         end do

! ## update H
         do i = 1 , 3
           H(i,i) = H(i,i) * exp( Vg(i,i) * deltat )
         end do

       else if( cBarostatMethod == 'AN' ) then

! ## update the particle positions
         Anaa   = exp( dt2 * Vg(1,1) )
         Anaa2  = Anaa * Anaa

         arg2 = ( Vg(1,1) * dt2 ) * ( Vg(1,1) * dt2 )
         poly = ((( e9 * arg2 + e7 ) * arg2 + e5 ) * arg2 + e3) * arg2 + 1.d0

         bb   = Anaa * poly * deltat

         do i = 1 , N

           Rnm(:,i,1) = Rnm(:,i,1) * Anaa2 + Vnm(:,i,1) * bb

         end do

! ## update H
         ClSc = exp(Vg(1,1) * deltat)
         do i = 1 , 3
           H(i,i) = H(i,i) * ClSc
         end do

       end if

       if( mod(istep,BookFreq) == 0 ) then
         call InversMatrix(H,InvH)
! ## periodic boundary condition
         call PBC_PI
       end if

       do i = 1, N
         R(:,i) = Rnm(:,i,1)
       end do

     end if

     call BcastRH
     Rcentroid = R

     Volume = det(H)

     call InversMatrix(H,InvH)
     call TransCellList

     if(QCorrectCutoff) then

       Virial_co = 0.d0

       Virial_co(1,1) = CorrectV / (3.d0*Volume)
       Virial_co(2,2) = Virial_co(1,1)
       Virial_co(3,3) = Virial_co(1,1)

       Ene_LJ_co = CorrectE / Volume

     end if

     call NormalModeTransform

     if(NHam/=0) call Rot_FixedPoint

! ## get the new force
!   ---------------------------------------------
! ## Book keeping 
     if( mod(istep,BookFreq) == 0 ) then
       R = Rcentroid
       call PairList
     end if

! ## dt_short
     E_fast = 0.d0
     V_fast = 0.d0

     do j = IniBead, FinBead

       do i = 1, N
         R(:,i) = Rpi(:,i,j)
       end do

       call GetForce(1,0)

       do i = 1, N
         F_fast(:,i,j) = ( Frc_Bond (:,i) &
       &                 + Frc_Angle(:,i) &
       &                 + Frc_UB   (:,i) &
       &                 + Frc_Dihed(:,i) &
       &                 + Frc_Impro(:,i) &
       &                 + Frc_OptC (:,i) ) * InvP

       end do

       E_fast = E_fast + Ene_Bond + Ene_Angle + Ene_UB &
       &      + Ene_Dihed + Ene_Impro + Ene_OptC

       V_fast = V_fast + Vir_Bond  + Vir_Angle + Vir_UB &
       &      + Vir_Dihed + Vir_Impro + Vir_OptC

     end do

     call GetNormalModeForce( F_fast, Fnm_fast )

     call SumFrcPI( Fnm_fast )

!   ---------------------------------------------
! Multiple time step

     if(mod(istep,lp) == 0) then

       E_mode = 0.d0
       V_mode = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call GetForce(2,0)

         do i = 1, N

           F_mode(:,i,j) = ( Frc_EAM  (:,i) &
         &                 + Frc_Ersp (:,i) ) * InvP

         end do

         E_mode = E_mode + Ene_LJ + Ene_EAM + Ene_Ersp

         V_mode = V_mode + Vir_Ersp + Vir_EAM

       end do

       call GetNormalModeForce( F_mode, Fnm_mode )

       call SumFrcPI( Fnm_mode )

     end if

! ------------------------------------------------------

     if(mod(istep,lk) == 0) then

       E_slow = 0.d0
       V_slow = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call GetForce(3,0)

         do i = 1, N

           F_slow(:,i,j) = Frc_Eksp(:,i) * InvP

         end do

         E_slow = E_slow + Ene_Eksp

         V_slow = V_slow + Vir_Eksp

       end do

       call GetNormalModeForce( F_slow, Fnm_slow )

       call SumFrcPI( Fnm_slow )

     end if

!   ---------------------------------------------

     if(QMasterPI) then
! ## update the particle velocities

       if(mod(istep,lk) == 0) then

         do j = IniBead, FinBead

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + ( dtk2 * Fnm_slow(:,i,j) &
             &                         + dtp2 * Fnm_mode(:,i,j) &
             &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

           end do

         end do

       else if(mod(istep,lp) == 0) then

         do j = IniBead, FinBead

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + ( dtp2 * Fnm_mode(:,i,j) &
             &                         + dt2  * Fnm_fast(:,i,j) ) * InvFictMass(i,j)

           end do

         end do

       else

         do j = IniBead, FinBead

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2  * Fnm_fast(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if

     end if

! ## update the virial
     if(mod(istep,lk) == 0) then

       Virial = V_fast + V_mode + V_slow

       call SumVir( Virial )

       if(QMaster) then

         Virial = Virial*InvP - Virial_co

!       --------------------------------
         call BaroThermostatPIcent(lk)
!       --------------------------------

       end if

     end if

!   - store parameters --------------------------

     if( mod(istep,lk) == 0 ) then
       EnePI  = E_fast + E_mode + E_slow
       call Print_Energy_NP_PI(istep)
     end if

     if( mod(istep,ixc) == 0 ) call Print_Config_PI
     if( mod(istep,ixv) == 0 ) call Print_Velocity_PI
     if( mod(istep,irs) == 0 ) call SaveParam
!   ---------------------------------------------

! ## check the cell strain
     if(mod(istep,lk) == 0) then

       call CheckCellShape

     end if

   end do

end subroutine IntegrEOM_NPT_PI


!######################################################################
!######################################################################


! ***********************************************************
! ** #######  Isolated System - NE or NT Ensemble  ######  **
! ** #######    <<<< Constraint Dynamics >>>>      ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine IntegrEOM_iso_PI

use Numbers, only : N
use CommonBlocks, only : QMaster, QHeatCap, Qdebug, QDelCellMove
use Configuration, only : R
use CommonPI
use OptConstraintParam, only : NHam, Rrot, RIni, Frc_OptC, Ene_OptC
use NonbondParam, only : Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : &
&   Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TimeParam, only : Nstep, ixc, ixv, lp, lk, isampleHC, Timeps, &
&   deltat, dt2, irs

implicit none

integer :: i, j, istep, iref_step
real(8) :: ww
real(8), dimension(3,N,Nbead) :: F_fast, Fnm_fast
real(8), dimension(3,N,Nbead) :: F_mode, Fnm_mode
real(8), dimension(3,N,Nbead) :: F_slow, Fnm_slow
real(8) :: E_fast, E_mode, E_slow
real(8) :: dt2ref

real(8), dimension(3,N,Nbead) :: FrcPI
real(8) :: CVtrace, Vtrace
real(8) :: HessCVTrace, HessTrace
real(8), dimension(3) :: Ftmp, Rtmp, Rtt
integer :: ioutput

   dt2ref = 0.5d0 * deltat_ref

   if(QHeatCap) then
     ioutput = 1000
     ioutput = ioutput/isampleHC
   end if

   if(NHam/=0) Rrot = RIni

!   if(QMaster.and.QOpFix) then
!     call ConstPrepare
!   end if

   call GetForce_Ref

   call NormalModeTransform

   F_fast = 0.d0
   F_mode = 0.d0
   F_slow = 0.d0

   do j = IniBead, FinBead

     do i = 1, N

       R(:,i) = Rpi(:,i,j)

     end do

     call GetForceIso(0)

     do i = 1, N

       F_fast(:,i,j) = ( Frc_Bond (:,i) &
       &               + Frc_Angle(:,i) &
       &               + Frc_UB   (:,i) &
       &               + Frc_OptC (:,i) ) * InvP

       F_mode(:,i,j) = ( Frc_Dihed(:,i) &
       &               + Frc_Impro(:,i) ) * InvP

       F_slow(:,i,j) = ( Frc_Ersp (:,i) ) * InvP

     end do

   end do

   call GetNormalModeForce( F_fast, Fnm_fast )
   call GetNormalModeForce( F_mode, Fnm_mode )
   call GetNormalModeForce( F_slow, Fnm_slow )

   call SumFrcPI( Fnm_fast )
   call SumFrcPI( Fnm_mode )
   call SumFrcPI( Fnm_slow )

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

! ## Multiple Time Scale
! ---------------------------------------------------
     if((QMaster).and.(mod(istep-1,lk) == 0)) then
!    ------------------------------
       call ThermostatPIcent(lk)
!    ------------------------------
     end if

! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)

     if(QMasterPI) then

       if(mod(istep-1,lk) == 0) then

         ww = dt2 * lk

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ww * Fnm_slow(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if

! ------------------------------------------------------

! ## Multiple Time Scale
! ------------------------------------------------------
       if(mod(istep-1,lp) == 0) then

! ## update the particle velocities
         ww = dt2 * lp

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ww * Fnm_mode(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if
! ------------------------------------------------------

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       do j = IniBead, FinBead

         do i = 1 , N

           Vnm(:,i,j) = Vnm(:,i,j) + dt2 * Fnm_fast(:,i,j) * InvFictMass(i,j)

         end do

       end do

! ## Multiple time step

       do iref_step = 1, Ncref

         call ThermostatPInonc

! -------------------------------
! update the particle positions
! -------------------------------
         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2ref * Fnm_ref(:,i,j) * InvFictMass(i,j)

           end do

         end do

         do j = IniBead, FinBead

           do i = 1, N

             Rnm(:,i,j) = Rnm(:,i,j) + deltat_ref*Vnm(:,i,j)

           end do

         end do

!         if(QOpFix) call AddConstR

         call GetForce_Ref

         do j = IniBead, FinBead

           if(j==1) cycle

           do i = 1, N

             Vnm(:,i,j) = Vnm(:,i,j) + dt2ref * Fnm_ref(:,i,j) * InvFictMass(i,j)

           end do

         end do

!         if(QOpFix) call AddConstV(istep)

         call ThermostatPInonc

       end do

     end if

     call NormalModeTransform

! ## get the new force
! ------------------------------------------------
     E_fast = 0.d0

     do j = IniBead, FinBead

       do i = 1, N
         R(:,i) = Rpi(:,i,j)
       end do

       call GetForceIso(1)

       do i = 1 , N

          F_fast(:,i,j) = ( Frc_Bond (:,i) &
        &                 + Frc_Angle(:,i) &
        &                 + Frc_UB   (:,i) &
        &                 + Frc_OptC (:,i) ) * InvP

       end do

       E_fast = E_fast + Ene_Bond  + Ene_Angle + Ene_UB &
       &      + Ene_OptC

     end do

     call GetNormalModeForce( F_fast, Fnm_fast )

     call SumFrcPI( Fnm_fast )

!   ---------------------------------------------

     if(QMasterPI) then
! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       do j = IniBead, FinBead

         do i = 1 , N

           Vnm(:,i,j) = Vnm(:,i,j) + dt2 * Fnm_fast(:,i,j) * InvFictMass(i,j)

         end do

       end do

     end if

! ## Multiple time step
! ------------------------------------------------------
     if(mod(istep,lp) == 0) then

       E_mode = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call GetForceIso(2)

         do i = 1 , N

           F_mode(:,i,j) = ( Frc_Dihed(:,i) &
           &               + Frc_Impro(:,i) ) * InvP

         end do

         E_mode = E_mode + Ene_Dihed + Ene_Impro

       end do

       call GetNormalModeForce( F_mode, Fnm_mode )

       call SumFrcPI( Fnm_mode )

       if(QMasterPI) then
! ## update the particle velocities
         ww = dt2 * lp

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ww * Fnm_mode(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if

     end if
! ------------------------------------------------------

! ## Multiple time step
! ------------------------------------------------------
     if(mod(istep,lk) == 0) then

       E_slow = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call GetForceIso(3)

         do i = 1 , N

           F_slow(:,i,j) = ( Frc_Ersp (:,i) ) * InvP

         end do

         E_slow = E_slow + Ene_LJ + Ene_Ersp

       end do

       call GetNormalModeForce( F_slow, Fnm_slow )

       call SumFrcPI( Fnm_slow )

       if(QMasterPI) then

         ww = dt2 * lk

         do j = IniBead, FinBead

           do i = 1 , N

             Vnm(:,i,j) = Vnm(:,i,j) + ww * Fnm_slow(:,i,j) * InvFictMass(i,j)

           end do

         end do

       end if

     end if
! ------------------------------------------------------

     if(QMaster) then
! ------------------
! ## bond constraint
! ------------------
!       if(QOpFix) call AddConstV(istep)

       if(mod(istep,lk) == 0) then

! ## thermostat velocities and thermostat positions
!    ------------------------------
         call ThermostatPIcent(lk)
!    ------------------------------

! ## remove the cell-momentum
         if(QDelCellMove) call Elim_CellMove

       end if

     end if

!   - save parameters ---------------------------------------

     if(QHeatCap.and.(mod(istep,isampleHC) == 0)) then

       FrcPI = F_fast + F_mode + F_slow

       call SumFrcPI( FrcPI )

       EnePI = E_fast + E_mode + E_slow
       EnePI = EnePI * InvP

       if(QMaster) then
         do i = 1, N
           R(:,i) = Rnm(:,i,1)
         end do
       end if
       call BcastR
       if(.not.QMaster) then
         do i = 1, N
           Rnm(:,i,1) = R(:,i)
         end do
       end if

       CVtrace = 0.d0
       Vtrace  = 0.d0

       if(QMasterPI) then

         do j = IniBead, FinBead

           do i = 1, N

             Ftmp(:) = - FrcPI(:,i,j)
             Rtmp(:) = Rpi(:,i,j) - Rnm(:,i,1)

             CVtrace = CVtrace + Rtmp(1)*Ftmp(1) + Rtmp(2)*Ftmp(2) + Rtmp(3)*Ftmp(3)

             Rtt(:) = Rpi(:,i,j)

             Vtrace = Vtrace + Rtt(1)*Ftmp(1) + Rtt(2)*Ftmp(2) + Rtt(3)*Ftmp(3)

           end do

         end do

       end if

       HessCVTrace = 0.d0
       HessTrace = 0.d0

       do j = IniBead, FinBead

         do i = 1, N
           R(:,i) = Rpi(:,i,j)
         end do

         call Calc_Hessian(HessCVTrace,HessTrace)

       end do

       HessTrace   = HessTrace   * InvP * 0.5d0
       HessCVTrace = HessCVTrace * InvP * 0.5d0

       call SumHC(EnePI, CVtrace, Vtrace, HessCVTrace, HessTrace)

       if(QMaster.and.Qdebug) then
       open(90,file='virial.check',status='unknown')

       write(90,'(d23.15)') EnePI
       write(90,'(d23.15)') CVtrace
       write(90,'(d23.15)') Vtrace
       write(90,'(d23.15)') HessCVTrace
       write(90,'(d23.15)') HessTrace

       close(90)
       end if

       if(QMaster) then
         call HeatCapacity(EnePI, CVtrace, Vtrace, HessCVTrace, HessTrace, istep/lk, ioutput)
       end if

     end if

     if( mod(istep,lk) == 0 ) then
       EnePI = E_fast + E_mode + E_slow
       call Print_Energy_iso_PI(istep)
     end if

     if( mod(istep,ixc) == 0 ) call Print_Config_PI
     if( mod(istep,ixv) == 0 ) call Print_Velocity_PI
     if( mod(istep,irs) == 0 ) call SaveParam
!   ----------------------------------------------------------

   end do

end subroutine IntegrEOM_iso_PI
