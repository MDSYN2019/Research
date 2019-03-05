

!######################################################################
!######################################################################


subroutine SetDynaLib(iflag)

use Numbers, only : N
use CommonBlocks, only : QPBC, QInitial, QBarostat, QThermostat, &
&   cThermostatMethod, cBarostatMethod
use BathParam
use Configuration
use QMDynamics
use UnitExParam, only : Avogadro
use CellParam, only : H
use AtomParam, only : Mass, InvMass
use TimeParam, only : Timeps

implicit none

integer :: iflag, i, j
real(8) :: Gauss, sv
external Gauss


if(iflag == 1) then

   open(1,file='mass.data',form='unformatted',status='old')

   read(1) N

   allocate( Mass(N) )
   allocate( InvMass(N) )

   allocate( R(3,N) )
   allocate( Vel(3,N) )
   allocate( FrcQM(3,N) )

   read(1) Mass

   close(1)

   Mass = Mass * 1.d-3 / Avogadro

   InvMass = 1.d0 / Mass

else if(iflag == 2) then

   open(7,file='restart.dat',form='unformatted',status='old')

   do i = 1, N
     read(7) R(:,i)
   end do

   if(QInitial) then

     call Gene_Velocity

   else

     do i = 1, N
       read(7) Vel(:,i)
     end do

   end if

   if(QPBC) then

     read(7) H

   end if

   if(QBarostat) then

     if(QInitial) then

     call Gene_Vg

     else

     read(7) Vg

     end if

   end if

   if(QThermostat) then

     if(QInitial) then

       if((cThermostatMethod == 'NHC').or. &
       &  (cThermostatMethod == 'NH')) then
         Rss = 0.d0

         do i = 1, NHchain
           sv = sqrt( kT * InvMts(i) )
           Vss(i) = sv * Gauss()
         end do

       else if(cThermostatMethod == 'MNHC') then
         RMNHC = 0.d0

         do i = 1, NumMNHC
           do j = 1, NHchain
             sv = sqrt( kT * InvMMNHC(j,i) )
             VMNHC(j,i) = sv * Gauss()
           end do
         end do

       end if

     else

       if((cThermostatMethod == 'NHC').or. &
       &  (cThermostatMethod == 'NH')) then

         do i = 1 , NHchain
           read(7) Rss(i), Vss(i)
         end do

       else if(cThermostatMethod == 'MNHC') then

         do i = 1 , NumMNHC

           do j = 1 , NHchain

             read(7) RMNHC(j,i), VMNHC(j,i)

           end do

         end do

       end if

     end if

   end if

   if(QInitial) then

     Timeps = 0.d0

   else

     read(7) timeps

   end if

   close(7)

end if

Contains

   subroutine Gene_Vg

   implicit none

   real(8) :: pref

     if( ( cBarostatMethod == 'PR' ).or. &
     &   ( cBarostatMethod == 'ST' ) ) then

       do i = 1 , 3
         do j = i , 3

           if(i==j) then
             pref = 1.d0
           else
             pref = 0.5d0
           end if

           Vg(i,j) = sqrt( pref * kT / Mp(i,j) ) * Gauss()

         end do
       end do

       Vg(2,1) = Vg(1,2)
       Vg(3,1) = Vg(1,3)
       Vg(3,2) = Vg(2,3)

     else if( cBarostatMethod == 'A3' ) then

       Vg = 0.d0
       do i = 1 , 3
         Vg(i,i) = sqrt( kT / Mp(i,i) ) * Gauss()
       end do

     else if( cBarostatMethod == 'A2' ) then

       Vg = 0.d0
       do i = 1 , 3
         Vg(i,i) = sqrt( kT / Mp(i,i) ) * Gauss()
       end do
       Vg(2,2) = Vg(1,1)

     else if( cBarostatMethod == 'AN' ) then

       Vg = 0.d0
       Vg(1,1) = sqrt( kT / Mp(1,1) ) * Gauss()

     end if

   end subroutine Gene_Vg

end subroutine SetDynaLib


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

subroutine DynaLib_NV

use Numbers, only : N
use CommonBlocks, only : QMaster, QSHAKE, QThermostat
use Configuration
use QMDynamics
use SHAKEparam, only : R_o
use OptConstraintParam, only : NHam
use AtomParam, only : InvMass
use TimeParam, only : Nstep, Timeps, deltat, dt2

implicit none

integer :: i

real(8), dimension(3,N) :: Acc


   if(Icurrent==0) then

     call CalcTemp

     if(NHam/=0) call Rot_FixedPoint

     do i = 1 , N

       Acc(:,i) = FrcQM(:,i) * InvMass(i)

     end do

!   - save parameters ---------------------------------------
     call Print_Energy_DynaLib_NV
!   ---------------------------------------------------------

     Icurrent = Icurrent + 1
     Timeps = Timeps + deltat

! ## Master - Slave
     if(QMaster) then

! ## thermostat velocities and thermostat positions
       if(QThermostat) call Thermostat(1,1)

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       do i = 1 , N
         Vel(:,i) = Vel(:,i) + dt2 * Acc(:,i)
       end do

! -------------------------------
! update the particle positions
! -------------------------------
       if(QSHAKE) R_o = R
       R   = R + Vel * deltat

!     -----------------------
       if(QSHAKE) call SHAKE
!     -----------------------

! ## periodic boundary condition
!     ----------
       call PBC
!     ----------

     end if

   else if(Icurrent==Nstep) then

     do i = 1 , N
       Acc(:,i) = FrcQM(:,i) * InvMass(i)
     end do

! ## Master - Slave
     if(QMaster) then

       do i = 1, N
         Vel(:,i) = Vel(:,i) + dt2 * Acc(:,i)
       end do

!     ------------------------
       if(QSHAKE) call RATTLE
!     ------------------------

! ## thermostat velocities and thermostat positions
       if(QThermostat) call Thermostat(1,2)

! ## remove the cell-momentum

     end if
! ------------------------------------------------------

!   - save parameters ---------------------------------------
     call Print_Energy_DynaLib_NV
!   ---------------------------------------------------------

   else

     do i = 1 , N

       Acc(:,i) = FrcQM(:,i) * InvMass(i)

     end do

! ## Master - Slave
     if(QMaster) then

       do i = 1, N
         Vel(:,i) = Vel(:,i) + dt2 * Acc(:,i)
       end do

!     ------------------------
       if(QSHAKE) call RATTLE
!     ------------------------

! ## thermostat velocities and thermostat positions
       if(QThermostat) call Thermostat(1,2)

! ## remove the cell-momentum

     end if
! ------------------------------------------------------

!   - save parameters ---------------------------------------
     call Print_Energy_DynaLib_NV
!   ---------------------------------------------------------

     Icurrent = Icurrent + 1
     Timeps = Timeps + deltat

! ## Master - Slave
     if(QMaster) then

! ## thermostat velocities and thermostat positions
       if(QThermostat) call Thermostat(1,1)

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       Vel = Vel + dt2 * Acc

! -------------------------------
! update the particle positions
! -------------------------------
       if(QSHAKE) R_o = R
       R   = R + Vel * deltat

!     -----------------------
       if(QSHAKE) call SHAKE
!     -----------------------

! ## periodic boundary condition
!     ----------
       call PBC
!     ----------

     end if

   end if

end subroutine DynaLib_NV


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

subroutine DynaLib_NP

use Numbers, only : N
use Configuration
use CommonBlocks, only : QMaster, QThermostat, cThermostatMethod, &
&   cBarostatMethod
use QMDynamics
use ThermoData, only : Ene_kin
use OptConstraintParam, only : NHam
use BathParam, only : Vg, gkT, Baro_kin
use CellParam, only : H, InvH, Volume
use AtomParam, only : InvMass
use TimeParam, only : Nstep, Timeps, deltat, dt2

implicit none

real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8), dimension(3,3) :: tempoH, Dummy33
real(8), dimension(3)   :: eigenValue, Dummy3
real(8), dimension(3)   :: fc1, fc2
real(8), dimension(3)   :: tempor, tempov
integer :: i
real(8) :: cf, arg2
real(8) :: poly, Dummy
real(8) :: bb
real(8) :: Anaa, Anaa2
real(8), dimension(3) :: aa, aa2, aa3
real(8), dimension(3) :: arg3, poly3
real(8), dimension(3) :: bb3

real(8), dimension(3,N) :: Acc

real(8) :: VelScale

real(8) :: det
External det


   if( Icurrent == 0 ) then

     call CalcTemp

     if(NHam/=0) call Rot_FixedPoint

     call Read_Force

     do i = 1 , N
       Acc(:,i) = FrcQM(:,i) * InvMass(i)
     end do

     call Print_Energy_DynaLib_NP

! ### start MD time evolution ###

     Timeps = Timeps + deltat
     Icurrent = Icurrent + 1

! ## Master - Slave
     if(QMaster) then

! ## baro & thermostat
       if( QThermostat ) then

         if( cThermostatMethod == 'NH' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_NH(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_NH(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_NH(Dummy,1)
           end if

         else if( cThermostatMethod == 'NHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN(Dummy,1)
           end if

         else if( cThermostatMethod == 'MNHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_MNHC(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_MNHC(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_MNHC(Dummy,1)
           end if

         else if( cThermostatMethod == 'VSCALE' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,1)
           end if

         end if

       else

         if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
           call BarostatPR(Dummy33,1)
         else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
           call BarostatA3(Dummy3,1)
         else if( cBarostatMethod == 'AN' ) then
           call BarostatAN(Dummy,1)
         end if

       end if

!       -----------------------------------------
! ## update the particle velocities
       Vel = Vel + dt2 * Acc

       if( cBarostatMethod == 'PR' ) then

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

           tempor = matmul( TeigenVec, R(:,i) )      ! cg^t*r
           tempov = matmul( TeigenVec, Vel(:,i) )    ! cg^t*v
           tempor = tempor * fc1 + tempov * fc2      ! Ie*cg^t*r + Is*cg^t*v*dt

           R(:,i) = matmul( eigenVec, tempor )       ! cg*(Ie*cg^t*r + Is*cg^t*v*dt)

         end do

! ## update H
         tempoH = matmul( TeigenVec, H )      ! cg^t*H

         do i = 1 , 3

           tempoH(:,i) = tempoH(:,i) * fc1    ! Ie*cg^t*H

         end do

         H = matmul( eigenVec, tempoH )       ! cg*Ie*cg^t*H

       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then

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
           R(:,i) = R(:,i) * aa2 + Vel(:,i) * bb3
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
           R(:,i) = R(:,i) * Anaa2 + Vel(:,i) * bb
         end do

! ## update H
         H(1,1) = H(1,1) * exp(Vg(1,1) * deltat)
         H(2,2) = H(1,1)
         H(3,3) = H(1,1)

       end if

       call InversMatrix(H,InvH)

! ## periodic boundary condition
!     ----------
       call PBC
!     ----------

     end if

     Volume = det(H)

     if(NHam/=0) call Rot_FixedPoint

   else if(Icurrent == Nstep) then

     call Read_Force

     do i = 1 , N

       Acc(:,i) = FrcQM(:,i) * InvMass(i)

     end do

!   ---------------------------------------------

     if(QMaster) then
! ## update the particle velocities

       do i = 1 , N
         Vel(:,i) = Vel(:,i) + dt2 * Acc(:,i)
       end do

!       ------------------------------------------
       if( QThermostat ) then

         if( cThermostatMethod == 'NH' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_NH(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_NH(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_NH(Dummy,1)
           end if

         else if( cThermostatMethod == 'NHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN(Dummy,1)
           end if

         else if( cThermostatMethod == 'MNHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_MNHC(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_MNHC(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_MNHC(Dummy,1)
           end if

         else if( cThermostatMethod == 'VSCALE' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,1)
           end if

           call CalcTemp
           call BathTemp
           VelScale = sqrt( gkT / (Ene_kin + Baro_kin) )

           do i = 1 , N
             Vel(:,i) = Vel(:,i) * VelScale
           end do

           Vg = Vg * VelScale

         end if

       else

         if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
           call BarostatPR(Dummy33,1)
         else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
           call BarostatA3(Dummy3,1)
         else if( cBarostatMethod == 'AN' ) then
           call BarostatAN(Dummy,1)
         end if

       end if
!       -------------------------------------------

     end if

!   - store parameters --------------------------
     call Print_Energy_DynaLib_NP
!   ---------------------------------------------

! ## check the cell strain
     if(mod(Icurrent,1) == 0) then

       call CheckCellShape

     end if

   else

     call Read_Force

     do i = 1 , N
       Acc(:,i) = FrcQM(:,i) * InvMass(i)
     end do

!   ---------------------------------------------

     if(QMaster) then
! ## update the particle velocities

       do i = 1 , N
         Vel(:,i) = Vel(:,i) + dt2 * Acc(:,i)
       end do

!       ------------------------------------------
       if( QThermostat ) then

         if( cThermostatMethod == 'NH' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_NH(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_NH(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_NH(Dummy,1)
           end if

         else if( cThermostatMethod == 'NHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN(Dummy,1)
           end if

         else if( cThermostatMethod == 'MNHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_MNHC(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_MNHC(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_MNHC(Dummy,1)
           end if

         else if( cThermostatMethod == 'VSCALE' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,1)
           end if

           call CalcTemp
           call BathTemp
           VelScale = sqrt( gkT / (Ene_kin + Baro_kin) )

           do i = 1 , N
             Vel(:,i) = Vel(:,i) * VelScale
           end do

           Vg = Vg * VelScale

         end if

       else

         if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
           call BarostatPR(Dummy33,1)
         else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
           call BarostatA3(Dummy3,1)
         else if( cBarostatMethod == 'AN' ) then
           call BarostatAN(Dummy,1)
         end if

       end if
!       -------------------------------------------

     end if

!   - store parameters --------------------------
     call Print_Energy_DynaLib_NP
!   ---------------------------------------------

! ## check the cell strain
     if(mod(Icurrent,1) == 0) then

       call CheckCellShape

     end if

     Timeps = Timeps + deltat
     Icurrent = Icurrent + 1

! ## Master - Slave
     if(QMaster) then

! ## baro & thermostat
       if( QThermostat ) then

         if( cThermostatMethod == 'NH' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_NH(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_NH(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_NH(Dummy,1)
           end if

         else if( cThermostatMethod == 'NHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN(Dummy,1)
           end if

         else if( cThermostatMethod == 'MNHC' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BaroThermostatPR_MNHC(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BaroThermostatA3_MNHC(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BaroThermostatAN_MNHC(Dummy,1)
           end if

         else if( cThermostatMethod == 'VSCALE' ) then

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,1)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,1)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,1)
           end if

         end if

       else

         if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
           call BarostatPR(Dummy33,1)
         else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
           call BarostatA3(Dummy3,1)
         else if( cBarostatMethod == 'AN' ) then
           call BarostatAN(Dummy,1)
         end if

       end if

! -----------------------------------------
! ## update the particle velocities
       Vel = Vel + dt2 * Acc

       if( cBarostatMethod == 'PR' ) then

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

           tempor = matmul( TeigenVec, R(:,i) )      ! cg^t*r
           tempov = matmul( TeigenVec, Vel(:,i) )    ! cg^t*v
           tempor = tempor * fc1 + tempov * fc2      ! Ie*cg^t*r + Is*cg^t*v*dt

           R(:,i) = matmul( eigenVec, tempor )       ! cg*(Ie*cg^t*r + Is*cg^t*v*dt)

         end do

! ## update H
         tempoH = matmul( TeigenVec, H )      ! cg^t*H

         do i = 1 , 3

           tempoH(:,i) = tempoH(:,i) * fc1    ! Ie*cg^t*H

         end do

         H = matmul( eigenVec, tempoH )       ! cg*Ie*cg^t*H

       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then

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
           R(:,i) = R(:,i) * aa2 + Vel(:,i) * bb3
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
           R(:,i) = R(:,i) * Anaa2 + Vel(:,i) * bb
         end do

! ## update H
         H(1,1) = H(1,1) * exp(Vg(1,1) * deltat)
         H(2,2) = H(1,1)
         H(3,3) = H(1,1)

       end if

       call InversMatrix(H,InvH)

! ## periodic boundary condition
!     ----------
       call PBC
!     ----------

     end if

     Volume = det(H)

     if(NHam/=0) call Rot_FixedPoint

   end if


end subroutine DynaLib_NP


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

subroutine DynaLib_iso

use Numbers, only : N
use Configuration
use CommonBlocks, only : QMaster, QSHAKE, QThermostat
use QMDynamics
use SHAKEparam, only : R_o
use OptConstraintParam, only : NHam, Rrot, RIni
use AtomParam, only : InvMass
use TimeParam, only : Nstep, Timeps, deltat, dt2

implicit none

integer :: i

real(8), dimension(3,N) :: Acc


   if(Icurrent == 0) then

     if(NHam/=0) Rrot = RIni

     call CalcTemp

! ---------------------
     call Read_Force
! ---------------------

     do i = 1 , N

       Acc(:,i) = FrcQM(:,i) * InvMass(i)

     end do

!   - save parameters ---------------------------------------
     call Print_Energy_DynaLib_iso
!   ----------------------------------------------------------

     Icurrent = Icurrent + 1
     Timeps = Timeps + deltat

     if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
!    ----------------------------------------
       if(QThermostat) call Thermostat(1,1)
!    ----------------------------------------

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       Vel = Vel + dt2 * Acc

! -------------------------------
! update the particle positions
! -------------------------------
       if(QSHAKE) R_o = R
       R   = R + Vel * deltat
! ------------------
! ## bond constraint
! ------------------
       if(QSHAKE) call SHAKE

     end if

   else if(Icurrent == Nstep) then

! ## get the new force
! ------------------------------------------------
      call Read_Force

      do i = 1 , N

        Acc(:,i) = FrcQM(:,i) * InvMass(i)

      end do

!   ---------------------------------------------

     if(QMaster) then
! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       Vel = Vel + dt2 * Acc

     end if

     if(QMaster) then
! ------------------
! ## bond constraint
! ------------------
       if(QSHAKE) call RATTLE

! ## thermostat velocities and thermostat positions
       if(QThermostat) call Thermostat(1,2)

     end if

!   - save parameters ---------------------------------------
     call Print_Energy_DynaLib_iso
!   ----------------------------------------------------------

   else

! ## get the new force
! ------------------------------------------------
      call Read_Force

      do i = 1 , N

        Acc(:,i) = FrcQM(:,i) * InvMass(i)

      end do

!   ---------------------------------------------

     if(QMaster) then
! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       Vel = Vel + dt2 * Acc

! ------------------
! ## bond constraint
! ------------------
       if(QSHAKE) call RATTLE

! ## thermostat velocities and thermostat positions
       if(QThermostat) call Thermostat(1,2)

     end if

!   - save parameters ---------------------------------------
     call Print_Energy_DynaLib_iso
!   ----------------------------------------------------------

     Icurrent = Icurrent + 1
     Timeps = Timeps + deltat

     if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
!    ----------------------------------------
       if(QThermostat) call Thermostat(1,1)
!    ----------------------------------------

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
       Vel = Vel + dt2 * Acc

! -------------------------------
! update the particle positions
! -------------------------------
       if(QSHAKE) R_o = R
       R   = R + Vel * deltat
! ------------------
! ## bond constraint
! ------------------
       if(QSHAKE) call SHAKE

     end if

   end if

end subroutine DynaLib_iso


!######################################################################
!######################################################################


subroutine Read_Force

use Numbers, only : N
use CommonBlocks, only : QBarostat
use QMDynamics
use ThermoData, only : Virial
use TimeParam, only : Timeps

implicit none

integer :: i
real(8) :: timeck

   open(1,file='force.dat',form='unformatted',status='old')

   do i = 1, N
     read(1) FrcQM(:,i)
   end do

   if(QBarostat) then
     read(1) Virial
   end if

   read(1) timeck

   close(1)

   if(timeck /= Timeps) then
     write(*,*) 'ERROR'
     write(*,*) 'two files "restart.dat" and "force.dat" are inconsistent'
     write(*,*) 'they have different time data'
     write(*,*) Timeps
     write(*,*) timeck
     call Finalize
   end if

end subroutine Read_Force


!######################################################################
!######################################################################


subroutine Fin_DynaLib

use Numbers, only : N
use Configuration
use IOparam, only : DirectoryName
use BathParam
use CommonBlocks, only : QPBC, QBarostat, QThermostat, cThermostatMethod
use CellParam, only : H
use TimeParam, only : Timeps

implicit none

integer :: i, j

   open(7,file=trim(DirectoryName)//'restart.dat',form='unformatted',status='old')

   do i = 1, N
     write(7) R(:,i)
   end do

   do i = 1, N
     write(7) Vel(:,i)
   end do

   if(QPBC) then

     write(7) H

   end if

   if(QBarostat) then

     write(7) Vg

   end if

   if(QThermostat) then

     if((cThermostatMethod == 'NHC').or. &
     &  (cThermostatMethod == 'NH')) then

       do i = 1 , NHchain
         write(7) Rss(i), Vss(i)
       end do

     else if(cThermostatMethod == 'MNHC') then

       do i = 1 , NumMNHC

         do j = 1 , NHchain

           write(7) RMNHC(j,i), VMNHC(j,i)

         end do

       end do

     end if

   end if

   write(7) Timeps

   close(7)


end subroutine Fin_DynaLib
