! ############################
! ## SUBROUTINE LIST 
! ## -- HMC_NV 
! ## -- HMC_NP 
! ## -- HMC_iso 
! ## -- HMC_NV_RIGID 
! ## -- HMC_NP_RIGID 
! ## -- HMC_iso_RIGID 
! ## -- RandomVelocity 
! ## -- RandomVelocityRB 
! ## -- RandomVelocityBarostat 
! ############################


!######################################################################
!######################################################################


! ***********************************************************
! ** #######        HMC Canonical Ensemble         ######  **
! ** #######    <<<< Flexible chain model >>>>     ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine HMC_NV

use Numbers, only : N
use CommonBlocks, only : QMaster, QCorrectCutoff
use Configuration
use CommonHMC
use EAM_param, only : Ene_EAM
use BathParam, only : Beta
use EwaldParam, only : Ene_Eksp
use OptConstraintParam, only : NHam, Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Elec, Ene_Ersp
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : Volume
use TailCorrect
use TimeParam, only : Nstep, itgn, ixc, ixv, lp, lk, BookFreq, &
&   deltat, dt2, irs
use ThermoData, only : Ene_kin, Virial

implicit none

integer :: i, istep, MCstep
real(8) :: ww

real(8), dimension(3,N) :: A_fast, A_mode, A_slow
real(8), dimension(3,N) :: A_fast_OLD, A_mode_OLD, A_slow_OLD
real(8), dimension(3,N) :: R_OLD, Vel_OLD
real(8) :: Ham, Ham_OLD, E_NEW, E_OLD
real(8), dimension(3,3) :: Virial_OLD
real(8) :: Ene_Bond_OLD, Ene_Angle_OLD, Ene_UB_OLD, Ene_Dihed_OLD, Ene_Impro_OLD
real(8) :: Ene_OptC_OLD, Ene_LJ_OLD, Ene_Ersp_OLD, Ene_Eksp_OLD, Ene_EAM_OLD

real(8) :: ranf
external ranf

   if(NHam/=0) call Rot_FixedPoint

! -------------------
   call GetForce(0,0)
! -------------------

   call GetAcc( A_fast, 1 )
   call GetAcc( A_mode, 2 )
   call GetAcc( A_slow, 3 )

   call SumFrc( A_fast )
   call SumFrc( A_mode )
   call SumFrc( A_slow )
                    
   if(QCorrectCutoff) then

     Virial_co = 0.d0

     Virial_co(1,1) = CorrectV / (3.d0*Volume)
     Virial_co(2,2) = Virial_co(1,1)
     Virial_co(3,3) = Virial_co(1,1)

     Ene_LJ_co = CorrectE / Volume

   end if

   call GetEnergy_NV(E_NEW)

   NumAccept = 0

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do MCstep = 1 , Nstep

     TimeMC = TimeMC + 1

     if(QMaster) then

       R_OLD   = R
       Vel_OLD = Vel

       E_OLD   = E_NEW

       A_fast_OLD = A_fast
       A_mode_OLD = A_mode
       A_slow_OLD = A_slow

       Ene_Bond_OLD  = Ene_Bond  
       Ene_Angle_OLD = Ene_Angle 
       Ene_UB_OLD    = Ene_UB    
       Ene_Dihed_OLD = Ene_Dihed 
       Ene_Impro_OLD = Ene_Impro 
       Ene_OptC_OLD  = Ene_OptC  
       Ene_LJ_OLD    = Ene_LJ    
       Ene_Ersp_OLD  = Ene_Ersp  
       Ene_Eksp_OLD  = Ene_Eksp  
       Ene_EAM_OLD   = Ene_EAM

       Virial_OLD = Virial

       call RandomVelocity

       Ham_OLD = E_OLD + Ene_kin * 0.5d0

     end if

     do istep = 1, MDstep

! ## Master - Slave
       if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
         if(mod(istep-1,lk) == 0) then

! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)
           ww = dt2 * lk
           do i = 1 , N
             Vel(:,i) = Vel(:,i) + ww * A_slow(:,i)
           end do

         end if
! ------------------------------------------------------

! ## Multiple Time Scale
! ------------------------------------------------------
         if(mod(istep-1,lp) == 0) then

! ## update the particle velocities
           ww = dt2 * lp
           do i = 1 , N
             Vel(:,i) = Vel(:,i) + ww * A_mode(:,i)
           end do

         end if
! ------------------------------------------------------

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
         do i = 1 , N
           Vel(:,i) = Vel(:,i) + dt2 * A_fast(:,i)
         end do

! -------------------------------
! update the particle positions
! -------------------------------
         R = R + Vel * deltat

       end if

       call BcastR

! ## get the new force
!   ---------------------------------------------
       call GetForce(1,0)

       call GetAcc( A_fast, 1 )
       call SumFrc( A_fast )

!   ---------------------------------------------

       if(QMaster) then
! ----------------------------------
! ## update the particle velocities
! ----------------------------------
         do i = 1 , N
           Vel(:,i) = Vel(:,i) + dt2 * A_fast(:,i)
         end do

       end if

! ## Multiple time step
! ------------------------------------------------------
       if( mod(istep-1,BookFreq) == 0 ) call PairList

       if(mod(istep,lp) == 0) then

         call GetForce(2,0)

         call GetAcc( A_mode, 2 )
         call SumFrc( A_mode )

         if(QMaster) then
! ## update the particle velocities
           ww = dt2 * lp

           do i = 1 , N
             Vel(:,i) = Vel(:,i) + ww * A_mode(:,i)
           end do

         end if

       end if
! ------------------------------------------------------

! ## Multiple time step
! ------------------------------------------------------
       if(mod(istep,lk) == 0) then

         call GetForce(3,0)

         call GetAcc( A_slow, 3 )
         call SumFrc( A_slow )

         if(QMaster) then

           ww = dt2 * lk

           do i = 1 , N
             Vel(:,i) = Vel(:,i) + ww * A_slow(:,i)
           end do

         end if

       end if
! ------------------------------------------------------

     end do

! ## acceptance

     call GetEnergy_NV(E_NEW)

     if(QMaster) then

       call CalcTemp

       Ham = E_new + Ene_kin * 0.5d0

       if( Ham < Ham_OLD ) then

! ## accepted
         NumAccept = NumAccept + 1

       else if( ranf() < exp( - Beta * (Ham - Ham_OLD ) ) ) then

! ## accepted
         NumAccept = NumAccept + 1

       else

! ## rejected
         Ene_Bond  = Ene_Bond_OLD  
         Ene_Angle = Ene_Angle_OLD 
         Ene_UB    = Ene_UB_OLD    
         Ene_Dihed = Ene_Dihed_OLD 
         Ene_Impro = Ene_Impro_OLD 
         Ene_OptC  = Ene_OptC_OLD  
         Ene_LJ    = Ene_LJ_OLD    
         Ene_Ersp  = Ene_Ersp_OLD  
         Ene_Eksp  = Ene_Eksp_OLD  
         Ene_EAM   = Ene_EAM_OLD  

         Virial = Virial_OLD

         E_NEW  = E_OLD

         R      = R_OLD
         Vel    = Vel_OLD
         A_fast = A_fast_OLD
         A_mode = A_mode_OLD
         A_slow = A_slow_OLD

       end if

       if( mod(MCstep,itgn) == 0 ) call Print_Energy_HMC_NV(MCstep)

     end if
!   - save parameters ---------------------------------------
     if( mod(MCstep,ixc) == 0 ) call Print_Config
     if( mod(MCstep,ixv) == 0 ) call Print_Velocity
     if( mod(MCstep,irs) == 0 ) call SaveParam
!   ---------------------------------------------------------

   end do

end subroutine HMC_NV


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

subroutine HMC_NP

use Numbers, only : N
use CommonBlocks, only : QMaster, QCorrectCutoff, ForceField, &
&   cBarostatMethod
use Configuration
use CommonHMC
use EAM_param, only : Vir_EAM, Ene_EAM
use BathParam, only: Beta, Baro_kin, Vg
use EwaldParam, only : Vir_Eksp, Ene_Eksp
use OptConstraintParam, only : NHam, Vir_OptC, Ene_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong, &
&   Ene_LJ, Ene_Elec, Ene_Ersp
use BondedParam, only : &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : H, InvH, Volume
use TailCorrect
use TimeParam, only : Nstep, itgn, ixc, ixv, lp, lk, BookFreq, &
&   deltat, dt2, irs
use ThermoData, only : Ene_kin, Virial

implicit none

integer :: i, istep, MCstep
real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8), dimension(3,3) :: tempoH, Dummy33
real(8), dimension(3)   :: eigenValue, Dummy3
real(8), dimension(3)   :: fc1, fc2
real(8), dimension(3)   :: tempor, tempov
real(8) :: ww, cf, arg2
real(8) :: poly, Dummy
real(8) :: bb
real(8) :: Anaa, Anaa2
real(8), dimension(3) :: aa, aa2, aa3
real(8), dimension(3) :: arg3, poly3
real(8), dimension(3) :: bb3

real(8), dimension(3,N) :: A_fast
real(8), dimension(3,N) :: A_mode
real(8), dimension(3,N) :: A_slow

real(8), dimension(3,N) :: A_fast_OLD, A_mode_OLD, A_slow_OLD
real(8), dimension(3,N) :: R_OLD, Vel_OLD
real(8) :: Ham, Ham_OLD, E_NEW, E_OLD
real(8), dimension(3,3) :: Virial_OLD, Vg_OLD, H_OLD
real(8) :: Ene_Bond_OLD, Ene_Angle_OLD, Ene_UB_OLD, Ene_Dihed_OLD, Ene_Impro_OLD
real(8) :: Ene_OptC_OLD, Ene_LJ_OLD, Ene_Ersp_OLD, Ene_Eksp_OLD, Ene_EAM_OLD
real(8) :: Ene_LJ_co_OLD
real(8) :: Pcompress, Volume_OLD, ExComp
real(8) :: det, ranf
External det
external ranf

if(ForceField(1:3) == 'EAM') then
open(13,file='CellMatrix.dat',form='unformatted')
end if

   if(NHam/=0) call Rot_FixedPoint

! ---------------
   call GetForce(0,0)
! ---------------

   call GetAcc( A_fast, 1 )
   call GetAcc( A_mode, 2 )
   call GetAcc( A_slow, 3 )

   call SumFrc( A_fast )
   call SumFrc( A_mode )
   call SumFrc( A_slow )

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

   Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
   &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
   &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

   call SumVir( Virial )

   Virial = Virial - Virial_co

   call GetEnergy_NP(E_NEW)

   NumAccept = 0

! ### start MD time evolution ###

   do MCstep = 1 , Nstep

     TimeMC = TimeMC + 1

     if(QMaster) then

       R_OLD   = R
       Vel_OLD = Vel
       H_OLD   = H
       Vg_OLD  = Vg

       Volume_OLD = det(H)

       A_fast_OLD = A_fast
       A_mode_OLD = A_mode
       A_slow_OLD = A_slow

       Virial_OLD = Virial

       E_OLD = E_NEW

       Ene_Bond_OLD  = Ene_Bond  
       Ene_Angle_OLD = Ene_Angle 
       Ene_UB_OLD    = Ene_UB    
       Ene_Dihed_OLD = Ene_Dihed 
       Ene_Impro_OLD = Ene_Impro 
       Ene_OptC_OLD  = Ene_OptC  
       Ene_LJ_OLD    = Ene_LJ    
       Ene_Ersp_OLD  = Ene_Ersp  
       Ene_Eksp_OLD  = Ene_Eksp  
       Ene_LJ_co_OLD = Ene_LJ_co 
       Ene_EAM_OLD   = Ene_EAM

       call RandomVelocity
       call RandomVelocityBarostat

       Ham_OLD = E_OLD + ( Ene_kin + Baro_kin ) * 0.5d0

     end if

     do istep = 1, MDstep

! ## Master - Slave
       if(QMaster) then

! ## baro & thermostat
         if( mod(istep-1,lk) == 0 ) then
!       -----------------------------------------

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,lk)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,lk)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,lk)
           end if

!       -----------------------------------------
! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)

           ww = dt2 * lk

           do i = 1 , N
             Vel(:,i) = Vel(:,i) + ww * A_slow(:,i)
           end do

         end if

         if(mod(istep-1,lp) == 0) then

! ## update the particle velocities
           ww = dt2 * lp

           do i = 1 , N
             Vel(:,i) = Vel(:,i) + ww * A_mode(:,i)
           end do

         end if

! ## update the particle velocities
         do i = 1 , N
           Vel(:,i) = Vel(:,i) + dt2 * A_fast(:,i)
         end do

         if( ( cBarostatMethod == 'PR' ) .or. &
         &   ( cBarostatMethod == 'ST' ) ) then

! ## update the particle positions
!         -------------------------------------
           call Jacobi(Vg,eigenVec,eigenValue) ! diagonalize Vg matrix
!         -------------------------------------

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

         else if( ( cBarostatMethod == 'A3' ).or.&
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

       end if

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

       if(NHam/=0) call Rot_FixedPoint

! ## get the new force
!   ---------------------------------------------
       call GetForce(1,0)

       call GetAcc( A_fast, 1 )
       call SumFrc( A_fast )

!   ---------------------------------------------

       if(QMaster) then
! ## update the particle velocities

         do i = 1 , N

           Vel(:,i) = Vel(:,i) + dt2 * A_fast(:,i)

         end do

       end if

!   ---------------------------------------------
! Multiple time step
       if( mod(istep-1,BookFreq) == 0 ) call PairList

       if(mod(istep,lp) == 0) then

         call GetForce(2,0)

         call GetAcc( A_mode, 2 )
         call SumFrc( A_mode )

         if(QMaster) then
! ## update the particle velocities
           ww = dt2 * lp

           do i = 1 , N

             Vel(:,i) = Vel(:,i) + ww * A_mode(:,i)

           end do

         end if

       end if

! ------------------------------------------------------

       if(mod(istep,lk) == 0) then

         call GetForce(3,0)

         call GetAcc( A_slow, 3 )
         call SumFrc( A_slow )

         if(QMaster) then

           ww = dt2 * lk

           do i = 1 , N

             Vel(:,i) = Vel(:,i) + ww * A_slow(:,i)

           end do

         end if

! ## update the virial
         Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
         &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
         &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

         call SumVir( Virial )

         if(QMaster) then

           Virial = Virial - Virial_co

!       ------------------------------------------
           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,lk)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,lk)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,lk)
           end if
!       -------------------------------------------

         end if

       end if

     end do

     call GetEnergy_NP(E_NEW)

! ## acceptance test
     if(QMaster) then

       call CalcTemp
       call BathTemp

       Ham = E_NEW + ( Ene_kin + Baro_kin ) * 0.5d0

       if( cBarostatMethod == 'AN' ) then

         ExComp = Beta * ( Ham - Ham_OLD )

       else

         Pcompress = - 2.d0 * log( Volume / Volume_OLD )

         ExComp = Beta * ( Ham - Ham_OLD) + Pcompress

       end if

       if( ExComp < 0.d0 ) then

! ## accepted
         NumAccept = NumAccept + 1

       else if( ranf() < exp( - ExComp ) ) then

! ## accepted
         NumAccept = NumAccept + 1

       else

! ## rejected
         Ene_Bond  = Ene_Bond_OLD  
         Ene_Angle = Ene_Angle_OLD 
         Ene_UB    = Ene_UB_OLD    
         Ene_Dihed = Ene_Dihed_OLD 
         Ene_Impro = Ene_Impro_OLD 
         Ene_OptC  = Ene_OptC_OLD  
         Ene_LJ    = Ene_LJ_OLD    
         Ene_Ersp  = Ene_Ersp_OLD  
         Ene_Eksp  = Ene_Eksp_OLD  
         Ene_LJ_co = Ene_LJ_co_OLD 
         Ene_EAM   = Ene_EAM_OLD

         Virial = Virial_OLD

         E_NEW  = E_OLD

         R      = R_OLD
         Vel    = Vel_OLD
         H      = H_OLD
         Vg     = Vg_OLD

         A_fast = A_fast_OLD
         A_mode = A_mode_OLD
         A_slow = A_slow_OLD

       end if

       if( mod(MCstep,itgn) == 0 ) call Print_Energy_HMC_NP(MCstep)

     end if

!   - store parameters --------------------------
     if( mod(MCstep,ixc) == 0 ) call Print_Config
     if( mod(MCstep,ixv) == 0 ) call Print_Velocity
     if( mod(MCstep,irs) == 0 ) call SaveParam
!   ---------------------------------------------

! ## check the cell strain
     call CheckCellShape

   end do

end subroutine HMC_NP


!######################################################################
!######################################################################


! ***********************************************************
! ** #######  Isolated System - NE or NT Ensemble  ######  **
! ** #######    <<<< Flexible chain model >>>>     ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine HMC_iso

implicit none

   write(*,*) '>>>> HMC in the isolated system <<<<'
   write(*,*) 'This has not been supported in this program yet!'
   call Finalize

end subroutine HMC_iso


!######################################################################
!######################################################################


! ***********************************************************
! ** ####### Microcanonical or Canonical Ensemble  ######  **
! ** #######       <<<< Rigid-body model >>>>      ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine HMC_NV_RIGID

use Numbers, only : N
use CommonBlocks, only : QMaster, QCorrectCutoff
use Configuration
use CommonHMC
use EAM_param, only : Ene_EAM
use RBparam, only : NumRB, NumRBType, NumRBAtom, InvInertiaRB, QSingle, &
&   Quaternion, R_RB, V_RB, Lmoment, InvMassRB, QLinear, RBType
use BathParam, only: Beta
use EwaldParam, only : Ene_Eksp
use OptConstraintParam, only : NHam, Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Elec, Ene_Ersp
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use TailCorrect
use AtomParam, only : InvMass
use TimeParam, only : Nstep, itgn, ixc, ixv, lp, lk, BookFreq, &
&   deltat, dt2, irs
use ThermoData, only : Ene_kin, Virial
use CellParam, only : Volume

implicit none

integer :: i, k, istep, MCstep
real(8) :: ww

real(8), dimension(3,N) :: A_fast, A_mode, A_slow
real(8), dimension(3,NumRB) :: G_fast, G_mode, G_slow
real(8), dimension(3,NumRB) :: T_fast, T_mode, T_slow
real(8), dimension(NumRBType) :: rizy, rizx
real(8) :: qn2, qn, th1, th2, Lmomentztemp
real(8) :: snt, cst
real(8), dimension(3) :: Lt, Omt
real(8) :: omg, th, diag, offd, omg2
real(8), dimension(4) :: qq
integer :: MyType

real(8) :: Ham, Ham_OLD, E_NEW, E_OLD
real(8), dimension(3,NumRB) :: RB_OLD, VB_OLD, L_OLD
real(8), dimension(4,NumRB) :: Q_OLD
real(8), dimension(3,NumRB) :: G_fast_OLD, G_mode_OLD, G_slow_OLD
real(8), dimension(3,NumRB) :: T_fast_OLD, T_mode_OLD, T_slow_OLD
real(8), dimension(3,3) :: Virial_OLD
real(8) :: Ene_Bond_OLD, Ene_Angle_OLD, Ene_UB_OLD, Ene_Dihed_OLD, Ene_Impro_OLD
real(8) :: Ene_OptC_OLD, Ene_LJ_OLD, Ene_Ersp_OLD, Ene_Eksp_OLD, Ene_EAM_OLD

real(8) :: ranf
external ranf

   do i = 1 , NumRBType

     if(NumRBAtom(i)==1) cycle

     rizy(i) = dt2 * ( InvInertiaRB(3,i) - InvInertiaRB(2,i) )
     rizx(i) = dt2 * ( InvInertiaRB(3,i) - InvInertiaRB(1,i) )

   end do

   do i = 1 , NumRB

     if(QSingle(i)) cycle

     qn2 = dot_product( Quaternion(:,i), Quaternion(:,i) )
     qn  = 1.d0 / sqrt( qn2 )
     Quaternion(:,i) = Quaternion(:,i) * qn

   end do

   if(NHam/=0) call Rot_FixedPoint

! -------------------
   call IntraMolVec

   call GetForce(0,0)
! -------------------

   call GetAcc( A_fast, 1 )
   call GetAcc( A_mode, 2 )
   call GetAcc( A_slow, 3 )

   call SumFrc( A_fast )
   call SumFrc( A_mode )
   call SumFrc( A_slow )

   if(QMaster) then

     call Force_Div_Component(A_fast,G_fast,T_fast)
     call Force_Div_Component(A_mode,G_mode,T_mode)
     call Force_Div_Component(A_slow,G_slow,T_slow)

   end if

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

   call GetEnergy_NV(E_NEW)

   NumAccept = 0

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do MCstep = 1 , Nstep

     TimeMC = TimeMC + 1

     if(QMaster) then

       RB_OLD = R_RB
       VB_OLD = V_RB
       Q_OLD  = Quaternion
       L_OLD  = Lmoment

       G_fast_OLD = G_fast
       G_mode_OLD = G_mode
       G_slow_OLD = G_slow

       T_fast_OLD = T_fast
       T_mode_OLD = T_mode
       T_slow_OLD = T_slow

       Virial_OLD = Virial

       E_OLD    = E_NEW

       Ene_Bond_OLD  = Ene_Bond  
       Ene_Angle_OLD = Ene_Angle 
       Ene_UB_OLD    = Ene_UB    
       Ene_Dihed_OLD = Ene_Dihed 
       Ene_Impro_OLD = Ene_Impro 
       Ene_OptC_OLD  = Ene_OptC  
       Ene_LJ_OLD    = Ene_LJ    
       Ene_Ersp_OLD  = Ene_Ersp  
       Ene_Eksp_OLD  = Ene_Eksp  
       Ene_EAM_OLD   = Ene_EAM

       call RandomVelocityRB

       Ham_OLD = E_OLD + Ene_kin * 0.5d0

     end if

     do istep = 1, MDstep

! ## Master - Slave
       if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
         if(mod(istep-1,lk) == 0) then

! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)
           ww = dt2 * lk
           k = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB   (:,i) = V_RB   (:,i) + ww * G_slow(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_slow(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_slow(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_slow(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if
! ------------------------------------------------------

! ## Multiple Time Scale
! ------------------------------------------------------
         if(mod(istep-1,lp) == 0) then

! ## update the particle velocities
           ww = dt2 * lp
           k = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB   (:,i) = V_RB   (:,i) + ww * G_mode(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_mode(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_mode(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_mode(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if
! ------------------------------------------------------

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
         k = 0

         do i = 1 , NumRB

           if(QSingle(i)) then

             k = k + 1
             V_RB(:,i) = V_RB(:,i) + dt2 * G_fast(:,i) * InvMass(k)

           else

             MyType = RBType(i)
             V_RB   (:,i) = V_RB   (:,i) + dt2 * G_fast(:,i) * InvMassRB(MyType)

             if(QLinear(i)) then

               Lmoment(2,i) = Lmoment(2,i) + dt2 * T_fast(2,i)
               Lmoment(3,i) = Lmoment(3,i) + dt2 * T_fast(3,i)

             else

               Lmoment(:,i) = Lmoment(:,i) + dt2 * T_fast(:,i)

               Lt = Lmoment(:,i)

               th1 = rizy(MyType) * Lt(2)
               cst = dcos(th1)
               snt = dsin(th1)
               Lmoment(1,i) =  cst * Lt(1) + snt * Lt(3)
               Lmomentztemp = -snt * Lt(1) + cst * Lt(3)

               th2 = rizx(MyType) * Lmoment(1,i)
               cst = dcos(th2)
               snt = dsin(th2)
               Lmoment(2,i) =  cst * Lt(2) - snt * Lmomentztemp
               Lmoment(3,i) =  snt * Lt(2) + cst * Lmomentztemp

             end if

             k = k + NumRBAtom(MyType)

           end if

         end do

! ##### Omg~(t+dt/2) = I^-1 L~(t+dt/2)
! ##### Q(t+dt) = exp(dt/2*A[Omg~(t+dt/2)]) Q(t)

         do i = 1 , NumRB

           if(QSingle(i)) cycle

           MyType = RBType(i)

           qq = Quaternion(:,i)

           if(QLinear(i)) then
             Omt(1) = 0.d0
             Omt(2) = Lmoment(2,i) * InvInertiaRB(2,MyType)
             Omt(3) = Lmoment(3,i) * InvInertiaRB(3,MyType)
           else
             Omt = Lmoment(:,i) * InvInertiaRB(:,MyType)
           end if

           omg2 = dot_product( Omt, Omt )
           omg  = sqrt( omg2 )
           th   = dt2 * omg
           diag = dcos(th)
           offd = dsin(th) / omg

           Quaternion(1,i) = diag * qq(1) + offd * ( - Omt(1) * qq(2) &
           &               - Omt(2) * qq(3) - Omt(3) * qq(4) )
           Quaternion(2,i) = diag * qq(2) + offd * (   Omt(1) * qq(1) &
           &               - Omt(2) * qq(4) + Omt(3) * qq(3) )
           Quaternion(3,i) = diag * qq(3) + offd * (   Omt(1) * qq(4) &
           &               + Omt(2) * qq(1) - Omt(3) * qq(2) )
           Quaternion(4,i) = diag * qq(4) + offd * ( - Omt(1) * qq(3) &
           &               + Omt(2) * qq(2) + Omt(3) * qq(1) )

         end do

! -------------------------------
! update the particle positions
! -------------------------------
         R_RB = R_RB + V_RB * deltat

       end if

       call BcastRgQuat

! ## update all interacting sites
!     ------------------
       call IntraMolVec 
!     ------------------

! ## get the new force
!     ---------------------------------------------
       call GetForce(1,0)

       call GetAcc( A_fast, 1 )
       call SumFrc( A_fast )

!     ---------------------------------------------

       if(QMaster) then

         call Force_Div_Component(A_fast,G_fast,T_fast)

! ----------------------------------
! ## update the particle velocities
! ----------------------------------
         k = 0

         do i = 1 , NumRB

           if(QSingle(i)) then

             k = k + 1
             V_RB(:,i) = V_RB(:,i) + dt2 * G_fast(:,i) * InvMass(k)

           else

             MyType = RBType(i)

             if(QLinear(i)) then

               Lmoment(2,i) = Lmoment(2,i) + dt2 * T_fast(2,i)
               Lmoment(3,i) = Lmoment(3,i) + dt2 * T_fast(3,i)

             else

               Lt = Lmoment(:,i)

               th2 = rizx(MyType) * Lt(1)
               cst = cos(th2)
               snt = sin(th2)
               Lmoment(2,i) =  cst * Lt(2) - snt * Lt(3)
               Lmomentztemp =  snt * Lt(2) + cst * Lt(3)

               th1 = rizy(MyType) * Lmoment(2,i)
               cst = cos(th1)
               snt = sin(th1)
               Lmoment(1,i) =  cst * Lt(1) + snt * Lmomentztemp
               Lmoment(3,i) = -snt * Lt(1) + cst * Lmomentztemp

               Lmoment(:,i) = Lmoment(:,i) + dt2 * T_fast(:,i)

             end if

             V_RB   (:,i) = V_RB   (:,i) + dt2 * G_fast(:,i) * InvMassRB(MyType)

             k = k + NumRBAtom(MyType)

           end if

         end do

       end if

! ## Multiple time step
! ------------------------------------------------------
       if( mod(istep-1,BookFreq) == 0 ) call PairList

       if(mod(istep,lp) == 0) then

         call GetForce(2,0)

         call GetAcc( A_mode, 2 )
         call SumFrc( A_mode )

         if(QMaster) then

           call Force_Div_Component(A_mode,G_mode,T_mode)

! ## update the particle velocities
           ww = dt2 * lp
           k  = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_mode(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_mode(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_mode(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if

       end if

! ------------------------------------------------------

! ## Multiple time step
! ------------------------------------------------------
       if(mod(istep,lk) == 0) then

         call GetForce(3,0)

         call GetAcc( A_slow, 3 )
         call SumFrc( A_slow )

         if(QMaster) then

           call Force_Div_Component(A_slow,G_slow,T_slow)

           ww = dt2 * lk
           k = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_slow(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_slow(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_slow(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if

       end if

! ------------------------------------------------------

     end do

     call GetEnergy_NV(E_NEW)

! ## acceptance test
     if(QMaster) then

       call CalcTemp

       Ham = E_new + Ene_kin * 0.5d0

       if( Ham < Ham_OLD ) then

! ## accepted
       NumAccept = NumAccept + 1

       else if( ranf() < exp( - Beta * (Ham - Ham_OLD ) ) ) then

! ## accepted
         NumAccept = NumAccept + 1

       else

! ## rejected
         Ene_Bond  = Ene_Bond_OLD  
         Ene_Angle = Ene_Angle_OLD 
         Ene_UB    = Ene_UB_OLD    
         Ene_Dihed = Ene_Dihed_OLD 
         Ene_Impro = Ene_Impro_OLD 
         Ene_OptC  = Ene_OptC_OLD  
         Ene_LJ    = Ene_LJ_OLD    
         Ene_Ersp  = Ene_Ersp_OLD  
         Ene_Eksp  = Ene_Eksp_OLD  
         Ene_EAM   = Ene_EAM_OLD  

         Virial = Virial_OLD

         E_NEW  = E_OLD

         R_RB       = RB_OLD
         V_RB       = VB_OLD
         Lmoment    = L_OLD
         Quaternion = Q_OLD

         G_fast = G_fast_OLD
         G_mode = G_mode_OLD
         G_slow = G_slow_OLD

         T_fast = T_fast_OLD
         T_mode = T_mode_OLD
         T_slow = T_slow_OLD

       end if

       if( mod(MCstep,itgn) == 0 ) call Print_Energy_HMC_NV(MCstep)

     end if


!   - save parameters ---------------------------------------
     if( mod(MCstep,ixc) == 0 ) call Print_Config
     if( mod(MCstep,ixv) == 0 ) call Print_Velocity
     if( mod(MCstep,irs) == 0 ) call SaveParam
!   ---------------------------------------------------------

   end do

end subroutine HMC_NV_RIGID


!######################################################################
!######################################################################


! ***********************************************************
! ** #######     Isothermal-Isobaric Ensemble      ######  **
! ** #######       <<<< Rigid-body model >>>>      ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine HMC_NP_RIGID

use Numbers, only : N
use CommonBlocks, only : QMaster, QCorrectCutoff, cBarostatMethod, &
&   ForceField
use Configuration
use CommonHMC
use EAM_param, only : Vir_EAM, Ene_EAM
use RBparam, only : NumRB, NumRBType, NumRBAtom, InvInertiaRB, QSingle, &
&   Quaternion, R_RB, V_RB, Lmoment, InvMassRB, QLinear, RBType
use BathParam, only: Beta, Baro_kin, Vg
use EwaldParam, only : Vir_Eksp, Ene_Eksp
use OptConstraintParam, only : NHam, Vir_OptC, Ene_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong, &
&   Ene_LJ, Ene_Elec, Ene_Ersp
use BondedParam, only : &
&   Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro, &
&   Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : H, InvH, Volume
use TailCorrect
use AtomParam, only : InvMass
use TimeParam, only : Nstep, itgn, ixc, ixv, lp, lk, BookFreq, &
&   deltat, dt2, irs
use ThermoData, only : Ene_kin, Virial

implicit none

real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

integer :: i, k, istep, MCstep

real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8), dimension(3,3) :: tempoH, Dummy33
real(8), dimension(3) :: eigenValue, Dummy3
real(8), dimension(3) :: fc1, fc2
real(8), dimension(3) :: tempor, tempov
real(8), dimension(3) :: aa, aa2, aa3
real(8), dimension(3) :: arg3, poly3
real(8), dimension(3) :: bb3

real(8) :: ww, cf, arg2
real(8) :: poly, Dummy
real(8) :: bb
real(8) :: Anaa, Anaa2

real(8), dimension(3,N) :: A_fast
real(8), dimension(3,N) :: A_mode
real(8), dimension(3,N) :: A_slow

real(8), dimension(3,NumRB) :: G_fast, G_mode, G_slow
real(8), dimension(3,NumRB) :: T_fast, T_mode, T_slow

real(8), dimension(NumRBType) :: rizy, rizx
real(8), dimension(3) :: Lt, Omt
real(8), dimension(4) :: qq
real(8) :: qn2, qn, th1, th2, Lmomentztemp
real(8) :: snt, cst
real(8) :: omg, th, diag, offd, omg2
integer :: MyType

real(8) :: Ham, Ham_OLD, E_NEW, E_OLD
real(8), dimension(3,NumRB) :: RB_OLD, VB_OLD, L_OLD
real(8), dimension(4,NumRB) :: Q_OLD
real(8), dimension(3,NumRB) :: G_fast_OLD, G_mode_OLD, G_slow_OLD
real(8), dimension(3,NumRB) :: T_fast_OLD, T_mode_OLD, T_slow_OLD
real(8), dimension(3,3) :: Virial_OLD, Vg_OLD, H_OLD
real(8) :: Ene_Bond_OLD, Ene_Angle_OLD, Ene_UB_OLD, Ene_Dihed_OLD, Ene_Impro_OLD
real(8) :: Ene_OptC_OLD, Ene_LJ_OLD, Ene_Ersp_OLD, Ene_Eksp_OLD, Ene_EAM_OLD
real(8) :: Ene_LJ_co_OLD
real(8) :: det, ranf
real(8) :: Pcompress, Volume_OLD, ExComp
External det
external ranf

if(ForceField(1:3) == 'EAM') then
open(13,file='CellMatrix.dat',form='unformatted')
end if

   do i = 1 , NumRBType

     if(NumRBAtom(i)==1) cycle

     rizy(i) = dt2 * ( InvInertiaRB(3,i) - InvInertiaRB(2,i) )
     rizx(i) = dt2 * ( InvInertiaRB(3,i) - InvInertiaRB(1,i) )

   end do

   do i = 1 , NumRB

     if(QSingle(i)) cycle

     qn2 = dot_product( Quaternion(:,i), Quaternion(:,i) )
     qn  = 1.d0 / sqrt( qn2 )
     Quaternion(:,i) = Quaternion(:,i) * qn

   end do

   if(NHam/=0) call Rot_FixedPoint

! --------------------
   call IntraMolVec

   call GetForce(0,0)
! --------------------

   call GetAcc( A_fast, 1 )
   call GetAcc( A_mode, 2 )
   call GetAcc( A_slow, 3 )

   call SumFrc( A_fast )
   call SumFrc( A_mode )
   call SumFrc( A_slow )

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

   Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
   &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
   &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

   call SumVir( Virial )

   Virial = Virial - Virial_co

   if(QMaster) then

     call Force_Div_Component(A_fast,G_fast,T_fast)
     call Force_Div_Component(A_mode,G_mode,T_mode)
     call Force_Div_Component(A_slow,G_slow,T_slow)

   end if

   call GetEnergy_NP(E_NEW)

   NumAccept = 0

! ### start MD time evolution ###

   do MCstep = 1 , Nstep

     TimeMC = TimeMC + 1

     if(QMaster) then

       RB_OLD = R_RB
       VB_OLD = V_RB
       Q_OLD  = Quaternion
       L_OLD  = Lmoment
       H_OLD  = H
       Vg_OLD = Vg

       Volume_OLD = det(H)

       G_fast_OLD = G_fast
       G_mode_OLD = G_mode
       G_slow_OLD = G_slow

       T_fast_OLD = T_fast
       T_mode_OLD = T_mode
       T_slow_OLD = T_slow

       Virial_OLD = Virial

       E_OLD      = E_NEW

       Ene_Bond_OLD  = Ene_Bond  
       Ene_Angle_OLD = Ene_Angle 
       Ene_UB_OLD    = Ene_UB    
       Ene_Dihed_OLD = Ene_Dihed 
       Ene_Impro_OLD = Ene_Impro 
       Ene_OptC_OLD  = Ene_OptC  
       Ene_LJ_OLD    = Ene_LJ    
       Ene_Ersp_OLD  = Ene_Ersp  
       Ene_Eksp_OLD  = Ene_Eksp  
       Ene_LJ_co_OLD = Ene_LJ_co
       Ene_EAM_OLD   = Ene_EAM

       call RandomVelocityRB
       call RandomVelocityBarostat

       Ham_OLD = E_OLD + ( Ene_kin + Baro_kin ) * 0.5d0

     end if

     do istep = 1, MDstep

! ## Master - Slave
       if(QMaster) then

! ## baro & thermostat
         if( mod(istep-1,lk) == 0 ) then
!       -----------------------------------------

           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,lk)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,lk)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,lk)
           end if

!       -----------------------------------------
! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)

           ww = dt2 * lk
           k = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_slow(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_slow(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_slow(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if

         if(mod(istep-1,lp) == 0) then

! ## update the particle velocities
           ww = dt2 * lp

           k = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_mode(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_mode(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_mode(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if

! ## update the particle velocities
         k = 0

         do i = 1 , NumRB

           if(QSingle(i)) then

             k = k + 1
             V_RB(:,i) = V_RB(:,i) + dt2 * G_fast(:,i) * InvMass(k)

           else

             MyType = RBType(i)
             V_RB(:,i) = V_RB(:,i) + dt2 * G_fast(:,i) * InvMassRB(MyType)

             if(QLinear(i)) then

               Lmoment(2,i) = Lmoment(2,i) + dt2 * T_fast(2,i)
               Lmoment(3,i) = Lmoment(3,i) + dt2 * T_fast(3,i)

             else

               Lmoment(:,i) = Lmoment(:,i) + dt2 * T_fast(:,i)

               Lt = Lmoment(:,i)

               th1 = rizy(MyType) * Lt(2)
               cst = dcos(th1)
               snt = dsin(th1)
               Lmoment(1,i) =  cst * Lt(1) + snt * Lt(3)
               Lmomentztemp = -snt * Lt(1) + cst * Lt(3)

               th2 = rizx(MyType) * Lmoment(1,i)
               cst = dcos(th2)
               snt = dsin(th2)
               Lmoment(2,i) =  cst * Lt(2) - snt * Lmomentztemp
               Lmoment(3,i) =  snt * Lt(2) + cst * Lmomentztemp

             end if

             k = k + NumRBAtom(MyType)

           end if

         end do

! ##### Omg~(t+dt/2) = I^-1 L~(t+dt/2)
! ##### Q(t+dt) = exp(dt/2*A[Omg~(t+dt/2)]) Q(t)

         do i = 1 , NumRB

           if(QSingle(i)) cycle

           MyType = RBType(i)

           qq = Quaternion(:,i)

           if(QLinear(i)) then
             Omt(1) = 0.d0
             Omt(2) = Lmoment(2,i) * InvInertiaRB(2,MyType)
             Omt(3) = Lmoment(3,i) * InvInertiaRB(3,MyType)
           else
             Omt = Lmoment(:,i) * InvInertiaRB(:,MyType)
           end if

           omg2 = dot_product( Omt, Omt )
           omg  = sqrt( omg2 )
           th   = dt2 * omg
           diag = dcos(th)
           offd = dsin(th) / omg

           Quaternion(1,i) = diag * qq(1) + offd * ( - Omt(1) * qq(2) &
           &               - Omt(2) * qq(3) - Omt(3) * qq(4) )
           Quaternion(2,i) = diag * qq(2) + offd * (   Omt(1) * qq(1) &
           &               - Omt(2) * qq(4) + Omt(3) * qq(3) )
           Quaternion(3,i) = diag * qq(3) + offd * (   Omt(1) * qq(4) &
           &               + Omt(2) * qq(1) - Omt(3) * qq(2) )
           Quaternion(4,i) = diag * qq(4) + offd * ( - Omt(1) * qq(3) &
           &               + Omt(2) * qq(2) + Omt(3) * qq(1) )

         end do

         if( ( cBarostatMethod == 'PR' ) .or. &
         &   ( cBarostatMethod == 'ST' ) ) then

! ## update the particle positions
!         -------------------------------------
           call Jacobi(Vg,eigenVec,eigenValue) ! diagonalize Vg matrix
!         -------------------------------------

           do i = 1 , 3

             cf     = exp( dt2 * eigenValue(i) )
             fc1(i) = cf * cf

             arg2   = ( eigenValue(i) * dt2 ) * ( eigenValue(i) * dt2 )
             poly   = ((( e9 * arg2 + e7 ) * arg2 + e5 ) * arg2 + e3) * arg2 + 1.d0

             fc2(i) = cf * poly * deltat

           end do

           TeigenVec = Transpose(eigenVec)

           do i = 1 , NumRB

             tempor = matmul( TeigenVec, R_RB(:,i) )   ! cg^t*r
             tempov = matmul( TeigenVec, V_RB(:,i) )   ! cg^t*v
             tempor = tempor * fc1 + tempov * fc2      ! Ie*cg^t*r + Is*cg^t*v*dt

             R_RB(:,i) = matmul( eigenVec, tempor )    ! cg*(Ie*cg^t*r + Is*cg^t*v*dt)

           end do

! ## update H
           tempoH = matmul( TeigenVec, H )      ! cg^t*H

           do i = 1 , 3

             tempoH(:,i) = tempoH(:,i) * fc1    ! Ie*cg^t*H

           end do

           H = matmul( eigenVec, tempoH )       ! cg*Ie*cg^t*H

         else if( ( cBarostatMethod == 'A3' ).or.&
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

           do i = 1, NumRB

             R_RB(:,i) = R_RB(:,i) * aa2 + V_RB(:,i) * bb3

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

           do i = 1 , NumRB

             R_RB(:,i) = R_RB(:,i) * Anaa2 + V_RB(:,i) * bb

           end do

! ## update H
           H(1,1) = H(1,1) * exp(Vg(1,1) * deltat)
           H(2,2) = H(1,1)
           H(3,3) = H(1,1)

         end if

       end if

       call BcastRgQuatH

       Volume = det(H)

       call InversMatrix(H,InvH)
       call TransCellList

! ## update all interacting sites
!     ------------------
       call IntraMolVec 
!     ------------------

       if(QCorrectCutoff) then

         Virial_co = 0.d0

         Virial_co(1,1) = CorrectV / (3.d0*Volume)
         Virial_co(2,2) = Virial_co(1,1)
         Virial_co(3,3) = Virial_co(1,1)

         Ene_LJ_co = CorrectE / Volume

       end if

! ## get the new force
!   ---------------------------------------------
       call GetForce(1,0)

       call GetAcc( A_fast, 1 )
       call SumFrc( A_fast )

!   ---------------------------------------------

       if(QMaster) then

         call Force_Div_Component(A_fast,G_fast,T_fast)

! ## update the particle velocities

         k = 0

         do i = 1 , NumRB

           if(QSingle(i)) then

             k = k + 1
             V_RB(:,i) = V_RB(:,i) + dt2 * G_fast(:,i) * InvMass(k)

           else

             MyType = RBType(i)

             if(QLinear(i)) then

               Lmoment(2,i) = Lmoment(2,i) + dt2 * T_fast(2,i)
               Lmoment(3,i) = Lmoment(3,i) + dt2 * T_fast(3,i)

             else

               Lt = Lmoment(:,i)

               th2 = rizx(MyType) * Lt(1)
               cst = cos(th2)
               snt = sin(th2)
               Lmoment(2,i) =  cst * Lt(2) - snt * Lt(3)
               Lmomentztemp =  snt * Lt(2) + cst * Lt(3)

               th1 = rizy(MyType) * Lmoment(2,i)
               cst = cos(th1)
               snt = sin(th1)
               Lmoment(1,i) =  cst * Lt(1) + snt * Lmomentztemp
               Lmoment(3,i) = -snt * Lt(1) + cst * Lmomentztemp

               Lmoment(:,i) = Lmoment(:,i) + dt2 * T_fast(:,i)

             end if

             V_RB(:,i) = V_RB(:,i) + dt2 * G_fast(:,i) * InvMassRB(MyType)

             k = k + NumRBAtom(MyType)

           end if

         end do

       end if

!   ---------------------------------------------
! Multiple time step
       if( mod(istep-1,BookFreq) == 0 ) call PairList

       if(mod(istep,lp) == 0) then

         call GetForce(2,0)

         call GetAcc( A_mode, 2 )
         call SumFrc( A_mode )

         if(QMaster) then

           call Force_Div_Component(A_mode,G_mode,T_mode)

! ## update the particle velocities
           ww = dt2 * lp
           k  = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB(:,i) = V_RB(:,i) + ww * G_mode(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_mode(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_mode(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_mode(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if

       end if

! ------------------------------------------------------

       if(mod(istep,lk) == 0) then

         call GetForce(3,0)

         call GetAcc( A_slow, 3 )
         call SumFrc( A_slow )

         if(QMaster) then

           call Force_Div_Component(A_slow,G_slow,T_slow)

           ww = dt2 * lk
           k = 0

           do i = 1 , NumRB

             if(QSingle(i)) then

               k = k + 1
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMass(k)

             else

               MyType = RBType(i)
               V_RB(:,i) = V_RB(:,i) + ww * G_slow(:,i) * InvMassRB(MyType)

               if(QLinear(i)) then
                 Lmoment(2,i) = Lmoment(2,i) + ww * T_slow(2,i)
                 Lmoment(3,i) = Lmoment(3,i) + ww * T_slow(3,i)
               else
                 Lmoment(:,i) = Lmoment(:,i) + ww * T_slow(:,i)
               end if

               k = k + NumRBAtom(MyType)

             end if

           end do

         end if

! ## update the virial
         Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
         &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
         &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

         call SumVir( Virial )

         if(QMaster) then

           Virial = Virial - Virial_co

!         ------------------------------------------
           if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
             call BarostatPR(Dummy33,lk)
           else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
             call BarostatA3(Dummy3,lk)
           else if( cBarostatMethod == 'AN' ) then
             call BarostatAN(Dummy,lk)
           end if
!       -------------------------------------------

         end if

       end if

     end do

     call GetEnergy_NP(E_NEW)

! ## acceptance test
     if(QMaster) then

       call CalcTemp
       call BathTemp

       Ham = E_NEW + ( Ene_kin + Baro_kin ) * 0.5d0

       if( cBarostatMethod == 'AN' ) then

         ExComp = Beta * ( Ham - Ham_OLD )

       else

         Pcompress = - 2.d0 * log( Volume / Volume_OLD )

         ExComp = Beta * ( Ham - Ham_OLD ) + Pcompress

       end if

       if( ExComp < 0.d0 ) then

! ## accepted
         NumAccept = NumAccept + 1

       else if( ranf() < exp( - ExComp ) ) then

! ## accepted
         NumAccept = NumAccept + 1

       else

! ## rejected
         Ene_Bond  = Ene_Bond_OLD  
         Ene_Angle = Ene_Angle_OLD 
         Ene_UB    = Ene_UB_OLD    
         Ene_Dihed = Ene_Dihed_OLD 
         Ene_Impro = Ene_Impro_OLD 
         Ene_OptC  = Ene_OptC_OLD  
         Ene_LJ    = Ene_LJ_OLD    
         Ene_Ersp  = Ene_Ersp_OLD  
         Ene_Eksp  = Ene_Eksp_OLD  
         Ene_LJ_co = Ene_LJ_co_OLD
         Ene_EAM   = Ene_EAM_OLD

         Virial = Virial_OLD

         E_NEW  = E_OLD

         R_RB       = RB_OLD
         V_RB       = VB_OLD
         Lmoment    = L_OLD
         Quaternion = Q_OLD
         H          = H_OLD
         Vg         = Vg_OLD

         G_fast = G_fast_OLD
         G_mode = G_mode_OLD
         G_slow = G_slow_OLD

         T_fast = T_fast_OLD
         T_mode = T_mode_OLD
         T_slow = T_slow_OLD

       end if

       if( mod(MCstep,itgn) == 0 ) call Print_Energy_HMC_NP(MCstep)

     end if

!   - store parameters --------------------------
     if( mod(MCstep,ixc) == 0 ) call Print_Config
     if( mod(MCstep,ixv) == 0 ) call Print_Velocity
     if( mod(MCstep,irs) == 0 ) call SaveParam
!   ---------------------------------------------

! ## check the cell strain
     call CheckCellShape

   end do

end subroutine HMC_NP_RIGID


!######################################################################
!######################################################################


! ***********************************************************
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine HMC_iso_RIGID

implicit none

   write(*,*) '>>>> HMC in the isolated system <<<<'
   write(*,*) 'This has not been supported in this program yet!'
   call Finalize

end subroutine HMC_iso_RIGID


!######################################################################
!######################################################################


subroutine RandomVelocity

use Numbers, only : N
use Configuration, only : Vel
use CommonHMC

implicit none

real(8) :: snt, cst
real(8), dimension(3,N) :: Vel_OLD
real(8) :: ranf, x
external ranf

   if(QPartial) then

     Vel_OLD = Vel

     call Gene_Velocity

     snt = sin(ThetaPartial)
     cst = cos(ThetaPartial)

     x = ranf()

     if(x <= 0.5) then

       Vel = - Vel * snt - Vel_OLD * cst

     else

       Vel = Vel * snt + Vel_OLD * cst

     end if

   else

     call Gene_Velocity

   end if

! ----------------
   call CalcTemp
! ----------------

end subroutine RandomVelocity


!######################################################################
!######################################################################


subroutine RandomVelocityRB

use CommonHMC
use RBparam, only : NumRB, V_RB, Lmoment

implicit none

real(8) :: snt, cst
real(8), dimension(3,NumRB) :: V_RB_OLD
real(8), dimension(3,NumRB) :: Lmoment_OLD
real(8) :: ranf, x
external ranf


   if(QPartial) then

     V_RB_OLD    = V_RB
     Lmoment_OLD = Lmoment

     call Gene_Velocity

     snt = sin(ThetaPartial)
     cst = cos(ThetaPartial)

     x = ranf()

     if(x <= 0.5) then

       V_RB    = - V_RB    * snt - V_RB_OLD    * cst
       Lmoment = - Lmoment * snt - Lmoment_OLD * cst

     else

       V_RB    =   V_RB    * snt + V_RB_OLD    * cst
       Lmoment =   Lmoment * snt + Lmoment_OLD * cst

     end if

   else

     call Gene_Velocity

   end if

! ----------------
   call CalcTemp
! ----------------

end subroutine RandomVelocityRB


!######################################################################
!######################################################################


subroutine RandomVelocityBarostat

use CommonBlocks, only : cBarostatMethod
use CommonHMC
use BathParam, only : kT, Vg, Mp

implicit none

real(8) :: snt, cst
real(8), dimension(3,3) :: Vg_OLD
real(8) :: ranf, x
external ranf

   if(QPartial) then

     Vg_OLD = Vg

     call Gene_Vg

     snt = sin(ThetaPartial)
     cst = cos(ThetaPartial)

     x = ranf()

     if(x <= 0.5) then

       Vg = - Vg * snt - Vg_OLD * cst

     else

       Vg =   Vg * snt + Vg_OLD * cst

     end if

   else

     call Gene_Vg

   end if

! ---------------
   call BathTemp
! ---------------

Contains

   subroutine Gene_Vg

   implicit none

   integer :: i, j
   real(8) :: Gauss, pref
   external Gauss

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

end subroutine RandomVelocityBarostat
