! ############################
! ## SUBROUTINE LIST 
! ## -- IntegrEOM_NV 
! ## -- IntegrEOM_NP_SHAKE 
! ## -- IntegrEOM_NP 
! ## -- IntegrEOM_iso 
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

subroutine IntegrEOM_NV

use Numbers, only : N
use CommonBlocks, only : QMaster, QHeatCap, QOpFix, QJarzynski, QSHAKE, &
&   QThermostat, QCorrectCutoff, QVolScale, ForceField, QMacro, Qwcformol, &
&   QDelCellMove, QInert, QCyl, QFSCyl, QEflux
use Configuration
use EAM_param, only : Ene_EAM
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use SimAnneal, only : QSimAnneal
use SHAKEparam, only : R_o
use EwaldParam, only : Vir_Eksp, Ene_Eksp
use OptConstraintParam, only : NHam, Vir_OptC, Ene_OptC
use NonbondParam, only : Ene_LJ, Ene_Ersp, &
&   Ene_NBlong, Ene_ELlong, Ene_NBshrt, Ene_ELshrt
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use CellParam, only : Volume, Vsc_Rate
use TailCorrect
use AtomParam, only : Mass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, isampleHC, BookFreq, Timeps, &
&   deltat, dt2, irs
use wcparam, only : Nsample_wc, Ronc
use CGball, only : NumSphere

implicit none

integer :: i, istep
real(8) :: ww1, ww2, Dummy

real(8), dimension(3,N) :: A_fast
real(8), dimension(3,N) :: A_mode
real(8), dimension(3,N) :: A_slow

real(8) :: Upot, det
real(8) :: Asx, Asy, Asz, Amx, Amy, Amz
real(8) :: Afx, Afy, Afz, Acx, Acy, Acz
integer :: ioutput
external det

   if(QHeatCap) then
     ioutput = 1000
     ioutput = ioutput / isampleHC
   end if

   ww2 = dt2*lk
   ww1 = dt2*lp

   if(QVolScale) Vsc_Rate(:) = deltat * Vsc_Rate(:)
   if(NHam/=0) call Rot_FixedPoint
   if(QMaster.and.QOpFix) call ConstPrepare
   if(QMaster.and.QJarzynski) call SMD_pre
   if(QSimAnneal) call PreAnneal
   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

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
   else
     Virial_co = 0.d0
     Ene_LJ_co = 0.d0
   end if

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

     if(QSimAnneal) call Annealing(istep)

! ## Master - Slave
     if(QMaster) then

! ## Multiple time step

       if( mod(istep-1,lk) == 0 ) then
! ## thermostat velocities and thermostat positions
         if(QThermostat) call Thermostat(lk,1)
! ## update the particle velocities
         do i = 1, N
           Asx = A_slow(1,i)*ww2
           Asy = A_slow(2,i)*ww2
           Asz = A_slow(3,i)*ww2
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Asx + Amx + Afx
           Acy = Asy + Amy + Afy
           Acz = Asz + Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else if(mod(istep-1,lp) == 0) then
         do i = 1, N
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Amx + Afx
           Acy = Amy + Afy
           Acz = Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else
         do i = 1, N
           Acx = A_fast(1,i)*dt2
           Acy = A_fast(2,i)*dt2
           Acz = A_fast(3,i)*dt2
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       end if

! ## update the particle positions
       if(QSHAKE) then
         do i = 1, N
           R_o(1,i) = R(1,i)
           R_o(2,i) = R(2,i)
           R_o(3,i) = R(3,i)
         end do
       end if

       do i = 1, N
         R(1,i) = R(1,i) + Vel(1,i) * deltat
         R(2,i) = R(2,i) + Vel(2,i) * deltat
         R(3,i) = R(3,i) + Vel(3,i) * deltat
       end do
! ## bond constraint
       if(QSHAKE) call SHAKE

       if(QOpFix) call AddConstR
       if(QVolScale) call VolScRupdate

     end if

     if(QJarzynski) call SMD_reference

     if(QVolScale) then
       call VolScVcorrect
     else
       call BcastR
     end if

! ## get the new force
!   ---------------------------------------------

     if( mod(istep,BookFreq) == 0 ) call PairList

     call GetForce(1,istep)
     call GetAcc( A_fast, 1 )
     call SumFrc( A_fast )

!   ---------------------------------------------

! Multiple time step
     if(mod(istep,lp) == 0) then
       call GetForce(2,istep)
       call GetAcc( A_mode, 2 )
       call SumFrc( A_mode )
     end if

! ------------------------------------------------------
! Multiple time step
     if(mod(istep,lk) == 0) then
       call GetForce(3,istep)
       call GetAcc( A_slow, 3 )
       call SumFrc( A_slow )
     end if

! ## Master - Slave
     if(QMaster) then

       if(mod(istep,lk) == 0) then
         do i = 1, N
           Asx = A_slow(1,i)*ww2 
           Asy = A_slow(2,i)*ww2 
           Asz = A_slow(3,i)*ww2 
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Asx + Amx + Afx
           Acy = Asy + Amy + Afy
           Acz = Asz + Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else if(mod(istep,lp) == 0) then
         do i = 1, N
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Amx + Afx
           Acy = Amy + Afy
           Acz = Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else
         do i = 1, N
           Acx = A_fast(1,i)*dt2
           Acy = A_fast(2,i)*dt2
           Acz = A_fast(3,i)*dt2
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       end if

!     ------------------------
       if(QSHAKE) call RATTLE
!     ------------------------

! ## thermostat velocities and thermostat positions
       if(mod(istep,lk)==0.and.QThermostat) call Thermostat(lk,2)

     end if

     if(QOpFix) call AddConstV(istep)

! ## remove the cell-momentum
     if(QMaster.and.QDelCellMove) call Elim_CellMove

! ------------------------------------------------------

! >> F monitor ##
     if( mod(istep,isampleF) == 0 ) then
       call Force_gA
       call SumFrcgA
       if(QMaster) then
         do i = NiniF + 1 , NfinF
           Fext(:,i) = ( A_fast(:,i) + A_mode(:,i) + A_slow(:,i) ) &
           &           * Mass(i) - Fint(:,i)
         end do
         call Monitor_Force
       end if
     end if
! << F monitor ##

! >> Heat Capacity
     if( QHeatCap .and. (mod(istep, isampleHC) == 0) ) then
       if(ForceField(1:2)=='CG') then
         Upot = Ene_Bond + Ene_Angle + Ene_UB   + Ene_Dihed  + Ene_Impro  + Ene_NBshrt &
         &    + Ene_Eksp + Ene_EAM   + Ene_OptC + Ene_NBlong + Ene_ELshrt + Ene_ELlong
       else
         Upot = Ene_Bond   + Ene_Angle + Ene_UB + Ene_Dihed + Ene_Impro + Ene_NBshrt &
         &    + Ene_Ersp   + Ene_Eksp  + Ene_LJ + Ene_EAM   + Ene_OptC  + Ene_NBlong &
         &    + Ene_ELshrt + Ene_ELlong
       end if
       call SumHC( Upot, Dummy, Dummy, Dummy, Dummy )
       if(QMaster) then
         call HeatCapacity( Upot, Dummy, Dummy, Dummy, Dummy, istep/isampleHC, ioutput )
       end if
     end if
! << Heat Capacity

! >> Steered MD ##
     if(QJarzynski.and.(mod(istep,lk)==0)) then
       call SMD_sample(istep)
     end if
! << Steered MD ##
! >> PMF of macrosphere
     if(QMacro.and.(mod(istep,lk)==0).and.NumSphere==1) then
       call MacroPMFsample(istep)
     end if
! << PMF of macrosphere
     if(QMaster.and.Qwcformol.and.(mod(istep,Nsample_wc)==0)) write(73,'(f12.7)') Ronc
     if(QMaster.and.QInert.and.(mod(istep,lk)==0)) call CntIn_Sample(istep)
     if((QCyl.or.QFSCyl).and.(mod(istep,lk)==0)) call CylPMFsample(istep)
! >> Eflux
     if(QEflux.and.QMaster.and.(mod(istep,lk)==0)) call Print_Eflux(istep)
! << Eflux

!   - save parameters ---------------------------------------
     if( mod(istep,lk)  == 0 ) call Print_Energy_NV(istep)
     if( mod(istep,ixc) == 0 ) call Print_Config
     if( mod(istep,ixv) == 0 ) call Print_Velocity
     if( mod(istep,irs) == 0 ) call SaveParam
!   ---------------------------------------------------------

#ifdef EnergyRep
     if(mod(istep,lk)==0) call enganal(istep,1)
#endif

   end do

   if(QMaster.and.QSimAnneal) close(50)
   if(QMaster.and.Qwcformol) close(73)

end subroutine IntegrEOM_NV


!######################################################################
!######################################################################


! ***********************************************************
! ** #######     Isothermal-Isobaric Ensemble      ######  **
! ** #######    <<<< Constraint Dynamics >>>>      ######  **
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! **  XI-RESPA for multiple time step integration          **
! ***********************************************************

subroutine IntegrEOM_NP_SHAKE

use Numbers, only : N
use CommonBlocks, only : QMaster, QOpFix, QJarzynski, &
&   cBarostatMethod, QCorrectCutoff, ForceField, QMacro, &
&   Qwcformol, QDelCellMove, QInert, QCyl, QFSCyl, QEflux
use Configuration
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use SimAnneal, only : QSimAnneal
use SHAKEparam, only : R_o, Frc_Const, Vir_Const
use EwaldParam, only : Vir_Eksp
use BathParam, only : NHchain, Rss, Vss, Vg
use OptConstraintParam, only : NHam, Vir_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong
use BondedParam, only : Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro
use CellParam, only : H, InvH, Volume
use TailCorrect
use EAM_param, only : Vir_EAM
use AtomParam, only : Mass, InvMass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, BookFreq, Timeps, &
&   deltat, dt2, irs
use ThermoData, only : Virial
use wcparam, only : Nsample_wc, Ronc
use CGball, only : NumSphere

implicit none

logical :: ROLL

real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

real(8), dimension(3,3) :: Rv, tempoRot, tempoH
real(8), dimension(3,3) :: Imatrix
real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8), dimension(3)   :: eigenValue
real(8), dimension(3)   :: fc1, fc2, fc3
real(8), dimension(3)   :: tempor, tempov
integer                 :: i, istep
real(8)                 :: cf, arg2, poly, ClSc
real(8), dimension(3)   :: aa, aa2, aa3, arg3, poly3, bb3
real(8)                 :: Anaa, Anaa2, Anaa3, bb

real(8), dimension(3,N) :: A_fast
real(8), dimension(3,N) :: A_mode
real(8), dimension(3,N) :: A_slow

real(8), dimension(3,3) :: VirTemp, Vir_fast, Vir_mode, Vir_slow

real(8), dimension(NHchain) :: Rss_o,Vss_o
real(8), dimension(3,3) :: Vg_o
real(8), dimension(3,N) :: V_o

real(8) :: Asx, Asy, Asz, Amx, Amy, Amz
real(8) :: Afx, Afy, Afz, Acx, Acy, Acz
real(8) :: det, ww1, ww2, xx
External det

   if(ForceField(1:3) == 'EAM') then
     open(13,file='CellMatrix.dat',form='unformatted')
   end if

   ww2 = dt2*lk
   ww1 = dt2*lp

   Imatrix = 0.d0
   do i = 1, 3
     Imatrix(i,i) = 1.d0
   end do

   if(NHam/=0) call Rot_FixedPoint

   if(QMaster.and.QOpFix) then
     call ConstPrepare
   end if

! ## Steered MD >>
   if(QMaster.and.QJarzynski) then
     call SMD_pre
   end if

! ## Simulated Annealing >>
   if(QSimAnneal) call PreAnneal

   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

! -------------------------------
   call GetForce(0,0)
! -------------------------------
   call GetAcc( A_fast, 1 )
   call GetAcc( A_mode, 2 )
   call GetAcc( A_slow, 3 )

   Vir_fast = Vir_Bond  + Vir_Angle + Vir_UB   &
   &        + Vir_Dihed + Vir_Impro + Vir_OptC
   Vir_mode = Vir_Ersp + Vir_EAM + Vir_NBshrt
   Vir_slow = Vir_Eksp + Vir_NBlong

   if(QCorrectCutoff) then
     Virial_co = 0.d0
     Virial_co(1,1) = CorrectV / (3.d0*Volume)
     Virial_co(2,2) = Virial_co(1,1)
     Virial_co(3,3) = Virial_co(1,1)
     Ene_LJ_co = CorrectE / Volume
   end if

! ## 
   call SumFrc( A_fast )
   call SumFrc( A_mode )
   call SumFrc( A_slow )
   call SumVir( Vir_fast )
   call SumVir( Vir_mode )
   call SumVir( Vir_slow )

   VirTemp = Vir_fast + Vir_mode * lp + Vir_slow * lk

   if(QCorrectCutoff) then
     VirTemp = VirTemp - Virial_co
   end if

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ### start MD time evolution ###

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

! ## Simulated Annealing >>
     if(QSimAnneal) call Annealing(istep)

! ## Master - Slave
     if(QMaster) then

! #### SHAKE/ROLL procedure ####

       ROLL=.True.

       Vg_o  = Vg
       Rss_o = Rss
       Vss_o = Vss
       R_o   = R
       V_o   = Vel

       do while (ROLL)

         Vg  = Vg_o
         Rss = Rss_o
         Vss = Vss_o
         R   = R_o
         Vel = V_o

! ## the Constraint Force
         call Force_Constraint

         Virial = VirTemp + Vir_Const

! ## baro & thermostat
         call Bath_NP_SHAKE(ROLL,1)

! ## update the particle velocities
         if(mod(istep-1,lk) == 0) then
           do i = 1, N
             xx  = InvMass(i)
             Asx = A_slow(1,i)*ww2 
             Asy = A_slow(2,i)*ww2 
             Asz = A_slow(3,i)*ww2 
             Amx = A_mode(1,i)*ww1
             Amy = A_mode(2,i)*ww1
             Amz = A_mode(3,i)*ww1
             Afx = (A_fast(1,i)+Frc_Const(1,i)*xx)*dt2
             Afy = (A_fast(2,i)+Frc_Const(2,i)*xx)*dt2
             Afz = (A_fast(3,i)+Frc_Const(3,i)*xx)*dt2
             Acx = Asx + Amx + Afx
             Acy = Asy + Amy + Afy
             Acz = Asz + Amz + Afz
             Vel(1,i) = Vel(1,i) + Acx
             Vel(2,i) = Vel(2,i) + Acy
             Vel(3,i) = Vel(3,i) + Acz
           end do
         else if(mod(istep-1,lp) == 0) then
           do i = 1, N
             xx  = InvMass(i)
             Amx = A_mode(1,i)*ww1
             Amy = A_mode(2,i)*ww1
             Amz = A_mode(3,i)*ww1
             Afx = (A_fast(1,i)+Frc_Const(1,i)*xx)*dt2
             Afy = (A_fast(2,i)+Frc_Const(2,i)*xx)*dt2
             Afz = (A_fast(3,i)+Frc_Const(3,i)*xx)*dt2
             Acx = Amx + Afx
             Acy = Amy + Afy
             Acz = Amz + Afz
             Vel(1,i) = Vel(1,i) + Acx
             Vel(2,i) = Vel(2,i) + Acy
             Vel(3,i) = Vel(3,i) + Acz
           end do
         else
           do i = 1, N
             Acx = (A_fast(1,i)+Frc_Const(1,i)*xx)*dt2
             Acy = (A_fast(2,i)+Frc_Const(2,i)*xx)*dt2
             Acz = (A_fast(3,i)+Frc_Const(3,i)*xx)*dt2
             Vel(1,i) = Vel(1,i) + Acx
             Vel(2,i) = Vel(2,i) + Acy
             Vel(3,i) = Vel(3,i) + Acz
           end do
         end if

! ## update the particle positions
! ----------------------------
! ## constrained <<<SHAKE>>>
! ----------------------------
         if( ( cBarostatMethod == 'PR' ) .or. & ! Parrinello-Rahman
         &   ( cBarostatMethod == 'ST' ) ) then

           call Jacobi(Vg,eigenVec,eigenValue) ! diagonalize Vg matrix

           do i = 1 , 3
             cf     = exp( dt2 * eigenValue(i) )
             fc1(i) = cf * cf
             arg2   = ( eigenValue(i) * dt2 ) * ( eigenValue(i) * dt2 )
             poly   = ( ( ( e9 * arg2 + e7 ) * arg2 + e5 ) * arg2 + e3 ) &
             &        * arg2 + 1.d0
             fc3(i) = cf * poly
             fc2(i) = fc3(i) * deltat
           end do

           TeigenVec = Transpose(eigenVec)

           do i = 1 , N
             tempor = matmul( TeigenVec, R(:,i) )      ! cg^t*r
             tempov = matmul( TeigenVec, Vel(:,i) )    ! cg^t*v
             tempor = tempor * fc1 + tempov * fc2      ! Ie*cg^t*r + Is*cg^t*v*dt
             R(:,i) = matmul( eigenVec, tempor )       ! cg*(Ie*cg^t*r + Is*cg^t*v*dt)
           end do

           tempoRot = matmul( TeigenVec, Imatrix )

           do i = 1 , 3
             tempoRot(:,i) = tempoRot(:,i) * fc3
           end do

           Rv = matmul( eigenVec, tempoRot )

           call SHAKEROLLPR(Rv,ROLL)

         else if( ( cBarostatMethod == 'A3'  ).or.& ! Anisotropic Andersen
         &        ( cBarostatMethod == 'A2' ) ) then

           do i = 1 , 3
             aa(i)   = exp( dt2 * Vg(i,i) )
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

           call SHAKEROLLA3(aa3, ROLL)

         else if( cBarostatMethod == 'AN' ) then ! Isotropic Andersen

           Anaa   = exp( dt2 * Vg(1,1) )
           Anaa2  = Anaa * Anaa

           arg2 = ( Vg(1,1) * dt2 ) * ( Vg(1,1) * dt2 )
           poly = ((( e9 * arg2 + e7 ) * arg2 + e5 ) * arg2 + e3) * arg2 + 1.d0
           Anaa3  = Anaa * poly
           bb   = Anaa3 * deltat

           do i = 1 , N
             R(:,i) = R(:,i) * Anaa2 + Vel(:,i) * bb
           end do

           call SHAKEROLL(Anaa3, ROLL)

         end if

       end do    ! ## end of SHAKE/ROLL

! ## update H ( Cell Matrix )
       if( ( cBarostatMethod == 'PR' ) .or. &
       &   ( cBarostatMethod == 'ST' ) ) then

         tempoH = matmul( TeigenVec, H )      ! cg^t*H

         do i = 1 , 3
           tempoH(:,i) = tempoH(:,i) * fc1    ! Ie*cg^t*H
         end do

         H = matmul( eigenVec, tempoH )       ! cg*Ie*cg^t*H

       else if( ( cBarostatMethod == 'A3'  ).or.& ! Anisotropic Andersen
       &        ( cBarostatMethod == 'A2' ) ) then

         do i = 1 , 3
           H(i,i) = H(i,i) * exp( Vg(i,i) * deltat )
         end do

       else if( cBarostatMethod == 'AN' ) then

         ClSc   = exp( Vg(1,1) * deltat )
         do i = 1 , 3
           H(i,i) = H(i,i) * ClSc
         end do

       end if

       if(QOpFix) call AddConstR

     end if

     if(QJarzynski) call SMD_reference

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
     call GetForce(1,istep)
     call GetAcc( A_fast, 1 )
     Vir_fast = Vir_Bond  + Vir_Angle + Vir_UB   &
     &        + Vir_Dihed + Vir_Impro + Vir_OptC
     call SumFrc( A_fast )
     call SumVir( Vir_fast )
     VirTemp = Vir_fast

!   ---------------------------------------------

     if( mod(istep,BookFreq) == 0 ) call PairList

! Multiple time step
     if(mod(istep,lp) == 0) then
       call GetForce(2,istep)
       call GetAcc( A_mode, 2 )
       Vir_mode = Vir_Ersp + Vir_EAM  + Vir_NBshrt
       call SumFrc( A_mode )
       call SumVir( Vir_mode )
       VirTemp = Vir_fast + Vir_mode * lp
     end if

! ------------------------------------------------------

     if(mod(istep,lk) == 0) then
       call GetForce(3,istep)
       call GetAcc( A_slow, 3 )
       Vir_slow = Vir_Eksp + Vir_NBlong
       call SumFrc( A_slow )
       call SumVir( Vir_slow )
       VirTemp = Vir_fast + Vir_mode * lp + Vir_slow * lk
     end if
!   ---------------------------------------------

     if(QCorrectCutoff) then
       VirTemp = VirTemp - Virial_co
     end if

! ## Master - Slave
     if(QMaster) then

! #### RATTLE/ROLL procedure ####

       ROLL=.True.

       Vg_o  = Vg
       Rss_o = Rss
       Vss_o = Vss
       V_o   = Vel

       do while (ROLL)

         Vg  = Vg_o
         Rss = Rss_o
         Vss = Vss_o
         Vel = V_o

! ## the Constraint Force
         call Force_Constraint

         if(mod(istep,lk) == 0) then
           do i = 1, N
             xx  = InvMass(i)
             Asx = A_slow(1,i)*ww2 
             Asy = A_slow(2,i)*ww2 
             Asz = A_slow(3,i)*ww2 
             Amx = A_mode(1,i)*ww1
             Amy = A_mode(2,i)*ww1
             Amz = A_mode(3,i)*ww1
             Afx = (A_fast(1,i)+Frc_Const(1,i)*xx)*dt2
             Afy = (A_fast(2,i)+Frc_Const(2,i)*xx)*dt2
             Afz = (A_fast(3,i)+Frc_Const(3,i)*xx)*dt2
             Acx = Asx + Amx + Afx
             Acy = Asy + Amy + Afy
             Acz = Asz + Amz + Afz
             Vel(1,i) = Vel(1,i) + Acx
             Vel(2,i) = Vel(2,i) + Acy
             Vel(3,i) = Vel(3,i) + Acz
           end do
         else if(mod(istep,lp) == 0) then
           do i = 1, N
             xx  = InvMass(i)
             Amx = A_mode(1,i)*ww1
             Amy = A_mode(2,i)*ww1
             Amz = A_mode(3,i)*ww1
             Afx = (A_fast(1,i)+Frc_Const(1,i)*xx)*dt2
             Afy = (A_fast(2,i)+Frc_Const(2,i)*xx)*dt2
             Afz = (A_fast(3,i)+Frc_Const(3,i)*xx)*dt2
             Acx = Amx + Afx
             Acy = Amy + Afy
             Acz = Amz + Afz
             Vel(1,i) = Vel(1,i) + Acx
             Vel(2,i) = Vel(2,i) + Acy
             Vel(3,i) = Vel(3,i) + Acz
           end do
         else
           do i = 1, N
             Acx = (A_fast(1,i)+Frc_Const(1,i)*xx)*dt2
             Acy = (A_fast(2,i)+Frc_Const(2,i)*xx)*dt2
             Acz = (A_fast(3,i)+Frc_Const(3,i)*xx)*dt2
             Vel(1,i) = Vel(1,i) + Acx
             Vel(2,i) = Vel(2,i) + Acy
             Vel(3,i) = Vel(3,i) + Acz
           end do
         end if

! ## update the virial
         Virial = VirTemp + Vir_Const

! ## barostat & thermostat velocities and thermostat positions
! ## constraint <<<RATTLE>>>
         call Bath_NP_SHAKE(ROLL,2)

       end do   ! ## end of RATTLE/ROLL

     end if

     if(QOpFix) call AddConstV(istep)

! ## remove the cell-momentum
     if(QMaster.and.QDelCellMove) call Elim_CellMove

! >> F monitor ##
     if( mod(istep,isampleF) == 0 ) then
       call Force_gA
       call SumFrcgA
       if(QMaster) then
         do i = NiniF + 1 , NfinF
           Fext(:,i) = ( A_fast(:,i) + A_mode(:,i) + A_slow(:,i) ) &
           &           * Mass(i) - Fint(:,i)
         end do
         do i = NiniF + 1 , NfinF
           Fint(:,i) = Fint(:,i) + Frc_Const(:,i)
         end do
         call Monitor_Force
       end if
     end if
! << F monitor ##

! >> Steered MD ##
     if(QJarzynski.and.(mod(istep,lk)==0)) then
       call SMD_sample(istep)
     end if
! << Steered MD ##
! >> PMF of macrosphere
     if(QMacro.and.(mod(istep,lk)==0).and.NumSphere==1) then
       call MacroPMFsample(istep)
     end if
! << PMF of macrosphere
     if(QMaster.and.Qwcformol.and.(mod(istep,Nsample_wc)==0)) write(73,'(f12.7)') Ronc
     if(QMaster.and.QInert.and.(mod(istep,lk)==0)) call CntIn_Sample(istep)
     if((QCyl.or.QFSCyl).and.(mod(istep,lk)==0)) call CylPMFsample(istep)
! >> Eflux
     if(QEflux.and.QMaster.and.(mod(istep,lk)==0)) call Print_Eflux(istep)
! << Eflux

!   - store parameters --------------------------
     if( mod(istep,lk)  == 0 ) call Print_Energy_NP(istep)
     if( mod(istep,ixc) == 0 ) call Print_Config
     if( mod(istep,ixv) == 0 ) call Print_Velocity
     if( mod(istep,irs) == 0 ) call SaveParam
!   ---------------------------------------------

! ## check the cell strain
     call CheckCellShape

#ifdef EnergyRep
     if(mod(istep,lk)==0) call enganal(istep,2)
#endif

   end do

   if(QMaster.and.QSimAnneal) close(50)
   if(QMaster.and.Qwcformol) close(73)

end subroutine IntegrEOM_NP_SHAKE


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

subroutine IntegrEOM_NP

use Numbers, only : N
use Configuration
use CommonBlocks, only : QMaster, QOpFix, QJarzynski, &
&   cBarostatMethod, ForceField, QCorrectCutoff, QMacro, &
&   Qwcformol, QDelCellMove, QInert, QCyl, QFSCyl, QEflux
use EAM_param, only : Vir_EAM
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use SimAnneal, only : QSimAnneal
use BathParam, only : Vg
use EwaldParam, only : Vir_Eksp
use OptConstraintParam, only : NHam, Vir_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong
use BondedParam, only : Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro
use CellParam, only : H, InvH, Volume
use TailCorrect
use AtomParam, only : Mass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, BookFreq, Timeps, &
&   deltat, dt2, irs
use ThermoData, only : Virial
use wcparam, only : Nsample_wc, Ronc
use CGball, only : NumSphere

implicit none

real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

integer :: i, istep
real(8), dimension(3,3) :: eigenVec
real(8), dimension(3)   :: eigenValue
real(8) :: ww1, ww2
real(8), dimension(3,N) :: A_fast
real(8), dimension(3,N) :: A_mode
real(8), dimension(3,N) :: A_slow
real(8) :: Asx, Asy, Asz, Amx, Amy, Amz
real(8) :: Afx, Afy, Afz, Acx, Acy, Acz
real(8) :: aux, auy, auz
real(8) :: fc1x, fc1y, fc1z, fc2x, fc2y, fc2z
real(8) :: cfx, cfy, cfz, argx, argy, argz
real(8) :: ar9x, ar9y, ar9z
real(8) :: ar7x, ar7y, ar7z
real(8) :: ar5x, ar5y, ar5z
real(8) :: ar3x, ar3y, ar3z
real(8) :: polx, poly, polz
real(8) :: Texx, Texy, Texz
real(8) :: Teyx, Teyy, Teyz
real(8) :: Tezx, Tezy, Tezz
real(8) :: tempRx, tempRy, tempRz
real(8) :: tempVx, tempVy, tempVz
real(8) :: tempHxx, tempHxy, tempHxz
real(8) :: tempHyx, tempHyy, tempHyz
real(8) :: tempHzx, tempHzy, tempHzz
real(8) :: Rx, Ry, Rz
real(8) :: Vx, Vy, Vz
real(8) :: Hxx, Hxy, Hxz
real(8) :: Hyx, Hyy, Hyz
real(8) :: Hzx, Hzy, Hzz
real(8) :: aax, aay, aaz
real(8) :: ClSc
real(8) :: det
External det

   if(ForceField(1:3) == 'EAM') then
     open(13,file='CellMatrix.dat',form='unformatted')
   end if

   ww2 = dt2*lk
   ww1 = dt2*lp

   if(NHam/=0) call Rot_FixedPoint
   if(QMaster.and.QOpFix) call ConstPrepare
   if(QMaster.and.QJarzynski) call SMD_pre
   if(QSimAnneal) call PreAnneal
   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

!  -------------------
   call GetForce(0,0)
!  -------------------

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

! ### start MD time evolution ###

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

     if(QSimAnneal) call Annealing(istep)

! ## Master - Slave
     if(QMaster) then

       if( mod(istep-1,lk) == 0 ) then
! ## baro & thermostat
         call Bath_NP(1)
! ## update the particle velocities
         do i = 1, N
           Asx = A_slow(1,i)*ww2 
           Asy = A_slow(2,i)*ww2 
           Asz = A_slow(3,i)*ww2 
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Asx + Amx + Afx
           Acy = Asy + Amy + Afy
           Acz = Asz + Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else if(mod(istep-1,lp) == 0) then
         do i = 1, N
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Amx + Afx
           Acy = Amy + Afy
           Acz = Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else
         do i = 1, N
           Acx = A_fast(1,i)*dt2
           Acy = A_fast(2,i)*dt2
           Acz = A_fast(3,i)*dt2
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       end if

       if( ( cBarostatMethod == 'PR' ) .or. &
       &   ( cBarostatMethod == 'ST' ) ) then

! ## update the particle positions
!       -------------------------------------
         call Jacobi(Vg,eigenVec,eigenValue) ! diagonalize Vg matrix
!       -------------------------------------

         aux = dt2 * eigenValue(1)
         auy = dt2 * eigenValue(2)
         auz = dt2 * eigenValue(3)
         cfx = exp(aux)
         cfy = exp(auy)
         cfz = exp(auz)
         fc1x = cfx * cfx
         fc1y = cfy * cfy
         fc1z = cfz * cfz
         argx = aux*aux
         argy = auy*auy
         argz = auz*auz
         ar9x = e9 * argx
         ar9y = e9 * argy
         ar9z = e9 * argz
         ar7x = (ar9x+e7)*argx
         ar7y = (ar9y+e7)*argy
         ar7z = (ar9z+e7)*argz
         ar5x = (ar7x+e5)*argx
         ar5y = (ar7y+e5)*argy
         ar5z = (ar7z+e5)*argz
         ar3x = (ar5x+e3)*argx
         ar3y = (ar5y+e3)*argy
         ar3z = (ar5z+e3)*argz
         polx = ar3x + 1.d0
         poly = ar3y + 1.d0
         polz = ar3z + 1.d0
         fc2x = cfx*polx*deltat
         fc2y = cfy*poly*deltat
         fc2z = cfz*polz*deltat

         Texx = eigenVec(1,1)
         Texy = eigenVec(2,1)
         Texz = eigenVec(3,1)
         Teyx = eigenVec(1,2)
         Teyy = eigenVec(2,2)
         Teyz = eigenVec(3,2)
         Tezx = eigenVec(1,3)
         Tezy = eigenVec(2,3)
         Tezz = eigenVec(3,3)

         do i = 1 , N
           Rx = R(1,i)
           Ry = R(2,i)
           Rz = R(3,i)
           Vx = Vel(1,i)
           Vy = Vel(2,i)
           Vz = Vel(3,i)
           tempRx = Texx*Rx + Texy*Ry + Texz*Rz
           tempRy = Teyx*Rx + Teyy*Ry + Teyz*Rz
           tempRz = Tezx*Rx + Tezy*Ry + Tezz*Rz
           tempVx = Texx*Vx + Texy*Vy + Texz*Vz
           tempVy = Teyx*Vx + Teyy*Vy + Teyz*Vz
           tempVz = Tezx*Vx + Tezy*Vy + Tezz*Vz
           tempRx = tempRx * fc1x + tempVx * fc2x
           tempRy = tempRy * fc1y + tempVy * fc2y
           tempRz = tempRz * fc1z + tempVz * fc2z
         end do

! ## update H
         Hxx = H(1,1)
         Hxy = H(1,2)
         Hxz = H(1,3)
         Hyx = H(2,1)
         Hyy = H(2,2)
         Hyz = H(2,3)
         Hzx = H(3,1)
         Hzy = H(3,2)
         Hzz = H(3,3)
         tempHxx = Texx*Hxx + Texy*Hyx + Texz*Hzx
         tempHxy = Texx*Hxy + Texy*Hyy + Texz*Hzy
         tempHxz = Texx*Hxz + Texy*Hyz + Texz*Hzz
         tempHyx = Teyx*Hxx + Teyy*Hyx + Teyz*Hzx
         tempHyy = Teyx*Hxy + Teyy*Hyy + Teyz*Hzy
         tempHyz = Teyx*Hxz + Teyy*Hyz + Teyz*Hzz
         tempHzx = Tezx*Hxx + Tezy*Hyx + Tezz*Hzx
         tempHzy = Tezx*Hxy + Tezy*Hyy + Tezz*Hzy
         tempHzz = Tezx*Hxz + Tezy*Hyz + Tezz*Hzz
         tempHxx = tempHxx * fc1x
         tempHxy = tempHxy * fc1x
         tempHxz = tempHxz * fc1x
         tempHyx = tempHyx * fc1y
         tempHyy = tempHyy * fc1y
         tempHyz = tempHyz * fc1y
         tempHzx = tempHzx * fc1z
         tempHzy = tempHzy * fc1z
         tempHzz = tempHzz * fc1z

         H(1,1) = Texx*tempHxx + Teyx*tempHyx + Tezx*tempHzx
         H(1,2) = Texx*tempHxy + Teyx*tempHyy + Tezx*tempHzy
         H(1,3) = Texx*tempHxz + Teyx*tempHyz + Tezx*tempHzz
         H(2,1) = Texy*tempHxx + Teyy*tempHyx + Tezy*tempHzx
         H(2,2) = Texy*tempHxy + Teyy*tempHyy + Tezy*tempHzy
         H(2,3) = Texy*tempHxz + Teyy*tempHyz + Tezy*tempHzz
         H(3,1) = Texz*tempHxx + Teyz*tempHyx + Tezz*tempHzx
         H(3,2) = Texz*tempHxy + Teyz*tempHyy + Tezz*tempHzy
         H(3,3) = Texz*tempHxz + Teyz*tempHyz + Tezz*tempHzz

       else if( ( cBarostatMethod == 'A3'  ).or. &
       &        ( cBarostatMethod == 'A2' ) ) then

         aux = dt2 * Vg(1,1)
         auy = dt2 * Vg(2,2)
         auz = dt2 * Vg(3,3)
         aax = exp(aux)
         aay = exp(auy)
         aaz = exp(auz)
         fc1x = aax * aax
         fc1y = aay * aay
         fc1z = aaz * aaz

         argx = aux * aux
         argy = auy * auy
         argz = auz * auz
         ar9x = e9 * argx
         ar9y = e9 * argy
         ar9z = e9 * argz
         ar7x = (ar9x+e7)*argx
         ar7y = (ar9y+e7)*argy
         ar7z = (ar9z+e7)*argz
         ar5x = (ar7x+e5)*argx
         ar5y = (ar7y+e5)*argy
         ar5z = (ar7z+e5)*argz
         ar3x = (ar5x+e3)*argx
         ar3y = (ar5y+e3)*argy
         ar3z = (ar5z+e3)*argz
         polx = ar3x + 1.d0
         poly = ar3y + 1.d0
         polz = ar3z + 1.d0
         fc2x = aax * polx * deltat
         fc2y = aay * poly * deltat
         fc2z = aaz * polz * deltat

         do i = 1, N
           R(1,i) = R(1,i) * fc1x + Vel(1,i) * fc2x
           R(2,i) = R(2,i) * fc1y + Vel(2,i) * fc2y
           R(3,i) = R(3,i) * fc1z + Vel(3,i) * fc2z
         end do

! ## update H
         H(1,1) = H(1,1) * exp( Vg(1,1) * deltat )
         H(2,2) = H(2,2) * exp( Vg(2,2) * deltat )
         H(3,3) = H(3,3) * exp( Vg(3,3) * deltat )

       else if( cBarostatMethod == 'AN' ) then

! ## update the particle positions
         aux = dt2 * Vg(1,1)
         aax = exp(aux)
         fc1x = aax * aax

         argx = aux * aux
         ar9x = e9 * argx
         ar7x = (ar9x+e7)*argx
         ar5x = (ar7x+e5)*argx
         ar3x = (ar5x+e3)*argx
         polx = ar3x + 1.d0

         fc2x = aax * polx * deltat

         do i = 1 , N
           R(1,i) = R(1,i) * fc1x + Vel(1,i) * fc2x
           R(2,i) = R(2,i) * fc1x + Vel(2,i) * fc2x
           R(3,i) = R(3,i) * fc1x + Vel(3,i) * fc2x
         end do

! ## update H
         ClSc = exp(Vg(1,1) * deltat)
         H(1,1) = H(1,1) * ClSc
         H(2,2) = H(2,2) * ClSc
         H(3,3) = H(3,3) * ClSc

       end if

       if(QOpFix) call AddConstR

     end if

     if(QJarzynski) call SMD_reference

     call BcastRH

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

     if(NHam/=0) call Rot_FixedPoint

! ## get the new force
!   ---------------------------------------------
     call GetForce(1,istep)
     call GetAcc( A_fast, 1 )
     call SumFrc( A_fast )
!   ---------------------------------------------

! ## update the particle velocities

     if( mod(istep,BookFreq) == 0 ) call PairList
!   ---------------------------------------------
! Multiple time step
     if(mod(istep,lp) == 0) then
       call GetForce(2,istep)
       call GetAcc( A_mode, 2 )
       call SumFrc( A_mode )
     end if
! ------------------------------------------------------

     if(mod(istep,lk) == 0) then
       call GetForce(3,istep)
       call GetAcc( A_slow, 3 )
       call SumFrc( A_slow )
! ## update the virial
       Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
       &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
       &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong
       call SumVir( Virial )
     end if

! ## update the particle velocities

     if(QMaster) then
       if(mod(istep,lk) == 0) then
         do i = 1, N
           Asx = A_slow(1,i)*ww2 
           Asy = A_slow(2,i)*ww2 
           Asz = A_slow(3,i)*ww2 
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Asx + Amx + Afx
           Acy = Asy + Amy + Afy
           Acz = Asz + Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
         Virial = Virial - Virial_co
         call Bath_NP(2)
       else if(mod(istep,lp) == 0) then
         do i = 1, N
           Amx = A_mode(1,i)*ww1
           Amy = A_mode(2,i)*ww1
           Amz = A_mode(3,i)*ww1
           Afx = A_fast(1,i)*dt2
           Afy = A_fast(2,i)*dt2
           Afz = A_fast(3,i)*dt2
           Acx = Amx + Afx
           Acy = Amy + Afy
           Acz = Amz + Afz
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       else
         do i = 1, N
           Acx = A_fast(1,i)*dt2
           Acy = A_fast(2,i)*dt2
           Acz = A_fast(3,i)*dt2
           Vel(1,i) = Vel(1,i) + Acx
           Vel(2,i) = Vel(2,i) + Acy
           Vel(3,i) = Vel(3,i) + Acz
         end do
       end if
     end if

     if(QOpFix) call AddConstV(istep)

! ## remove the cell-momentum
     if(QMaster.and.QDelCellMove) call Elim_CellMove

! >> F monitor ##
     if( mod(istep,isampleF) == 0 ) then
       call Force_gA
       call SumFrcgA
       if(QMaster) then
         do i = NiniF + 1 , NfinF
           Fext(:,i) = ( A_fast(:,i) + A_mode(:,i) + A_slow(:,i) ) &
           &           * Mass(i) - Fint(:,i)
         end do
         call Monitor_Force
       end if
     end if
! << F monitor ##

! >> Steered MD ##
     if(QJarzynski.and.(mod(istep,lk)==0)) then
       call SMD_sample(istep)
     end if
! << Steered MD ##
! >> PMF of macrosphere
     if(QMacro.and.(mod(istep,lk)==0).and.NumSphere==1) then
       call MacroPMFsample(istep)
     end if
! << PMF of macrosphere
     if(QMaster.and.Qwcformol.and.(mod(istep,Nsample_wc)==0)) write(73,'(f12.7)') Ronc
     if(QMaster.and.QInert.and.(mod(istep,lk)==0)) call CntIn_Sample(istep)
     if((QCyl.or.QFSCyl).and.(mod(istep,lk)==0)) call CylPMFsample(istep)
! >> Eflux
     if(QEflux.and.QMaster.and.(mod(istep,lk)==0)) call Print_Eflux(istep)
! << Eflux

!   - store parameters --------------------------
     if( mod(istep,lk)  == 0 ) call Print_Energy_NP(istep)
     if( mod(istep,ixc) == 0 ) call Print_Config
     if( mod(istep,ixv) == 0 ) call Print_Velocity
     if( mod(istep,irs) == 0 ) call SaveParam
!   ---------------------------------------------

! ## check the cell strain
     if(mod(istep,lk) == 0) then
       call CheckCellShape
     end if

#ifdef EnergyRep
     if(mod(istep,lk)==0) call enganal(istep,2)
#endif

   end do

   if(QMaster.and.QSimAnneal) close(50)
   if(QMaster.and.Qwcformol) close(73)

end subroutine IntegrEOM_NP


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

subroutine IntegrEOM_iso

use Numbers, only : N
use CommonBlocks, only : QMaster, QHeatCap, QOpFix, QJarzynski, QSHAKE, &
&   QThermostat, Qdebug, Qwcformol, QDelCellMove, QInert, QCyl, QFSCyl
use Configuration
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use SimAnneal, only : QSimAnneal
use CGdata
use SHAKEparam, only : R_o
use OptConstraintParam, only : NHam, Rrot, RIni, Ene_OptC
use NonbondParam, only : Ene_Elec, Ene_LJ, Ene_NBshrt, Ene_NBlong, &
&   Ene_ELshrt, Ene_ELlong
use BondedParam, only : Ene_Bond, Ene_Angle, Ene_UB, Ene_Dihed, Ene_Impro
use AtomParam, only : Mass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, isampleHC, Timeps, &
&   deltat, dt2, irs
use wcparam, only : Nsample_wc, Ronc

implicit none

integer :: i, istep

real(8), dimension(3,N) :: A_fast
real(8), dimension(3,N) :: A_mode
real(8), dimension(3,N) :: A_slow

real(8) :: Upot, Dummy
integer :: ioutput

   if(QHeatCap) then
     ioutput = 1000
     ioutput = ioutput / isampleHC
   end if

   if(NHam/=0) Rrot = RIni

   if(QMaster.and.QOpFix) then
     call ConstPrepare
   end if
   if(QMaster.and.QJarzynski) then
     call SMD_pre
   end if

   if(QSimAnneal) call PreAnneal

   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

   if(Qdebug) then
     call checkf
     print *, 'Rrespa=', Rrespa
     print *, 'Rheal=',Rheal
   end if

! ---------------------
   call GetForceIso(0)
! ---------------------
   call GetAcc( A_fast, 1 )
   call GetAcc( A_mode, 2 )
   call GetAcc( A_slow, 3 )

   call SumFrc( A_fast )
   call SumFrc( A_mode )
   call SumFrc( A_slow )

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

     if(QSimAnneal) call Annealing(istep)

     if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
       if(mod(istep-1,lk) == 0) then
! ## thermostat
         if(QThermostat) call Thermostat(lk,1)
! ## update the particle velocities
         call VelUpdate(A_slow,dt2*lk)
       end if

! ## update the particle velocities
       if(mod(istep-1,lp) == 0) call VelUpdate(A_mode,dt2*lp)
       call VelUpdate(A_fast,dt2)

! update the particle positions
       if(QSHAKE) R_o = R
       R = R + Vel * deltat
! ## bond constraint
       if(QSHAKE) call SHAKE

       if(QOpFix) call AddConstR

     end if

     if(QJarzynski) call SMD_reference

     call BcastR

! ## get the new force
! ------------------------------------------------
     call GetForceIso(1)
     call GetAcc( A_fast, 1 )
     call SumFrc( A_fast )
!   ---------------------------------------------
! ## update the particle velocities
     if(QMaster) call VelUpdate(A_fast,dt2)

! ## Multiple time step
! ------------------------------------------------------
     if(mod(istep,lp) == 0) then
       call GetForceIso(2)
       call GetAcc( A_mode, 2 )
       call SumFrc( A_mode )
       if(QMaster) call VelUpdate(A_mode,dt2*lp)
     end if
! ------------------------------------------------------

! ## Multiple time step
! ------------------------------------------------------
     if(mod(istep,lk) == 0) then
       call GetForceIso(3)
       call GetAcc( A_slow, 3 )
       call SumFrc( A_slow )
       if(QMaster) call VelUpdate(A_slow,dt2*lk)
     end if
! ------------------------------------------------------

     if(QMaster) then
! ## bond constraint
       if(QSHAKE) call RATTLE
! ## thermostat velocities and thermostat positions
       if((mod(istep,lk)==0).and.QThermostat) call Thermostat(lk,2)
     end if

     if(QOpFix) call AddConstV(istep)

! ## remove the cell-momentum
     if(QMaster.and.QDelCellMove) call Elim_CellMove

! >> F monitor ##
     if( mod(istep,isampleF) == 0 ) then
       call Force_gA
       call SumFrcgA
       if(QMaster) then
         do i = NiniF + 1 , NfinF
           Fext(:,i) = ( A_fast(:,i) + A_mode(:,i) + A_slow(:,i) ) &
           &           * Mass(i) - Fint(:,i)
         end do
         call Monitor_Force
       end if
     end if
! << F monitor ##

! >> Heat Capacity
     if( QHeatCap .and. (mod(istep, isampleHC) == 0) ) then
       Upot = Ene_Bond   + Ene_Angle + Ene_UB   + Ene_Dihed  + Ene_Impro &
       &    + Ene_Elec   + Ene_LJ    + Ene_OptC + Ene_NBshrt + Ene_NBlong&
       &    + Ene_ELshrt + Ene_ELlong
       call SumHC( Upot, Dummy, Dummy, Dummy, Dummy )
       if(QMaster) then
         call HeatCapacity( Upot, Dummy, Dummy, Dummy, Dummy, istep/isampleHC, ioutput )
       end if
     end if
! << Heat Capacity

! >> Steered MD ##
     if(QJarzynski.and.(mod(istep,lk)==0)) then
       call SMD_sample(istep)
     end if
! << Steered MD ##
     if(QMaster.and.Qwcformol.and.(mod(istep,Nsample_wc)==0)) write(73,'(f12.7)') Ronc
     if(QMaster.and.QInert.and.(mod(istep,lk)==0)) call CntIn_Sample(istep)
     if((QCyl.or.QFSCyl).and.(mod(istep,lk)==0)) call CylPMFsample(istep)

!   - save parameters ---------------------------------------
     if( mod(istep,lk)  == 0 ) call Print_Energy_iso(istep)
     if( mod(istep,ixc) == 0 ) call Print_Config
     if( mod(istep,ixv) == 0 ) call Print_Velocity
     if( mod(istep,irs) == 0 ) call SaveParam
!   ----------------------------------------------------------

#ifdef EnergyRep
     if(mod(istep,lk)==0) call enganal(istep,1)
#endif

   end do

   if(QMaster.and.QSimAnneal) close(50)
   if(QMaster.and.Qwcformol) close(73)

contains

subroutine checkf

   use CGdata, only: Rrespa, Rheal

   integer :: ia, ib
   real(8) :: xx, x2, yy, Fswitch

   open(111,file='checkf.dat')

   ia = nint( Rrespa*10.)
   ib = nint( (Rrespa+Rheal)*10.)

   do i = ia, ib

     xx = dble(i)*0.1d0
     x2 = xx * xx

     yy = Fswitch(x2)

     write(111,*) xx, yy

   end do

   close(111)

end subroutine checkf

end subroutine IntegrEOM_iso


!######################################################################
!######################################################################


subroutine VelUpdate(Acc,ww)

use Configuration, only : Vel
use Numbers, only : N

implicit none

integer :: i
real(8), dimension(3,N) :: Acc
real(8) :: ww

   do i = 1, N
     Vel(:,i) = Vel(:,i) + ww * Acc(:,i)
   end do

end subroutine VelUpdate
