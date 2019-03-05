! ############################
! ## SUBROUTINE LIST 
! ## -- IntegrEOM_NV_RIGID 
! ## -- IntegrEOM_NP_RIGID 
! ## -- IntegrEOM_iso_RIGID 
! ############################


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

subroutine IntegrEOM_NV_RIGID

use Numbers, only : N
use CommonBlocks, only : QMaster, QVolScale, QOpFix, QJarzynski, QCyl, QFSCyl, &
&   QTICR, QCorrectCutoff, QThermostat, QMacro, Qwcformol, QDelCellMove, QInert,&
&   QEflux
use Configuration, only : R
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use RBparam, only : NumRB, NumRBType, NumRBAtom, InvInertiaRB, QSingle, &
&   Quaternion, V_RB, Lmoment, InvMassRB, QLinear, RBType, R_RB
use FEparam
use SimAnneal, only :QSimAnneal
use OptConstraintParam, only : NHam
use CellParam, only : Volume, Vsc_Rate
use TailCorrect
use AtomParam, only : InvMass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, BookFreq, Timeps, &
&   deltat, dt2, irs
use ThermoData, only : Virial
use wcparam, only : Nsample_wc, Ronc
use CGball, only : NumSphere

implicit none

integer :: i, k, istep
real(8) :: ww1, ww2

real(8), dimension(3,N) :: A_fast, A_mode, A_slow
real(8), dimension(3,NumRB) :: G_fast, G_mode, G_slow
real(8), dimension(3,NumRB) :: T_fast, T_mode, T_slow
real(8), dimension(NumRBType) :: rizy, rizx
real(8) :: qn2, qn, th1, th2, Lmomentztemp
real(8) :: snt, cst, det
real(8), dimension(3) :: Lt, Omt
real(8) :: omg, th, diag, offd, omg2
real(8), dimension(4) :: qq
real(8) :: Gsx, Gsy, Gsz, Gmx, Gmy, Gmz, Gfx, Gfy, Gfz, Gcx, Gcy, Gcz
real(8) :: Tsx, Tsy, Tsz, Tmx, Tmy, Tmz, Tfx, Tfy, Tfz, Tcx, Tcy, Tcz
integer :: MyType
! >> TICR ##
integer :: npoint, niter
! << TICR ##
external det

   ww2 = dt2*lk
   ww1 = dt2*lp

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

   if(QVolScale) then
     Vsc_Rate(:) = deltat * Vsc_Rate(:)
   end if

   if(NHam/=0) call Rot_FixedPoint

   if(QMaster.and.QOpFix) then
     call ConstPrepare
   end if
   if(QMaster.and.QJarzynski) then
     call SMD_pre
   end if

! >> TICR ##
   if(QTICR) then
     npoint = 1
     niter  = 0
     call SetTICR(0,npoint)
   end if
! << TICR ##

   if(QSimAnneal) call PreAnneal

   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

! ----------------------------
   call IntraMolVec

   call GetForce(0,0)

   if(QTICR.and.QLJ) then
     call Force_TICR(npoint)
   end if
! ----------------------------

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

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

     if(QSimAnneal) call Annealing(istep)

! ## Master - Slave
     if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
       if(mod(istep-1,lk) == 0) then
! ## thermostat velocities and thermostat positions
         if(QThermostat) call Thermostat(lk,1)
! ##  V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)
         k = 0
         do i = 1 , NumRB
           Gsx = G_slow(1,i)*ww2
           Gsy = G_slow(2,i)*ww2
           Gsz = G_slow(3,i)*ww2
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gsx + Gmx + Gfx
           Gcy = Gsy + Gmy + Gfy
           Gcz = Gsz + Gmz + Gfz
           Tsx = T_slow(1,i)*ww2
           Tsy = T_slow(2,i)*ww2
           Tsz = T_slow(3,i)*ww2
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tsx + Tmx + Tfx
           Tcy = Tsy + Tmy + Tfy
           Tcz = Tsz + Tmz + Tfz

           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
       else if(mod(istep-1,lp) == 0) then
         k = 0
         do i = 1 , NumRB
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gmx + Gfx
           Gcy = Gmy + Gfy
           Gcz = Gmz + Gfz
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tmx + Tfx
           Tcy = Tmy + Tfy
           Tcz = Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
       else
         k = 0
         do i = 1 , NumRB
           Gcx = G_fast(1,i)*dt2
           Gcy = G_fast(2,i)*dt2
           Gcz = G_fast(3,i)*dt2
           Tcx = T_fast(1,i)*dt2
           Tcy = T_fast(2,i)*dt2
           Tcz = T_fast(3,i)*dt2
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
       end if

! ##### Omg~(t+dt/2) = I^-1 L~(t+dt/2)
! ##### Q(t+dt) = exp(dt/2*A[Omg~(t+dt/2)]) Q(t)

       do i = 1 , NumRB

         if(QSingle(i)) cycle

         MyType = RBType(i)

         qq(:) = Quaternion(:,i)

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
       do i = 1 , NumRB
         R_RB(1,i) = R_RB(1,i) + V_RB(1,i) * deltat
         R_RB(2,i) = R_RB(2,i) + V_RB(2,i) * deltat
         R_RB(3,i) = R_RB(3,i) + V_RB(3,i) * deltat
       end do

       if(QVolScale) call VolScRupdateRB

     end if

     if(QJarzynski) call SMD_reference

     if(QVolScale) then
       call VolScVcorrectRB
     else
       call BcastRgQuat
     end if

! ## update all interacting sites
!   ------------------
     call IntraMolVec 
!   ------------------

     if( mod(istep,BookFreq) == 0 ) then
       if(QTICR) then
         call PairListTICR
       else
         call PairList
       end if
     end if
! ## get the new force
!   ---------------------------------------------
     call GetForce(1,istep)
     call GetAcc( A_fast, 1 )
     call SumFrc( A_fast )
     if(QMaster) call Force_Div_Component(A_fast,G_fast,T_fast)
!   ---------------------------------------------
     if(mod(istep,lp) == 0) then
       call GetForce(2,istep)
       if(QTICR.and.QLJ) then
         call Force_TICR(npoint)
       end if
       call GetAcc( A_mode, 2 )
       call SumFrc( A_mode )
       if(QMaster) call Force_Div_Component(A_mode,G_mode,T_mode)
     end if
!   ---------------------------------------------
     if(mod(istep,lk) == 0) then
       call GetForce(3,istep)
       call GetAcc( A_slow, 3 )
       call SumFrc( A_slow )
       if(QMaster) call Force_Div_Component(A_slow,G_slow,T_slow)
     end if
!   ---------------------------------------------
! ----------------------------------
! ## update the particle velocities
! ----------------------------------
     if(QMaster) then
       if(mod(istep,lk) == 0) then
         k = 0
         do i = 1 , NumRB
           Gsx = G_slow(1,i)*ww2
           Gsy = G_slow(2,i)*ww2
           Gsz = G_slow(3,i)*ww2
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gsx + Gmx + Gfx
           Gcy = Gsy + Gmy + Gfy
           Gcz = Gsz + Gmz + Gfz
           Tsx = T_slow(1,i)*ww2
           Tsy = T_slow(2,i)*ww2
           Tsz = T_slow(3,i)*ww2
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tsx + Tmx + Tfx
           Tcy = Tsy + Tmy + Tfy
           Tcz = Tsz + Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             end if
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             k = k + NumRBAtom(MyType)
           end if
         end do
! ## thermostat velocities and thermostat positions
         if(QThermostat) call Thermostat(lk,2)
       else if(mod(istep,lp) == 0) then
         k  = 0
         do i = 1 , NumRB
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gmx + Gfx
           Gcy = Gmy + Gfy
           Gcz = Gmz + Gfz
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tmx + Tfx
           Tcy = Tmy + Tfy
           Tcz = Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             end if
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             k = k + NumRBAtom(MyType)
           end if
         end do
       else
         k = 0
         do i = 1 , NumRB
           Gcx = G_fast(1,i)*dt2
           Gcy = G_fast(2,i)*dt2
           Gcz = G_fast(3,i)*dt2
           Tcx = T_fast(1,i)*dt2
           Tcy = T_fast(2,i)*dt2
           Tcz = T_fast(3,i)*dt2
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             end if
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             k = k + NumRBAtom(MyType)
           end if
         end do
       end if
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
           Fext(:,i) = ( A_fast(:,i) + A_mode(:,i) + A_slow(:,i) ) - Fint(:,i)
         end do
         call Monitor_Force
       end if
     end if
! << F monitor ##

! >> TICR ##
     if( mod(istep,isampleTI) == 0 ) then
       call TICRSampling(niter,npoint)
     end if
! << TICR ##

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

end subroutine IntegrEOM_NV_RIGID


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

subroutine IntegrEOM_NP_RIGID

use Numbers, only : N
use CommonBlocks, only : QMaster, QOpFix, QJarzynski, &
&   QTICR, QCorrectCutoff, QThermostat, Qwcformol, &
&   cThermostatMethod, cBarostatMethod, ForceField, &
&   QMacro, QDelCellMove, QInert, QCyl, QFSCyl, QEflux
use Configuration, only : R
use EAM_param, only : Vir_EAM
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use RBparam, only : NumRB, NumRBType, NumRBAtom, InvInertiaRB, QSingle, &
&   Quaternion, V_RB, Lmoment, InvMassRB, QLinear, RBType, R_RB
use FEparam
use SimAnneal, only : QSimAnneal
use BathParam, only: gkT, Baro_kin, Vg
use EwaldParam, only : Vir_Eksp
use OptConstraintParam, only : NHam, Vir_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong
use BondedParam, only : Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro
use CellParam, only : H, InvH, Volume
use TailCorrect
use AtomParam, only : InvMass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, BookFreq, Timeps, &
&   deltat, dt2, irs
use ThermoData, only : Ene_kin, Virial
use wcparam, only : Nsample_wc, Ronc
use CGball, only : NumSphere

implicit none

real(8), parameter :: e3 = 1.d0 / 6.d0
real(8), parameter :: e5 = e3   / 20.d0
real(8), parameter :: e7 = e5   / 42.d0
real(8), parameter :: e9 = e7   / 72.d0

integer :: i, k, istep

real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8), dimension(3,3) :: tempoH
real(8), dimension(3) :: eigenValue
real(8), dimension(3) :: fc1, fc2
real(8), dimension(3) :: tempor, tempov
real(8), dimension(3) :: aa, aa2, aa3
real(8), dimension(3) :: arg3, poly3
real(8), dimension(3) :: bb3

real(8) :: Gsx, Gsy, Gsz, Gmx, Gmy, Gmz, Gfx, Gfy, Gfz, Gcx, Gcy, Gcz
real(8) :: Tsx, Tsy, Tsz, Tmx, Tmy, Tmz, Tfx, Tfy, Tfz, Tcx, Tcy, Tcz
real(8) :: ww1, ww2, cf, arg2
real(8) :: poly
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

real(8) :: ClSc

! ## TICR
integer :: npoint, niter

real(8) :: det
External det

if(ForceField(1:3) == 'EAM') then
open(13,file='CellMatrix.dat',form='unformatted')
end if

   ww2 = dt2*lk
   ww1 = dt2*lp

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

   if(QMaster.and.QOpFix) call ConstPrepare

   if(QMaster.and.QJarzynski) call SMD_pre

! >> TICR ##
   if(QTICR) then
     npoint = 1
     niter  = 0
     call SetTICR(0,npoint)
   end if
! << TICR ##

   if(QSimAnneal) call PreAnneal

   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

   call IntraMolVec

   call GetForce(0,0)

   if(QTICR.and.QLJ) then
     call Force_TICR(npoint)
   end if

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

   if(QMaster) then
     call Force_Div_Component(A_fast,G_fast,T_fast)
     call Force_Div_Component(A_mode,G_mode,T_mode)
     call Force_Div_Component(A_slow,G_slow,T_slow)
   end if

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
! ## V(t+l*dt/2)=V(t)+(l*dt/2)*F(t)
         k = 0
         do i = 1 , NumRB
           Gsx = G_slow(1,i)*ww2
           Gsy = G_slow(2,i)*ww2
           Gsz = G_slow(3,i)*ww2
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gsx + Gmx + Gfx
           Gcy = Gsy + Gmy + Gfy
           Gcz = Gsz + Gmz + Gfz
           Tsx = T_slow(1,i)*ww2
           Tsy = T_slow(2,i)*ww2
           Tsz = T_slow(3,i)*ww2
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tsx + Tmx + Tfx
           Tcy = Tsy + Tmy + Tfy
           Tcz = Tsz + Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
       else if(mod(istep-1,lp) == 0) then
         k = 0
         do i = 1 , NumRB
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gmx + Gfx
           Gcy = Gmy + Gfy
           Gcz = Gmz + Gfz
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tmx + Tfx
           Tcy = Tmy + Tfy
           Tcz = Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
       else
! ## update the particle velocities
         k = 0
         do i = 1 , NumRB
           Gcx = G_fast(1,i)*dt2
           Gcy = G_fast(2,i)*dt2
           Gcz = G_fast(3,i)*dt2
           Tcx = T_fast(1,i)*dt2
           Tcy = T_fast(2,i)*dt2
           Tcz = T_fast(3,i)*dt2
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
       end if

! ##### Omg~(t+dt/2) = I^-1 L~(t+dt/2)
! ##### Q(t+dt) = exp(dt/2*A[Omg~(t+dt/2)]) Q(t)

       do i = 1 , NumRB

         if(QSingle(i)) cycle

         MyType = RBType(i)

         qq(:) = Quaternion(:,i)

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

       else if( ( cBarostatMethod == 'A3' ).or. &
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
         ClSc = exp(Vg(1,1) * deltat)
         do i = 1 , 3
           H(i,i) = H(i,i) * ClSc
         end do

       end if

       if(QOpFix) call AddConstR

     end if

     if(QJarzynski) call SMD_reference

     call BcastRgQuatH

! ## update all interacting sites
!   ------------------
     call IntraMolVec 
!   ------------------

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

     if( mod(istep,BookFreq) == 0 ) then
       if(QTICR) then
         call PairListTICR
       else
         call PairList
       end if
     end if
! ## get the new force
!   ---------------------------------------------
     call GetForce(1,istep)
     call GetAcc( A_fast, 1 )
     call SumFrc( A_fast )
     if(QMaster) call Force_Div_Component(A_fast,G_fast,T_fast)
!   ---------------------------------------------
     if(mod(istep,lp) == 0) then
       call GetForce(2,istep)
       if(QTICR.and.QLJ) then
         call Force_TICR(npoint)
       end if
       call GetAcc( A_mode, 2 )
       call SumFrc( A_mode )
       if(QMaster) call Force_Div_Component(A_mode,G_mode,T_mode)
     end if
!   ---------------------------------------------
     if(mod(istep,lk) == 0) then
       call GetForce(3,istep)
       call GetAcc( A_slow, 3 )
       call SumFrc( A_slow )
       if(QMaster) call Force_Div_Component(A_slow,G_slow,T_slow)
     end if
!   ---------------------------------------------

! ## update the particle velocities
     if(QMaster) then
       if(mod(istep,lk) == 0) then
         k = 0
         do i = 1 , NumRB
           Gsx = G_slow(1,i)*ww2
           Gsy = G_slow(2,i)*ww2
           Gsz = G_slow(3,i)*ww2
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gsx + Gmx + Gfx
           Gcy = Gsy + Gmy + Gfy
           Gcz = Gsz + Gmz + Gfz
           Tsx = T_slow(1,i)*ww2
           Tsy = T_slow(2,i)*ww2
           Tsz = T_slow(3,i)*ww2
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tsx + Tmx + Tfx
           Tcy = Tsy + Tmy + Tfy
           Tcz = Tsz + Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             end if
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             k = k + NumRBAtom(MyType)
           end if
         end do

! ## update the virial
         Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
         &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
         &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong
         call SumVir( Virial )
         if(QMaster) then
           Virial = Virial - Virial_co
           call Bath_NP(2)
         end if

       else if(mod(istep,lp) == 0) then
         k  = 0
         do i = 1 , NumRB
           Gmx = G_mode(1,i)*ww1
           Gmy = G_mode(2,i)*ww1
           Gmz = G_mode(3,i)*ww1
           Gfx = G_fast(1,i)*dt2
           Gfy = G_fast(2,i)*dt2
           Gfz = G_fast(3,i)*dt2
           Gcx = Gmx + Gfx
           Gcy = Gmy + Gfy
           Gcz = Gmz + Gfz
           Tmx = T_mode(1,i)*ww1
           Tmy = T_mode(2,i)*ww1
           Tmz = T_mode(3,i)*ww1
           Tfx = T_fast(1,i)*dt2
           Tfy = T_fast(2,i)*dt2
           Tfz = T_fast(3,i)*dt2
           Tcx = Tmx + Tfx
           Tcy = Tmy + Tfy
           Tcz = Tmz + Tfz
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             end if
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             k = k + NumRBAtom(MyType)
           end if
         end do
       else
         k = 0
         do i = 1 , NumRB
           Gcx = G_fast(1,i)*dt2
           Gcy = G_fast(2,i)*dt2
           Gcz = G_fast(3,i)*dt2
           Tcx = T_fast(1,i)*dt2
           Tcy = T_fast(2,i)*dt2
           Tcz = T_fast(3,i)*dt2
           if(QSingle(i)) then
             k = k + 1
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMass(k)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMass(k)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMass(k)
           else
             MyType = RBType(i)
             if(QLinear(i)) then
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             else
               Lt(1) = Lmoment(1,i)
               Lt(2) = Lmoment(2,i)
               Lt(3) = Lmoment(3,i)
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
               Lmoment(1,i) = Lmoment(1,i) + Tcx
               Lmoment(2,i) = Lmoment(2,i) + Tcy
               Lmoment(3,i) = Lmoment(3,i) + Tcz
             end if
             V_RB(1,i) = V_RB(1,i) + Gcx * InvMassRB(MyType)
             V_RB(2,i) = V_RB(2,i) + Gcy * InvMassRB(MyType)
             V_RB(3,i) = V_RB(3,i) + Gcz * InvMassRB(MyType)
             k = k + NumRBAtom(MyType)
           end if
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
           &           - Fint(:,i)
         end do
         call Monitor_Force
       end if
     end if
! << F monitor ##

! >> TICR ##
     if( mod(istep,isampleTI) == 0 ) then
       call TICRSampling(niter,npoint)
     end if
! << TICR ##

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

!   - store parameters ------------------------------
     if( mod(istep,lk)  == 0 ) call Print_Energy_NP(istep)
     if( mod(istep,ixc) == 0 ) call Print_Config
     if( mod(istep,ixv) == 0 ) call Print_Velocity
     if( mod(istep,irs) == 0 ) call SaveParam
!   -------------------------------------------------

! ## check the cell strain
     if(mod(istep,lk) == 0) call CheckCellShape

#ifdef EnergyRep
     if(mod(istep,lk)==0) call enganal(istep,2)
#endif

   end do

   if(QMaster.and.QSimAnneal) close(50)
   if(QMaster.and.Qwcformol) close(73)

end subroutine IntegrEOM_NP_RIGID


!######################################################################
!######################################################################


! ***********************************************************
! ** MD main part : integration of the equations of motion **
! ** time evolution of the particle coordinate             **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine IntegrEOM_iso_RIGID

use Numbers, only : N
use CommonBlocks, only : QMaster, QOpFix, QJarzynski, &
&   QTICR, QThermostat, Qwcformol, QDelCellMove, QInert, QCyl, QFSCyl
use F_monitor, only : isampleF, NiniF, NfinF, Fint, Fext
use RBparam, only : NumRB, NumRBType, NumRBAtom, InvInertiaRB, QSingle, &
&   Quaternion, V_RB, Lmoment, InvMassRB, QLinear, RBType, R_RB
use FEparam
use SimAnneal, only : QSimAnneal
use OptConstraintParam, only : NHam, Rrot, RIni
use AtomParam, only : InvMass
use TimeParam, only : Nstep, ixc, ixv, lp, lk, Timeps, &
&   deltat, dt2, irs
use wcparam, only : Nsample_wc, Ronc

implicit none

integer :: i, k, istep
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

   if(NHam/=0) Rrot = RIni

   if(QTICR) then
     write(*,*) 'ERROR : TICR cannot be used for isolated systems.'
     call Finalize
   end if

   if(QMaster.and.QOpFix) then
     call ConstPrepare
   end if
   if(QMaster.and.QJarzynski) then
     call SMD_pre
   end if

   if(QSimAnneal) call PreAnneal

   if(QMaster.and.Qwcformol) open(73,file='Position_constmol.dat',status='unknown')

! -------------------
   call IntraMolVec

   call GetForceIso(0)
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

! ----------------------------------------------------------------------
!                  ### start MD time evolution ###
! ----------------------------------------------------------------------

   do istep = 1 , Nstep

     Timeps = Timeps + deltat

     if(QSimAnneal) call Annealing(istep)

! ## Master - Slave
     if(QMaster) then

! ## Multiple Time Scale
! ------------------------------------------------------
       if(mod(istep-1,lk) == 0) then

! ## thermostat velocities and thermostat positions
         if(QThermostat) call Thermostat(lk,1)

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

! ---------------------------------
! ## update the particle positions
! ---------------------------------
       R_RB = R_RB + V_RB * deltat

       if(QOpFix) call AddConstR

     end if

     if(QJarzynski) call SMD_reference

     call BcastRgQuat

! ## update all interacting sites
!   ------------------
     call IntraMolVec 
!   ------------------

! ## get the new force
!   ---------------------------------------------
     call GetForceIso(1)

     call GetAcc( A_fast, 1 )
     call SumFrc( A_fast )
!   ---------------------------------------------

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
     if(mod(istep,lp) == 0) then

       call GetForceIso(2)

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

     end if
! ------------------------------------------------------

! ## Multiple time step
! ------------------------------------------------------
     if(mod(istep,lk) == 0) then

       call GetForceIso(3)

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

! ## thermostat velocities and thermostat positions
         if(QThermostat) call Thermostat(lk,2)

       end if

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
           &           - Fint(:,i)
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
     if(QMaster.and.Qwcformol.and.(mod(istep,Nsample_wc)==0)) write(73,'(f12.7)') Ronc
     if(QMaster.and.QInert.and.(mod(istep,lk)==0)) call CntIn_Sample(istep)
     if((QCyl.or.QFSCyl).and.(mod(istep,lk)==0)) call CylPMFsample(istep)

!   - save parameters ------------------------------------------
     if( mod(istep,lk)  == 0 ) call Print_Energy_iso(istep)
     if( mod(istep,ixc) == 0 ) call Print_Config
     if( mod(istep,ixv) == 0 ) call Print_Velocity
     if( mod(istep,irs) == 0 ) call SaveParam
!   ------------------------------------------------------------

#ifdef EnergyRep
     if(mod(istep,lk)==0) call enganal(istep,1)
#endif

   end do

   if(QMaster.and.QSimAnneal) close(50)
   if(QMaster.and.Qwcformol) close(73)

end subroutine IntegrEOM_iso_RIGID
