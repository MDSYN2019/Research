! ############################
! ## SUBROUTINE LIST 
! ## -- Thermostat
! ## -- Thermostat_NH 
! ## -- Thermostat_NHC 
! ## -- Thermostat_MNHC 
! ## -- Bath_NP_SHAKE 
! ## -- Bath_NP 
! ## -- BaroThermostatAN 
! ## -- BaroThermostatA3 
! ## -- BaroThermostatPR 
! ## -- BarostatAN 
! ## -- BarostatA3 
! ## -- BarostatPR 
! ## -- BaroThermostatAN_MNHC 
! ## -- BaroThermostatA3_MNHC 
! ## -- BaroThermostatPR_MNHC 
! ## -- BaroThermostatAN_NH 
! ## -- BaroThermostatA3_NH 
! ## -- BaroThermostatPR_NH 
! ############################


!######################################################################
!######################################################################


subroutine Thermostat(ll,m)

use CommonBlocks, only : QRigidBody, cThermostatMethod
use RBparam, only : V_RB, Lmoment
use Configuration, only : Vel
use BathParam, only : gkT
use ThermoData, only : Ene_kin

implicit none

integer :: ll, m
real(8) :: VelScale

   if( cThermostatMethod == 'NH' ) then

     call Thermostat_NH(ll)

   else if( cThermostatMethod == 'NHC' ) then

     call Thermostat_NHC(ll)

   else if( cThermostatMethod == 'MNHC' ) then

     call Thermostat_MNHC(ll)

   else if(( cThermostatMethod == 'VSCALE' ).and.(m == 2)) then

     call CalcTemp

     VelScale = sqrt( gkT / Ene_kin )
     if(QRigidBody) then
       V_RB    = V_RB    * VelScale
       Lmoment = Lmoment * VelScale
     else
       Vel = Vel * VelScale
     end if

   end if

end subroutine Thermostat


!######################################################################
!######################################################################


! ***********************************************************
! ** integration of the equations of motion                **
! ** time evolution of thermostat                          **
! ** Nose-Hoover chain ( length = NHchain )                **
! **                                                       **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine Thermostat_NH(ll)

use CommonBlocks, only : QRigidBody
use BathParam
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use ThermoData, only : Ene_kin

implicit none

integer :: i, j, ll
real(8) :: Scale, aa, w2, w4
real(8) :: Gms

! ------------------------------------
! Nose-Hoover chain length [NHchain]
! ------------------------------------
   Scale=1.d0

! ----------------
   call CalcTemp
! ----------------

   Gms = ( Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

   do i = 1 , Nsc

     do j = 1 , NYoshid

       w2 = wdti2(j) * dble(ll)
       w4 = wdti4(j) * dble(ll)

! ## update the themostat velocities coupling to particles

       Vss(1) = Vss(1) + Gms * w4

! --------------------------------
! update the particle velocities
! --------------------------------
       aa      = exp( -w2 * Vss(1) )
       Scale   = Scale * aa
       Ene_kin = Ene_kin * aa * aa

! --------------------------------
! update the thermostat positions
! --------------------------------
       Rss = Rss + Vss * w2

! ------------------
! update the forces
! ------------------
       Gms = ( Ene_kin - gkT ) * InvMts(1)

! -----------------------------------
! update the thermostat velocities
! -----------------------------------
       Vss(1) = Vss(1) + Gms * w4

     end do

   end do

! ---------------------------------
! update the particle velocities
! ---------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     do i = 1 , NumRB

       V_RB(:,i)    = V_RB(:,i)    * Scale

       if(QSingle(i)) cycle

       if(QLinear(i)) then

         Lmoment(2,i) = Lmoment(2,i) * Scale
         Lmoment(3,i) = Lmoment(3,i) * Scale

       else

         Lmoment(:,i) = Lmoment(:,i) * Scale

       end if

     end do

! ## flexible ## -------------------------------------------------------

   else

     Vel = Vel * Scale

   end if

!   write(*,*) 'Scale=',Scale

end subroutine Thermostat_NH


!######################################################################
!######################################################################


! ***********************************************************
! ** integration of the equations of motion                **
! ** time evolution of thermostat                          **
! ** Nose-Hoover chain ( length = NHchain )                **
! **                                                       **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine Thermostat_NHC(ll)

use CommonBlocks, only : QRigidBody
use BathParam
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use ThermoData, only : Ene_kin

implicit none

integer :: i, j, k, ll
real(8) :: Scale, aa, w2, w4, w8
real(8), dimension(NHchain) :: Gms

! ------------------------------------
! Nose-Hoover chain length [NHchain]
! ------------------------------------
   Scale=1.d0

! ----------------
   call CalcTemp
! ----------------

   Gms(1) = ( Ene_kin - gkT ) * InvMts(1)

   do i = 1 , NHchain-1

     Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)

   end do


! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

   do i = 1 , Nsc

     do j = 1 , NYoshid

       w2 = wdti2(j) * dble(ll)
       w4 = wdti4(j) * dble(ll)
       w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

       Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

       do k = NHchain-1 , 1 , -1

         aa     = exp( -w8 * Vss(k+1) )
         Vss(k) = ( Vss(k) * aa + Gms(k) * w4 ) * aa

       end do

! --------------------------------
! update the particle velocities
! --------------------------------
       aa      = exp( -w2 * Vss(1) )
       Scale   = Scale * aa
       Ene_kin = Ene_kin * aa * aa

! --------------------------------
! update the thermostat positions
! --------------------------------
       Rss = Rss + Vss * w2

! ------------------
! update the forces
! ------------------
       Gms(1) = ( Ene_kin - gkT ) * InvMts(1)

! -----------------------------------
! update the thermostat velocities
! -----------------------------------
       do k = 1 , NHchain-1

         aa       = exp( -w8 * Vss(k+1) )
         Vss(k)   = ( Vss(k) * aa + w4 * Gms(k) ) * aa
         Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)

       end do

       Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

     end do

   end do

! ---------------------------------
! update the particle velocities
! ---------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     do i = 1 , NumRB

       V_RB(:,i)    = V_RB(:,i)    * Scale

       if(QSingle(i)) cycle

       if(QLinear(i)) then

         Lmoment(2,i) = Lmoment(2,i) * Scale
         Lmoment(3,i) = Lmoment(3,i) * Scale

       else

         Lmoment(:,i) = Lmoment(:,i) * Scale

       end if

     end do

! ## flexible ## -------------------------------------------------------

   else

     Vel = Vel * Scale

   end if

!   write(*,*) 'Scale=',Scale

end subroutine Thermostat_NHC


!######################################################################
!######################################################################


! ***********************************************************
! ** integration of the equations of motion of thermostat  **
! ** Massive Nose-Hoover chain ( length = NHchain )        **
! **                                                       **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine Thermostat_MNHC(ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody
use BathParam
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear, RBType, NumRBAtom, &
&   MassRB, InvInertiaRB
use Configuration, only : Vel
use AtomParam, only : Mass

implicit none

integer :: i, j, k, ll, l, m, MyType, Num
real(8) :: aa, w2, w4, w8
real(8), dimension(NHchain,NumMNHC) :: Gms

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

! ## translocation

     k = 0
     do i = 1, NumRB
       j = (i-1) * 3
       if(QSingle(i)) then
         k = k + 1
         Gms(1,j+1) = ( Mass(k) * V_RB(1,i) * V_RB(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(k) * V_RB(2,i) * V_RB(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(k) * V_RB(3,i) * V_RB(3,i) - kT ) * InvMMNHC(1,j+3)
       else
         MyType = RBType(i)
         k = k + NumRBAtom(MyType)
         Gms(1,j+1) = ( MassRB(MyType) * V_RB(1,i) * V_RB(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( MassRB(MyType) * V_RB(2,i) * V_RB(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( MassRB(MyType) * V_RB(3,i) * V_RB(3,i) - kT ) * InvMMNHC(1,j+3)
       end if
     end do

! ## rotation

     Num = NfT
     do i = 1, NumRB
       if(QSingle(i)) cycle
       MyType = RBType(i)
       if(QLinear(i)) then
         Gms(1,Num+1) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) - kT ) &
         &              * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) - kT ) &
         &              * InvMMNHC(1,Num+2)
         Num = Num + 2
       else
         Gms(1,Num+1) = ( Lmoment(1,i) * Lmoment(1,i) * InvInertiaRB(1,MyType) - kT ) &
         &              * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) - kT ) &
         &              * InvMMNHC(1,Num+2)
         Gms(1,Num+3) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) - kT ) &
         &              * InvMMNHC(1,Num+3)
         Num = Num + 3
       end if
     end do

     do i = 1 , NHchain-1
       do j = 1, NumMNHC
         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)
       end do
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1
           do l = 1, NumMNHC
             aa = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + Gms(k,l) * w4 ) * aa
           end do
         end do

! --------------------------------
! update the particle velocities
! --------------------------------
         l = 0
         do k = 1, NumRB
           V_RB(1,k) = V_RB(1,k) * exp( -w2 * VMNHC(1,l+1) )
           V_RB(2,k) = V_RB(2,k) * exp( -w2 * VMNHC(1,l+2) )
           V_RB(3,k) = V_RB(3,k) * exp( -w2 * VMNHC(1,l+3) )
           l = l + 3
         end do

         do k = 1 , NumRB
           if(QSingle(k)) cycle
           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+2) )
             l = l + 2
           else
             Lmoment(1,k) = Lmoment(1,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+2) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+3) )
             l = l + 3
           end if
         end do

! --------------------------------
! update the thermostat positions
! --------------------------------
         RMNHC = RMNHC + VMNHC * w2

! ------------------
! update the forces
! ------------------
         k = 0

         do l = 1, NumRB
           m = (l-1) * 3
           if(QSingle(l)) then
             k = k + 1
             Gms(1,m+1) = ( Mass(k) * V_RB(1,l) * V_RB(1,l) - kT ) * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( Mass(k) * V_RB(2,l) * V_RB(2,l) - kT ) * InvMMNHC(1,m+2)
             Gms(1,m+3) = ( Mass(k) * V_RB(3,l) * V_RB(3,l) - kT ) * InvMMNHC(1,m+3)
           else
             MyType = RBType(l)
             k = k + NumRBAtom(MyType)
             Gms(1,m+1) = ( MassRB(MyType) * V_RB(1,l) * V_RB(1,l) - kT ) * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( MassRB(MyType) * V_RB(2,l) * V_RB(2,l) - kT ) * InvMMNHC(1,m+2)
             Gms(1,m+3) = ( MassRB(MyType) * V_RB(3,l) * V_RB(3,l) - kT ) * InvMMNHC(1,m+3)
           end if
         end do

         Num = NfT

         do l = 1, NumRB
           if(QSingle(l)) cycle
           MyType = RBType(l)
           if(QLinear(l)) then
             Gms(1,Num+1) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Num = Num + 2
           else
             Gms(1,Num+1) = ( Lmoment(1,l) * Lmoment(1,l) * InvInertiaRB(1,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Gms(1,Num+3) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+3)
             Num = Num + 3
           end if
         end do

! -----------------------------------
! update the thermostat velocities
! -----------------------------------
         do k = 1 , NHchain-1
           do l = 1, NumMNHC
             aa         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l)   = ( VMNHC(k,l) * aa + w4 * Gms(k,l) ) * aa
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) * InvMMNHC(k+1,l)
           end do
         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

! ## flexible ## ------------------------------------------------------

   else

     j = 0
     do i = 1, N
       Gms(1,j+1) = ( Mass(i) * Vel(1,i) * Vel(1,i) - kT ) * InvMMNHC(1,j+1)
       Gms(1,j+2) = ( Mass(i) * Vel(2,i) * Vel(2,i) - kT ) * InvMMNHC(1,j+2)
       Gms(1,j+3) = ( Mass(i) * Vel(3,i) * Vel(3,i) - kT ) * InvMMNHC(1,j+3)
       j = j + 3
     end do

     do i = 1 , NHchain-1
       do j = 1, NumMNHC
         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)
       end do
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1
           do l = 1, NumMNHC
             aa = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + Gms(k,l) * w4 ) * aa
           end do
         end do

! --------------------------------
! update the particle velocities
! --------------------------------
         l = 0
         do k = 1, N
           Vel(1,k) = Vel(1,k) * exp( -w2 * VMNHC(1,l+1) )
           Vel(2,k) = Vel(2,k) * exp( -w2 * VMNHC(1,l+2) )
           Vel(3,k) = Vel(3,k) * exp( -w2 * VMNHC(1,l+3) )
           l = l + 3
         end do

! --------------------------------
! update the thermostat positions
! --------------------------------
         RMNHC = RMNHC + VMNHC * w2

! ------------------
! update the forces
! ------------------
         l = 0
         do k = 1, N
           Gms(1,l+1) = ( Mass(k) * Vel(1,k) * Vel(1,k) - kT ) * InvMMNHC(1,l+1)
           Gms(1,l+2) = ( Mass(k) * Vel(2,k) * Vel(2,k) - kT ) * InvMMNHC(1,l+2)
           Gms(1,l+3) = ( Mass(K) * Vel(3,k) * Vel(3,k) - kT ) * InvMMNHC(1,l+3)
           l = l + 3
         end do

! -----------------------------------
! update the thermostat velocities
! -----------------------------------
         do k = 1 , NHchain-1
           do l = 1, NumMNHC
             aa         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + w4 * Gms(k,l) ) * aa
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) * InvMMNHC(k+1,l)
           end do
         end do

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

       end do

     end do

   end if

end subroutine Thermostat_MNHC


!#####################################################################
!#####################################################################


subroutine Bath_NP_SHAKE(ROLL,iflag)

use CommonBlocks, only : QThermostat, cThermostatMethod, cBarostatMethod

logical :: ROLL
integer :: iflag
real(8), dimension(3,3) :: VelRotation
real(8), dimension(3)   :: Vscale3
real(8)                 :: Vscale

if(iflag==1) then

   if( QThermostat ) then

     if( cThermostatMethod == 'NH' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR_NH(VelRotation,1)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3_NH(Vscale3,1)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN_NH(Vscale,1)
       end if

     else if( cThermostatMethod == 'NHC' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR(VelRotation,1)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3(Vscale3,1)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN(Vscale,1)
       end if

     end if

   else

     if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
       call BarostatPR(VelRotation,1)
     else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
       call BarostatA3(Vscale3,1)
     else if( cBarostatMethod == 'AN' ) then
       call BarostatAN(Vscale,1)
     end if

   end if

else if(iflag==2) then

   if( QThermostat ) then

     if( cThermostatMethod == 'NH' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR_NH(VelRotation,1)
         call RATTLEROLLPR(VelRotation,ROLL)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3_NH(Vscale3,1)
         call RATTLEROLLA3(Vscale3,ROLL)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN_NH(Vscale,1)
         call RATTLEROLL(Vscale,ROLL)
       end if

     else if( cThermostatMethod == 'NHC' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR(VelRotation,1)
         call RATTLEROLLPR(VelRotation,ROLL)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3(Vscale3,1)
         call RATTLEROLLA3(Vscale3,ROLL)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN(Vscale,1)
         call RATTLEROLL(Vscale,ROLL)
       end if

     end if

   else

     if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
       call BarostatPR(VelRotation,1)
       call RATTLEROLLPR(VelRotation,ROLL)
     else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
       call BarostatA3(Vscale3,1)
       call RATTLEROLLA3(Vscale3,ROLL)
     else if( cBarostatMethod == 'AN' ) then
       call BarostatAN(Vscale,1)
       call RATTLEROLL(Vscale,ROLL)
     end if

   end if

end if

end subroutine Bath_NP_SHAKE


!#####################################################################
!#####################################################################


subroutine Bath_NP(iflag)

use CommonBlocks, only : QThermostat, cThermostatMethod, cBarostatMethod, QRigidBody
use TimeParam, only : lk
use BathParam, only : gkT, Baro_kin, Vg
use ThermoData, only : Ene_kin
use Configuration, only : Vel
use RBparam, only : V_RB, Lmoment

integer :: iflag
real(8), dimension(3,3) :: Dummy33
real(8), dimension(3)   :: Dummy3
real(8)                 :: Dummy, VelScale

if(iflag==1) then

   if( QThermostat ) then

     if( cThermostatMethod == 'NH' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR_NH(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3_NH(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN_NH(Dummy,lk)
       end if

     else if( cThermostatMethod == 'NHC' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN(Dummy,lk)
       end if

     else if( cThermostatMethod == 'MNHC' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR_MNHC(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BaroThermostatA3_MNHC(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN_MNHC(Dummy,lk)
       end if

     else if( cThermostatMethod == 'VSCALE' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BarostatPR(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' )) then
         call BarostatA3(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BarostatAN(Dummy,lk)
       end if

     end if

   else

     if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
       call BarostatPR(Dummy33,lk)
     else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
       call BarostatA3(Dummy3,lk)
     else if( cBarostatMethod == 'AN' ) then
       call BarostatAN(Dummy,lk)
     end if

   end if

else if(iflag==2) then

   if( QThermostat ) then

     if( cThermostatMethod == 'NH' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR_NH(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
         call BaroThermostatA3_NH(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN_NH(Dummy,lk)
       end if

     else if( cThermostatMethod == 'NHC' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
         call BaroThermostatA3(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN(Dummy,lk)
       end if

     else if( cThermostatMethod == 'MNHC' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BaroThermostatPR_MNHC(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
         call BaroThermostatA3_MNHC(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BaroThermostatAN_MNHC(Dummy,lk)
       end if

     else if( cThermostatMethod == 'VSCALE' ) then

       if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
         call BarostatPR(Dummy33,lk)
       else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
         call BarostatA3(Dummy3,lk)
       else if( cBarostatMethod == 'AN' ) then
         call BarostatAN(Dummy,lk)
       end if

       call CalcTemp
       call BathTemp

       VelScale = sqrt( gkT / (Ene_kin + Baro_kin) )

       if(QRigidBody) then
         V_RB    = V_RB    * VelScale
         Lmoment = Lmoment * VelScale
       else
         Vel = Vel * VelScale
       end if

       Vg = Vg * VelScale

     end if

   else

     if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
       call BarostatPR(Dummy33,lk)
     else if( ( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
       call BarostatA3(Dummy3,lk)
     else if( cBarostatMethod == 'AN' ) then
       call BarostatAN(Dummy,lk)
     end if

   end if

end if

end subroutine Bath_NP


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** The Andersen barostat                                 **
! ** Nose-Hoover chain ( length = NHchain )                **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatAN(Scale,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, QPathInt
use BathParam
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use CellParam, only : Volume
use ThermoData, only : Ene_kin, Ene_kinT, Virial

implicit NONE

integer :: i, j, k, ll
real(8) :: Scale, aa, bb, w2, w4, w8, odnf, Ekin
real(8), dimension(NHchain) :: Gms
real(8)                     :: Gmp
real(8) :: TraceVirial, P_o

   TraceVirial = Virial(1,1) + Virial(2,2) + Virial(3,3)
   P_o         = Pressure_o * 3.d0 * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Scale=1.d0

! ------------------------------
   if(QPathInt) then

     Ekin = 0.d0

     do i = 1, N

       Ekin = Ekin + FictMass(i,1) *             &
       &      dot_product( Vnm(:,i,1), Vnm(:,i,1) )

     end do

   else

     call CalcTemp          ! total kinetic energy

     Ekin = Ene_kin

   end if
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 / dble(NfT)
     Gmp = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

     Gms(1) = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

     do i = 1 , NHchain-1

       Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)

     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

         do k = NHchain-1 , 1 , -1

           aa     = exp( -w8 * Vss(k+1) )
           Vss(k) = ( Vss(k) * aa + Gms(k) * w4 ) * aa

         end do

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

! ## update the particle velocities
         aa      = exp( -w2 * Vss(1) )
         bb      = exp( -w2 * odnf * Vg(1,1) )
         V_RB    = V_RB * aa *  bb

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * aa
             Lmoment(3,k) = Lmoment(3,k) * aa
           else
             Lmoment(:,k) = Lmoment(:,k) * aa
           end if

         end do

!   ---------------------
         call CalcTemp
!   ---------------------
         Gmp     = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

         Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

! ## update the forces
         Gms(1) = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           aa       = exp( -w8 * Vss(k+1) )
           Vss(k)   = ( Vss(k) * aa + w4 * Gms(k) ) * aa
           Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)

         end do

         Vss(NHchain)=Vss(NHchain)+Gms(NHchain)*w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 * InvNf
     Gmp = ( odnf * Ekin + TraceVirial - P_o ) * InvMp(1,1)

     Gms(1) = ( Baro_kin + Ekin - gkT ) * InvMts(1)

     do i = 1 , NHchain-1

       Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)

     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

         do k = NHchain-1 , 1 , -1

           aa     = exp( -w8 * Vss(k+1) )
           Vss(k) = ( Vss(k) * aa + Gms(k) * w4 ) * aa

         end do

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

! ## update the particle velocities
         aa     = exp( -w2 * ( Vss(1) + odnf * Vg(1,1) ) )
         Scale  = Scale * aa
         Ekin   = Ekin * aa * aa
         Gmp    = ( odnf * Ekin + TraceVirial - P_o ) * InvMp(1,1)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

         Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

! ## update the forces
         Gms(1) = ( Baro_kin + Ekin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           aa       = exp( -w8 * Vss(k+1) )
           Vss(k)   = ( Vss(k) * aa + w4 * Gms(k) ) * aa
           Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)

         end do

         Vss(NHchain)=Vss(NHchain)+Gms(NHchain)*w4

       end do

     end do

! ---------------------------------
! update the particle velocities
! ---------------------------------
     if(QPathInt) then

       Vnm(:,:,1) = Vnm(:,:,1) * Scale

     else

       Vel = Vel * Scale

     end if

   end if

end subroutine BaroThermostatAN


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** Parrinello-Rahman barostat ( rectangular cell )       **
! ** Nose-Hoover chain ( length = NHchain )                **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatA3(Vscale,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, QPathInt, cBarostatMethod
use BathParam
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use CellParam, only : Volume
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, ll, i1, i2
real(8) :: aa, w2, w4, w8, odnf
real(8), dimension(NHchain) :: Gms
real(8), dimension(3)     :: Gmp, cf
real(8) :: P_o
real(8), dimension(3) :: Vscale
real(8) :: TrVg, dterm
real(8), dimension(3) :: Sckin
real(8) :: Ekin
logical :: QIxy

   if(cBarostatMethod=='A2') then
     QIxy = .True.
     i1   = CoupleEdge(1)
     i2   = CoupleEdge(2)
   else
     QIxy = .False.
   end if

   P_o = Pressure_o * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Vscale = 1.d0

! ------------------------------
   if(QPathInt) then

     Sckin = 0.d0
     do i = 1, N
       Sckin(:) = Sckin(:) + FictMass(i,1) * Vnm(:,i,1) * Vnm(:,i,1)
     end do

     Ekin = Sckin(1) + Sckin(2) + Sckin(3)

   else

     call CalcTemp          ! total kinetic energy

     do i = 1 , 3
       Sckin(i) = Pkinp(i,i)
     end do

     Ekin = Ene_kin

   end if
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = 0.d0
     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)
     end do

     odnf = 1.d0 / dble(NfT)

     do i = 1 , 3
       Gmp(i) = ( Pkinp(i,i) + odnf * Ene_kinT + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

     Gms(1) = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

     do i = 1 , NHchain-1
       Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

         do k = NHchain-1 , 1 , -1
           aa     = exp( -w8 * Vss(k+1) )
           Vss(k) = ( Vss(k) * aa + Gms(k) * w4 ) * aa
         end do

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         dterm = TrVg + Vss(1)
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( dterm + Vg(k,k) ) )
         end do

         do k = 1, NumRB
           V_RB(:,k) = V_RB(:,k) * cf(:)
         end do

         aa = exp( -w2 * Vss(1) )

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * aa
             Lmoment(3,k) = Lmoment(3,k) * aa
           else
             Lmoment(:,k) = Lmoment(:,k) * aa
           end if

         end do

!       ---------------
         call CalcTemp
!       ---------------

         do k = 1, 3
           Gmp(k) = ( Pkinp(k,k) + odnf * Ene_kinT + Virial(k,k) - P_o ) * InvMp(k,k)
         end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         Baro_kin = 0.d0
         do k = 1 , 3
           Baro_kin = Baro_kin + Mp(k,k) * Vg(k,k) * Vg(k,k)
         end do

! ## update the forces
         Gms(1) = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         do k = 1 , NHchain-1
           aa       = exp( -w8 * Vss(k+1) )
           Vss(k)   = ( Vss(k) * aa + w4 * Gms(k) ) * aa
           Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)
         end do

         Vss(NHchain)=Vss(NHchain)+Gms(NHchain)*w4

       end do

     end do

! ## flexible ## -------------------------------------------------------
   else

     Baro_kin = 0.d0

     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)
     end do

     odnf = InvNf

     do i = 1 , 3
       Gmp(i) = ( Sckin(i) + odnf * Ekin + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

     Gms(1) = ( Baro_kin + Ekin - gkT ) * InvMts(1)

     do i = 1 , NHchain-1
       Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

         do k = NHchain-1 , 1 , -1
           aa     = exp( -w8 * Vss(k+1) )
           Vss(k) = ( Vss(k) * aa + Gms(k) * w4 ) * aa
         end do

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         dterm = TrVg + Vss(1)
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( dterm + Vg(k,k) ) )
         end do

         Vscale   = Vscale * cf

         if(QPathInt) then

           do k = 1, N
             Vnm(:,k,1) = Vnm(:,k,1) * cf
           end do

           Sckin = 0.d0
           do k = 1, N
             Sckin(:) = Sckin(:) + FictMass(k,1) * Vnm(:,k,1) * Vnm(:,k,1)
           end do
           Ekin = Sckin(1) + Sckin(2) + Sckin(3)

         else

           do k = 1, N
             Vel(:,k) = Vel(:,k) * cf(:)
           end do

!         ---------------
           call CalcTemp
!         ---------------

           do k = 1 , 3
             Sckin(k) = Pkinp(k,k)
           end do

           Ekin = Ene_kin

         end if

         do k = 1, 3
           Gmp(k) = ( Sckin(k) + odnf * Ekin + Virial(k,k) - P_o ) * InvMp(k,k)
         end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         Baro_kin = 0.d0
         do k = 1 , 3
           Baro_kin = Baro_kin + Mp(k,k) * Vg(k,k) * Vg(k,k)
         end do

! ## update the forces
         Gms(1) = ( Baro_kin + Ekin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         do k = 1 , NHchain-1
           aa       = exp( -w8 * Vss(k+1) )
           Vss(k)   = ( Vss(k) * aa + w4 * Gms(k) ) * aa
           Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)
         end do

         Vss(NHchain)=Vss(NHchain)+Gms(NHchain)*w4

       end do

     end do

   end if

   if(QIxy) Vg(i2,i2) = Vg(i1,i1) ! to reduce the numerical error

end subroutine BaroThermostatA3


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** Parrinello-Rahman barostat ( fully flexible cell )    **
! ** Nose-Hoover chain ( length = NHchain )                **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatPR(VelRotation,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, QPathInt, cBarostatMethod
use BathParam
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use CellParam, only : H, Volume
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, l, ll, ii
real(8) :: cf, w2, w4, w8
real(8), dimension(3) :: Scale, tempov, eigenValue
real(8), dimension(NHchain) :: Gms
real(8), dimension(3,3) :: Gmp
real(8), dimension(3,3) :: RotV, Imatrix
real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8) :: TrVg, diagTerm, Pext
real(8), dimension(3,3) :: VelRotation, tempoRot
! ## Stress >>
real(8), dimension(3,3) :: StressTerm, Htrans, tempM
! ## << Stress
real(8), dimension(3,3) :: Sckin
real(8) :: Ekin

   Imatrix = 0.d0

   do i = 1 , 3

     Imatrix(i,i) = 1.d0

   end do

   VelRotation = Imatrix

   Pext = Pressure_o * Volume

   if(cBarostatMethod == 'ST') then

     Htrans = transpose( H )

     tempM = matmul( H, SigmaS )

     StressTerm = - matmul( tempM, Htrans )

   else

     StressTerm = 0.d0

   end if

! ------------------------------
   if(QPathInt) then

     Sckin = 0.d0

     do i = 1, N

       do j = 1, 3
         do k = j, 3

           Sckin(j,k) = Sckin(j,k) + FictMass(i,1) * Vnm(j,i,1) * Vnm(k,i,1)

         end do
       end do

     end do

     Sckin(2,1) = Sckin(1,2)
     Sckin(3,1) = Sckin(1,3)
     Sckin(3,2) = Sckin(2,3)

     Ekin = Sckin(1,1) + Sckin(2,2) + Sckin(3,3)

   else

     call CalcTemp          ! total kinetic energy

     do i = 1 , 3
       do j = 1 , 3

         Sckin(i,j) = Pkinp(i,j)

       end do
     end do

     Ekin = Ene_kin

   end if
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = 0.d0

     do j = 1 , 3

       do i = 1 , 3

         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)

       end do

     end do

     diagTerm = Ene_kinT / dble(NfT) - Pext

     Gmp = (Pkinp + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

     Gms(1) = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

     do i = 1 , NHchain-1

       Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)

     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

         do k = NHchain-1 , 1 , -1

           cf     = exp( -w8 * Vss(k+1) )
           Vss(k) = ( Vss(k) * cf + Gms(k) * w4 ) * cf

         end do

! ## update dlog(v)/dt

         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

! ## update the particle velocities
         TrVg     = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         diagTerm = TrVg + Vss(1)
         RotV   = Vg + diagTerm * Imatrix
!       ----------------------------------------
         call Jacobi(RotV,eigenVec,eigenValue)
!       ----------------------------------------

         TeigenVec = Transpose( eigenVec )

         Scale = exp( -w2 * eigenValue )

         do k = 1 , NumRB

           tempov   = matmul( TeigenVec, V_RB(:,k) )
           tempov   = tempov * Scale

           V_RB(:,k) = matmul( eigenVec, tempov )

         end do

         cf = exp( -w2 * Vss(1) )

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * cf
             Lmoment(3,k) = Lmoment(3,k) * cf
           else
             Lmoment(:,k) = Lmoment(:,k) * cf
           end if

         end do

!       ---------------
         call CalcTemp
!       ---------------

         diagTerm = Ene_kinT / dble(NfT) - Pext

         Gmp = (Pkinp + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

         Baro_kin = 0.d0

         do k = 1 , 3

           do l = 1 , 3

             Baro_kin = Baro_kin + Mp(k,l) * Vg(k,l) * Vg(k,l)

           end do

         end do

! ## update the forces
         Gms(1) = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           cf       = exp( -w8 * Vss(k+1) )
           Vss(k)   = ( Vss(k) * cf + w4 * Gms(k) ) * cf
           Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)

         end do

         Vss(NHchain)=Vss(NHchain)+Gms(NHchain)*w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = 0.d0

     do j = 1 , 3

       do i = 1 , 3

         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)

       end do

     end do

     diagTerm = Ekin * InvNf - Pext

     Gmp = (Sckin + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

     Gms(1) = ( Baro_kin + Ekin - gkT ) * InvMts(1)

     do i = 1 , NHchain-1

       Gms(i+1) = ( Mts(i) * Vss(i) * Vss(i) - kT ) * InvMts(i+1)

     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(NHchain) = Vss(NHchain) + Gms(NHchain) * w4

         do k = NHchain-1 , 1 , -1

           cf     = exp( -w8 * Vss(k+1) )
           Vss(k) = ( Vss(k) * cf + Gms(k) * w4 ) * cf

         end do

! ## update dlog(v)/dt

         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

! ## update the particle velocities
         TrVg     = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         diagTerm = TrVg + Vss(1)
         RotV   = Vg + diagTerm * Imatrix
!       ----------------------------------------
         call Jacobi(RotV,eigenVec,eigenValue)
!       ----------------------------------------

         TeigenVec = Transpose( eigenVec )

         Scale = exp( -w2 * eigenValue )

         if(QPathInt) then

           do k = 1 , N

             tempov   = matmul( TeigenVec, Vnm(:,k,1) )
             tempov   = tempov * Scale

             Vnm(:,k,1) = matmul( eigenVec, tempov )

           end do

           tempoRot = matmul( TeigenVec, VelRotation )

           do k = 1 , 3

             tempoRot(:,k) = tempoRot(:,k) * Scale

           end do

           VelRotation = matmul( eigenVec, tempoRot )

           Sckin = 0.d0

           do ii = 1, N

             do k = 1, 3
               do l = k, 3

                 Sckin(k,l) = Sckin(k,l) + FictMass(ii,1) * Vnm(k,ii,1) * Vnm(l,ii,1)

               end do
             end do

           end do

           Sckin(2,1) = Sckin(1,2)
           Sckin(3,1) = Sckin(1,3)
           Sckin(3,2) = Sckin(2,3)

           Ekin = Sckin(1,1) + Sckin(2,2) + Sckin(3,3)

         else

           do k = 1 , N

             tempov   = matmul( TeigenVec, Vel(:,k) )
             tempov   = tempov * Scale

             Vel(:,k) = matmul( eigenVec, tempov )

           end do

           tempoRot = matmul( TeigenVec, VelRotation )

           do k = 1 , 3

             tempoRot(:,k) = tempoRot(:,k) * Scale

           end do

           VelRotation = matmul( eigenVec, tempoRot )

!         ---------------
           call CalcTemp
!         ---------------

           Sckin = Pkinp
           Ekin = Ene_kin

         end if

         diagTerm = Ekin * InvNf - Pext

         Gmp = (Sckin + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

         Baro_kin = 0.d0

         do k = 1 , 3

           do l = 1 , 3

             Baro_kin = Baro_kin + Mp(k,l) * Vg(k,l) * Vg(k,l)

           end do

         end do

! ## update the forces
         Gms(1) = ( Baro_kin + Ekin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           cf       = exp( -w8 * Vss(k+1) )
           Vss(k)   = ( Vss(k) * cf + w4 * Gms(k) ) * cf
           Gms(k+1) = ( Mts(k) * Vss(k) * Vss(k) - kT ) * InvMts(k+1)

         end do

         Vss(NHchain)=Vss(NHchain)+Gms(NHchain)*w4

       end do

     end do

   end if

end subroutine BaroThermostatPR


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion of bath        **
! ** time evolution of barostat                            **
! ** the Andersen barostat ( cubic cell )                  **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BarostatAN(Scale,ll)

use Numbers, only : NfT
use CommonBlocks, only : QRigidBody
use BathParam
use RBparam, only : V_RB
use Configuration, only : Vel
use CellParam, only : Volume
use ThermoData, only : Ene_kin, Ene_kinT, Virial

implicit NONE

integer :: i, j, ll
real(8) :: Scale, aa, bb, w2, w4, odnf
real(8) :: Gmp
real(8) :: TraceVirial, P_o

   Scale = 1.d0

   TraceVirial = Virial(1,1) + Virial(2,2) + Virial(3,3)
   P_o         = Pressure_o * 3.d0 * Volume

! ------------------------------
   call CalcTemp          ! total kinetic energy
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 / dble(NfT)
     Gmp = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)

! ## update dlog(v)/dt

         Vg(1,1) = Vg(1,1) + Gmp * w4

! ## update the particle velocities
         bb      = exp( -w2 * odnf * Vg(1,1) )
         V_RB    = V_RB * bb

!   ---------------------
         call CalcTemp
!   ---------------------
         Gmp     = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

! ## update dlog(v)/dt
         Vg(1,1) = Vg(1,1) + Gmp * w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     odnf=1.d0 + 3.d0 * InvNf
     Gmp = ( odnf * Ene_kin + TraceVirial - P_o ) * InvMp(1,1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)

! ## update dlog(v)/dt

         Vg(1,1) = Vg(1,1) + Gmp * w4

! ## update the particle velocities
         aa      = exp( -w2 * odnf * Vg(1,1) )
         Scale   = Scale * aa
         Ene_kin = Ene_kin * aa * aa
         Gmp     = ( odnf * Ene_kin + TraceVirial - P_o ) * InvMp(1,1)

! ## update dlog(v)/dt
         Vg(1,1) = Vg(1,1) + Gmp * w4

       end do

     end do

! ---------------------------------
! update the particle velocities
! ---------------------------------
     Vel = Vel * Scale

   end if

end subroutine BarostatAN


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion for bath       **
! ** time evolution of the barostat                        **
! ** Parrinello-Rahman barostat ( rectangular cell )       **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BarostatA3(Vscale,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, cBarostatMethod
use BathParam
use RBparam, only : NumRB, V_RB
use Configuration, only : Vel
use CellParam, only : Volume
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, ll, i1, i2
real(8) :: w2, w4, odnf
real(8), dimension(3) :: Gmp, cf
real(8) :: P_o
real(8), dimension(3) :: Vscale
real(8) :: TrVg, dterm
logical :: QIxy

   if(cBarostatMethod=='A2') then
     QIxy = .True.
     i1   = CoupleEdge(1)
     i2   = CoupleEdge(2)
   else
     QIxy = .False.
   end if

   P_o = Pressure_o * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Vscale = 1.d0

! ------------------------------
   call CalcTemp          ! total kinetic energy
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     odnf = 1.d0 / dble(NfT)

     do i = 1 , 3
       Gmp(i) = ( Pkinp(i,i) + odnf * Ene_kinT + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)

! ## update dlog(v)/dt

         do k = 1 , 3
           Vg(k,k) = Vg(k,k) + Gmp(k) * w4
         end do

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         dterm = TrVg
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( dterm + Vg(k,k) ) )
         end do

         do k = 1, NumRB
           V_RB(:,k) = V_RB(:,k) * cf(:)
         end do

!       ---------------
         call CalcTemp
!       ---------------

         do k = 1, 3
           Gmp(k) = ( Pkinp(k,k) + odnf * Ene_kinT + Virial(k,k) - P_o ) * InvMp(k,k)
         end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update dlog(v)/dt

         do k = 1 , 3
           Vg(k,k) = Vg(k,k) + Gmp(k) * w4
         end do

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     odnf = InvNf

     do i = 1 , 3
       Gmp(i) = ( Pkinp(i,i) + odnf * Ene_kin + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if
! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)

! ## update dlog(v)/dt

         do k = 1 , 3
           Vg(k,k) = Vg(k,k) + Gmp(k) * w4
         end do

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         dterm = TrVg
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( dterm + Vg(k,k) ) )
         end do

         Vscale   = Vscale * cf

         do k = 1, N
           Vel(:,k) = Vel(:,k) * cf(:)
         end do

!       ---------------
         call CalcTemp
!       ---------------

         do k = 1, 3
           Gmp(k) = ( Pkinp(k,k) + odnf * Ene_kin + Virial(k,k) - P_o ) * InvMp(k,k)
         end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update dlog(v)/dt

         do k = 1 , 3
           Vg(k,k) = Vg(k,k) + Gmp(k) * w4
         end do

       end do

     end do

   end if

   if(QIxy) Vg(i2,i2) = Vg(i1,i1)

end subroutine BarostatA3


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion                **
! ** time evolution of the barostat                        **
! ** the Parrinello-Rahman barostat ( fully flexible cell )**
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BarostatPR(VelRotation,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, cBarostatMethod
use BathParam
use RBparam, only : NumRB, V_RB
use Configuration, only : Vel
use CellParam, only : H, Volume
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, ll
real(8) :: w2, w4
real(8), dimension(3) :: Scale, tempov, eigenValue
real(8), dimension(3,3) :: Gmp
real(8), dimension(3,3) :: RotV, Imatrix
real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8) :: TrVg, diagTerm, Pext
real(8), dimension(3,3) :: VelRotation, tempoRot
! ## Stress >>
real(8), dimension(3,3) :: StressTerm, Htrans, tempM

   Imatrix = 0.d0

   do i = 1 , 3

     Imatrix(i,i) = 1.d0

   end do

   VelRotation = Imatrix

   Pext = Pressure_o * Volume

   if(cBarostatMethod == 'ST') then

     Htrans = transpose( H )

     tempM = matmul( H, SigmaS )

     StressTerm = - matmul( tempM, Htrans )

   else

     StressTerm = 0.d0

   end if

! ------------------------------
   call CalcTemp          ! total kinetic energy
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     diagTerm = Ene_kinT / dble(NfT) - Pext

     Gmp = (Pkinp + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)

! ## update dlog(v)/dt

         Vg  = Vg + Gmp * w4

! ## update the particle velocities
         TrVg     = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         diagTerm = TrVg
         RotV   = Vg + diagTerm * Imatrix
!       ----------------------------------------
         call Jacobi(RotV,eigenVec,eigenValue)
!       ----------------------------------------

         TeigenVec = Transpose( eigenVec )

         Scale = exp( -w2 * eigenValue )

         do k = 1 , NumRB

           tempov   = matmul( TeigenVec, V_RB(:,k) )
           tempov   = tempov * Scale

           V_RB(:,k) = matmul( eigenVec, tempov )

         end do

!       ---------------
         call CalcTemp
!       ---------------

         diagTerm = Ene_kinT / dble(NfT) - Pext

         Gmp = (Pkinp + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update dlog(v)/dt
         Vg  = Vg + Gmp * w4

       end do

     end do

   else

     diagTerm = Ene_kin * InvNf - Pext

     Gmp = (Pkinp + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)

! ## update dlog(v)/dt

         Vg  = Vg + Gmp * w4

! ## update the particle velocities
         TrVg     = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         diagTerm = TrVg
         RotV   = Vg + diagTerm * Imatrix
!       ----------------------------------------
         call Jacobi(RotV,eigenVec,eigenValue)
!       ----------------------------------------

         TeigenVec = Transpose( eigenVec )

         Scale = exp( -w2 * eigenValue )

         do k = 1 , N

           tempov   = matmul( TeigenVec, Vel(:,k) )
           tempov   = tempov * Scale

           Vel(:,k) = matmul( eigenVec, tempov )

         end do

         tempoRot = matmul( TeigenVec, VelRotation )

         do k = 1 , 3

           tempoRot(:,k) = tempoRot(:,k) * Scale

         end do

         VelRotation = matmul( eigenVec, tempoRot )

!       ---------------
         call CalcTemp
!       ---------------

         diagTerm = Ene_kin * InvNf - Pext

         Gmp = (Pkinp + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update dlog(v)/dt
         Vg  = Vg + Gmp * w4

       end do

     end do

   end if

end subroutine BarostatPR


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** The Andersen barostat                                 **
! ** Massive Nose-Hoover chain ( length = NHchain )        **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatAN_MNHC(Scale,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, QPathInt
use BathParam
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear, RBType, &
&   NumRBAtom, MassRB, InvInertiaRB
use Configuration, only : Vel
use CellParam, only : Volume
use AtomParam, only : Mass
use ThermoData, only : Ene_kin, Ene_kinT, Virial

implicit NONE

integer :: i, j, k, l, ll, Num, MyType, m
real(8) :: Scale, aa, bb, w2, w4, w8, odnf, Ekin
real(8), dimension(NHchain,NumMNHC) :: Gms
real(8) :: Gmp
real(8) :: TraceVirial, P_o

   TraceVirial = Virial(1,1) + Virial(2,2) + Virial(3,3)
   P_o         = Pressure_o * 3.d0 * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Scale=1.d0 ! Dummy

! ------------------------------
   if(QPathInt) then

     Ekin = 0.d0
     do i = 1, N
       Ekin = Ekin + FictMass(i,1) *             &
       &      dot_product( Vnm(:,i,1), Vnm(:,i,1) )
     end do

   else

     call CalcTemp          ! total kinetic energy
     Ekin = Ene_kin

   end if
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 / dble(NfT)
     Gmp = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

     k = 0
     do i = 1, NumRB
       j = (i-1) * 3
       if(QSingle(i)) then
         k = k + 1
         Gms(1,j+1) = ( Mass(k) * V_RB(1,i) * V_RB(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(k) * V_RB(2,i) * V_RB(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(k) * V_RB(3,i) * V_RB(3,i) - kT ) * InvMMNHC(1,j+3)
       else
         MyType = RBType(i)
         k = k + NumRBAtom(MyType)
         Gms(1,j+1) = ( MassRB(MyType) * V_RB(1,i) * V_RB(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( MassRB(MyType) * V_RB(2,i) * V_RB(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( MassRB(MyType) * V_RB(3,i) * V_RB(3,i) - kT ) * InvMMNHC(1,j+3)
       end if
     end do

     Num = NfT

     do i = 1, NumRB
       if(QSingle(i)) cycle
       MyType = RBType(i)
       if(QLinear(i)) then
         Gms(1,Num+1) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) &
         &              - kT ) * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) &
         &              - kT ) * InvMMNHC(1,Num+2)
         Num = Num + 2
       else
         Gms(1,Num+1) = ( Lmoment(1,i) * Lmoment(1,i) * InvInertiaRB(1,MyType) &
         &              - kT ) * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) &
         &              - kT ) * InvMMNHC(1,Num+2)
         Gms(1,Num+3) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) &
         &              - kT ) * InvMMNHC(1,Num+3)
         Num = Num + 3
       end if
     end do

     Gms(1,NumMNHC) = ( Baro_kin - kT ) * InvMMNHC(1,NumMNHC)

     do i = 1 , NHchain-1
       do j = 1, NumMNHC
         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)
       end do
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1
           do l = 1, NumMNHC
             aa = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + Gms(k,l) * w4 ) * aa
           end do
         end do

! ## update dlog(v)/dt

         aa  = exp( -w8 * VMNHC(1,NumMNHC) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

! ## update the particle velocities
         bb      = exp( -w2 * odnf * Vg(1,1) )
         V_RB    = V_RB *  bb

         do k = 1, NumRB
           l = (k-1) * 3
           V_RB(1,k) = V_RB(1,k) * exp( -w2 * VMNHC(1,l+1) )
           V_RB(2,k) = V_RB(2,k) * exp( -w2 * VMNHC(1,l+2) )
           V_RB(3,k) = V_RB(3,k) * exp( -w2 * VMNHC(1,l+3) )
         end do

         l = NfT

         do k = 1 , NumRB
           if(QSingle(k)) cycle
           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+2) )
             l = l + 2
           else
             Lmoment(1,k) = Lmoment(1,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+2) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+3) )
             l = l + 3
           end if
         end do

!   ---------------------
         call CalcTemp
!   ---------------------
         Gmp     = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

! ## update the thermostat positions
         RMNHC = RMNHC + VMNHC * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * VMNHC(1,NumMNHC) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

         Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

! ## update the forces
         k = 0
         do l = 1, NumRB
           m = (l-1) * 3
           if(QSingle(l)) then
             k = k + 1
             Gms(1,m+1) = ( Mass(k) * V_RB(1,l) * V_RB(1,l) - kT ) * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( Mass(k) * V_RB(2,l) * V_RB(2,l) - kT ) * InvMMNHC(1,m+1)
             Gms(1,m+3) = ( Mass(k) * V_RB(3,l) * V_RB(3,l) - kT ) * InvMMNHC(1,m+1)
           else
             MyType = RBType(l)
             k = k + NumRBAtom(MyType)
             Gms(1,m+1) = ( MassRB(MyType) * V_RB(1,l) * V_RB(1,l) - kT ) &
             &          * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( MassRB(MyType) * V_RB(2,l) * V_RB(2,l) - kT ) &
             &          * InvMMNHC(1,m+1)
             Gms(1,m+3) = ( MassRB(MyType) * V_RB(3,l) * V_RB(3,l) - kT ) &
             &          * InvMMNHC(1,m+1)
           end if
         end do

         Num = NfT
         do l = 1, NumRB
           if(QSingle(l)) cycle
           MyType = RBType(l)
           if(QLinear(l)) then
             Gms(1,Num+1) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Num = Num + 2
           else
             Gms(1,Num+1) = ( Lmoment(1,l) * Lmoment(1,l) * InvInertiaRB(1,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Gms(1,Num+3) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+3)
             Num = Num + 3
           end if
         end do

         Gms(1,NumMNHC) = ( Baro_kin - kT ) * InvMMNHC(1,NumMNHC)

! ## update the thermostat velocities
         do k = 1 , NHchain-1
           do l = 1, NumMNHC
             aa       = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + w4 * Gms(k,l) ) * aa
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) &
             &          * InvMMNHC(k+1,l)
           end do
         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 * InvNf
     Gmp = ( odnf * Ekin + TraceVirial - P_o ) * InvMp(1,1)

     j = 0
     if(QPathInt) then
       do i = 1, N
         Gms(1,j+1) = ( FictMass(i,1) * Vnm(1,i,1) * Vnm(1,i,1) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( FictMass(i,1) * Vnm(2,i,1) * Vnm(2,i,1) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( FictMass(i,1) * Vnm(3,i,1) * Vnm(3,i,1) - kT ) * InvMMNHC(1,j+3)
         j = j + 3
       end do
     else
       do i = 1, N
         Gms(1,j+1) = ( Mass(i) * Vel(1,i) * Vel(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(i) * Vel(2,i) * Vel(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(i) * Vel(3,i) * Vel(3,i) - kT ) * InvMMNHC(1,j+3)
         j = j + 3
       end do
     end if

     Gms(1,NumMNHC) = ( Baro_kin - kT ) * InvMMNHC(1,NumMNHC)

     do i = 1 , NHchain-1
       do j = 1, NumMNHC
         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)
       end do
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

         do k = NHchain-1 , 1 , -1
           do l = 1, NumMNHC
             aa = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + Gms(k,l) * w4 ) * aa
           end do
         end do

! ## update dlog(v)/dt

         aa  = exp( -w8 * VMNHC(1,NumMNHC) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

! ## update the particle velocities

         aa = exp( -w2 * odnf * Vg(1,1) )

         if(QPathInt) then
           l = 0
           do k = 1 , N
             Vnm(1,k,1) = Vnm(1,k,1) * exp( -w2 * VMNHC(1,l+1) ) * aa
             Vnm(2,k,1) = Vnm(2,k,1) * exp( -w2 * VMNHC(1,l+2) ) * aa
             Vnm(3,k,1) = Vnm(3,k,1) * exp( -w2 * VMNHC(1,l+3) ) * aa
             l = l + 3
           end do
           Ekin = 0.d0
           do k = 1, N
             Ekin = Ekin + FictMass(k,1) *             &
             &      dot_product( Vnm(:,k,1), Vnm(:,k,1) )
           end do
         else
           l = 0
           do k = 1 , N
             Vel(1,k) = Vel(1,k) * exp( -w2 * VMNHC(1,l+1) ) * aa
             Vel(2,k) = Vel(2,k) * exp( -w2 * VMNHC(1,l+2) ) * aa
             Vel(3,k) = Vel(3,k) * exp( -w2 * VMNHC(1,l+3) ) * aa
             l = l + 3
           end do
!         ---------------
           call CalcTemp
!         ---------------
           Ekin = Ene_kin
         end if

         Gmp     = ( odnf * Ekin + TraceVirial - P_o ) * InvMp(1,1)

! ## update the thermostat positions
         RMNHC = RMNHC + VMNHC * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * VMNHC(1,NumMNHC) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

         Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

! ## update the forces
         if(QPathInt) then

           do k = 1, N

             l = (k-1) * 3

             Gms(1,l+1) = ( FictMass(k,1) * Vnm(1,k,1) * Vnm(1,k,1) - kT ) &
             &          * InvMMNHC(1,l+1)
             Gms(1,l+2) = ( FictMass(k,1) * Vnm(2,k,1) * Vnm(2,k,1) - kT ) &
             &          * InvMMNHC(1,l+2)
             Gms(1,l+3) = ( FictMass(K,1) * Vnm(3,k,1) * Vnm(3,k,1) - kT ) &
             &          * InvMMNHC(1,l+3)

           end do

         else

           do k = 1, N

             l = (k-1) * 3

             Gms(1,l+1) = ( Mass(k) * Vel(1,k) * Vel(1,k) - kT ) * InvMMNHC(1,l+1)
             Gms(1,l+2) = ( Mass(k) * Vel(2,k) * Vel(2,k) - kT ) * InvMMNHC(1,l+2)
             Gms(1,l+3) = ( Mass(K) * Vel(3,k) * Vel(3,k) - kT ) * InvMMNHC(1,l+3)

           end do

         end if

         Gms(1,NumMNHC) = ( Baro_kin - kT ) * InvMMNHC(1,NumMNHC)

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           do l = 1, NumMNHC

             aa         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + w4 * Gms(k,l) ) * aa
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) &
             &          * InvMMNHC(k+1,l)

           end do

         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

   end if

end subroutine BaroThermostatAN_MNHC


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** Parrinello-Rahman barostat ( rectangular cell )       **
! ** Massive Nose-Hoover chain ( length = NHchain )        **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatA3_MNHC(Vscale,ll)

use Numbers, only : N, NfT, NfR
use CommonBlocks, only : QRigidBody, QPathInt, cBarostatMethod
use BathParam
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear, RBType, &
&   NumRBAtom, MassRB, InvInertiaRB
use Configuration, only : Vel
use CellParam, only : Volume
use AtomParam, only : Mass
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, ll, l, MyType, Num, m, i1, i2, i3
real(8) :: aa, w2, w4, w8, odnf
real(8), dimension(NHchain,NumMNHC) :: Gms
real(8), dimension(3)     :: Gmp, cf
real(8) :: P_o
real(8), dimension(3) :: Vscale
real(8) :: TrVg
real(8), dimension(3) :: Sckin
real(8) :: Ekin
logical :: QIxy

   if(cBarostatMethod=='A2') then
     QIxy = .True.
     i1   = CoupleEdge(1)
     i2   = CoupleEdge(2)
     do i = 1, 3
       if(i==i1.or.i==i2) cycle
       i3 = i
     end do
   else
     QIxy = .False.
   end if

   P_o = Pressure_o * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Vscale = 1.d0

! ------------------------------
   if(QPathInt) then

     Sckin = 0.d0
     do i = 1, N
       Sckin(:) = Sckin(:) + FictMass(i,1) * Vnm(:,i,1) * Vnm(:,i,1)
     end do

     Ekin = Sckin(1) + Sckin(2) + Sckin(3)

   else

     call CalcTemp          ! total kinetic energy

     do i = 1 , 3
       Sckin(i) = Pkinp(i,i)
     end do

     Ekin = Ene_kin

   end if
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = 0.d0
     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)
     end do

     odnf = 1.d0 / dble(NfT)

     do i = 1 , 3
       Gmp(i) = ( Pkinp(i,i) + odnf * Ene_kinT + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

     k = 0
     do i = 1, NumRB

       j = (i-1) * 3

       if(QSingle(i)) then

         k = k + 1
         Gms(1,j+1) = ( Mass(k) * V_RB(1,i) * V_RB(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(k) * V_RB(2,i) * V_RB(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(k) * V_RB(3,i) * V_RB(3,i) - kT ) * InvMMNHC(1,j+3)

       else

         MyType = RBType(i)
         k = k + NumRBAtom(MyType)

         Gms(1,j+1) = ( MassRB(MyType) * V_RB(1,i) * V_RB(1,i) - kT ) &
         &          * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( MassRB(MyType) * V_RB(2,i) * V_RB(2,i) - kT ) &
         &          * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( MassRB(MyType) * V_RB(3,i) * V_RB(3,i) - kT ) &
         &          * InvMMNHC(1,j+3)

       end if

     end do

     Num = NfT

     do i = 1, NumRB

       if(QSingle(i)) cycle

       MyType = RBType(i)

       if(QLinear(i)) then

         Gms(1,Num+1) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) - kT ) &
         &              * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) - kT ) &
         &              * InvMMNHC(1,Num+2)
         Num = Num + 2

       else

         Gms(1,Num+1) = ( Lmoment(1,i) * Lmoment(1,i) * InvInertiaRB(1,MyType) - kT ) &
         &              * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) - kT ) &
         &              * InvMMNHC(1,Num+2)
         Gms(1,Num+3) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) - kT ) &
         &              * InvMMNHC(1,Num+3)
         Num = Num + 3

       end if

     end do

     Num = NfT + NfR

     if(QIxy) then
       Gms(1,Num+1) = ( Mp(i1,i1) * Vg(i1,i1) * Vg(i1,i1) - kT ) * InvMMNHC(1,Num+1)
       Gms(1,Num+2) = ( Mp(i3,i3) * Vg(i3,i3) * Vg(i3,i3) - kT ) * InvMMNHC(1,Num+2)
     else
       do i = 1, 3
         Num = Num + 1
         Gms(1,Num) = ( Mp(i,i) * Vg(i,i) * Vg(i,i) - kT ) * InvMMNHC(1,Num)
       end do
     end if

     do i = 1 , NHchain-1
       do j = 1, NumMNHC
         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)
       end do
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1
           do l = 1, NumMNHC
             aa = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + Gms(k,l) * w4 ) * aa
           end do
         end do

! ## update dlog(v)/dt

         Num = NfT + NfR

         if(QIxy) then
           aa  = exp( -w8 * VMNHC(1,Num+1) )
           Vg(i1,i1) = (Vg(i1,i1) * aa + Gmp(i1) * w4 ) * aa
           Vg(i2,i2) = Vg(i1,i1)
           aa  = exp( -w8 * VMNHC(1,Num+2) )
           Vg(i3,i3) = (Vg(i3,i3) * aa + Gmp(i3) * w4 ) * aa
         else
           do k = 1, 3
             Num = Num + 1
             aa  = exp( -w8 * VMNHC(1,Num) )
             Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
           end do
         end if

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( TrVg + Vg(k,k) ) )
         end do

         do k = 1, NumRB
           l = 3 * (k-1)
           V_RB(1,k) = V_RB(1,k) * cf(1) * exp( -w2 * VMNHC(1,l+1) )
           V_RB(2,k) = V_RB(2,k) * cf(2) * exp( -w2 * VMNHC(1,l+2) )
           V_RB(3,k) = V_RB(3,k) * cf(3) * exp( -w2 * VMNHC(1,l+3) )
         end do

         l = NfT

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+2) )
             l = l + 2
           else
             Lmoment(1,k) = Lmoment(1,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+2) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+3) )
             l = l + 3
           end if

         end do

!       ---------------
         call CalcTemp
!       ---------------

          do k = 1, 3
            Gmp(k) = ( Pkinp(k,k) + odnf * Ene_kinT + Virial(k,k) - P_o ) * InvMp(k,k)
          end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update the thermostat positions
         RMNHC = RMNHC + VMNHC * w2

! ## update dlog(v)/dt
         Num = NfT + NfR

         if(QIxy) then
           aa  = exp( -w8 * VMNHC(1,Num+1) )
           Vg(i1,i1) = (Vg(i1,i1) * aa + Gmp(i1) * w4 ) * aa
           Vg(i2,i2) = Vg(i1,i1)
           aa  = exp( -w8 * VMNHC(1,Num+2) )
           Vg(i3,i3) = (Vg(i3,i3) * aa + Gmp(i3) * w4 ) * aa
         else
           do k = 1 , 3
             Num = Num + 1
             aa  = exp( -w8 * VMNHC(1,Num) )
             Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
           end do
         end if

         Baro_kin = 0.d0

         do k = 1 , 3
           Baro_kin = Baro_kin + Mp(k,k) * Vg(k,k) * Vg(k,k)
         end do

! ## update the forces
         k = 0

         do l = 1, NumRB

           m = (l-1) * 3

           if(QSingle(l)) then
             k = k + 1
             Gms(1,m+1) = ( Mass(k) * V_RB(1,l) * V_RB(1,l) - kT ) * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( Mass(k) * V_RB(2,l) * V_RB(2,l) - kT ) * InvMMNHC(1,m+2)
             Gms(1,m+3) = ( Mass(k) * V_RB(3,l) * V_RB(3,l) - kT ) * InvMMNHC(1,m+3)
           else
             MyType = RBType(l)
             k = k + NumRBAtom(MyType)
             Gms(1,m+1) = ( MassRB(MyType) * V_RB(1,l) * V_RB(1,l) - kT ) &
             &          * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( MassRB(MyType) * V_RB(2,l) * V_RB(2,l) - kT ) &
             &          * InvMMNHC(1,m+2)
             Gms(1,m+3) = ( MassRB(MyType) * V_RB(3,l) * V_RB(3,l) - kT ) &
             &          * InvMMNHC(1,m+3)
           end if

         end do

         Num = NfT

         do l = 1, NumRB

           if(QSingle(l)) cycle

           MyType = RBType(l)

           if(QLinear(l)) then
             Gms(1,Num+1) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Num = Num + 2
           else
             Gms(1,Num+1) = ( Lmoment(1,l) * Lmoment(1,l) * InvInertiaRB(1,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Gms(1,Num+3) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+3)
             Num = Num + 3
           end if

         end do

         Num = NfT + NfR

         if(QIxy) then
           Gms(1,Num+1) = ( Mp(i1,i1) * Vg(i1,i1) * Vg(i1,i1) - kT ) * InvMMNHC(1,Num+1)
           Gms(1,Num+2) = ( Mp(i3,i3) * Vg(i3,i3) * Vg(i3,i3) - kT ) * InvMMNHC(1,Num+2)
         else
           do k = 1, 3
             Num = Num + 1
             Gms(1,Num) = ( Mp(k,k) * Vg(k,k) * Vg(k,k) - kT ) * InvMMNHC(1,Num)
           end do
         end if

! ## update the thermostat velocities
         do k = 1 , NHchain-1
           do l = 1, NumMNHC
             aa         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l)   = ( VMNHC(k,l) * aa + w4 * Gms(k,l) ) * aa
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) * InvMMNHC(k+1,l)
           end do
         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = 0.d0
     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)
     end do

     odnf = InvNf

     do i = 1 , 3
       Gmp(i) = ( Sckin(i) + odnf * Ekin + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

     if(QPathInt) then

       do i = 1, N
         j = (i-1) * 3
         Gms(1,j+1) = ( FictMass(i,1) * Vnm(1,i,1) * Vnm(1,i,1) - kT ) &
         &          * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( FictMass(i,1) * Vnm(2,i,1) * Vnm(2,i,1) - kT ) &
         &          * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( FictMass(i,1) * Vnm(3,i,1) * Vnm(3,i,1) - kT ) &
         &          * InvMMNHC(1,j+3)
       end do

     else

       do i = 1, N
         j = (i-1) * 3
         Gms(1,j+1) = ( Mass(i) * Vel(1,i) * Vel(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(i) * Vel(2,i) * Vel(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(i) * Vel(3,i) * Vel(3,i) - kT ) * InvMMNHC(1,j+3)
       end do

     end if

     Num = 3 * N

     if(QIxy) then
       Gms(1,Num+1) = ( Mp(i1,i1) * Vg(i1,i1) * Vg(i1,i1) - kT ) * InvMMNHC(1,Num+1)
       Gms(1,Num+2) = ( Mp(i3,i3) * Vg(i3,i3) * Vg(i3,i3) - kT ) * InvMMNHC(1,Num+2)
     else
       do i = 1, 3
         Num = Num + 1
         Gms(1,Num) = ( Mp(i,i) * Vg(i,i) * Vg(i,i) - kT ) * InvMMNHC(1,Num)
       end do
     end if

     do i = 1 , NHchain-1
       do j = 1, NumMNHC
         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)
       end do
     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1
           do l = 1, NumMNHC
             aa = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + Gms(k,l) * w4 ) * aa
           end do
         end do

! ## update dlog(v)/dt

         Num = 3 * N

         if(QIxy) then
           aa  = exp( -w8 * VMNHC(1,Num+1) )
           Vg(i1,i1) = (Vg(i1,i1) * aa + Gmp(i1) * w4 ) * aa
           Vg(i2,i2) = Vg(i1,i1)
           aa  = exp( -w8 * VMNHC(1,Num+2) )
           Vg(i3,i3) = (Vg(i3,i3) * aa + Gmp(i3) * w4 ) * aa
         else
           do k = 1 , 3
             Num = Num + 1
             aa      = exp( -w8 * VMNHC(1,Num) )
             Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
           end do
         end if

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( TrVg + Vg(k,k) ) )
         end do

         if(QPathInt) then

           do k = 1, N
             l = 3 * (k-1)
             Vnm(1,k,1) = Vnm(1,k,1) * cf(1) * exp( -w2 * VMNHC(1,l+1) )
             Vnm(2,k,1) = Vnm(2,k,1) * cf(2) * exp( -w2 * VMNHC(1,l+2) )
             Vnm(3,k,1) = Vnm(3,k,1) * cf(3) * exp( -w2 * VMNHC(1,l+3) )
           end do

           Sckin = 0.d0
           do k = 1, N
             Sckin(:) = Sckin(:) + FictMass(k,1) * Vnm(:,k,1) * Vnm(:,k,1)
           end do

           Ekin = Sckin(1) + Sckin(2) + Sckin(3)

         else

           do k = 1, N
             l = 3 * (k-1)
             Vel(1,k) = Vel(1,k) * cf(1) * exp( -w2 * VMNHC(1,l+1) )
             Vel(2,k) = Vel(2,k) * cf(2) * exp( -w2 * VMNHC(1,l+2) )
             Vel(3,k) = Vel(3,k) * cf(3) * exp( -w2 * VMNHC(1,l+3) )
           end do

!         ---------------
           call CalcTemp
!         ---------------

           do k = 1 , 3
             Sckin(k) = Pkinp(k,k)
           end do

           Ekin = Ene_kin

         end if

          do k = 1, 3
            Gmp(k) = ( Sckin(k) + odnf * Ekin + Virial(k,k) - P_o ) * InvMp(k,k)
          end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update the thermostat positions
         RMNHC = RMNHC + VMNHC * w2

! ## update dlog(v)/dt
         Num = 3 * N

         if(QIxy) then
           aa  = exp( -w8 * VMNHC(1,Num+1) )
           Vg(i1,i1) = (Vg(i1,i1) * aa + Gmp(i1) * w4 ) * aa
           Vg(i2,i2) = Vg(i1,i1)
           aa  = exp( -w8 * VMNHC(1,Num+2) )
           Vg(i3,i3) = (Vg(i3,i3) * aa + Gmp(i3) * w4 ) * aa
         else
           do k = 1 , 3
             Num = Num + 1
             aa  = exp( -w8 * VMNHC(1,Num) )
             Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
           end do
         end if

         Baro_kin = 0.d0

         do k = 1 , 3
           Baro_kin = Baro_kin + Mp(k,k) * Vg(k,k) * Vg(k,k)
         end do

! ## update the forces
         if(QPathInt) then

           do k = 1, N
             l = (k-1) * 3
             Gms(1,l+1) = ( FictMass(k,1) * Vnm(1,k,1) * Vnm(1,k,1) - kT ) &
             &          * InvMMNHC(1,l+1)
             Gms(1,l+2) = ( FictMass(k,1) * Vnm(2,k,1) * Vnm(2,k,1) - kT ) &
             &          * InvMMNHC(1,l+2)
             Gms(1,l+3) = ( FictMass(K,1) * Vnm(3,k,1) * Vnm(3,k,1) - kT ) &
             &          * InvMMNHC(1,l+3)
           end do

         else

           do k = 1, N
             l = (k-1) * 3
             Gms(1,l+1) = ( Mass(k) * Vel(1,k) * Vel(1,k) - kT ) * InvMMNHC(1,l+1)
             Gms(1,l+2) = ( Mass(k) * Vel(2,k) * Vel(2,k) - kT ) * InvMMNHC(1,l+2)
             Gms(1,l+3) = ( Mass(K) * Vel(3,k) * Vel(3,k) - kT ) * InvMMNHC(1,l+3)
           end do

         end if

         Num = 3 * N

         if(QIxy) then
           Gms(1,Num+1) = ( Mp(i1,i1) * Vg(i1,i1) * Vg(i1,i1) - kT ) * InvMMNHC(1,Num+1)
           Gms(1,Num+2) = ( Mp(i3,i3) * Vg(i3,i3) * Vg(i3,i3) - kT ) * InvMMNHC(1,Num+2)
         else
           do k = 1, 3
             Num = Num + 1
             Gms(1,Num) = ( Mp(k,k) * Vg(k,k) * Vg(k,k) - kT ) * InvMMNHC(1,Num)
           end do
         end if

! ## update the thermostat velocities
         do k = 1 , NHchain-1
           do l = 1, NumMNHC
             aa         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * aa + w4 * Gms(k,l) ) * aa
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) &
             &          * InvMMNHC(k+1,l)
           end do
         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

   end if

end subroutine BaroThermostatA3_MNHC


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** Parrinello-Rahman barostat ( fully flexible cell )    **
! ** Massive Nose-Hoover chain ( length = NHchain )        **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatPR_MNHC(VelRotation,ll)

use Numbers, only : N, NfT, NfR
use CommonBlocks, only : QRigidBody, QPathInt, cBarostatMethod
use BathParam
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear, RBType, &
&   NumRBAtom, MassRB, InvInertiaRB
use Configuration, only : Vel
use CellParam, only : H, Volume
use AtomParam, only : Mass
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, l, ll, Num, MyType, m, ii
real(8) :: cf, w2, w4, w8
real(8), dimension(3) :: Scale, tempov, eigenValue
real(8), dimension(NHchain,NumMNHC) :: Gms
real(8), dimension(3,3) :: Gmp
real(8), dimension(3,3) :: RotV, RotVpre, Imatrix
real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8) :: TrVg, diagTerm, Pext, pref
real(8), dimension(3,3) :: VelRotation, tempoRot
! ## Stress >>
real(8), dimension(3,3) :: StressTerm, Htrans, tempM
! ## << Stress
real(8), dimension(3,3) :: Sckin
real(8) :: Ekin

   Imatrix = 0.d0

   do i = 1 , 3

     Imatrix(i,i) = 1.d0

   end do

   VelRotation = Imatrix

   Pext = Pressure_o * Volume

   if(cBarostatMethod == 'ST') then

     Htrans = transpose( H )

     tempM = matmul( H, SigmaS )

     StressTerm = - matmul( tempM, Htrans )

   else

     StressTerm = 0.d0

   end if

! ------------------------------
   if(QPathInt) then

     Sckin = 0.d0

     do i = 1, N

       do j = 1, 3
         do k = j, 3

           Sckin(j,k) = Sckin(j,k) + FictMass(i,1) * Vnm(j,i,1) * Vnm(k,i,1)

         end do
       end do

     end do

     Sckin(2,1) = Sckin(1,2)
     Sckin(3,1) = Sckin(1,3)
     Sckin(3,2) = Sckin(2,3)

     Ekin = Sckin(1,1) + Sckin(2,2) + Sckin(3,3)

   else

     call CalcTemp          ! total kinetic energy

     do i = 1 , 3
       do j = 1 , 3

         Sckin(i,j) = Pkinp(i,j)

       end do
     end do

     Ekin = Ene_kin

   end if
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = 0.d0

     do i = 1 , 3
       do j = 1 , 3
         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)
       end do
     end do

     diagTerm = Ene_kinT / dble(NfT) - Pext

     Gmp = (Pkinp + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

     k = 0

     do i = 1, NumRB

       j = (i-1) * 3

       if(QSingle(i)) then

         k = k + 1

         Gms(1,j+1) = ( Mass(k) * V_RB(1,i) * V_RB(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(k) * V_RB(2,i) * V_RB(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(k) * V_RB(3,i) * V_RB(3,i) - kT ) * InvMMNHC(1,j+3)

       else

         MyType = RBType(i)
         k = k + NumRBAtom(MyType)

         Gms(1,j+1) = ( MassRB(MyType) * V_RB(1,i) * V_RB(1,i) - kT ) &
         &          * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( MassRB(MyType) * V_RB(2,i) * V_RB(2,i) - kT ) &
         &          * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( MassRB(MyType) * V_RB(3,i) * V_RB(3,i) - kT ) &
         &          * InvMMNHC(1,j+3)

       end if

     end do

     Num = NfT

     do i = 1, NumRB

       if(QSingle(i)) cycle

       MyType = RBType(i)

       if(QLinear(i)) then

         Gms(1,Num+1) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) - kT ) &
         &              * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) - kT ) &
         &              * InvMMNHC(1,Num+2)
         Num = Num + 2

       else

         Gms(1,Num+1) = ( Lmoment(1,i) * Lmoment(1,i) * InvInertiaRB(1,MyType) - kT ) &
         &              * InvMMNHC(1,Num+1)
         Gms(1,Num+2) = ( Lmoment(2,i) * Lmoment(2,i) * InvInertiaRB(2,MyType) - kT ) &
         &              * InvMMNHC(1,Num+2)
         Gms(1,Num+3) = ( Lmoment(3,i) * Lmoment(3,i) * InvInertiaRB(3,MyType) - kT ) &
         &              * InvMMNHC(1,Num+3)
         Num = Num + 3

       end if

     end do

     Num = NfT + NfR

     do i = 1, 3

       do j = i, 3

         if(i==j) then
           pref = 1.d0
         else
           pref = 2.d0
         end if

         Num = Num + 1
         Gms(1,Num) = ( pref * Mp(i,j) * Vg(i,j) * Vg(i,j) - kT ) * InvMMNHC(1,Num)

       end do

     end do

     do i = 1 , NHchain-1

       do j = 1, NumMNHC

         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) * InvMMNHC(i+1,j)

       end do

     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1

           do l = 1, NumMNHC

             cf = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * cf + Gms(k,l) * w4 ) * cf

           end do

         end do

! ## update dlog(v)/dt

         Num = NfT + NfR

         do k = 1, 3

           do l = k, 3

             Num = Num + 1
             cf  = exp( -w8 * VMNHC(1,Num) )
             Vg(k,l) = (Vg(k,l) * cf + Gmp(k,l) * w4 ) * cf

           end do

         end do

         Vg(2,1) = Vg(1,2)
         Vg(3,1) = Vg(1,3)
         Vg(3,2) = Vg(2,3)

! ## update the particle velocities
         TrVg    = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         RotVpre = Vg + TrVg * Imatrix
         RotV    = RotVpre

         l = 0

         do k = 1, NumRB

           RotV(1,1) = RotVpre(1,1) + VMNHC(1,l+1)
           RotV(2,2) = RotVpre(2,2) + VMNHC(1,l+2)
           RotV(3,3) = RotVpre(3,3) + VMNHC(1,l+3)

           l = l + 3

!         ----------------------------------------
           call Jacobi(RotV,eigenVec,eigenValue)
!         ----------------------------------------

           TeigenVec = Transpose( eigenVec )

           Scale = exp( -w2 * eigenValue )

           tempov   = matmul( TeigenVec, V_RB(:,k) )
           tempov   = tempov * Scale

           V_RB(:,k) = matmul( eigenVec, tempov )

         end do

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then

             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+2) )

             l = l + 2

           else

             Lmoment(1,k) = Lmoment(1,k) * exp( -w2 * VMNHC(1,l+1) )
             Lmoment(2,k) = Lmoment(2,k) * exp( -w2 * VMNHC(1,l+2) )
             Lmoment(3,k) = Lmoment(3,k) * exp( -w2 * VMNHC(1,l+3) )

             l = l + 3

           end if

         end do

!       ---------------
         call CalcTemp
!       ---------------

         diagTerm = Ene_kinT / dble(NfT) - Pext

         Gmp = (Pkinp + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update the thermostat positions
         RMNHC = RMNHC + VMNHC * w2

! ## update dlog(v)/dt

         Num = NfT + NfR

         do k = 1, 3

           do l = k, 3

             Num = Num + 1
             cf  = exp( -w8 * VMNHC(1,Num) )
             Vg(k,l)  = (Vg(k,l) * cf + Gmp(k,l) * w4 ) * cf

           end do

         end do

         Vg(2,1) = Vg(1,2)
         Vg(3,1) = Vg(1,3)
         Vg(3,2) = Vg(2,3)

         Baro_kin = 0.d0

         do k = 1 , 3

           do l = 1 , 3

             Baro_kin = Baro_kin + Mp(k,l) * Vg(k,l) * Vg(k,l)

           end do

         end do

! ## update the forces
         k = 0

         do l = 1, NumRB

           m = (l-1) * 3

           if(QSingle(l)) then

             k = k + 1

             Gms(1,m+1) = ( Mass(k) * V_RB(1,l) * V_RB(1,l) - kT ) * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( Mass(k) * V_RB(2,l) * V_RB(2,l) - kT ) * InvMMNHC(1,m+2)
             Gms(1,m+3) = ( Mass(k) * V_RB(3,l) * V_RB(3,l) - kT ) * InvMMNHC(1,m+3)

           else

             MyType = RBType(l)
             k = k + NumRBAtom(MyType)

             Gms(1,m+1) = ( MassRB(MyType) * V_RB(1,l) * V_RB(1,l) - kT ) &
             &          * InvMMNHC(1,m+1)
             Gms(1,m+2) = ( MassRB(MyType) * V_RB(2,l) * V_RB(2,l) - kT ) &
             &          * InvMMNHC(1,m+2)
             Gms(1,m+3) = ( MassRB(MyType) * V_RB(3,l) * V_RB(3,l) - kT ) &
             &          * InvMMNHC(1,m+3)

           end if

         end do

         Num = NfT

         do l = 1, NumRB

           if(QSingle(l)) cycle

           MyType = RBType(l)

           if(QLinear(l)) then

             Gms(1,Num+1) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Num = Num + 2

           else

             Gms(1,Num+1) = ( Lmoment(1,l) * Lmoment(1,l) * InvInertiaRB(1,MyType) &
             &                - kT ) * InvMMNHC(1,Num+1)
             Gms(1,Num+2) = ( Lmoment(2,l) * Lmoment(2,l) * InvInertiaRB(2,MyType) &
             &                - kT ) * InvMMNHC(1,Num+2)
             Gms(1,Num+3) = ( Lmoment(3,l) * Lmoment(3,l) * InvInertiaRB(3,MyType) &
             &                - kT ) * InvMMNHC(1,Num+3)
             Num = Num + 3

           end if

         end do

         Num = NfT + NfR

         do k = 1, 3

           do l = k, 3

             if(k==l) then
               pref = 1.d0
             else
               pref = 2.d0
             end if

             Num = Num + 1
             Gms(1,Num) = ( pref * Mp(k,l) * Vg(k,l) * Vg(k,l) - kT ) &
             &          * InvMMNHC(1,Num)

           end do

         end do

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           do l = 1, NumMNHC

             cf         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * cf + w4 * Gms(k,l) ) * cf
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) &
             &          * InvMMNHC(k+1,l)

           end do

         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = 0.d0

     do i = 1 , 3

       do j = 1 , 3

         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)

       end do

     end do

     diagTerm = Ekin * InvNf - Pext

     Gmp = (Sckin + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

     if(QPathInt) then

       do i = 1, N

         j = (i-1) * 3
         Gms(1,j+1) = ( FictMass(i,1) * Vnm(1,i,1) * Vnm(1,i,1) - kT ) &
         &          * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( FictMass(i,1) * Vnm(2,i,1) * Vnm(2,i,1) - kT ) &
         &          * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( FictMass(i,1) * Vnm(3,i,1) * Vnm(3,i,1) - kT ) &
         &          * InvMMNHC(1,j+3)

       end do

     else

       do i = 1, N

         j = (i-1) * 3
         Gms(1,j+1) = ( Mass(i) * Vel(1,i) * Vel(1,i) - kT ) * InvMMNHC(1,j+1)
         Gms(1,j+2) = ( Mass(i) * Vel(2,i) * Vel(2,i) - kT ) * InvMMNHC(1,j+2)
         Gms(1,j+3) = ( Mass(i) * Vel(3,i) * Vel(3,i) - kT ) * InvMMNHC(1,j+3)

       end do

     end if

     Num = 3 * N

     do i = 1, 3

       do j = i, 3

         if(i==j) then
           pref = 1.d0
         else
           pref = 2.d0
         end if

         Num = Num + 1
         Gms(1,Num) = ( pref * Mp(i,j) * Vg(i,j) * Vg(i,j) - kT ) &
         &          * InvMMNHC(1,Num)

       end do

     end do

     do i = 1 , NHchain-1

       do j = 1, NumMNHC

         Gms(i+1,j) = ( MMNHC(i,j) * VMNHC(i,j) * VMNHC(i,j) - kT ) &
         &          * InvMMNHC(i+1,j)

       end do

     end do

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         VMNHC(NHchain,:) = VMNHC(NHchain,:) + Gms(NHchain,:) * w4

         do k = NHchain-1 , 1 , -1

           do l = 1, NumMNHC

             cf = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * cf + Gms(k,l) * w4 ) * cf

           end do

         end do

! ## update dlog(v)/dt

         Num = 3 * N

         do k = 1 , 3

           do l = k, 3

             Num = Num + 1
             cf      = exp( -w8 * VMNHC(1,Num) )
             Vg(k,l) = ( Vg(k,l) * cf + Gmp(k,l) * w4 ) * cf

           end do

         end do

         Vg(2,1) = Vg(1,2)
         Vg(3,1) = Vg(1,3)
         Vg(3,2) = Vg(2,3)

! ## update the particle velocities
         TrVg    = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         RotVpre = Vg + TrVg * Imatrix
         RotV    = RotVpre

         l = 0

         if(QPathInt) then

           do k = 1, N

             RotV(1,1) = RotVpre(1,1) + VMNHC(1,l+1)
             RotV(2,2) = RotVpre(2,2) + VMNHC(1,l+2)
             RotV(3,3) = RotVpre(3,3) + VMNHC(1,l+3)
!
             l = l + 3
!           ----------------------------------------
             call Jacobi(RotV,eigenVec,eigenValue)
!           ----------------------------------------

             TeigenVec = Transpose( eigenVec )

             Scale = exp( -w2 * eigenValue )

             tempov   = matmul( TeigenVec, Vnm(:,k,1) )
             tempov   = tempov * Scale

             Vnm(:,k,1) = matmul( eigenVec, tempov )

           end do

         else

           do k = 1, N

             RotV(1,1) = RotVpre(1,1) + VMNHC(1,l+1)
             RotV(2,2) = RotVpre(2,2) + VMNHC(1,l+2)
             RotV(3,3) = RotVpre(3,3) + VMNHC(1,l+3)
!
             l = l + 3
!           ----------------------------------------
             call Jacobi(RotV,eigenVec,eigenValue)
!           ----------------------------------------

             TeigenVec = Transpose( eigenVec )

             Scale = exp( -w2 * eigenValue )

             tempov   = matmul( TeigenVec, Vel(:,k) )
             tempov   = tempov * Scale

             Vel(:,k) = matmul( eigenVec, tempov )

           end do

         end if

         tempoRot = matmul( TeigenVec, VelRotation )

         do k = 1 , 3

           tempoRot(:,k) = tempoRot(:,k) * Scale

         end do

         VelRotation = matmul( eigenVec, tempoRot )

!       ---------------
         if(QPathInt) then

           Sckin = 0.d0

           do ii = 1, N

             do k = 1, 3
               do l = j, 3

                 Sckin(k,l) = Sckin(k,l) + FictMass(ii,1) * Vnm(k,ii,1) * Vnm(l,ii,1)

               end do
             end do

           end do

           Sckin(2,1) = Sckin(1,2)
           Sckin(3,1) = Sckin(1,3)
           Sckin(3,2) = Sckin(2,3)

           Ekin = Sckin(1,1) + Sckin(2,2) + Sckin(3,3)

         else

           call CalcTemp
           Sckin = Pkinp
           Ekin = Ene_kin

         end if
!       ---------------

         diagTerm = Ekin * InvNf - Pext

         Gmp = (Sckin + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update the thermostat positions
         RMNHC = RMNHC + VMNHC * w2

! ## update dlog(v)/dt

         Num = 3 * N

         do k = 1, 3

           do l = k, 3

             Num = Num + 1
             cf  = exp( -w8 * VMNHC(1,Num) )
             Vg(k,l)  = (Vg(k,l) * cf + Gmp(k,l) * w4 ) * cf

           end do

         end do

         Vg(2,1) = Vg(1,2)
         Vg(3,1) = Vg(1,3)
         Vg(3,2) = Vg(2,3)

         Baro_kin = 0.d0

         do k = 1 , 3

           do l = 1 , 3

             Baro_kin = Baro_kin + Mp(k,l) * Vg(k,l) * Vg(k,l)

           end do

         end do

! ## update the forces
         if(QPathInt) then

           do k = 1, N

             l = (k-1) * 3
             Gms(1,l+1) = ( FictMass(k,1) * Vnm(1,k,1) * Vnm(1,k,1) - kT ) &
             &          * InvMMNHC(1,l+1)
             Gms(1,l+2) = ( FictMass(k,1) * Vnm(2,k,1) * Vnm(2,k,1) - kT ) &
             &          * InvMMNHC(1,l+2)
             Gms(1,l+3) = ( FictMass(K,1) * Vnm(3,k,1) * Vnm(3,k,1) - kT ) &
             &          * InvMMNHC(1,l+3)

           end do

         else

           do k = 1, N

             l = (k-1) * 3
             Gms(1,l+1) = ( Mass(k) * Vel(1,k) * Vel(1,k) - kT ) * InvMMNHC(1,l+1)
             Gms(1,l+2) = ( Mass(k) * Vel(2,k) * Vel(2,k) - kT ) * InvMMNHC(1,l+2)
             Gms(1,l+3) = ( Mass(K) * Vel(3,k) * Vel(3,k) - kT ) * InvMMNHC(1,l+3)

           end do

         end if

         Num = 3 * N

         do k = 1, 3

           do l = k, 3

             if(k==l) then
               pref = 1.d0
             else
               pref = 2.d0
             end if

             Num = Num + 1
             Gms(1,Num) = ( pref * Mp(k,l) * Vg(k,l) * Vg(k,l) - kT ) &
             &          * InvMMNHC(1,Num)

           end do

         end do

! ## update the thermostat velocities
         do k = 1 , NHchain-1

           do l = 1, NumMNHC

             cf         = exp( -w8 * VMNHC(k+1,l) )
             VMNHC(k,l) = ( VMNHC(k,l) * cf + w4 * Gms(k,l) ) * cf
             Gms(k+1,l) = ( MMNHC(k,l) * VMNHC(k,l) * VMNHC(k,l) - kT ) * InvMMNHC(k+1,l)

           end do

         end do

         VMNHC(NHchain,:)=VMNHC(NHchain,:)+Gms(NHchain,:)*w4

       end do

     end do

   end if

end subroutine BaroThermostatPR_MNHC


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** The Andersen barostat                                 **
! ** Nose-Hoover chain ( length = NHchain )                **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatAN_NH(Scale,ll)

use Numbers, only : NfT
use CommonBlocks, only : QRigidBody
use BathParam
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use CellParam, only : Volume
use ThermoData, only : Ene_kin, Ene_kinT, Virial

implicit NONE

integer :: i, j, k, ll
real(8) :: Scale, aa, bb, w2, w4, w8, odnf
real(8) :: Gms
real(8) :: Gmp
real(8) :: TraceVirial, P_o

   TraceVirial = Virial(1,1) + Virial(2,2) + Virial(3,3)
   P_o         = Pressure_o * 3.d0 * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Scale=1.d0

! ------------------------------
   call CalcTemp          ! total kinetic energy
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 / dble(NfT)
     Gmp = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

     Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(1) = Vss(1) + Gms * w4

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

! ## update the particle velocities
         aa      = exp( -w2 * Vss(1) )
         bb      = exp( -w2 * odnf * Vg(1,1) )
         V_RB    = V_RB * aa *  bb

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * aa
             Lmoment(3,k) = Lmoment(3,k) * aa
           else
             Lmoment(:,k) = Lmoment(:,k) * aa
           end if

         end do

!   ---------------------
         call CalcTemp
!   ---------------------
         Gmp     = ( odnf * Ene_kinT + TraceVirial - P_o ) * InvMp(1,1)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

         Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

! ## update the forces
         Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         Vss(1) = Vss(1) + Gms * w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

     odnf=1.d0 + 3.d0 * InvNf
     Gmp = ( odnf * Ene_kin + TraceVirial - P_o ) * InvMp(1,1)

     Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(1) = Vss(1) + Gms * w4

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

! ## update the particle velocities
         aa      = exp( -w2 * ( Vss(1) + odnf * Vg(1,1) ) )
         Scale   = Scale * aa
         Ene_kin = Ene_kin * aa * aa
         Gmp     = ( odnf * Ene_kin + TraceVirial - P_o ) * InvMp(1,1)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )
         Vg(1,1) = (Vg(1,1) * aa + Gmp * w4 ) * aa

         Baro_kin = Mp(1,1) * Vg(1,1) * Vg(1,1)

! ## update the forces
         Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         Vss(1) = Vss(1) + Gms * w4

       end do

     end do

! ---------------------------------
! update the particle velocities
! ---------------------------------
     Vel = Vel * Scale

   end if

end subroutine BaroThermostatAN_NH


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** Parrinello-Rahman barostat ( rectangular cell )       **
! ** Nose-Hoover chain ( length = NHchain )                **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatA3_NH(Vscale,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, cBarostatMethod
use BathParam
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use CellParam, only : Volume
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, ll, i1, i2
real(8) :: aa, w2, w4, w8, odnf
real(8) :: Gms
real(8), dimension(3) :: Gmp, cf
real(8) :: P_o
real(8), dimension(3) :: Vscale
real(8) :: TrVg, dterm
logical :: QIxy

   if(cBarostatMethod=='A2') then
     QIxy = .True.
     i1   = CoupleEdge(1)
     i2   = CoupleEdge(2)
   else
     QIxy = .False.
   end if

   P_o = Pressure_o * Volume

! ------------------------------------
! Nose-Hoover chain length [NHchain=5]
! ------------------------------------
   Vscale = 1.d0

! ------------------------------
   call CalcTemp          ! total kinetic energy
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = 0.d0
     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)
     end do

     odnf = 1.d0 / dble(NfT)

     do i = 1 , 3
       Gmp(i) = ( Pkinp(i,i) + odnf * Ene_kinT + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

     Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(1) = Vss(1) + Gms * w4

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         dterm = TrVg + Vss(1)
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( dterm + Vg(k,k) ) )
         end do

         do k = 1, NumRB
           V_RB(:,k) = V_RB(:,k) * cf
         end do

         aa = exp( -w2 * Vss(1) )

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then
             Lmoment(2,k) = Lmoment(2,k) * aa
             Lmoment(3,k) = Lmoment(3,k) * aa
           else
             Lmoment(:,k) = Lmoment(:,k) * aa
           end if

         end do

!       ---------------
         call CalcTemp
!       ---------------

         do k = 1, 3
           Gmp(k) = ( Pkinp(k,k) + odnf * Ene_kinT + Virial(k,k) - P_o ) * InvMp(k,k)
         end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         Baro_kin = 0.d0
         do k = 1 , 3
           Baro_kin = Baro_kin + Mp(k,k) * Vg(k,k) * Vg(k,k)
         end do

! ## update the forces
         Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         Vss(1) = Vss(1) + Gms * w4

       end do

     end do

   else

     Baro_kin = 0.d0
     do i = 1 , 3
       Baro_kin = Baro_kin + Mp(i,i) * Vg(i,i) * Vg(i,i)
     end do

     odnf = InvNf

     do i = 1 , 3
       Gmp(i) = ( Pkinp(i,i) + odnf * Ene_kin + Virial(i,i) - P_o ) * InvMp(i,i)
     end do

     if(QIxy) then
       Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
       Gmp(i2) = Gmp(i1)
     end if

     Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(1) = Vss(1) + Gms * w4

! ## update dlog(v)/dt

         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         TrVg  = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         dterm = TrVg + Vss(1)
! ## update the particle velocities

         do k = 1 , 3
           cf(k) = exp( -w2 * ( dterm + Vg(k,k) ) )
         end do

         Vscale   = Vscale * cf

         do k = 1, N
           Vel(:,k) = Vel(:,k) * cf(:)
         end do

!       ---------------
         call CalcTemp
!       ---------------

         do k = 1, 3
           Gmp(k) = ( Pkinp(k,k) + odnf * Ene_kin + Virial(k,k) - P_o ) * InvMp(k,k)
         end do

         if(QIxy) then
           Gmp(i1) = (Gmp(i1)+Gmp(i2)) * 0.5d0
           Gmp(i2) = Gmp(i1)
         end if

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         aa  = exp( -w8 * Vss(1) )

         do k = 1 , 3
           Vg(k,k) = (Vg(k,k) * aa + Gmp(k) * w4 ) * aa
         end do

         Baro_kin = 0.d0
         do k = 1 , 3
           Baro_kin = Baro_kin + Mp(k,k) * Vg(k,k) * Vg(k,k)
         end do

! ## update the forces
         Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         Vss(1) = Vss(1) + Gms * w4

       end do

     end do

   end if

   if(QIxy) Vg(i2,i2) = Vg(i1,i1)

end subroutine BaroThermostatA3_NH


!#####################################################################
!#####################################################################


! ***********************************************************
! ** integration of the equations of motion ( bath )       **
! ** time evolution of the barostat and thermostats        **
! ** using SHAKE/ROLL and RATTLE/ROLL procedure            **
! ** Parrinello-Rahman barostat ( fully flexible cell )    **
! ** Nose-Hoover chain ( length = NHchain )                **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine BaroThermostatPR_NH(VelRotation,ll)

use Numbers, only : N, NfT
use CommonBlocks, only : QRigidBody, cBarostatMethod
use BathParam
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear
use Configuration, only : Vel
use CellParam, only : H, Volume
use ThermoData, only : Ene_kin, Ene_kinT, Pkinp, Virial

implicit NONE

integer :: i, j, k, l, ll
real(8) :: cf, w2, w4, w8
real(8), dimension(3) :: Scale, tempov, eigenValue
real(8) :: Gms
real(8), dimension(3,3) :: Gmp
real(8), dimension(3,3) :: RotV, Imatrix
real(8), dimension(3,3) :: eigenVec, TeigenVec
real(8) :: TrVg, diagTerm, Pext
real(8), dimension(3,3) :: VelRotation, tempoRot
real(8), dimension(3,3) :: StressTerm, Htrans, tempM

   Imatrix = 0.d0

   do i = 1 , 3

     Imatrix(i,i) = 1.d0

   end do

   VelRotation = Imatrix

   Pext = Pressure_o * Volume

   if(cBarostatMethod == 'ST') then

     Htrans = transpose( H )

     tempM = matmul( H, SigmaS )

     StressTerm = - matmul( tempM, Htrans )

   else

     StressTerm = 0.d0

   end if

! ------------------------------
   call CalcTemp          ! total kinetic energy
! ------------------------------

! ## rigid-body ## -----------------------------------------------------

   if(QRigidBody) then

     Baro_kin = 0.d0

     do j = 1 , 3
       do i = 1 , 3
         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)
       end do
     end do

     diagTerm = Ene_kinT / dble(NfT) - Pext

     Gmp = (Pkinp + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

     Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(1) = Vss(1) + Gms * w4

! ## update dlog(v)/dt

         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

! ## update the particle velocities
         TrVg     = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) / dble(NfT)
         diagTerm = TrVg + Vss(1)
         RotV   = Vg + diagTerm * Imatrix
!       ----------------------------------------
         call Jacobi(RotV,eigenVec,eigenValue)
!       ----------------------------------------

         TeigenVec = Transpose( eigenVec )

         Scale = exp( -w2 * eigenValue )

         do k = 1 , NumRB

           tempov   = matmul( TeigenVec, V_RB(:,k) )
           tempov   = tempov * Scale

           V_RB(:,k) = matmul( eigenVec, tempov )

         end do

         cf = exp( -w2 * Vss(1) )

         do k = 1 , NumRB

           if(QSingle(k)) cycle

           if(QLinear(k)) then

             Lmoment(2,k) = Lmoment(2,k) * cf
             Lmoment(3,k) = Lmoment(3,k) * cf

           else

             Lmoment(:,k) = Lmoment(:,k) * cf

           end if

         end do

!       ---------------
         call CalcTemp
!       ---------------

         diagTerm = Ene_kinT / dble(NfT) - Pext

         Gmp = (Pkinp + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

         Baro_kin = 0.d0

         do k = 1 , 3
           do l = 1 , 3
             Baro_kin = Baro_kin + Mp(k,l) * Vg(k,l) * Vg(k,l)
           end do
         end do

! ## update the forces
         Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         Vss(1) = Vss(1) + Gms * w4

       end do

     end do

! ## flexible ## -------------------------------------------------------

   else

     Baro_kin = 0.d0

     do j = 1 , 3
       do i = 1 , 3
         Baro_kin = Baro_kin + Mp(i,j) * Vg(i,j) * Vg(i,j)  ! Wg*Tr(Vg^t*Vg)
       end do
     end do

     diagTerm = Ene_kin * InvNf - Pext

     Gmp = (Pkinp + Virial + StressTerm + diagTerm*Imatrix ) * InvMp

     Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
     Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
     Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
     Gmp(2,1) = Gmp(1,2)
     Gmp(3,1) = Gmp(1,3)
     Gmp(3,2) = Gmp(2,3)

     Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ---------------------------------------
! start the multiple time step procedure
! ---------------------------------------

     do i = 1 , Nsc

       do j = 1 , NYoshid

         w2 = wdti2(j) * dble(ll)
         w4 = wdti4(j) * dble(ll)
         w8 = wdti8(j) * dble(ll)

! ## update the themostat velocities coupling to particles

         Vss(1) = Vss(1) + Gms * w4

! ## update dlog(v)/dt

         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

! ## update the particle velocities
         TrVg     = ( Vg(1,1) + Vg(2,2) + Vg(3,3) ) * InvNf
         diagTerm = TrVg + Vss(1)
         RotV   = Vg + diagTerm * Imatrix
!       ----------------------------------------
         call Jacobi(RotV,eigenVec,eigenValue)
!       ----------------------------------------

         TeigenVec = Transpose( eigenVec )

         Scale = exp( -w2 * eigenValue )

         do k = 1 , N

           tempov   = matmul( TeigenVec, Vel(:,k) )
           tempov   = tempov * Scale

           Vel(:,k) = matmul( eigenVec, tempov )

         end do

         tempoRot = matmul( TeigenVec, VelRotation )

         do k = 1 , 3

           tempoRot(:,k) = tempoRot(:,k) * Scale

         end do

         VelRotation = matmul( eigenVec, tempoRot )

!       ---------------
         call CalcTemp
!       ---------------

         diagTerm = Ene_kin * InvNf - Pext

         Gmp = (Pkinp + Virial + StressTerm + diagTerm * Imatrix ) * InvMp

         Gmp(1,2) = ( Gmp(1,2) + Gmp(2,1) ) * 0.5d0
         Gmp(1,3) = ( Gmp(1,3) + Gmp(3,1) ) * 0.5d0
         Gmp(2,3) = ( Gmp(2,3) + Gmp(3,2) ) * 0.5d0
         Gmp(2,1) = Gmp(1,2)
         Gmp(3,1) = Gmp(1,3)
         Gmp(3,2) = Gmp(2,3)

! ## update the thermostat positions
         Rss = Rss + Vss * w2

! ## update dlog(v)/dt
         cf  = exp( -w8 * Vss(1) )
         Vg  = (Vg * cf + Gmp * w4 ) * cf

         Baro_kin = 0.d0

         do k = 1 , 3

           do l = 1 , 3

             Baro_kin = Baro_kin + Mp(k,l) * Vg(k,l) * Vg(k,l)

           end do

         end do

! ## update the forces
         Gms = ( Baro_kin + Ene_kin - gkT ) * InvMts(1)

! ## update the thermostat velocities
         Vss(1) = Vss(1) + Gms * w4

       end do

     end do

   end if

end subroutine BaroThermostatPR_NH
