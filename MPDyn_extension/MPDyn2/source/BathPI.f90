! ############################
! ## SUBROUTINE LIST 
! ## -- ThermostatPIcent
! ## -- ThermostatPIcent_NHC
! ## -- ThermostatPIcent_MNHC
! ## -- ThermostatPInonc
! ## -- BaroThermostatPIcent
! ############################


!######################################################################
!######################################################################


subroutine ThermostatPIcent(ll)

use CommonBlocks,only : cThermostatMethod

implicit none

integer :: ll

   if( cThermostatMethod == 'NHC' ) then

     call ThermostatPIcent_NHC(ll)

   else if( cThermostatMethod == 'MNHC' ) then

     call ThermostatPIcent_MNHC(ll)

   end if

end subroutine ThermostatPIcent


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

subroutine ThermostatPIcent_NHC(ll)

use Numbers, only : N
use BathParam
use CommonPI

implicit none

integer :: i, j, k, ll
real(8) :: Scale, aa, w2, w4, w8, Sckin
real(8), dimension(NHchain) :: Gms

! ------------------------------------
! Nose-Hoover chain length [NHchain]
! ------------------------------------
   Scale=1.d0

! ----------------
   Sckin = 0.d0

   do i = 1, N

     Sckin = Sckin + FictMass(i,1) *             &
     &       dot_product( Vnm(:,i,1), Vnm(:,i,1) )

   end do

! ----------------

   Gms(1) = ( Sckin - gkT ) * InvMts(1)

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
       aa    = exp( -w2 * Vss(1) )
       Scale = Scale * aa
       Sckin = Sckin * aa * aa

! --------------------------------
! update the thermostat positions
! --------------------------------
       Rss = Rss + Vss * w2

! ------------------
! update the forces
! ------------------
       Gms(1) = ( Sckin - gkT ) * InvMts(1)

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

   do i = 1 , N
     Vnm(:,i,1) = Vnm(:,i,1) * Scale
   end do

end subroutine ThermostatPIcent_NHC


!######################################################################
!######################################################################


! ***********************************************************
! ** integration of the equations of motion of thermostat  **
! ** Massive Nose-Hoover chain ( length = NHchain )        **
! **                                                       **
! ** integrated by                                         **
! **  <reversible REfference System Propagator Algorithm>  **
! ***********************************************************

subroutine ThermostatPIcent_MNHC(ll)

use Numbers, only : N
use BathParam
use CommonPI

implicit none

integer :: i, j, k, ll, l
real(8) :: aa, w2, w4, w8
real(8), dimension(NHchain,NumMNHC) :: Gms

   do i = 1, N

     j = (i-1) * 3

     Gms(1,j+1) = ( FictMass(i,1) * Vnm(1,i,1) * Vnm(1,i,1) - kT ) * InvMMNHC(1,j+1)
     Gms(1,j+2) = ( FictMass(i,1) * Vnm(2,i,1) * Vnm(2,i,1) - kT ) * InvMMNHC(1,j+2)
     Gms(1,j+3) = ( FictMass(i,1) * Vnm(3,i,1) * Vnm(3,i,1) - kT ) * InvMMNHC(1,j+3)

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
       do k = 1, N

         l = 3 * (k-1)

         Vnm(1,k,1) = Vnm(1,k,1) * exp( -w2 * VMNHC(1,l+1) )
         Vnm(2,k,1) = Vnm(2,k,1) * exp( -w2 * VMNHC(1,l+2) )
         Vnm(3,k,1) = Vnm(3,k,1) * exp( -w2 * VMNHC(1,l+3) )

       end do

! --------------------------------
! update the thermostat positions
! --------------------------------
       RMNHC = RMNHC + VMNHC * w2

! ------------------
! update the forces
! ------------------
       do k = 1, N

         l = (k-1) * 3

         Gms(1,l+1) = ( FictMass(k,1) * Vnm(1,k,1) * Vnm(1,k,1) - kT ) * InvMMNHC(1,l+1)
         Gms(1,l+2) = ( FictMass(k,1) * Vnm(2,k,1) * Vnm(2,k,1) - kT ) * InvMMNHC(1,l+2)
         Gms(1,l+3) = ( FictMass(K,1) * Vnm(3,k,1) * Vnm(3,k,1) - kT ) * InvMMNHC(1,l+3)

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

end subroutine ThermostatPIcent_MNHC


!######################################################################
!######################################################################


subroutine ThermostatPInonc

use Numbers, only : N
use BathParam
use CommonPI

implicit none

integer :: iatom, i, j, mode
real(8) :: w2, w4, w8
real(8), dimension(3) :: Dkin, Scale, pfact, vfact
real(8) :: pdt

   pdt = 1.d0 / dble(Ncref)

   do j = 1 , NYoshid

     w2 = wdti2(j) * pdt
     w4 = wdti4(j) * pdt
     w8 = wdti8(j) * pdt

! ## update the themostat velocities coupling to particles

     do mode = IniBead, FinBead

       if(mode == 1) cycle

       do iatom = 1, N

! ## calculate kinetic energy of normal mode coordinate for each atom
         Dkin(:) = FictMass(iatom,mode) * Vnm(:,iatom,mode) * Vnm(:,iatom,mode)
         Scale(:) = 1.d0

         Fbath(:,iatom,1,mode) = (Dkin(:) - kT) * InvQmass(mode)

         do i = 1, NHchain - 1
           Fbath(:,iatom,i+1,mode) = &
           &    (Qmass(mode)*Vbath(:,iatom,i,mode)*Vbath(:,iatom,i,mode) &
           &     - kT) * InvQmass(mode)
         end do

!     /* update the thermostat velocities */
         Vbath(:,iatom,NHchain,mode) = Vbath(:,iatom,NHchain,mode) &
         &                           + Fbath(:,iatom,NHchain,mode) * w4

         do i = NHchain-1, 1, -1

           vfact(:) = exp(- Vbath(:,iatom,i+1,mode) * w8 )

           Vbath(:,iatom,i,mode) = ( Vbath(:,iatom,i,mode) * vfact(:)  &
           &                       + Fbath(:,iatom,i,mode) * w4 ) * vfact(:)

         end do

! --------------------------------
! update the particle velocities
! --------------------------------
         pfact(:) = exp( - w2 * Vbath(:,iatom,1,mode) )
         Scale(:) = Scale(:) * pfact(:)

!     /* update the force */
         Fbath(:,iatom,1,mode) = ( Scale(:) * Scale(:) * Dkin(:) - kT) * InvQmass(mode)

!     /* update the thermostat position */
         do i = 1, NHchain
           Rbath(:,iatom,i,mode) = Rbath(:,iatom,i,mode) &
           &                     + Vbath(:,iatom,i,mode) * w2
         end do

!     /* update the tehrmostat velocities */
         do i = 1, NHchain-1

           vfact(:) = exp( - w8 * Vbath(:,iatom,i+1,mode) )

           Vbath(:,iatom,i,mode) = ( Vbath(:,iatom,i,mode) * vfact(:)  &
           &                       + Fbath(:,iatom,i,mode) * w4 ) * vfact(:)

           Fbath(:,iatom,i+1,mode) = &
           &    (Qmass(mode) * Vbath(:,iatom,i,mode) * Vbath(:,iatom,i,mode) &
           &    - kT) * InvQmass(mode)

         end do

         Vbath(:,iatom,NHchain,mode) = Vbath(:,iatom,NHchain,mode) &
         &                           + Fbath(:,iatom,NHchain,mode) * w4

!     /* update the paricle velocities */
         Vnm(:,iatom,mode) = Vnm(:,iatom,mode) * Scale(:)

       end do

     end do

   end do

end subroutine ThermostatPInonc


!######################################################################
!######################################################################


subroutine BaroThermostatPIcent(ll)

use CommonBlocks, only : cThermostatMethod, cBarostatMethod

implicit none

integer :: ll
real(8) :: Scale
real(8), dimension(3) :: VScale
real(8), dimension(3,3) :: VelRotation


   if( cThermostatMethod == 'NHC' ) then

     if( ( cBarostatMethod == 'PR' ).or.( cBarostatMethod == 'ST') ) then
       call BaroThermostatPR(VelRotation,ll)
     else if(( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
       call BaroThermostatA3(VScale,ll)
     else if( cBarostatMethod == 'AN' ) then
       call BaroThermostatAN(Scale,ll)
     end if

   else if( cThermostatMethod == 'MNHC' ) then

     if( ( cBarostatMethod == 'PR' ).or.(cBarostatMethod == 'ST') ) then
       call BaroThermostatPR_MNHC(VelRotation,ll)
     else if(( cBarostatMethod == 'A3' ).or.( cBarostatMethod == 'A2' ) ) then
       call BaroThermostatA3_MNHC(VScale,ll)
     else if( cBarostatMethod == 'AN' ) then
       call BaroThermostatAN_MNHC(Scale,ll)
     end if

   end if

end subroutine BaroThermostatPIcent
