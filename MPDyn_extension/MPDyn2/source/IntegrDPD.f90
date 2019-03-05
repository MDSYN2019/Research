! ############################
! ## SUBROUTINE LIST 
! ## -- IntegrDPD_VVerlet 
! ## -- IntegrDPD_VVerletSC 
! ## -- IntegrDPD_LeapFrogSC 
! ## -- IntegrDPD_Lowe 
! ############################


!######################################################################
!######################################################################


subroutine IntegrDPD_VVerlet

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration
use CommonDPD
use RBparam
use CellParam, only : CellL
use TimeParam, only : Nstep, deltat, dt2, Timeps

implicit none

!
! #
! # time evolution of particles
! # reversible RESPA (Velocity-Verlet) like scheme
! #
! # The algorithm used here is based upon Trotter expansion of
! # Liouville equation in order to obtain stable dynamics with
! # larger time step.
! #
! # Ref: R.D.Groot and P.B.Warren, JCP 107 4423(1997).
! #      J.M.V.A.Koelman and P.J.Hoogerbrugge, Europhys.Lett.
! #                                             21 363(1993).
! #      N.Matubayasi and M.Nakahara, JCP 110 3291(1999).
! #      ( time reversible scheme of integration for quarternion )
! #
! # NOTE!
! # If lamda is not equal to 1/2, the algorithm for rotational motion
! # has to be changed.
! #
!

real(8) :: rizy, rizx
real(8) :: Lamdt, deltath, ShearD
integer :: ll, i, j, k, l
real(8) :: Lmomentztemp, th1, th2
real(8), dimension(3) :: Lt
real(8), dimension(3) :: Omt
real(8) :: cst, snt
real(8) :: omg, omg2
real(8) :: th, diag, offd
real(8), dimension(4) :: qq

!************************************
!*  main part of the DPD simulation *
!************************************

   allocate( Velt(3,N) )
   allocate( dR(3,N)   )
   allocate( SdR(3,N)  )

   dt2     = 0.5d0  * deltat * deltat
   Lamdt   = Lambda * deltat
   deltath = 0.5d0  * deltat

   if(ShearRate == 0.) then

     QSheared = .False.

   else

     QSheared = .True.
     ShearD  = ShearRate * deltat * CellL(2)

   end if

   SdR = 0.d0

   if(QColloid) then

     allocate( SdRc(3,NumColloid) )
     allocate(  dRc(3,NumColloid) )

     SdRc = 0.d0

   end if

   l = 0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       rizy = InvInertiaRB(3,1) - InvInertiaRB(2,1)
       rizx = InvInertiaRB(3,1) - InvInertiaRB(1,1)

       do j = 1 , NumColloid
         Omega(:,j) = Lmoment(:,j) * InvInertiaRB(:,1)
       end do

       V_RBt  = V_RB
       Omegat = Omega

       l = l + NumColloid * NumCollAtm

     else

       do j = 1 , NumMol(i)

         do k = 1 , NumAtm(i)

           l = l + 1
           Velt(:,l) = Vel(:,l)

         end do

       end do

     end if

   end do

!
! space-fixed coordinates and velocities of the colloid components
! ---------------------------------
   if(QColloid) then

     call RotationMatrix
     call SiteInfo_Colloid

   end if
! ---------------------------------

!
! Periodic Boundary Condition
! ---------------------------------
   call PBC_DPD
! ---------------------------------

! ## Force
! ----------------------------------------------------------
   call GetForceVV

   FrcDP = FrcDPt

   if(QColloid) then

     call ForceTorque_COM_Colloid

     FrcCo = FrcCOt

   end if
! ----------------------------------------------------------
!
! #########################################################
!                      MAIN ROUTINE
! #########################################################
!
   do ll = 1 , Nstep

     timeps = timeps + deltat
!
! #### v~(t+dt)=v(t)+(lamda*dt)*F(t)/m
!
! ######################
     if(QColloid) then
! ######################

       l = 0

       do i = 1 , NumSpec

         if(ColloidFlag(i)) then

           do j = 1 , NumColloid

             V_RBt(:,j) = V_RB(:,j) + Lamdt * FrcCo(:,j) * InvMassRB(1)

           end do

! ##### L~(t+lamda*dt)=L(t)+(lamda*dt)*N(t)

           Lmoment  = Lmoment + deltath * Torque
           Lmomentt = Lmoment + Lamdt   * Torque

           do j = 1 , NumColloid

             Lt = Lmoment(:,j)

             th1 = deltath * rizy * Lt(2)
             cst = cos(th1)
             snt = sin(th1)
             Lmoment(1,j) =   cst * Lt(1) + snt * Lt(3)
             Lmomentztemp = - snt * Lt(1) + cst * Lt(3)

             th2 = deltath * rizx * Lmoment(1,j)
             cst = cos(th2)
             snt = sin(th2)
             Lmoment(2,j) =   cst * Lt(2) - snt * Lmomentztemp
             Lmoment(3,j) =   snt * Lt(2) + cst * Lmomentztemp

           end do

           do j = 1 , NumColloid

             Lt = Lmomentt(:,j)

             th1 = Lamdt * rizy * Lt(2)
             cst = cos(th1)
             snt = sin(th1)
             Lmomentt(1,j) =   cst * Lt(1) + snt * Lt(3)
             Lmomentztemp  = - snt * Lt(1) + cst * Lt(3)

             th2 = Lamdt * rizx * Lmomentt(1,j)
             cst = cos(th2)
             snt = sin(th2)
             Lmomentt(2,j) =   cst * Lt(2) - snt * Lmomentztemp
             Lmomentt(3,j) =   snt * Lt(2) + cst * Lmomentztemp

           end do

! ##### Omg~(t+lamda*dt) = I^-1 L~(t+lamda*dt)

           do j = 1 , NumColloid

             Omega(:,j)  = Lmoment(:,j)  * InvInertiaRB(:,1)
             Omegat(:,j) = Lmomentt(:,j) * InvInertiaRB(:,1)

           end do

! ##### Q(t+dt)=exp(dt/2*A[Omg~(t+lamda*dt)]) Q(t)

           do j = 1 , NumColloid

             qq = Quaternion(:,j)

             Omt = Omega(:,j)

             omg2 = dot_product( Omt, Omt )
             omg  = sqrt( omg2 )
             th   = deltath * omg
             diag = dcos(th)
             offd = dsin(th) / omg

             Quaternion(1,j) = diag * qq(1) + offd * ( - Omt(1) * qq(2) &
             &               - Omt(2) * qq(3) - Omt(3) * qq(4) )
             Quaternion(2,j) = diag * qq(2) + offd * (   Omt(1) * qq(1) &
             &               - Omt(2) * qq(4) + Omt(3) * qq(3) )
             Quaternion(3,j) = diag * qq(3) + offd * (   Omt(1) * qq(4) &
             &               + Omt(2) * qq(1) - Omt(3) * qq(2) )
             Quaternion(4,j) = diag * qq(4) + offd * ( - Omt(1) * qq(3) &
             &               + Omt(2) * qq(2) + Omt(3) * qq(1) )

           end do
!
! ##### A(t+dt)  from  Q(t+dt)
! ----------------------------------------------------------
           call RotationMatrix
! ----------------------------------------------------------

           do j = 1 , NumColloid

             dRc(:,j) = deltat * V_RB(:,j) + dt2 * FrcCo(:,j) * InvMassRB(1)

           end do

           R_RB = R_RB + dRc


           l = l + NumColloid * NumCollAtm

         else

           do j = 1 , NumMol(i)

             do k = 1 , NumAtm(i)

               l = l + 1
               Velt(:,l) = Vel(:,l) + Lamdt * FrcDP(:,l)
               dR  (:,l) = deltat * Vel(:,l) + dt2 * FrcDP(:,l)
               R   (:,l) = R(:,l) + dR(:,l)

             end do

           end do

         end if

       end do

! ######################
     else
! ######################

       Velt = Vel + FrcDP * Lamdt

       dR = Vel * deltat + FrcDP * dt2

       R  = R + dR

! ######################
     end if
! ######################

! Lees-Edwards

     if(QSheared) then

       SlideGap = SlideGap + ShearD
       if( SlideGap > CellL(1) ) SlideGap = SlideGap - CellL(1)

     end if

! ##
! ## NOTE!  Cell List Method
! ## The SiteInfo_Colloid subroutine must be called before the PBC one.
! ##
! space-fixed coordinates and velocities of the colloid components
! --------------------------------------
     if(QColloid) call SiteInfo_Colloid
! --------------------------------------

! Periodic Boundary Condition (Sheared)
! ---------------------------------
     call PBC_DPD
! ---------------------------------

! ------------------------------------------------------------------
     call GetForceVV

     if(QColloid) call ForceTorque_COM_Colloid
! ------------------------------------------------------------------
!
! ##### Free Rotation on the Lx axis, succeedingly on the Ly axis
!
! ######################
     if(QColloid) then
! ######################

       l = 0
       do i = 1 , NumSpec

         if(ColloidFlag(i)) then

           do j = 1 , NumColloid

             Lt = Lmoment(:,j)

             th2 = deltath * rizx * Lt(1)
             cst = cos(th2)
             snt = sin(th2)
             Lmoment(2,j) =   cst * Lt(2) - snt * Lt(3)
             Lmomentztemp =   snt * Lt(2) + cst * Lt(3)

             th1 = deltath * rizy * Lmoment(2,j)
             cst = cos(th1)
             snt = sin(th1)
             Lmoment(1,j) =   cst * Lt(1) + snt * Lmomentztemp
             Lmoment(3,j) = - snt * Lt(1) + cst * Lmomentztemp

           end do

! ##### L(t+dt)=L~(t+dt)+(dt/2)*N(t+dt)

           Lmoment = Lmoment + deltath * Torque

! #### v(t+dt)=v(t)+(dt/2)*(F(t)+F(t+dt))/m

           do j = 1 , NumColloid

             V_RB(:,j) = V_RB(:,j) + deltath * (FrcCo(:,j) + FrcCot(:,j)) * InvMassRB(1)

           end do

           l = l + NumColloid * NumCollAtm

         else

           do j = 1 , NumMol(i)

             do k = 1 , NumAtm(i)

               l = l + 1
               Vel(:,l) = Vel(:,l) + deltath * (FrcDP(:,l) + FrcDPt(:,l))

             end do

           end do

         end if

       end do

! ######################
     else
! ######################

       Vel = Vel + deltath * (FrcDP + FrcDPt)

! ######################
     end if
! ######################

     FrcDP  = FrcDPt

     if(QColloid) FrcCo = FrcCot

!     call Elim_CellMoveDP

     call DPD_monitor(ll)

   end do

end subroutine IntegrDPD_VVerlet


!######################################################################
!######################################################################


subroutine IntegrDPD_VVerletSC

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration
use CommonDPD
use RBparam
use BookParam, only : Npair, ListIJ
use CellParam, only : CellL
use TimeParam, only : Nstep, deltat, dt2, Timeps
use ThermoData, only : Virial

implicit none

!
! # 
! # Self-consistent velocity Verlet 
! # I. Vattulainen et al., J. Chem. Phys. 116 3967 (2002).
! # 
!

real(8) :: rizy, rizx
real(8) :: deltath, ShearD
integer :: ll, i, j, k, l, ii
real(8) :: Lmomentztemp, th1, th2
real(8), dimension(3) :: Lt
real(8), dimension(3) :: Omt, Vij
real(8) :: cst, snt
real(8) :: omg, omg2
real(8) :: th, diag, offd
real(8), dimension(4) :: qq
real(8) :: RV
real(8), dimension(3) :: Fij

!************************************
!*  main part of the DPD simulation *
!************************************

   allocate( Velt(3,N) )
   allocate( dR(3,N) )
   allocate( SdR(3,N) )

   allocate( VelP(3,N) )

   dt2     = 0.5d0  * deltat * deltat
   deltath = 0.5d0  * deltat

! -----------------------------
   if(ShearRate == 0.) then
     QSheared = .False.
   else
     QSheared = .True.
     ShearD  = ShearRate * deltat * CellL(2)
   end if

   SdR = 0.d0

   if(QColloid) then

     allocate( SdRc(3,NumColloid) )
     allocate(  dRc(3,NumColloid) )
     SdRc = 0.d0

   end if

!
! space-fixed coordinates and velocities of the colloid components
! ---------------------------------
   if(QColloid) then
     rizy = InvInertiaRB(3,1) - InvInertiaRB(2,1)
     rizx = InvInertiaRB(3,1) - InvInertiaRB(1,1)

     do j = 1 , NumColloid
       Omega(:,j) = Lmoment(:,j) * InvInertiaRB(:,1)
     end do

     call RotationMatrix
     call SiteInfo_ColloidSC
   end if
! ---------------------------------

!
! Periodic Boundary Condition (Sheared)
! ---------------------------------
   call PBC_DPD
! ---------------------------------

! ## Force
! ----------------------------------------------------------
   call GetForceVVSC

! ## Dissipatvie Force Calculation

   FrcDPd = 0.d0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)

     Vij    = Vel(:,i) - Vel(:,j)
     Vij(1) = Vij(1) + SLList(l)
     RV     = dot_product( Vij, dRList(:,l) )
     Fij    = pfList(l) * RV * dRList(:,l)

     FrcDPd(:,i) = FrcDPd(:,i) + Fij
     FrcDPd(:,j) = FrcDPd(:,j) - Fij

   end do

! ----------------------------------------------------------

   if(QColloid) then

     call ForceTorque_COM_Colloid
     call ForceTorque_COM_Dissipative

   end if

! ----------------------------------------------------------

!
! #########################################################
!                      MAIN ROUTINE
! #########################################################
!
   do ll = 1 , Nstep

     timeps = timeps + deltat
!
! #### v~(t+0.5*dt)=v(t)+(0.5*dt)*F(t)/m
!
     if(QColloid) then

       l = 0

       do i = 1 , NumSpec

         if(ColloidFlag(i)) then

           do j = 1 , NumColloid

             V_RB(:,j) = V_RB(:,j) + deltath * ( FrcCot(:,j) + FrcCod(:,j) ) * InvMassRB(1)

           end do

! ##### L~(t+0.5*dt)=L(t)+(0.5*dt)*N(t)

           Lmoment  = Lmoment + deltath * ( Torque + Torqued )

           do j = 1 , NumColloid

             Lt = Lmoment(:,j)

             th1 = deltath * rizy * Lt(2)
             cst = cos(th1)
             snt = sin(th1)
             Lmoment(1,j) =   cst * Lt(1) + snt * Lt(3)
             Lmomentztemp = - snt * Lt(1) + cst * Lt(3)

             th2 = deltath * rizx * Lmoment(1,j)
             cst = cos(th2)
             snt = sin(th2)
             Lmoment(2,j) =   cst * Lt(2) - snt * Lmomentztemp
             Lmoment(3,j) =   snt * Lt(2) + cst * Lmomentztemp

           end do

! ##### Omg~(t+0.5*dt) = I^-1 L~(t+0.5*dt)

           do j = 1 , NumColloid

             Omega(:,j)  = Lmoment(:,j)  * InvInertiaRB(:,1)

           end do

! ##### Q(t+dt)=exp(dt/2*A[Omg~(t+lamda*dt)]) Q(t)

           do j = 1 , NumColloid

             qq = Quaternion(:,j)

             Omt = Omega(:,j)

             omg2 = dot_product( Omt, Omt )
             omg  = sqrt( omg2 )
             th   = deltath * omg
             diag = dcos(th)
             offd = dsin(th) / omg

             Quaternion(1,j) = diag * qq(1) + offd * ( - Omt(1) * qq(2) &
             &               - Omt(2) * qq(3) - Omt(3) * qq(4) )
             Quaternion(2,j) = diag * qq(2) + offd * (   Omt(1) * qq(1) &
             &               - Omt(2) * qq(4) + Omt(3) * qq(3) )
             Quaternion(3,j) = diag * qq(3) + offd * (   Omt(1) * qq(4) &
             &               + Omt(2) * qq(1) - Omt(3) * qq(2) )
             Quaternion(4,j) = diag * qq(4) + offd * ( - Omt(1) * qq(3) &
             &               + Omt(2) * qq(2) + Omt(3) * qq(1) )

           end do
!
! ##### A(t+dt)  from  Q(t+dt)
! ----------------------------------------------------------
           call RotationMatrix
! ----------------------------------------------------------

           do j = 1 , NumColloid

             dRc(:,j) = deltat * V_RB(:,j)

           end do

           R_RB = R_RB + dRc


           l = l + NumColloid * NumCollAtm

         else

           do j = 1 , NumMol(i)

             do k = 1 , NumAtm(i)

               l = l + 1
               Vel(:,l) = Vel(:,l) + deltath * ( FrcDPt(:,l) + FrcDPd(:,l) )
               dR (:,l) = deltat * Vel(:,l)
               R  (:,l) = R(:,l) + dR(:,l)

             end do

           end do

         end if

       end do

     else

       Vel = Vel + deltath * ( FrcDPt + FrcDPd )

       dR = deltat * Vel
       R  = R + dR

     end if

! Lees-Edwards

     if(QSheared) then

       SlideGap = SlideGap + ShearD
       if( SlideGap > CellL(1) ) SlideGap = SlideGap - CellL(1)

     end if

! ##
! ## NOTE!  Cell List Method
! ## The SiteInfo_Colloid subroutine must be called before the PBC one.
! ##
! space-fixed coordinates and velocities of the colloid components
! --------------------------------------
     if(QColloid) call SiteInfo_Colloid
! --------------------------------------

! Periodic Boundary Condition (Sheared)
! ---------------------------------
     call PBC_DPD
! ---------------------------------

! ------------------------------------------------------------------
     call GetForceVVSC

     if(QColloid) call ForceTorque_COM_Colloid
! ------------------------------------------------------------------

     Velt = Vel

     if(QColloid) then
       V_RBt    = V_RB
       Lmomentt = Lmoment
     end if

! ## Start Self-Consistent Force Calculation

     do ii = 1 , Iterate

! ## Dissipatvie Force Calculation

       if(QColloid) call SiteVeloc_Colloid

       FrcDPd = 0.d0

       if(ii /= Iterate) then

         do l = 1 , Npair

           i = ListIJ(1,l)
           j = ListIJ(2,l)

           Vij    = Vel(:,i) - Vel(:,j)
           Vij(1) = Vij(1) + SLList(l)
           RV     = dot_product( Vij, dRList(:,l) )
           Fij    = pfList(l) * RV * dRList(:,l)

           FrcDPd(:,i) = FrcDPd(:,i) + Fij
           FrcDPd(:,j) = FrcDPd(:,j) - Fij

         end do

       else

         VelP = Vel

         VirialD = 0.d0

         do l = 1 , Npair

           i = ListIJ(1,l)
           j = ListIJ(2,l)

           Vij    = Vel(:,i) - Vel(:,j)
           Vij(1) = Vij(1) + SLList(l)
           RV     = dot_product( Vij, dRList(:,l) )
           Fij    = pfList(l) * RV * dRList(:,l)

           FrcDPd(:,i) = FrcDPd(:,i) + Fij
           FrcDPd(:,j) = FrcDPd(:,j) - Fij

           VirialD(1,1) = VirialD(1,1) + Fij(1) * dRList(1,l) * R1List(l)
           VirialD(1,2) = VirialD(1,2) + Fij(1) * dRList(2,l) * R1List(l)
           VirialD(1,3) = VirialD(1,3) + Fij(1) * dRList(3,l) * R1List(l)
           VirialD(2,2) = VirialD(2,2) + Fij(2) * dRList(2,l) * R1List(l)
           VirialD(2,3) = VirialD(2,3) + Fij(2) * dRList(3,l) * R1List(l)
           VirialD(3,3) = VirialD(3,3) + Fij(3) * dRList(3,l) * R1List(l)

         end do

         VirialD(2,1) = VirialD(1,2)
         VirialD(3,1) = VirialD(1,3)
         VirialD(3,2) = VirialD(2,3)

         if(QSheared) Virial = Virial + VirialD

       end if

!
! ##### Free Rotation on the Lx axis, succeedingly on the Ly axis
!
       if(QColloid) then

         call ForceTorque_COM_Dissipative

         l = 0
         do i = 1 , NumSpec

           if(ColloidFlag(i)) then

             do j = 1 , NumColloid

               Lt = Lmomentt(:,j)

               th2 = deltath * rizx * Lt(1)
               cst = cos(th2)
               snt = sin(th2)
               Lmoment(2,j) =   cst * Lt(2) - snt * Lt(3)
               Lmomentztemp  =   snt * Lt(2) + cst * Lt(3)

               th1 = deltath * rizy * Lmoment(2,j)
               cst = cos(th1)
               snt = sin(th1)
               Lmoment(1,j) =   cst * Lt(1) + snt * Lmomentztemp
               Lmoment(3,j) = - snt * Lt(1) + cst * Lmomentztemp

             end do

! ##### L(t+dt)=L~(t+dt)+(dt/2)*N(t+dt)

             Lmoment = Lmoment + deltath * ( Torque + Torqued )

! #### v(t+dt)=v(t)+(dt/2)*(F(t)+F(t+dt))/m

             do j = 1 , NumColloid

               V_RB(:,j) = V_RBt(:,j) + deltath * ( FrcCot(:,j) + FrcDPd(:,j) ) * InvMassRB(1)

             end do

             l = l + NumColloid * NumCollAtm

           else

             do j = 1 , NumMol(i)

               do k = 1 , NumAtm(i)

                 l = l + 1
                 Vel(:,l) = Velt(:,l) + deltath * ( FrcDPt(:,l) + FrcDPd(:,l) )

               end do

             end do

           end if

         end do

       else

         Vel = Velt + deltath * ( FrcDPt + FrcDPd )

       end if

     end do

     if(QColloid) call SiteVeloc_Colloid

! -----------------------------------------------

     call Elim_CellMoveDP

     call DPD_monitor(ll)

   end do

end subroutine IntegrDPD_VVerletSC


!######################################################################
!######################################################################


subroutine IntegrDPD_LeapFrogSC

use Configuration
use CommonDPD

implicit none

   write(*,*) 'this scheme have not been implemented yet!'
   call Finalize

end subroutine IntegrDPD_LeapFrogSC


!######################################################################
!######################################################################


subroutine IntegrDPD_Lowe

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration
use CommonDPD
use RBparam
use BookParam, only : Npair, ListIJ
use CellParam, only : CellL
use TimeParam, only : Nstep, deltat, dt2, Timeps
use ThermoData, only : Virial

implicit none

!
! #
! # 1. Lowe-Andersen thermostat:
! #    C. P. Lowe, Europhys. Lett. 47 145 (1999).
! #
! # 2. Peters thermostat:
! #    E. A. J. F. Peters, Europhys. Lett. 66 311 (2004).
! #
!

real(8) :: rizy, rizx
real(8) :: deltath, ShearD
integer :: ll, i, j, k, l, ii, jj
real(8) :: Lmomentztemp, th1, th2
real(8), dimension(3) :: Lt
real(8), dimension(3) :: Omt
real(8) :: cst, snt
real(8) :: omg, omg2
real(8) :: th, diag, offd, Wr, Wd
real(8), dimension(4) :: qq
real(8), dimension(3) :: Vij, Sij, deltaV
real(8) :: RV, theta, ranf, gauss
real(8) :: ex, ex2, uij, aij, bsq, bij
! ## for check
real(8), dimension(3) :: Fdis, Rij
#ifdef DPDcheck
real(8), dimension(3) :: Fran
#endif
real(8) :: Invdt
! ##

external ranf, gauss

!************************************
!*  main part of the DPD simulation *
!************************************

   allocate( Velt(3,N) )
   allocate( dR(3,N)   )
   allocate( SdR(3,N)  )

   if(QColloid) then
     write(*,*) 'ERROR : Lowe integrator cannot be used for Colloidal system'
     call Finalize
   end if

   dt2     = 0.5d0  * deltat * deltat
   deltath = 0.5d0  * deltat
! ## for check
   Invdt   = 1.d0 / deltat
! ##

   if(ShearRate == 0.) then

     QSheared = .False.

   else

     QSheared = .True.
     ShearD  = ShearRate * deltat * CellL(2)

   end if

   SdR = 0.d0

   if(QColloid) then

     allocate( SdRc(3,NumColloid) )
     allocate(  dRc(3,NumColloid) )

     SdRc = 0.d0

   end if

   l = 0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       rizy = InvInertiaRB(3,1) - InvInertiaRB(2,1)
       rizx = InvInertiaRB(3,1) - InvInertiaRB(1,1)

       do j = 1 , NumColloid
         Omega(:,j) = Lmoment(:,j) * InvInertiaRB(:,1)
       end do

       V_RBt  = V_RB
       Omegat = Omega

       l = l + NumColloid * NumCollAtm

     else

       do j = 1 , NumMol(i)

         do k = 1 , NumAtm(i)

           l = l + 1
           Velt(:,l) = Vel(:,l)

         end do

       end do

     end if

   end do

!
! space-fixed coordinates and velocities of the colloid components
! ---------------------------------
   if(QColloid) then

     call RotationMatrix
     call SiteInfo_Colloid

   end if
! ---------------------------------

!
! Periodic Boundary Condition
! ---------------------------------
   call PBC_DPD
! ---------------------------------

! ## Force
! ----------------------------------------------------------
   call GetForceLowePeters

   FrcDP = FrcDPt

   if(QColloid) then

     call ForceTorque_COM_Colloid

     FrcCo = FrcCOt

   end if

! ----------------------------------------------------------
!
! #########################################################
!                      MAIN ROUTINE
! #########################################################
!
   do ll = 1 , Nstep

     timeps = timeps + deltat
!
! #### v~(t+dt)=v(t)+(lamda*dt)*F(t)/m
!
! ######################
     if(QColloid) then
! ######################

       l = 0

       do i = 1 , NumSpec

         if(ColloidFlag(i)) then

           do j = 1 , NumColloid

             V_RBt(:,j) = V_RB(:,j) + deltath * FrcCo(:,j) * InvMassRB(1)

           end do

! ##### L~(t+lamda*dt)=L(t)+(lamda*dt)*N(t)

           Lmoment  = Lmoment + deltath * Torque
           Lmomentt = Lmoment + deltath * Torque

           do j = 1 , NumColloid

             Lt = Lmoment(:,j)

             th1 = deltath * rizy * Lt(2)
             cst = cos(th1)
             snt = sin(th1)
             Lmoment(1,j) =   cst * Lt(1) + snt * Lt(3)
             Lmomentztemp = - snt * Lt(1) + cst * Lt(3)

             th2 = deltath * rizx * Lmoment(1,j)
             cst = cos(th2)
             snt = sin(th2)
             Lmoment(2,j) =   cst * Lt(2) - snt * Lmomentztemp
             Lmoment(3,j) =   snt * Lt(2) + cst * Lmomentztemp

           end do

           do j = 1 , NumColloid

             Lt = Lmomentt(:,j)

             th1 = deltath * rizy * Lt(2)
             cst = cos(th1)
             snt = sin(th1)
             Lmomentt(1,j) =   cst * Lt(1) + snt * Lt(3)
             Lmomentztemp  = - snt * Lt(1) + cst * Lt(3)

             th2 = deltath * rizx * Lmomentt(1,j)
             cst = cos(th2)
             snt = sin(th2)
             Lmomentt(2,j) =   cst * Lt(2) - snt * Lmomentztemp
             Lmomentt(3,j) =   snt * Lt(2) + cst * Lmomentztemp

           end do

! ##### Omg~(t+lamda*dt) = I^-1 L~(t+lamda*dt)

           do j = 1 , NumColloid

             Omega(:,j)  = Lmoment(:,j)  * InvInertiaRB(:,1)
             Omegat(:,j) = Lmomentt(:,j) * InvInertiaRB(:,1)

           end do

! ##### Q(t+dt)=exp(dt/2*A[Omg~(t+lamda*dt)]) Q(t)

           do j = 1 , NumColloid

             qq = Quaternion(:,j)

             Omt = Omega(:,j)

             omg2 = dot_product( Omt, Omt )
             omg  = sqrt( omg2 )
             th   = deltath * omg
             diag = dcos(th)
             offd = dsin(th) / omg

             Quaternion(1,j) = diag * qq(1) + offd * ( - Omt(1) * qq(2) &
             &               - Omt(2) * qq(3) - Omt(3) * qq(4) )
             Quaternion(2,j) = diag * qq(2) + offd * (   Omt(1) * qq(1) &
             &               - Omt(2) * qq(4) + Omt(3) * qq(3) )
             Quaternion(3,j) = diag * qq(3) + offd * (   Omt(1) * qq(4) &
             &               + Omt(2) * qq(1) - Omt(3) * qq(2) )
             Quaternion(4,j) = diag * qq(4) + offd * ( - Omt(1) * qq(3) &
             &               + Omt(2) * qq(2) + Omt(3) * qq(1) )

           end do
!
! ##### A(t+dt)  from  Q(t+dt)
! ----------------------------------------------------------
           call RotationMatrix
! ----------------------------------------------------------

           do j = 1 , NumColloid

             dRc(:,j) = deltat * V_RB(:,j) + dt2 * FrcCo(:,j) * InvMassRB(1)

           end do

           R_RB = R_RB + dRc


           l = l + NumColloid * NumCollAtm

         else

           do j = 1 , NumMol(i)

             do k = 1 , NumAtm(i)

               l = l + 1
               Velt(:,l) = Vel(:,l) + deltath * FrcDP(:,l)
               dR  (:,l) = deltat * Vel(:,l) + dt2 * FrcDP(:,l)
               R   (:,l) = R(:,l) + dR(:,l)

             end do

           end do

         end if

       end do

! ######################
     else
! ######################

       Velt = Vel + FrcDP * deltath

       dR = Velt * deltat

       R  = R + dR

! ######################
     end if
! ######################

! Lees-Edwards

     if(QSheared) then

       SlideGap = SlideGap + ShearD
       if( SlideGap > CellL(1) ) SlideGap = SlideGap - CellL(1)

     end if

! ##
! ## NOTE!  Cell List Method
! ## The SiteInfo_Colloid subroutine must be called before the PBC one.
! ##
! space-fixed coordinates and velocities of the colloid components
! --------------------------------------
     if(QColloid) call SiteInfo_Colloid
! --------------------------------------

! Periodic Boundary Condition (Sheared)
! ---------------------------------
     call PBC_DPD
! ---------------------------------

! ------------------------------------------------------------------
   call GetForceLowePeters

   if(QColloid) call ForceTorque_COM_Colloid
! ------------------------------------------------------------------
!
! ##### Free Rotation on the Lx axis, succeedingly on the Ly axis
!
! ######################
     if(QColloid) then
! ######################

       l = 0
       do i = 1 , NumSpec

         if(ColloidFlag(i)) then

           do j = 1 , NumColloid

             Lt = Lmoment(:,j)

             th2 = deltath * rizx * Lt(1)
             cst = cos(th2)
             snt = sin(th2)
             Lmoment(2,j) =   cst * Lt(2) - snt * Lt(3)
             Lmomentztemp =   snt * Lt(2) + cst * Lt(3)

             th1 = deltath * rizy * Lmoment(2,j)
             cst = cos(th1)
             snt = sin(th1)
             Lmoment(1,j) =   cst * Lt(1) + snt * Lmomentztemp
             Lmoment(3,j) = - snt * Lt(1) + cst * Lmomentztemp

           end do

! ##### L(t+dt)=L~(t+dt)+(dt/2)*N(t+dt)

           Lmoment = Lmoment + deltath * Torque

! #### v(t+dt)=v(t)+(dt/2)*(F(t)+F(t+dt))/m

           do j = 1 , NumColloid

             V_RB(:,j) = V_RB(:,j) + deltath * (FrcCo(:,j) + FrcCot(:,j)) * InvMassRB(1)

           end do

           l = l + NumColloid * NumCollAtm

         else

           do j = 1 , NumMol(i)

             do k = 1 , NumAtm(i)

               l = l + 1
               Vel(:,l) = Vel(:,l) + deltath * (FrcDP(:,l) + FrcDPt(:,l))

             end do

           end do

         end if

       end do

! ######################
     else
! ######################

       Vel = Vel + deltath * (FrcDP + FrcDPt)

! ######################
     end if
! ######################

     FrcDP  = FrcDPt

     if(QColloid) FrcCo = FrcCot

! ## thermalization

     if(QColloid) call SiteVeloc_Colloid

     if(IntegrMethod=='Lowe') then

       do l = 1, Npair

         if(ranf() < gammt) then

           i = ListIJ(1,l)
           j = ListIJ(2,l)

           Vij    = Vel(:,i) - Vel(:,j)
           Vij(1) = Vij(1) + SLList(l)
           Sij    = dRList(:,l)

           RV     = dot_product( Vij, Sij )

           theta  = gauss() * sigmt

           deltaV = 0.5d0 * ( theta - RV ) * Sij

           Vel(:,i) = Vel(:,i) + deltaV
           Vel(:,j) = Vel(:,j) - deltaV

         end if

       end do

     else if(IntegrMethod=='Peters') then

! ## for check
       VirialD = 0.d0
#ifdef DPDcheck
       VirialR = 0.d0
#endif
! ##

       do l = 1, Npair

         ii = Npair - l + 1
         jj = int(ranf()*ii) + 1

         i = ListIJ(1,jj)
         j = ListIJ(2,jj)

         Vij    = Vel(:,i) - Vel(:,j)
         Vij(1) = Vij(1) + SLList(jj)
         Sij    = dRList(:,jj)

         RV  = dot_product( Vij, Sij )

         gammt = GammDP(TypeNum(i),TypeNum(j))

         Wr = 1.d0 - R1List(jj)
         Wd = Wr * Wr

         ex  = exp( - 2.d0 * gammt * Wd )
         ex2 = ex * ex

         uij = 0.5d0
         aij = - uij * (1.d0 - ex) * RV

         bsq   = uij * (1.d0 - ex2)
         theta = ranf() - 0.5d0
         bij   = sigmt * sqrt(bsq) * theta

! ## for check
         Fdis  = aij * Sij * Invdt
#ifdef DPDcheck
         Fran  = bij * Sij * Invdt
#endif

         Rij = Sij * R1List(jj)

         VirialD(1,1) = VirialD(1,1) + Fdis(1) * Rij(1)
         VirialD(1,2) = VirialD(1,2) + Fdis(1) * Rij(2)
         VirialD(1,3) = VirialD(1,3) + Fdis(1) * Rij(3)
         VirialD(2,2) = VirialD(2,2) + Fdis(2) * Rij(2)
         VirialD(2,3) = VirialD(2,3) + Fdis(2) * Rij(3)
         VirialD(3,3) = VirialD(3,3) + Fdis(3) * Rij(3)
#ifdef DPDcheck
         VirialR(1,1) = VirialR(1,1) + Fran(1) * Rij(1)
         VirialR(1,2) = VirialR(1,2) + Fran(1) * Rij(2)
         VirialR(1,3) = VirialR(1,3) + Fran(1) * Rij(3)
         VirialR(2,2) = VirialR(2,2) + Fran(2) * Rij(2)
         VirialR(2,3) = VirialR(2,3) + Fran(2) * Rij(3)
         VirialR(3,3) = VirialR(3,3) + Fran(3) * Rij(3)
#endif
         deltaV = (aij + bij) * Sij

         Vel(:,i) = Vel(:,i) + deltaV
         Vel(:,j) = Vel(:,j) - deltaV

         ListIJ(1,jj) = ListIJ(1,ii)
         ListIJ(2,jj) = ListIJ(2,ii)

         dRList(:,jj) = dRList(:,ii)

         SLList(jj) = SLList(ii)
         R1List(jj) = R1List(ii)

       end do

       VirialD(2,1) = VirialD(1,2)
       VirialD(3,1) = VirialD(1,3)
       VirialD(3,2) = VirialD(2,3)
#ifdef DPDcheck
       VirialR(2,1) = VirialR(1,2)
       VirialR(3,1) = VirialR(1,3)
       VirialR(3,2) = VirialR(2,3)
#endif

       if(QSheared) Virial = Virial + VirialD

     end if

!     call Elim_CellMoveDP

     call DPD_monitor(ll)

   end do

end subroutine IntegrDPD_Lowe
