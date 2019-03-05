! ############################
! ## SUBROUTINE LIST 
! ## -- SHAKE 
! ## -- RATTLE 
! ## -- SHAKEROLL 
! ## -- RATTLEROLL 
! ## -- SHAKEROLLA3 
! ## -- RATTLEROLLA3 
! ## -- SHAKEROLLPR 
! ## -- RATTLEROLLPR 
! ############################


!#####################################################################
!#####################################################################


! *******************************************************************
! ** Constraint dynamics of a chain of atoms using SHAKE.          **
! ** Reference:                                                    **
! ** HC Andersen, J. Comput. Phys. 52, 24, 1983.                   **
! *******************************************************************
! ***************
! **   SHAKE   **
! ***************

subroutine SHAKE

use CommonBlocks, only : QPBC
use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat, dt2

implicit NONE

logical done
logical, dimension(MaxHcBond+1) :: Moving, Moved
real(8), dimension(MaxHcBond+1) :: RMass
real(8), dimension(3,MaxHcBond+1) :: Ri, Rf, Vf
real(8), dimension(3) :: Rfij, dr, drs
real(8), dimension(3,3) :: VircA
integer :: Iteration

real(8), parameter :: RpTol=1.0d-28
real(8) :: Tolerance2, InvDt, d2, dq, R2, drf2
real(8) :: InvMa, InvMb, gab

integer :: i, j, k, kk

!*******************************************************************

   Tolerance2  = 2.d0 * TolA
   InvDt = 1.d0 / deltat

   if(QPBC) VircA = 0.d0

!  ** loop over molecules **

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Ri(:,kk)   = R_o(:,i)
       Rf(:,kk)   = R  (:,i)
       Vf(:,kk)   = Vel(:,i)
       RMass (kk) = InvMass(i)
       Moving(kk) = .False.
       Moved (kk) = .True.

     end do

     Iteration = 0
     done = .False.

!  ** start of iterative loop **

150  if((.not.done).and.(Iteration <= MaxIteration)) then

       done = .True.

       do kk = 1 , NCoupleBond(k)

         i = CouplePair(k,kk,1)
         j = CouplePair(k,kk,2)

         if( Moved(i) .or. Moved(j) ) then

           Rfij = Rf(:,i) - Rf(:,j)

           R2 = dot_product(Rfij,Rfij)
           d2 = rSHAKE(k,kk)
           dq = d2 - R2

           if(abs(dq) > (d2*Tolerance2)) then

           dr = Ri(:,i) - Ri(:,j)
           drf2 = dot_product(dr,Rfij)

           if( drf2 < (d2*RpTol) ) then

             write(*,'(a)') ' Error : constraint '
             write(*,*) 'Group=',k,'Bond=',kk
             call Finalize

           endif

           InvMa = RMass(i)
           InvMb = RMass(j)
           gab = dq / ( 2.d0 * ( InvMa + InvMb ) * drf2 )
           drs = dr * gab

           if(QPBC) then

             VircA(1,1) = VircA(1,1) + gab * dr(1) * dr(1)
             VircA(1,2) = VircA(1,2) + gab * dr(1) * dr(2)
             VircA(1,3) = VircA(1,3) + gab * dr(1) * dr(3)
             VircA(2,2) = VircA(2,2) + gab * dr(2) * dr(2)
             VircA(2,3) = VircA(2,3) + gab * dr(2) * dr(3)
             VircA(3,3) = VircA(3,3) + gab * dr(3) * dr(3)

           end if

           Rf(:,i) = Rf(:,i) + InvMa * drs
           Rf(:,j) = Rf(:,j) - InvMb * drs

           drs = drs * InvDt

           Vf(:,i) = Vf(:,i) + InvMa * drs
           Vf(:,j) = Vf(:,j) - InvMb * drs

           Moving(i) = .True.
           Moving(j) = .True.
           done      = .False.

           endif

         endif

       end do

       do kk = 1 , NCoupleAtom(k)

         Moved(kk)  = Moving(kk)
         Moving(kk) = .False.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!   ** end of iterative loop **

     if(.not.done) then

       write(*,'(a)') ' Error : Too many iterations in SHAKE '
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

!  ** store away new values    **

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       R  (:,i) = Rf(:,kk)
       Vel(:,i) = Vf(:,kk)

     end do

!    write(*,*) 'it=',it

   end do

   if(QPBC) then

     VircA(2,1) = VircA(1,2)
     VircA(3,1) = VircA(1,3)
     VircA(3,2) = VircA(2,3)

     Vir_SHAKE = VircA / dt2

   end if

!   ** end of loop over molecules **

end subroutine SHAKE


!#####################################################################
!#####################################################################


! ****************
! **   RATTLE   **
! ****************

subroutine RATTLE

use CommonBlocks, only : QPBC
use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : dt2

implicit NONE

logical :: done
logical, dimension(MaxHcBond+1) :: Moving,Moved
real(8), dimension(MaxHcBond+1) :: RMass
real(8), dimension(3,MaxHcBond+1) :: Ri,Vi
real(8), dimension(3,3) :: VircB
real(8) :: RV, InvMa, InvMb, gab
real(8), dimension(3) :: dv, dr, drs

integer :: i, j, k, kk, Iteration

!*******************************************************************
!** second part of velocity verlet with constraints               **
!*******************************************************************

   if(QPBC) VircB = 0.d0

!  ** loop over all molecules **

   do k = 1 , NSHAKEGroup

!  ** velocity verlet algorithm part b **

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Ri(:,kk)   = R  (:,i)
       Vi(:,kk)   = Vel(:,i)
       RMass (kk) = InvMass(i)
       Moving(kk) = .False.
       Moved (kk) = .True.

     end do

!  ** start of iterative loop **

     Iteration = 0
     done = .false.

150  if((.not.done) .and. (Iteration < MaxIteration)) then

       done=.True.

       do kk = 1 , NCoupleBond(k)

         i = CouplePair(k,kk,1)
         j = CouplePair(k,kk,2)

         if( Moved(i) .or. Moved(j) ) then

           dv = Vi(:,i) - Vi(:,j)
           dr = Ri(:,i) - Ri(:,j)
           RV = dot_product(dr,dv)
           InvMa = RMass(i)
           InvMb = RMass(j)
           gab = -RV / ( ( InvMa + InvMb ) * rSHAKE(k,kk) )

           if( abs(gab) > TolB ) then

             if(QPBC) then

               VircB(1,1) = VircB(1,1) + gab * dr(1) * dr(1)
               VircB(1,2) = VircB(1,2) + gab * dr(1) * dr(2)
               VircB(1,3) = VircB(1,3) + gab * dr(1) * dr(3)
               VircB(2,2) = VircB(2,2) + gab * dr(2) * dr(2)
               VircB(2,3) = VircB(2,3) + gab * dr(2) * dr(3)
               VircB(3,3) = VircB(3,3) + gab * dr(3) * dr(3)

             end if

             drs = dr * gab

             Vi(:,i) = Vi(:,i) + InvMa * drs
             Vi(:,j) = Vi(:,j) - InvMb * drs

             Moving(i)=.True.
             Moving(j)=.True.
             done=.False.

           endif

         endif

       end do

       do kk = 1 , NCoupleAtom(k)

         Moved(kk)  = Moving(kk)
         Moving(kk) = .False.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!  ** end of iterative loop **

     if (.not. done) then

       write(*,'(a)') ' Error : Too many constraint iterations in Rattle '
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Vel(:,i) = Vi(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration

   end do

!  ** end of loop over molecules **

   if(QPBC) then

     VircB(2,1) = VircB(1,2)
     VircB(3,1) = VircB(1,3)
     VircB(3,2) = VircB(2,3)

     Vir_RATTL = VircB / dt2

   end if

end subroutine RATTLE


!#####################################################################
!#####################################################################


! *******************************************************************
! ** constraint dynamics of a chain of atoms using SHAKE.          **
! ** Reference:                                                    **
! ** HC Andersen, J. Comput. Phys. 52, 24, 1983.                   **
! *******************************************************************
! *********************
! **   SHAKE/ROLL    **
! *********************

subroutine SHAKEROLL(Rv,ROLL)

use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat

implicit NONE

logical :: done, ROLL
logical, dimension(MaxHcBond+1) :: moving, moved
real(8), dimension(MaxHcBond+1) :: RMass
real(8), dimension(3,MaxHcBond+1)   :: Ro, Rguess, Vguess
real(8), dimension(3)      :: Rlamdaini, Roij, Rorot, dR, dV
real(8) :: Rv
integer :: Iteration

real(8), parameter :: RpTol=1.0d-28

real(8) :: Tol2, InvDt, InvDt2, Rlamdaini2, d2, Rgap, RlRor
real(8) :: InvMa, InvMb, dlamda
integer :: i, j, k
integer :: kk

!*******************************************************************

   Tol2  = 2.d0 * TolA
   InvDt = 1.d0 / deltat
   InvDt2= InvDt * InvDt

   ROLL = .False.

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Ro    (:,kk) = R_o(:,i)
       Rguess(:,kk) = R  (:,i)
       Vguess(:,kk) = Vel(:,i)
       RMass (kk)   = InvMass(i)
       moving(kk)   = .False.
       moved (kk)   = .True.

     end do

     Iteration = 0
     done = .False.

!  ** start of iterative loop **

150  if((.not.done).and.(Iteration <= MaxIteration)) then

       done = .true.

       do kk = 1 , NCoupleBond(k)

         i = CouplePair(k,kk,1)
         j = CouplePair(k,kk,2)

         if( moved(i) .or. moved(j) ) then

           Rlamdaini = Rguess(:,i) - Rguess(:,j)

           Rlamdaini2 = dot_product(Rlamdaini,Rlamdaini)
           d2         = rSHAKE(k,kk)
           Rgap       = d2 - Rlamdaini2

           if( abs(Rgap) > (d2*Tol2) ) then

             ROLL = .True.

             Roij  = Ro(:,i) - Ro(:,j)
             Rorot = Rv * Roij
             RlRor = dot_product( Rorot , Rlamdaini )

             if( RlRor < ( d2 * RpTol ) ) then

               write(*,'(a)') ' Error : constraint '
               write(*,*) 'Group=',k,'Bond=',kk
               call Finalize

             endif

             InvMa  = RMass(i)
             InvMb  = RMass(j)
             dlamda = Rgap / ( ( InvMa + InvMb ) * RlRor )

             Lagmultip(k,kk) = Lagmultip(k,kk) + dlamda * InvDt2

             dR = Rorot * dlamda * 0.5d0

             Rguess(:,i) = Rguess(:,i) + InvMa * dR
             Rguess(:,j) = Rguess(:,j) - InvMb * dR

             dV = Roij * dlamda * InvDt * 0.5d0

             Vguess(:,i) = Vguess(:,i) + InvMa * dV
             Vguess(:,j) = Vguess(:,j) - InvMb * dV

             moving(i) = .True.
             moving(j) = .True.
             done      = .False.

           endif

         endif

       end do

       do kk = 1 , NCoupleAtom(k)

         moved(kk)  = moving(kk)
         moving(kk) = .false.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!   ** end of iterative loop **

     if(.not.done) then

       write(*,'(a)') ' Error : Too many constraint iterations in SHAKE '
       write(*,*) Iteration
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

!  ** store away new values    **

     do kk = 1, NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       R  (:,i) = Rguess(:,kk)
       Vel(:,i) = Vguess(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration

   end do

end subroutine SHAKEROLL


!#####################################################################
!#####################################################################


! **********************
! **   RATTLE/ROLL    **
! **********************

subroutine RATTLEROLL(Vscale,ROLL)

use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat

implicit NONE

logical :: done, ROLL
logical, dimension(MaxHcBond+1) :: moving, moved
real(8), dimension(MaxHcBond+1) :: RMass
real(8), dimension(3,MaxHcBond+1) :: Rfixed, Vguess
real(8) :: RV, InvMa, InvMb, dlamda, Vscale
real(8), dimension(3) :: dV, Vlamdaini, Rdeltat

real(8) :: InvDt
integer :: i, j, k, Iteration
integer :: kk

   InvDt = 1.d0 / deltat

   ROLL = .False.

   do k = 1 , NSHAKEGroup

!  ** velocity verlet algorithm part b **

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Rfixed(:,kk) = R  (:,i)
       Vguess(:,kk) = Vel(:,i)
       RMass (kk)   = InvMass(i)
       moving(kk)   = .false.
       moved (kk)   = .true.

     end do

!  ** start of iterative loop **

     Iteration = 0
     done = .false.

150 if((.not.done) .and. (Iteration < MaxIteration)) then

       done=.true.

       do kk = 1 , NCoupleBond(k)

         i  = CouplePair(k,kk,1)
         j  = CouplePair(k,kk,2)

         if( moved(i) .or. moved(j) ) then

           Vlamdaini = Vguess(:,i) - Vguess(:,j)
           Rdeltat   = Rfixed(:,i) - Rfixed(:,j)
           RV        = dot_product( Rdeltat , Vlamdaini )
           InvMa = RMass(i)
           InvMb = RMass(j)
           dlamda = -RV / ( Vscale * ( InvMa + InvMb ) &
           &            * dot_product( Rdeltat, Rdeltat ) )

           if( abs(dlamda) > TolB ) then

             ROLL = .True.

             Lagmultip(k,kk) = Lagmultip(k,kk) + dlamda * 2.d0 * InvDt

             dV = Vscale * dlamda * Rdeltat

             Vguess(:,i) = Vguess(:,i) + InvMa * dV
             Vguess(:,j) = Vguess(:,j) - InvMb * dV

             moving(i)=.true.
             moving(j)=.true.
             done=.false.

           endif

         endif

       end do

       do kk = 1 , NCoupleBond(k)

         moved(kk)  = moving(kk)
         moving(kk) = .false.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!  ** end of iterative loop **

     if (.not. done) then

       write(*,'(a)') ' Error : Too many constraint iterations in Rattle '
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Vel(:,i) = Vguess(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration
   end do

end subroutine RATTLEROLL


!#####################################################################
!#####################################################################


! *******************************************************************
! ** constraint dynamics of a chain of atoms using SHAKE.          **
! ** Reference:                                                    **
! ** HC Andersen, J. Comput. Phys. 52, 24, 1983.                   **
! *******************************************************************
! *********************
! **   SHAKE/ROLL    **
! *********************

subroutine SHAKEROLLA3(Rv,ROLL)

use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat

implicit NONE

logical :: done, ROLL
logical, dimension(MaxHcBond+1) :: moving, moved
real(8), dimension(MaxHcBond+1) :: RMass
real(8), dimension(3,MaxHcBond+1)   :: Ro, Rguess, Vguess
real(8), dimension(3) :: Rlamdaini, Roij, Rorot, dR, dV
real(8), dimension(3) :: Rv
integer :: Iteration

real(8), parameter :: RpTol=1.0d-28

real(8) :: Tol2, InvDt, InvDt2, Rlamdaini2, d2, Rgap, RlRor
real(8) :: InvMa, InvMb, dlamda
integer :: i, j, k
integer :: kk

!*******************************************************************

   Tol2  = 2.d0 * TolA
   InvDt = 1.d0 / deltat
   InvDt2= InvDt * InvDt

   ROLL = .False.

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Ro    (:,kk) = R_o(:,i)
       Rguess(:,kk) = R  (:,i)
       Vguess(:,kk) = Vel(:,i)
       RMass (kk)   = InvMass(i)
       moving(kk)   = .False.
       moved (kk)   = .True.

     end do

     Iteration = 0
     done = .False.

!  ** start of iterative loop **

150  if((.not.done).and.(Iteration <= MaxIteration)) then

       done = .true.

       do kk = 1 , NCoupleBond(k)

         i  = CouplePair(k,kk,1)
         j  = CouplePair(k,kk,2)

         if( moved(i) .or. moved(j) ) then

           Rlamdaini = Rguess(:,i) - Rguess(:,j)

           Rlamdaini2 = dot_product(Rlamdaini,Rlamdaini)
           d2         = rSHAKE(k,kk)
           Rgap       = d2 - Rlamdaini2

           if( abs(Rgap) > (d2*Tol2) ) then

             ROLL = .True.

             Roij  = Ro(:,i) - Ro(:,j)
             Rorot = Rv * Roij
             RlRor = dot_product( Rorot , Rlamdaini )

             if( RlRor < ( d2 * RpTol ) ) then

               write(*,'(a)') ' Error : constraint '
               write(*,*) 'Group=',k,'Bond=',kk
               call Finalize

             endif

             InvMa  = RMass(i)
             InvMb  = RMass(j)
             dlamda = Rgap / ( ( InvMa + InvMb ) * RlRor )

             Lagmultip(k,kk) = Lagmultip(k,kk) + dlamda * InvDt2

             dR = Rorot * dlamda * 0.5d0

             Rguess(:,i) = Rguess(:,i) + InvMa * dR
             Rguess(:,j) = Rguess(:,j) - InvMb * dR

             dV = Roij * dlamda * InvDt * 0.5d0

             Vguess(:,i) = Vguess(:,i) + InvMa * dV
             Vguess(:,j) = Vguess(:,j) - InvMb * dV

             moving(i) = .True.
             moving(j) = .True.
             done      = .False.

           endif

         endif

       end do

       do kk = 1 , NCoupleAtom(k)

         moved(kk)  = moving(kk)
         moving(kk) = .false.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!   ** end of iterative loop **

     if(.not.done) then

       write(*,'(a)') ' Error : Too many constraint iterations in SHAKE '
       write(*,*) Iteration
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

!  ** store away new values    **

     do kk = 1, NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       R  (:,i) = Rguess(:,kk)
       Vel(:,i) = Vguess(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration

   end do

end subroutine SHAKEROLLA3


!#####################################################################
!#####################################################################


! **********************
! **   RATTLE/ROLL    **
! **********************

subroutine RATTLEROLLA3(Vscale,ROLL)

use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat

implicit NONE

logical :: done, ROLL
logical, dimension(MaxHcBond+1) :: moving, moved
real(8), dimension(MaxHcBond+1) :: RMass
real(8), dimension(3,MaxHcBond+1) :: Rfixed, Vguess
real(8) :: RV, InvMa, InvMb, dlamda
real(8), dimension(3) :: Vscale
real(8), dimension(3) :: dV, Vlamdaini, Rdeltat, RotRdeltat

real(8) :: InvDt
integer :: i, j, k, Iteration
integer :: kk

   InvDt = 1.d0 / deltat

   ROLL = .False.

   do k = 1 , NSHAKEGroup

!  ** velocity verlet algorithm part b **

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Rfixed(:,kk) = R  (:,i)
       Vguess(:,kk) = Vel(:,i)
       RMass (kk)   = InvMass(i)
       moving(kk)   = .false.
       moved (kk)   = .true.

     end do

!  ** start of iterative loop **

     Iteration = 0
     done = .false.

150 if((.not.done) .and. (Iteration < MaxIteration)) then

       done = .true.

       do kk = 1 , NCoupleBond(k)

         i  = CouplePair(k,kk,1)
         j  = CouplePair(k,kk,2)

         if( moved(i) .or. moved(j) ) then

           Vlamdaini = Vguess(:,i) - Vguess(:,j)
           Rdeltat   = Rfixed(:,i) - Rfixed(:,j)
           RV        = dot_product( Rdeltat , Vlamdaini )
           InvMa = RMass(i)
           InvMb = RMass(j)
           RotRdeltat = Vscale * Rdeltat
           dlamda = -RV / ( ( InvMa + InvMb ) &
           &            * dot_product( RotRdeltat, Rdeltat ) )

           if( abs(dlamda) > TolB ) then

             ROLL = .True.

             Lagmultip(k,kk) = Lagmultip(k,kk) + dlamda * 2.d0 * InvDt

             dV = Vscale * dlamda * Rdeltat

             Vguess(:,i) = Vguess(:,i) + InvMa * dV
             Vguess(:,j) = Vguess(:,j) - InvMb * dV

             moving(i)=.true.
             moving(j)=.true.
             done=.false.

           endif

         endif

       end do

       do kk = 1 , NCoupleBond(k)

         moved(kk)  = moving(kk)
         moving(kk) = .false.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!  ** end of iterative loop **

     if (.not. done) then

       write(*,'(a)') ' Error : Too many constraint iterations in Rattle '
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Vel(:,i) = Vguess(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration

   end do

end subroutine RATTLEROLLA3


!#####################################################################
!#####################################################################


! *******************************************************************
! ** constraint dynamics of a chain of atoms using SHAKE.          **
! ** Reference:                                                    **
! ** HC Andersen, J. Comput. Phys. 52, 24, 1983.                   **
! *******************************************************************
! *********************
! **   SHAKE/ROLL    **
! *********************

subroutine SHAKEROLLPR(Rv,ROLL)

use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat

implicit NONE

logical :: done, ROLL
logical, dimension(MaxHcBond+1) :: moving, moved
real(8), dimension(MaxHcBond+1) :: RMass

real(8), dimension(3,MaxHcBond+1) :: Ro, Rguess, Vguess
real(8), dimension(3)   :: Rlamdaini, Roij, Rorot, dR, dV
real(8), dimension(3,3) :: Rv
integer :: Iteration

real(8), parameter :: RpTol=1.0d-28

real(8) :: Tol2, InvDt, InvDt2, Rlamdaini2, d2, Rgap, RlRor
real(8) :: InvMa, InvMb, dlamda
integer :: i, j, k
integer :: kk

!*******************************************************************

   Tol2  = 2.d0 * TolA
   InvDt = 1.d0 / deltat
   InvDt2= InvDt * InvDt

   ROLL = .False.

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Ro    (:,kk) = R_o(:,i)
       Rguess(:,kk) = R  (:,i)
       Vguess(:,kk) = Vel(:,i)
       RMass (kk)   = InvMass(i)
       moving(kk)   = .False.
       moved (kk)   = .True.

     end do

     Iteration = 0
     done = .False.

!  ** start of iterative loop **

150  if((.not.done) .and. (Iteration <= MaxIteration)) then

       done = .true.

       do kk = 1 , NCoupleBond(k)

         i  = CouplePair(k,kk,1)
         j  = CouplePair(k,kk,2)

         if( moved(i) .or. moved(j) ) then

           Rlamdaini = Rguess(:,i) - Rguess(:,j)

           Rlamdaini2 = dot_product(Rlamdaini,Rlamdaini)
           d2         = rSHAKE(k,kk)
           Rgap       = d2 - Rlamdaini2

           if( abs(Rgap) > (d2*Tol2) ) then

             ROLL = .True.

             Roij  = Ro(:,i) - Ro(:,j)
             Rorot = matmul( Rv , Roij )
             RlRor = dot_product( Rorot , Rlamdaini )

             if( RlRor < ( d2 * RpTol ) ) then

               write(*,'(a)') ' Error : constraint '
               write(*,*) 'Group=',k,'Bond=',kk
               call Finalize

             endif

             InvMa  = RMass(i)
             InvMb  = RMass(j)
             dlamda = Rgap / ( ( InvMa + InvMb ) * RlRor )

             Lagmultip(k,kk) = Lagmultip(k,kk) + dlamda * InvDt2

             dR = Rorot * dlamda * 0.5d0

             Rguess(:,i) = Rguess(:,i) + InvMa * dR
             Rguess(:,j) = Rguess(:,j) - InvMb * dR

             dV = Roij * dlamda * InvDt * 0.5d0

             Vguess(:,i) = Vguess(:,i) + InvMa * dV
             Vguess(:,j) = Vguess(:,j) - InvMb * dV

             moving(i) = .True.
             moving(j) = .True.
             done      = .False.

           endif

         endif

       end do

       do kk = 1 , NCoupleAtom(k)

         moved(kk)  = moving(kk)
         moving(kk) = .false.

       end do

       Iteration = Iteration + 1
       goto 150

     endif

!   ** end of iterative loop **

     if(.not.done) then

       write(*,'(a)') ' Error : Too many constraint iterations in SHAKE '
       write(*,*) Iteration
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

!  ** store away new values    **

     do kk = 1, NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       R  (:,i) = Rguess(:,kk)
       Vel(:,i) = Vguess(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration

   end do

end subroutine SHAKEROLLPR


!#####################################################################
!#####################################################################


! **********************
! **   RATTLE/ROLL    **
! **********************

subroutine RATTLEROLLPR(VelRotation,ROLL)

use Configuration
use SHAKEparam
use AtomParam, only : InvMass
use TimeParam, only : deltat

implicit NONE

logical :: done, ROLL
logical, dimension(MaxHcBond+1) :: moving, moved
real(8), dimension(MaxHcBond+1) :: RMass

real(8), dimension(3,MaxHcBond+1) :: Rfixed, Vguess
real(8) :: RV, InvMa, InvMb, dlamda
real(8), dimension(3,3) :: VelRotation
real(8), dimension(3) :: dV, Vlamdaini, Rdeltat, RotRdeltat

real(8) :: InvDt
integer :: i, j, k, Iteration
integer :: kk

   InvDt = 1.d0 / deltat

   ROLL = .False.

   do k = 1 , NSHAKEGroup

!  ** velocity verlet algorithm part b **

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Rfixed(:,kk) = R  (:,i)
       Vguess(:,kk) = Vel(:,i)
       RMass (kk)   = InvMass(i)
       moving(kk)   = .false.
       moved (kk)   = .true.

     end do

!  ** start of iterative loop **

     Iteration = 0
     done = .false.

150  if((.not.done) .and. (Iteration < MaxIteration)) then

       done=.true.

       do kk = 1 , NCoupleBond(k)

         i = CouplePair(k,kk,1)
         j = CouplePair(k,kk,2)

         if( moved(i) .or. moved(j) ) then

           Vlamdaini = Vguess(:,i) - Vguess(:,j)
           Rdeltat   = Rfixed(:,i) - Rfixed(:,j)
           RV        = dot_product( Rdeltat , Vlamdaini )

           InvMa = RMass(i)
           InvMb = RMass(j)
           RotRdeltat = matmul( VelRotation, Rdeltat )
           dlamda = -RV / ( ( InvMa + InvMb ) &
           &            * dot_product( RotRdeltat, Rdeltat ) )

           if( abs(dlamda) > TolB ) then

             ROLL = .True.

             Lagmultip(k,kk) = Lagmultip(k,kk) + dlamda * 2.d0 * InvDt

             dV = dlamda * RotRdeltat

             Vguess(:,i) = Vguess(:,i) + InvMa * dV
             Vguess(:,j) = Vguess(:,j) - InvMb * dV

             moving(i)=.true.
             moving(j)=.true.
             done=.false.

           end if

         end if

       end do

       do kk = 1 , NCoupleBond(k)

         moved(kk)  = moving(kk)
         moving(kk) = .false.

       end do

       Iteration = Iteration + 1
       go to 150

     endif

!  ** end of iterative loop **

     if (.not. done) then

       write(*,'(a)') ' Error : Too many constraint iterations in Rattle '
       write(*,'(a,i5)') ' molecule or subunit ',k
       call Finalize

     endif

     do kk = 1 , NCoupleAtom(k)

       i = CoupleAtom(k,kk)
       Vel(:,i) = Vguess(:,kk)

     end do

!       write(*,*) 'Iteration=',Iteration

   end do

end subroutine RATTLEROLLPR
