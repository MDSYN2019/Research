! ############################
! ## SUBROUTINE LIST 
! ## -- Elim_CellMove 
! ## -- Elim_CellMoveDP 
! ## -- Elim_CellRot 
! ############################


!#####################################################################
!#####################################################################


! *********************************************************
! ** In this subroutine average momentum of all particle **
! ** is calculated. This calculation is necessary for    **
! ** eliminating the cell(system) translational motion.  **
! *********************************************************

subroutine Elim_CellMove

use Numbers, only : N
use CommonBlocks, only : QRigidBody
use Configuration, only : Vel
use RBparam, only : NumRB, V_RB, QSingle, RBType, MassRB, NumRBAtom
use AtomParam, only : Mass
use ThermoData, only : CellMomentum

implicit NONE

integer :: i, k, MyType
real(8) :: TotalMass

   CellMomentum = 0.d0
   TotalMass    = 0.d0

   if(QRigidBody) then

     k = 0

     do i = 1 , NumRB

       if(QSingle(i)) then

         k = k + 1
         CellMomentum = CellMomentum + Mass(k) * V_RB(:,i)
         TotalMass    = TotalMass    + Mass(k)

       else

         MyType = RBType(i)
         CellMomentum = CellMomentum + MassRB(MyType) * V_RB(:,i)
         TotalMass    = TotalMass    + MassRB(MyType)
         k = k + NumRBAtom(MyType)

       end if


     end do

     CellMomentum = CellMomentum / TotalMass

     do i = 1 , NumRB

       V_RB(:,i) = V_RB(:,i) - CellMomentum

     end do

   else

     do i = 1 , N

       CellMomentum = CellMomentum + Mass(i) * Vel(:,i)
       TotalMass    = TotalMass    + Mass(i)

     end do

     CellMomentum = CellMomentum / TotalMass

     do i = 1 , N

       Vel(:,i) = Vel(:,i) - CellMomentum

     end do

   end if

end subroutine Elim_CellMove


!#####################################################################
!#####################################################################


! *********************************************************
! ** In this subroutine average momentum of all particle **
! ** is calculated. This calculation is necessary for    **
! ** eliminating the cell(system) translational motion.  **
! *********************************************************

subroutine Elim_CellMoveDP

use Numbers, only : N, NumSpec, NumMol, NumAtm
use Configuration, only : Vel
use CommonDPD
use RBparam, only : MassRB, V_RB
use AtomParam, only : Mass
use ThermoData, only : CellMomentum

implicit NONE

integer :: i , j, k, l
real(8) :: TotalMass

   Cellmomentum = 0.d0
   TotalMass    = 0.d0

! ########################
   if(QColloid) then
! ########################

     l = 0

     do i = 1 , NumSpec

       if(ColloidFlag(i)) then

         do j = 1 , NumColloid

           Cellmomentum = Cellmomentum + MassRB(1) * V_RB(:,j)
           TotalMass    = TotalMass    + MassRB(1)

         end do

         l = l + NumColloid * NumCollAtm

       else

         do j = 1 , NumMol(i)

           do k = 1 , NumAtm(i)

             l = l + 1
             Cellmomentum = Cellmomentum + Vel(:,l)
             TotalMass    = TotalMass    + 1.d0

           end do

         end do

       end if

     end do

! ########################
   else
! ########################

     do i = 1 , N
       Cellmomentum = Cellmomentum + Vel(:,i)
     end do

     TotalMass = dble(N)

! ########################
   end if
! ########################

   Cellmomentum = Cellmomentum / TotalMass

! ########################
   if(QColloid) then
! ########################

     l = 0

     do i = 1 , NumSpec

       if(ColloidFlag(i)) then

         do j = 1 , NumColloid

           V_RB(:,j) = V_RB(:,j) - Cellmomentum

         end do

         l = l + NumColloid * NumCollAtm

       else

         do j = 1 , NumMol(i)

           do k = 1 , NumAtm(i)

             l = l + 1
             Vel(:,l) = Vel(:,l) - Cellmomentum

           end do

         end do

       end if

     end do

! ########################
   else
! ########################

     do i = 1 , N

       Vel(:,i) = Vel(:,i) - Cellmomentum

     end do

! ########################
   end if
! ########################

end subroutine Elim_CellMoveDP


!#####################################################################
!#####################################################################


! *********************************************************
! ** In this subroutine average momentum of all particle **
! ** is calculated. This calculation is necessary for    **
! ** eliminating the cell(system) translational motion.  **
! *********************************************************

subroutine Elim_CellRot

use Numbers, only : N
use CommonBlocks, only : QMaster, QRigidBody
use Configuration
use RBparam, only : NumRB, R_RB, V_RB, QSingle, RBType, MassRB, NumRBAtom
use AtomParam, only : Mass

implicit NONE

integer :: i, k, MyType
real(8) :: M, TotalMass, rgmax
real(8), dimension(3,3) :: Isys, InvIsys
real(8), dimension(3) :: Om_sys, Lsys, RVtemp, Rtemp, Vtemp, Rg
intrinsic max

   Lsys = 0.d0
   Isys = 0.d0

   if(QRigidBody) then

     k = 0

     Rg        = 0.d0
     TotalMass = 0.d0

     do i = 1, NumRB

       if(QSingle(i)) then

         k = k + 1
         M = Mass(k)

       else

         MyType = RBType(i)
         M = MassRB(MyType)
         k = k + NumRBAtom(MyType)

       end if

       Rg(:)     = Rg(:)     + M * R_RB(:,i)
       TotalMass = TotalMass + M

     end do

     Rg(:) = Rg(:) / TotalMass

     do i = 1, NumRB

       R_RB(:,i) = R_RB(:,i) - Rg(:)

     end do

     k = 0

     do i = 1 , NumRB

       if(QSingle(i)) then

         k = k + 1
         M = Mass(k)

       else

         MyType = RBType(i)
         M = MassRB(MyType)
         k = k + NumRBAtom(MyType)

       end if

       Rtemp(:) = R_RB(:,i)
       Vtemp(:) = V_RB(:,i)
       call Prodv(Rtemp,Vtemp,RVtemp)
       Lsys(:) = Lsys(:) + M * RVtemp(:)

       Isys(1,1) = Isys(1,1) + M * ( R_RB(2,i)**2 + R_RB(3,i)**2 )
       Isys(2,2) = Isys(2,2) + M * ( R_RB(3,i)**2 + R_RB(1,i)**2 )
       Isys(3,3) = Isys(3,3) + M * ( R_RB(1,i)**2 + R_RB(2,i)**2 )
       Isys(1,2) = Isys(1,2) - M * R_RB(1,i) * R_RB(2,i)
       Isys(1,3) = Isys(1,3) - M * R_RB(1,i) * R_RB(3,i)
       Isys(2,3) = Isys(2,3) - M * R_RB(2,i) * R_RB(3,i)

     end do

     call InversMatrix( Isys, InvIsys )

     Om_sys = matmul( InvIsys, Lsys )

     do i = 1 , NumRB

       Rtemp(:) = R_RB(:,i)
       call Prodv(Om_sys, Rtemp, Vtemp)

       V_RB(:,i) = V_RB(:,i) - Vtemp(:)

     end do

   else

     Rg        = 0.d0
     TotalMass = 0.d0

     do i = 1, N

       Rg(:)     = Rg(:)     + Mass(i) * R(:,i)
       TotalMass = TotalMass + Mass(i)

     end do

     Rg(:) = Rg(:) / TotalMass

     do i = 1, N

       R(:,i) = R(:,i) - Rg(:)

     end do

     do i = 1 , N

       M = Mass(i)

       Rtemp(:) = R(:,i)
       Vtemp(:) = Vel(:,i)
       call Prodv(Rtemp,Vtemp,RVtemp)
       Lsys(:) = Lsys(:) + M * RVtemp(:)

       Isys(1,1) = Isys(1,1) + M * ( R(2,i)**2 + R(3,i)**2 )
       Isys(2,2) = Isys(2,2) + M * ( R(3,i)**2 + R(1,i)**2 )
       Isys(3,3) = Isys(3,3) + M * ( R(1,i)**2 + R(2,i)**2 )
       Isys(1,2) = Isys(1,2) - M * R(1,i) * R(2,i)
       Isys(1,3) = Isys(1,3) - M * R(1,i) * R(3,i)
       Isys(2,3) = Isys(2,3) - M * R(2,i) * R(3,i)

     end do

     call InversMatrix( Isys, InvIsys )

     Om_sys = matmul( InvIsys, Lsys )

     do i = 1 , N

       Rtemp(:) = R(:,i)
       call Prodv(Om_sys, Rtemp, Vtemp)

       Vel(:,i) = Vel(:,i) - Vtemp(:)

     end do

   end if

   rgmax = max( Rg(1), Rg(2), Rg(3) )

   if(rgmax > 1.d-8 .and. QMaster ) then
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     write(*,*) 'WARNING : the system is traslocated to have its COM at zero!'
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
   end if


Contains

   subroutine Prodv(x,y,z)

   implicit none

   real(8), dimension(3) :: x, y, z

      z(1) = x(2) * y(3) - y(2) * x(3)
      z(2) = x(3) * y(1) - y(3) * x(1)
      z(3) = x(1) * y(2) - y(1) * x(2)

   end subroutine Prodv

end subroutine Elim_CellRot
