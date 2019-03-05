! ############################
! ## SUBROUTINE LIST 
! ## -- PBC 
! ## -- PBC_PI 
! ## -- PBC_DPD 
! ############################


!######################################################################
!######################################################################


! *******************************
! * Periodic Boundary Condition *
! *******************************

subroutine PBC

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QRigidBody
use Configuration, only : R
use RBparam, only : NumRBinMol, QSingle, R_RB, MassRB, NumRBAtom, RBType
use CellParam, only : H, InvH
use AtomParam, only : Mass

implicit none

integer :: i, j, Num, moltype
integer :: Ix, Iy, Iz
real(8) :: Rx, Ry, Rz
real(8) :: Sx, Sy, Sz
real(8) :: Tx, Ty, Tz
real(8) :: SumM, ms, xx
integer :: k, Ng, MyType

   if(QRigidBody) then

     Ng = 0
     k  = 0

     do moltype = 1 , NumSpec

       do i = 1 , NumMol(moltype)

         Rx   = 0.d0
         Ry   = 0.d0
         Rz   = 0.d0
         SumM = 0.d0

         do j = 1 , NumRBinMol(moltype)

           Ng = Ng + 1

           if(QSingle(Ng)) then

             k  = k + 1
             ms = Mass(k)
             Rx   = Rx   + ms * R_RB(1,Ng)
             Ry   = Ry   + ms * R_RB(2,Ng)
             Rz   = Rz   + ms * R_RB(3,Ng)
             SumM = SumM + ms

           else

             MyType = RBType(Ng)
             ms     = MassRB(MyType)
             Rx   = Rx   + ms * R_RB(1,Ng)
             Ry   = Ry   + ms * R_RB(2,Ng)
             Rz   = Rz   + ms * R_RB(3,Ng)
             SumM = SumM + ms

             k = k + NumRBAtom(MyType)

           end if

         end do

         xx = 1.d0 / SumM

         Rx = Rx * xx
         Ry = Ry * xx
         Rz = Rz * xx

         Sx = InvH(1,1)*Rx + InvH(1,2)*Ry + InvH(1,3)*Rz
         Sy = InvH(2,1)*Rx + InvH(2,2)*Ry + InvH(2,3)*Rz
         Sz = InvH(3,1)*Rx + InvH(3,2)*Ry + InvH(3,3)*Rz
         if(Sx>0.5) then
           Ix = -1
         else if(Sx<-0.5) then
           Ix =  1
         else
           Ix =  0
         end if
         if(Sy>0.5) then
           Iy = -1
         else if(Sy<-0.5) then
           Iy =  1
         else
           Iy =  0
         end if
         if(Sz>0.5) then
           Iz = -1
         else if(Sz<-0.5) then
           Iz =  1
         else
           Iz =  0
         end if
         Tx = H(1,1)*Ix + H(1,2)*Iy + H(1,3)*Iz
         Ty = H(2,1)*Ix + H(2,2)*Iy + H(2,3)*Iz
         Tz = H(3,1)*Ix + H(3,2)*Iy + H(3,3)*Iz

         do j = Ng - NumRBinMol(moltype) + 1, Ng

           R_RB(1,j) = R_RB(1,j) + Tx
           R_RB(2,j) = R_RB(2,j) + Ty
           R_RB(3,j) = R_RB(3,j) + Tz

         end do

       end do

     end do

   else

     Num = 0

     do moltype = 1 , NumSpec

       do i = 1 , NumMol(moltype)

         Rx   = 0.d0
         Ry   = 0.d0
         Rz   = 0.d0
         SumM = 0.d0

         do j = 1 , NumAtm(moltype)

           Num  = Num  + 1
           ms   = Mass(Num)
           Rx   = Rx   + ms * R(1,Num)
           Ry   = Ry   + ms * R(2,Num)
           Rz   = Rz   + ms * R(3,Num)
           SumM = SumM + ms

         end do

         xx = 1.d0 / SumM

         Rx = Rx * xx
         Ry = Ry * xx
         Rz = Rz * xx

         Sx = InvH(1,1)*Rx + InvH(1,2)*Ry + InvH(1,3)*Rz
         Sy = InvH(2,1)*Rx + InvH(2,2)*Ry + InvH(2,3)*Rz
         Sz = InvH(3,1)*Rx + InvH(3,2)*Ry + InvH(3,3)*Rz
         if(Sx>0.5) then
           Ix = -1
         else if(Sx<-0.5) then
           Ix =  1
         else
           Ix =  0
         end if
         if(Sy>0.5) then
           Iy = -1
         else if(Sy<-0.5) then
           Iy =  1
         else
           Iy =  0
         end if
         if(Sz>0.5) then
           Iz = -1
         else if(Sz<-0.5) then
           Iz =  1
         else
           Iz =  0
         end if
         Tx = H(1,1)*Ix + H(1,2)*Iy + H(1,3)*Iz
         Ty = H(2,1)*Ix + H(2,2)*Iy + H(2,3)*Iz
         Tz = H(3,1)*Ix + H(3,2)*Iy + H(3,3)*Iz

         do j = Num - NumAtm(moltype) + 1 , Num

           R(1,j) = R(1,j) + Tx
           R(2,j) = R(2,j) + Ty
           R(3,j) = R(3,j) + Tz

         end do

       end do

     end do

   end if

end subroutine PBC


!######################################################################
!######################################################################


! *******************************
! * Periodic Boundary Condition *
! *******************************

subroutine PBC_PI

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonPI
use CellParam, only : H, InvH

implicit none

integer :: i, j, Num, moltype
integer, dimension(3) :: Ic
real(8), dimension(3) :: Rg
real(8), dimension(3) :: Si
real(8), dimension(3) :: Tc
real(8) :: SumM

   Num = 0

   do moltype = 1 , NumSpec

     do i = 1 , NumMol(moltype)

       Rg   = 0.d0
       SumM = 0.d0

       do j = 1 , NumAtm(moltype)

         Num  = Num  + 1
         Rg   = Rg   + FictMass(Num,1) * Rnm(:,Num,1)
         SumM = SumM + FictMass(Num,1)

       end do

       Rg = Rg / SumM

       Si = matmul( InvH, Rg )
       Ic = -nint( Si )
       Tc = matmul( H , Ic )

       do j = Num - NumAtm(moltype) + 1 , Num

         Rnm(:,j,1) = Rnm(:,j,1) + Tc

       end do

     end do

   end do

end subroutine PBC_PI


!#####################################################################
!#####################################################################


subroutine PBC_DPD

use Numbers, only : N
use Configuration
use CommonDPD
use RBparam, only : R_RB, V_RB
use CellParam, only : CellL, InvCL

implicit none

integer :: i, j
real(8), dimension(3) :: sr
real(8) :: vv

! Periodic Boundary Condition

! ###################
   if(QSheared) then
! ###################

     if(QColloid) then   ! NOTE!! only for colloidal center of mass

       do j = 1 , NumColloid

         vv    = -nint(R_RB(2,j) * InvCL(2))
         sr(2) = vv * CellL(2)
         sr(3) = -nint(R_RB(3,j) * InvCL(3)) * CellL(3)

         R_RB(1,j) = R_RB(1,j) + vv * SlideGap
         sr(1) = -nint(R_RB(1,j) * InvCL(1)) * CellL(1)

         R_RB(:,j) = R_RB(:,j) + sr(:)

         V_RB (1,j) = V_RB (1,j) + sr(2) * ShearRate
         V_RBt(1,j) = V_RBt(1,j) + sr(2) * ShearRate

       end do

     end if

     do j = 1 , N

       vv    = -nint(R(2,j) * InvCL(2))
       sr(2) = vv * CellL(2)
       sr(3) = -nint(R(3,j) * InvCL(3)) * CellL(3)

       R(1,j) = R(1,j) + vv * SlideGap
       sr(1)  = -nint( R(1,j) * InvCL(1)) * CellL(1)

       R(:,j) = R(:,j) + sr(:)

       Vel (1,j) = Vel (1,j) + sr(2) * ShearRate
       Velt(1,j) = Velt(1,j) + sr(2) * ShearRate

     end do

! ###################
   else
! ###################

     if(QColloid) then   ! NOTE!! only for colloidal center of mass

       do j = 1 , NumColloid

         R_RB(:,j) = R_RB(:,j) - nint( R_RB(:,j) * InvCL(:) ) * CellL(:)

       end do

     end if

     do i = 1 , N

       R(:,i) = R(:,i) - nint( R(:,i) * InvCL(:) ) * CellL(:)

     end do

! ###################
   end if
! ###################

end subroutine PBC_DPD
