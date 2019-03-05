

!######################################################################
!######################################################################


subroutine ContInertia

use ParamInertia
use OptConstraintParam, only : Ene_OptC
use CommonMPI, only : NProcs
use Configuration, only : R
use Numbers, only : N

implicit none

real(8), dimension(3,3) :: RMat, Rot
real(8), dimension(3) :: Evalue, v1, v2, v3
integer :: ie3, i, ii
integer, dimension(2) :: neg
real(8) :: EneInert, dd
real(8) :: e1, e2, e3, xx, pref0, pref1, pref3
real(8), dimension(3) :: Rg

   call CalcInertia(RMat,Rg)

   do i = 1, N
     R(:,i) = R(:,i) - Rg(:)
   end do

   call Jacobi(RMat,Rot,Evalue)
   call Sort3(Evalue,ie3)

   e3 = Evalue(ie3)
   v3(:) = Rot(:,ie3)
   ii = 0
   do i = 1, 3
     if(i==ie3) cycle
     ii = ii + 1
     neg(ii) = i
   end do

   e1 = Evalue(neg(1))
   e2 = Evalue(neg(2))
   v1(:) = Rot(:,neg(1))
   v2(:) = Rot(:,neg(2))

   xx = 1.d0 / (e1 + e2)

   Dop = 1.d0 - 2.d0 * e3 * xx  ! order parameter
   dd = Dop - TargetD
   pref0 = k_D * dd

   EneInert = 0.5d0 * pref0 * dd

   Ene_OptC = Ene_OptC + EneInert / dble(NProcs)

! ## e3 term
   pref3 = - 2.d0 * xx * pref0
   call Calc_DLDR(v3,pref3)

   pref1 = 2.d0 * e3 * xx * xx * pref0
   call Calc_DLDR(v1,pref1)
   call Calc_DLDR(v2,pref1)

end subroutine ContInertia


!######################################################################
!######################################################################


subroutine CalcInertia(RMat,Rg)

use AtomParam, only : Mass
use Configuration, only : R
use ParamInertia

implicit none

integer :: i, j
real(8), dimension(3,3) :: RMat
real(8), dimension(3) :: Rg
real(8) :: SumM, xmass, x, y, z

   Rg   = 0.d0
   SumM = 0.d0

   do j = 1, NumInert

     i = ListInert(j)
     Rsel(:,j) = R(:,i)
     Msel(j) = Mass(i)
     Rg   = Rg   + Mass(i) * R(:,i)
     SumM = SumM + Mass(i)

   end do

   Rg = Rg / SumM

   do i = 1, NumInert
     Rsel(:,i) = Rsel(:,i) - Rg(:)
   end do


   RMat(:,:) = 0.d0

   do i = 1, NumInert
     xmass = Msel(i)
     x = Rsel(1,i)
     y = Rsel(2,i)
     z = Rsel(3,i)
     RMat(1,1) = RMat(1,1) + xmass * ( y * y + z * z )
     RMat(2,2) = RMat(2,2) + xmass * ( z * z + x * x )
     RMat(3,3) = RMat(3,3) + xmass * ( x * x + y * y )
     RMat(1,2) = RMat(1,2) + xmass *   x * y
     RMat(2,3) = RMat(2,3) + xmass *   y * z
     RMat(3,1) = RMat(3,1) + xmass *   z * x
   end do
   RMat(2,1) = RMat(1,2)
   RMat(3,2) = RMat(2,3)
   RMat(1,3) = RMat(3,1)

end subroutine CalcInertia


!######################################################################
!######################################################################


subroutine Sort3(Evalue,ie3)

implicit none

real(8), dimension(3) :: Evalue
integer :: ie3
real(8) :: d12, d23, d13

   d12 = abs( Evalue(1) - Evalue(2) )
   d23 = abs( Evalue(2) - Evalue(3) )
   d13 = abs( Evalue(1) - Evalue(3) )

   if(d12 < d23) then
     if(d13 < d12) then
       ie3 = 2
     else
       ie3 = 3
     end if
   else if(d13 < d23) then
     ie3 = 2
   else
     ie3 = 1
   end if

end subroutine Sort3


!######################################################################
!######################################################################


subroutine Calc_DLDR(EigenVec,prefinert)

use ParamInertia
use CommonMPI
use OptConstraintParam, only : Frc_OptC

implicit none

integer :: i, Nas
real(8) :: xi, yi, zi, xmass
real(8) :: dX1,dX2,dX3
real(8) :: dY1,dY2,dY3
real(8) :: dZ1,dZ2,dZ3
real(8), dimension(3) :: EigenVec
real(8) :: Ux, Uy, Uz
real(8) :: dLdX, dLdY, dLdZ
real(8) :: prefinert


   Nas = NProcs - MyRank

   Ux = EigenVec(1)
   Uy = EigenVec(2)
   Uz = EigenVec(3)

   do i = Nas, NumInert, NProcs
     xmass = prefinert * Msel(i)
     xi = Rsel(1,i)
     yi = Rsel(2,i)
     zi = Rsel(3,i)

     dX1 = - yi * Uy - zi * Uz
     dX2 = - yi * Ux + 2.d0 * xi * Uy
     dX3 = - zi * Ux + 2.d0 * xi * Uz
     dLdX= Ux * dX1 + Uy * dX2 + Uz * dX3

     dY1 = 2.d0 * yi * Ux - xi * Uy
     dY2 = - xi * Ux - zi * Uz
     dY3 = - zi * Uy + 2.d0 * yi * Uz
     dLdY = Ux * dY1 + Uy * dY2 + Uz * dY3

     dZ1 = 2.d0 * zi * Ux - xi * Uz
     dZ2 = 2.d0 * zi * Uy - yi * Uz
     dZ3 = - xi * Ux - yi * Uy
     dLdZ = Ux * dZ1 + Uy * dZ2 + Uz * dZ3

     dLdX = xmass * dLdX
     dLdY = xmass * dLdY
     dLdZ = xmass * dLdZ

     Frc_OptC(1,i) = Frc_OptC(1,i) + dLdX
     Frc_OptC(2,i) = Frc_OptC(2,i) + dLdY
     Frc_OptC(3,i) = Frc_OptC(3,i) + dLdZ
   end do

end subroutine Calc_DLDR


!######################################################################
!######################################################################


subroutine CntIn_Sample(step)

use TimeParam, only : lk, Timeps
use ParamInertia
use IOparam, only : DirectoryName
use CommonBlocks, only : Job_name

integer :: step
real(8) :: dd, ave

   if(step==lk) then
     open(27,file=trim(DirectoryName)//trim(adjustl(Job_name))//'PMF_INERT.dat',status='unknown')
     DPMF = 0.d0
   end if

   dd = Dop - TargetD
   DPMF = DPMF + k_D * dd

   ave = DPMF / dble(step/lk)

   write(27,'(f12.4,2e16.8)') Timeps, ave, Dop

end subroutine CntIn_Sample
