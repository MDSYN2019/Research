! ############################
! ## SUBROUTINE LIST 
! ## -- ErrorFuncList 
! ## -- ErrorFunction 
! ## -- RecLatticeList 
! ## -- Ewald_SelfTerm 
! ## -- Ewald_Reciprocal 
! ############################


!#####################################################################
!#####################################################################


! ****************************************************
! ** subroutine for making error function data list **
! ** EFList                                         **
! ****************************************************

subroutine ErrorFuncList

use EwaldParam, only : msh, EFList

implicit NONE

integer, parameter :: max = 5
integer :: i,EFdim
real(8) :: x, ErrorFunction
external ErrorFunction

   EFdim = max * msh

   allocate( EFList(0:Efdim) )

   do i = 0 , msh*max

     x = dble( i ) / msh
     EFList(i) = ErrorFunction(x)

   end do

end subroutine ErrorFuncList


!********************************************************************
!********************************************************************


! **************************************
! *     error function                 *
! **************************************

Function ErrorFunction(arg)

implicit NONE

real(8), parameter, dimension(5) :: a = (/0.254829592d0,   &
   &                                     -0.284496736d0,   &
   &                                      1.421413740d0,   &
   &                                     -1.453152027d0,   &
   &                                      1.061405429d0/)
real(8), parameter :: p = 0.3275911d0
real(8)            :: t, arg, poly, tn, ErrorFunction
integer :: i

   t    = 1.d0 / ( 1.d0 + p * arg )
   tn   = t
   poly = a(1) * tn

   do i = 2 , 5

     tn   = tn * t
     poly = poly + a(i) * tn

   end do

   ErrorFunction = poly * dexp( - arg * arg )

end Function ErrorFunction


!#####################################################################
!#####################################################################


! ***************************************************
! ** subroutine for making reciprocal lattice list **
! ** NOTE: Cell Shape is not taken into account !! **
! ***************************************************

subroutine RecLatticeList

use CommonBlocks, only : QMaster, QPathInt, Qstdout
use CommonMPI
use CommonPI
use EwaldParam, only : ih2mx, Nh, TmpNh, ih, kmaxx, kmaxy, kmaxz
use CellParam, only : H, InvH

implicit NONE

integer :: NAddition, Count
integer :: i, ihb, k, ix, iy, iz, ii, jj
real(8) :: a2, b2, c2, l2, hr2mx, hr2
real(8), dimension(3) :: a, b, c, hr
real(8), dimension(3,3) :: InvHt
integer, dimension(:,:), allocatable :: ih_temp
integer, dimension(:,:), allocatable :: ih_temporary
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   a(:) = H(:,1)
   b(:) = H(:,2)
   c(:) = H(:,3)
   a2 = dot_product( a , a )
   b2 = dot_product( b , b )
   c2 = dot_product( c , c )
   l2 = (a2*b2*c2)**(1./3.d0)

   ihb    = int( sqrt( dble( ih2mx ) ) ) / 2 * 3

   ii = ihb * 2 / 3 + 1
   jj = int(4.2d0 * ii**3) + 1
   allocate( ih_temp(3,jj) )

   hr2mx  = dble(ih2mx) / l2

   InvHt = Transpose(InvH)

   k = 0

   do ix =    0 , ihb

     do iy = -ihb , ihb

       do iz = -ihb , ihb

         hr(1) = InvHt(1,1)*ix + InvHt(1,2)*iy + InvHt(1,3)*iz
         hr(2) = InvHt(2,1)*ix + InvHt(2,2)*iy + InvHt(2,3)*iz
         hr(3) = InvHt(3,1)*ix + InvHt(3,2)*iy + InvHt(3,3)*iz

         hr2 = dot_product( hr , hr )

         if( hr2 > hr2mx ) cycle
         if( ( ix == 0 ) .and. ( iy <  0 ) ) cycle
         if( ( ix == 0 ) .and. ( iy == 0 ) .and. ( iz <= 0 ) ) cycle

         k = k + 1
         if(k==jj) then
           allocate( ih_temporary(3,jj) )
           do ii = 1, jj
             ih_temporary(:,ii) = ih_temp(:,ii)
           end do
           deallocate( ih_temp )
           allocate( ih_temp(3,jj*2) )
           do ii = 1, jj
             ih_temp(:,ii) = ih_temporary(:,ii)
           end do
           deallocate( ih_temporary )
           jj = jj*2
         end if

         ih_temp(1,k) = ix
         ih_temp(2,k) = iy
         ih_temp(3,k) = iz

       end do

     end do

   end do

   TmpNh = k

   do ix = 1, ihb

     hr2 = (ix * InvHt(1,1))**2

     if( hr2 > hr2mx ) cycle

     kmaxx = ix

   end do

   do iy = 1, ihb

     hr2 = (iy * InvHt(2,2))**2

     if( hr2 > hr2mx ) cycle

     kmaxy = iy

   end do

   do iz = 1, ihb

     hr2 = (iz * InvHt(3,3))**2

     if( hr2 > hr2mx ) cycle

     kmaxz = iz

   end do

   if(QMaster.and.Qstdout) then

     write(*,'(a,3i5)') 'kmax = ', kmaxx, kmaxy, kmaxz

   end if

   Nh        = TmpNh / NProcsTemp
   NAddition = mod(TmpNh,NProcsTemp)

   if(MyRankTemp < NAddition) Nh = Nh + 1

   allocate( ih(3,Nh) )

   Count = 0

   do i = MyRankTemp+1 , TmpNh , NProcsTemp

     Count = Count + 1
     ih(:,Count) = ih_temp(:,i)

   end do

   deallocate( ih_temp )


end subroutine RecLatticeList


!#####################################################################
!#####################################################################


! *****************************
! **  Ewald Self Energy term **
! *****************************

subroutine Ewald_SelfTerm

use Numbers, only : N
use NonbondParam, only : Charge
use UnitExParam, only : pi
use EwaldParam, only : Alpha, Ene_Eslf

implicit NONE

integer :: i

   Ene_Eslf = 0.d0
   do i = 1 , N
     Ene_Eslf = Ene_Eslf + Charge(i) * Charge(i)
   end do
   Ene_Eslf = - Ene_Eslf * Alpha / sqrt( pi )

end subroutine Ewald_SelfTerm


!#####################################################################
!#####################################################################


subroutine Ewald_Reciprocal

use Numbers, only : N
use CommonBlocks, only : QRigidBody
use Configuration, only : R
use CommonMPI
use RBparam, only : NumRB, RBType, NumRBAtom, QSingle, Rmolec
use UnitExParam, only : InvPi, sqpi, pi2
use EwaldParam, only : Nel, Nh, Nelist, ih, alp2, Frc_Eksp, Vir_Eksp, Ene_Eksp, PCh
use CellParam, only : InvH, Volume

implicit NONE

integer :: i, j, k, l

real(8) :: Wkx, Wky, Wkz
real(8) :: CsSum, SnSum, tyo, pref
real(8) :: et, et1, pa2, arkn2
real(8) :: cf, cf1, ep, er4
integer :: Nc, MyType, ii
real(8) :: zz
real(8) :: Saxx, Saxy, Saxz, Sayx, Sayy, Sayz, Sazx, Sazy, Sazz
real(8) :: kn2, Epkn, Trp, Invkn2, InvVOL
real(8) :: IHxx, IHxy, IHxz, IHyx, IHyy, IHyz, IHzx, IHzy, IHzz
real(8) :: knx, kny, knz, Fx, Fy, Fz, Rx, Ry, Rz
integer :: ihx, ihy, ihz
real(8), dimension(:), allocatable :: Csl, Snl
real(8), dimension(:,:), allocatable :: Rel, Fel
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz

   pref = -sqpi * alp2
   InvVOL = 1.d0 / Volume
   er4 = 4.d0 * InvVOL

   et  = InvPi * InvVOL
   pa2 = -sqpi * alp2 * 2.d0

   IHxx = InvH(1,1)
   IHxy = InvH(1,2)
   IHxz = InvH(1,3)
   IHyx = InvH(2,1)
   IHyy = InvH(2,2)
   IHyz = InvH(2,3)
   IHzx = InvH(3,1)
   IHzy = InvH(3,2)
   IHzz = InvH(3,3)

   allocate( Rel(3,Nel) )
   allocate( Fel(3,Nel) )
   allocate( Snl(Nel) )
   allocate( Csl(Nel) )

   do i = 1, Nel
     j = Nelist(i)
     Rel(1,i) = R(1,j)
     Rel(2,i) = R(2,j)
     Rel(3,i) = R(3,j)
   end do
   Fel = 0.d0

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

!--------------------------------------------------------------------------

   do l = 1 , Nh

     ihx = ih(1,l)
     ihy = ih(2,l)
     ihz = ih(3,l)

     knx = IHxx*ihx + IHyx*ihy + IHzx*ihz
     kny = IHxy*ihx + IHyy*ihy + IHzy*ihz
     knz = IHxz*ihx + IHyz*ihy + IHzz*ihz
     kn2 = knx*knx + kny*kny + knz*knz
     Invkn2 = 1.d0 / kn2
     Epkn  = exp( pref * kn2 ) * Invkn2

     Wkx = knx * pi2
     Wky = kny * pi2
     Wkz = knz * pi2
     CsSum = 0.d0
     SnSum = 0.d0

     do i = 1, Nel

       zz = Pch(i)

       tyo = Wkx*Rel(1,i) + Wky*Rel(2,i) + Wkz*Rel(3,i)

       Csl(i) = zz * cos(tyo)
       Snl(i) = zz * sin(tyo)

       CsSum = CsSum + Csl(i)
       SnSum = SnSum + Snl(i)

     end do

     Trp = CsSum * CsSum + SnSum * SnSum

     ep = er4 * Epkn

     do i = 1 , Nel

       cf1 = Snl(i) * CsSum - Csl(i) * SnSum
       cf  = cf1 * ep

       Fel(1,i) = Fel(1,i) + cf * knx
       Fel(2,i) = Fel(2,i) + cf * kny
       Fel(3,i) = Fel(3,i) + cf * knz

     end do

     et1   =  et * Epkn * Trp
     arkn2 = ( -2.d0 * Invkn2 + pa2 ) * et1

     Vxx = Vxx + et1 + arkn2 * knx * knx
     Vxy = Vxy       + arkn2 * knx * kny
     Vxz = Vxz       + arkn2 * knx * knz
     Vyy = Vyy + et1 + arkn2 * kny * kny
     Vyz = Vyz       + arkn2 * kny * knz
     Vzz = Vzz + et1 + arkn2 * knz * knz

     Ene_Eksp = Ene_Eksp + Epkn * Trp

   end do

   Vir_Eksp(1,1) = Vxx
   Vir_Eksp(1,2) = Vxy
   Vir_Eksp(1,3) = Vxz
   Vir_Eksp(2,1) = Vxy
   Vir_Eksp(2,2) = Vyy
   Vir_Eksp(2,3) = Vyz
   Vir_Eksp(3,1) = Vxz
   Vir_Eksp(3,2) = Vyz
   Vir_Eksp(3,3) = Vzz

   do i = 1, Nel
     j = Nelist(i)
     Frc_Eksp(1,j) = Fel(1,i)
     Frc_Eksp(2,j) = Fel(2,i)
     Frc_Eksp(3,j) = Fel(3,i)
   end do

   if(QRigidBody) then

     k = 0
     Saxx = 0.d0
     Saxy = 0.d0
     Saxz = 0.d0
     Sayx = 0.d0
     Sayy = 0.d0
     Sayz = 0.d0
     Sazx = 0.d0
     Sazy = 0.d0
     Sazz = 0.d0

     do i = 1 , NumRB

       if(QSingle(i)) then

         k = k + 1

       else

         MyType = RBType(i)
         Nc     = NumRBAtom(MyType)

         do j = 1 , Nc

           ii = k + j

           Fx = Frc_Eksp(1,ii)
           Fy = Frc_Eksp(2,ii)
           Fz = Frc_Eksp(3,ii)
           Rx = Rmolec(1,j,i)
           Ry = Rmolec(2,j,i)
           Rz = Rmolec(3,j,i)

           Saxx = Saxx + Fx * Rx
           Saxy = Saxy + Fx * Ry
           Saxz = Saxz + Fx * Rz
           Sayx = Sayx + Fy * Rx
           Sayy = Sayy + Fy * Ry
           Sayz = Sayz + Fy * Rz
           Sazx = Sazx + Fz * Rx
           Sazy = Sazy + Fz * Ry
           Sazz = Sazz + Fz * Rz

         end do

         k = k + Nc

       end if

     end do

     Vir_Eksp(1,1) = Vir_Eksp(1,1) - Saxx
     Vir_Eksp(1,2) = Vir_Eksp(1,2) - Saxy
     Vir_Eksp(1,3) = Vir_Eksp(1,3) - Saxz
     Vir_Eksp(2,1) = Vir_Eksp(2,1) - Sayx
     Vir_Eksp(2,2) = Vir_Eksp(2,2) - Sayy
     Vir_Eksp(2,3) = Vir_Eksp(2,3) - Sayz
     Vir_Eksp(3,1) = Vir_Eksp(3,1) - Sazx
     Vir_Eksp(3,2) = Vir_Eksp(3,2) - Sazy
     Vir_Eksp(3,3) = Vir_Eksp(3,3) - Sazz

   end if


!----------------------------------------------------------------------

   Ene_Eksp = Ene_Eksp * InvVOL * InvPi

!----------------------------------------------------------------------
   deallocate( Rel, Fel, Snl, Csl )

end subroutine Ewald_Reciprocal
