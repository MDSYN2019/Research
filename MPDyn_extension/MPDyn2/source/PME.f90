! ############################
! ## SUBROUTINE LIST 
! ## -- PMEwald_pre
! ## -- PME_GetSize
! ## -- PME_setup
! ## -- PME_Reciprocal
! ## -- Scaled_coord 
! ## -- Load_Bspline_moduli
! ## -- dFTModulus
! ## -- Fill_charge_grid 
! ## -- PME_energy
! ## -- PME_force
! ## -- Get_Bspline_coeffs
! ## -- fill_bspline
! ## -- init
! ## -- one_pass
! ## -- diff
! ## -- Get_FFTdimension
! ## -- FFTsetup
! ## -- FFT_forward
! ## -- FFT_back
! ## -- pubz3di
! ## -- sgifft3di
! ## -- fftw_ini
! ## -- pubz3d
! ## -- sgifft3d
! ## -- fftw3d 
! ############################


!######################################################################
!######################################################################


subroutine PMEwald_pre

use EwaldParam, only : Nel
use PMEparam

implicit none

   allocate( ScRs(3,Nel) )
   allocate( BthetaX(Bsp_order,Nel) )
   allocate( BthetaY(Bsp_order,Nel) )
   allocate( BthetaZ(Bsp_order,Nel) )
   allocate( dBthetaX(Bsp_order,Nel) )
   allocate( dBthetaY(Bsp_order,Nel) )
   allocate( dBthetaZ(Bsp_order,Nel) )
   allocate( BsplineModuleX(Nfft(1)) )
   allocate( BsplineModuleY(Nfft(2)) )
   allocate( BsplineModuleZ(Nfft(3)) )

   call PME_GetSize
   call PME_setup

end subroutine PMEwald_pre


!######################################################################
!######################################################################


subroutine PME_GetSize

use EwaldParam, only : Nel
use CommonBlocks, only : QMaster, QPathInt
use PMEparam, only : Nfft, Nfftdim, Bsp_order, SizeBtheta, SizeGridQ, MaxGrid
use CommonMPI, only : NProcs
use CommonPI, only : NumProcess

implicit none

integer :: Sizeheap,Sizestack
integer :: SizeFFTtable, SizeFFTwork
integer :: ii, jj
integer :: NProcsTemp

   if(QPathInt) then
     NProcsTemp = NumProcess
   else
     NProcsTemp = NProcs
   end if

   call Get_FFTdimension(SizeFFTtable, SizeFFTwork)

   SizeBtheta = Nel * Bsp_order
   SizeGridQ  = 2 * Nfftdim(1) * Nfftdim(2) * Nfftdim(3)
   Sizeheap   = Nfft(1) + Nfft(2) + Nfft(3) + SizeFFTtable
   Sizestack  = SizeGridQ + 6 * SizeBtheta + SizeFFTwork + 3 * Nel

   if(NProcsTemp /= 1) then
     if(Nfftdim(1) >= Nfftdim(3)) then
       ii = Nfftdim(1) / NProcsTemp
       jj = mod( Nfftdim(1), NProcsTemp )
       if(jj /= 0) ii = ii + 1
       MaxGrid = ii * Nfftdim(2) * Nfftdim(3) * 2
     else
       ii = Nfftdim(3) / NProcsTemp
       jj = mod( Nfftdim(3), NProcsTemp )
       if(jj /= 0) ii = ii + 1
       MaxGrid = ii * Nfftdim(2) * Nfftdim(1) * 2
     end if

   else

     MaxGrid = SizeGridQ

   end if

end subroutine PME_GetSize


!######################################################################
!######################################################################


subroutine PME_setup

implicit none

   call Load_Bspline_moduli
   call FFTsetup
   call Prep_Atom_to_Mesh ! parallel 

end subroutine PME_setup


!######################################################################
!######################################################################


subroutine PME_Reciprocal

use CellParam, only : InvH
use PMEparam, only : SizeGridQ

implicit none

real(8), dimension(3,3) :: InvHt
real(8), dimension(SizeGridQ) :: gridQ

   InvHt = Transpose(InvH)

   call Scaled_coord

   call Get_Bspline_coeffs

   call Fill_charge_grid(gridQ)

   call FFT_back(gridQ)

   call PME_energy(gridQ,InvHt)

   call FFT_forward(gridQ)

   call PME_force(gridQ,InvHt)

end subroutine PME_Reciprocal


!######################################################################
!######################################################################


subroutine Scaled_coord

use EwaldParam, only : Nel, Nelist
use CommonBlocks, only : QPathInt
use CellParam, only : InvH
use Configuration, only : R
use PMEparam, only : ScRs, Nfft
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess

implicit none

! INPUT:
!      N: number of atoms
!      R: arrays of cartesian coords
!      InvH: Inverse of the cell matrix
! OUTPUT:
!     ScRs the scaled and shifted fractional coords

integer :: i, j, Nas
integer :: NProcsTemp, MyRankTemp
real(8) :: Sx, Sy, Sz, Rx, Ry, Rz
real(8) :: IHxx, IHxy, IHxz, IHyx, IHyy, IHyz, IHzx, IHzy, IHzz

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   IHxx = InvH(1,1)
   IHxy = InvH(1,2)
   IHxz = InvH(1,3)
   IHyx = InvH(2,1)
   IHyy = InvH(2,2)
   IHyz = InvH(2,3)
   IHzx = InvH(3,1)
   IHzy = InvH(3,2)
   IHzz = InvH(3,3)

   do i = Nas, Nel, NProcsTemp
     j  = Nelist(i)
     Rx = R(1,j)
     Ry = R(2,j)
     Rz = R(3,j)
     Sx = IHxx*Rx + IHxy*Ry + IHxz*Rz
     Sy = IHyx*Rx + IHyy*Ry + IHyz*Rz
     Sz = IHzx*Rx + IHzy*Ry + IHzz*Rz
     if(Sx >  0.5d0) Sx = Sx - 1.0d0
     if(Sx < -0.5d0) Sx = Sx + 1.0d0
     if(Sy >  0.5d0) Sy = Sy - 1.0d0
     if(Sy < -0.5d0) Sy = Sy + 1.0d0
     if(Sz >  0.5d0) Sz = Sz - 1.0d0
     if(Sz < -0.5d0) Sz = Sz + 1.0d0
     Sx = Sx + 0.5d0
     Sy = Sy + 0.5d0
     Sz = Sz + 0.5d0
     ScRs(1,i) = Nfft(1) * Sx
     ScRs(2,i) = Nfft(2) * Sy
     ScRs(3,i) = Nfft(3) * Sz
   end do

end subroutine Scaled_coord


!######################################################################
!######################################################################


subroutine Load_Bspline_moduli

use CommonBlocks, only : QMaster
use PMEparam, only : Nfft, Bsp_order, BsplineModuleX, BsplineModuleY, BsplineModuleZ

implicit none

integer, parameter :: MAXORDER=25
integer, parameter :: MAXNFFT=1000

real(8), dimension(MAXORDER) :: array,darray
real(8) :: w
real(8), dimension(MAXNFFT) :: bsp_arr
integer :: i,maxn

   if ( Bsp_order > MAXORDER )then
     if(QMaster) then
       write(6,*) 'order too large! check on MAXORDER'
     end if
     call Finalize
   endif

   maxn = max(Nfft(1),Nfft(2),Nfft(3))

   if ( maxn > MAXNFFT )then 
     if(QMaster) then
       write(6,*)'nfft1-3 too large! check on MAXNFFT'
     end if
     call Finalize
   endif

   w = 0.d0

   call fill_bspline(w,Bsp_order,array,darray)

   do i = 1,maxn
     bsp_arr(i) = 0.d0
   end do

   do i = 2, (Bsp_order + 1)
     bsp_arr(i) = array(i-1)
   end do

   call dFTModulus(BsplineModuleX,bsp_arr,Nfft(1))
   call dFTModulus(BsplineModuleY,bsp_arr,Nfft(2))
   call dFTModulus(BsplineModuleZ,bsp_arr,Nfft(3))

end subroutine Load_Bspline_moduli


!######################################################################
!######################################################################


subroutine dFTModulus(bsp_mod,bsp_arr,Numfft)

use UnitExParam, only : pi2

implicit none

integer :: Numfft
real(8), dimension(Numfft) :: bsp_mod,bsp_arr

! Computes the modulus of the discrete fourier transform of bsp_arr,
!  storing it into bsp_mod 

integer :: i,j
real(8) :: cst,snt,arg,tiny

   tiny = 1.d-7

   do i = 1,Numfft

     cst = 0.d0
     snt = 0.d0

     do j = 1,Numfft

       arg = pi2 * (i-1) * (j-1) / Numfft
       cst = cst + bsp_arr(j)*cos(arg)
       snt = snt + bsp_arr(j)*sin(arg)

     end do

     bsp_mod(i) = cst*cst + snt*snt

   end do

   do i = 1 , Numfft

       if ( bsp_mod(i) < tiny ) then
         bsp_mod(i) = 0.5d0*(bsp_mod(i-1) + bsp_mod(i+1))
       end if

   end do

end subroutine dFTModulus


!######################################################################
!######################################################################


subroutine Fill_charge_grid(gridQ)

use EwaldParam, only : Nel, PCh
use CommonBlocks, only : QPathInt
use PMEparam, only : Bsp_order, Nfft, ScRs, BthetaX, BthetaY, &
                   & BthetaZ, NfftDim
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess

! INPUT:
!      N:  number of atoms
!      Charge: the array of atomic charges
!      BthetaX,BthetaY,BthetaZ: the spline coeff arrays
!      ScRs : the scaled and shifted fractional coords
!      Nfft: the charge grid dimensions
!      Nfftdim: physical charge grid dims
!      Bsp_order: the order of spline interpolation
! OUTPUT:
!      gridQ the charge grid

implicit none

integer :: l, ith1, ith2, ith3
integer :: i, j, k, ii, jj, kk, Nas
real(8), dimension(2,NfftDim(1),NfftDim(2),NfftDim(3)) :: gridQ
real(8), dimension(NfftDim(1),NfftDim(2),NfftDim(3)) :: Q
real(8) :: fct0, fct1, fct2
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   Q = 0.d0

   do l = Nas, Nel, NProcsTemp

     fct0 = PCh(l)

     kk = int(ScRs(3,l)) - Bsp_order

     do ith3 = 1 , Bsp_order

       kk = kk + 1
       k  = kk + 1 + (Nfft(3) - isign(Nfft(3),kk))/2

       fct1 = fct0 * BthetaZ(ith3,l)

       jj = int(ScRs(2,l)) - Bsp_order

       do ith2 = 1 , Bsp_order

         jj = jj + 1
         j  = jj + 1 + (Nfft(2) - isign(Nfft(2),jj))/2

         fct2 = BthetaY(ith2,l) * fct1

         ii = int(ScRs(1,l)) - Bsp_order

         do ith1 = 1 , Bsp_order

           ii = ii + 1
           i  = ii + 1 + (Nfft(1) - isign(Nfft(1),ii))/2

           Q(i,j,k) = Q(i,j,k) + BthetaX(ith1,l) * fct2

         end do

       end do

     end do

   end do

   if(NProcsTemp/=1) call SumChargeDens(Q)

   do i = 1, NfftDim(1)

     do j = 1, NfftDim(2)

       do k = Nas, NfftDim(3), NProcsTemp

         gridQ(1,i,j,k) = Q(i,j,k)
         gridQ(2,i,j,k) = 0.d0

       end do

     end do

   end do

!   print *, 'GRIDQ',gridQ(1,1,1,1),gridQ(2,1,1,1)

end subroutine Fill_charge_grid


!######################################################################
!######################################################################


subroutine PME_energy(gridQ)


use CommonBlocks, only : QPathInt, Qdebug, QMaster
use CellParam, only : Volume, InvH
use EwaldParam, only : alp2, Ene_Eksp, Vir_Eksp
use PMEparam, only : Nfft, BsplineModuleX, BsplineModuleY, BsplineModuleZ, &
                   & Nfftdim
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess
use UnitExParam, only : pi, sqpi

implicit none


real(8) :: pref, denom, eterm, vterm, energy
integer :: k1, k2, k3, Nff, ini
integer :: ihx, ihy, ihz
integer :: Nfx, Nfy, Nfz, Nas
real(8) :: struc2
real(8), dimension(2,NfftDim(1),NfftDim(2),NfftDim(3)) :: gridQ
real(8) :: kn2, Bspx, Bspy
integer :: NProcsTemp, MyRankTemp
real(8) :: IHxx, IHxy, IHxz, IHyx, IHyy, IHyz, IHzx, IHzy, IHzz
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
real(8) :: knx, kny, knz, est2, vtmx, vtmy, vtmz

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   IHxx = InvH(1,1)
   IHxy = InvH(1,2)
   IHxz = InvH(1,3)
   IHyx = InvH(2,1)
   IHyy = InvH(2,2)
   IHyz = InvH(2,3)
   IHzx = InvH(3,1)
   IHzy = InvH(3,2)
   IHzz = InvH(3,3)

   Nff = Nfft(1) * Nfft(2)
   pref = - sqpi * alp2

   Nfx = Nfft(1) / 2
   if ( 2*Nfx < Nfft(1) ) Nfx = Nfx + 1

   Nfy = Nfft(2) / 2
   if ( 2*Nfy < Nfft(2) ) Nfy = Nfy + 1

   Nfz = Nfft(3) / 2
   if ( 2*Nfz < Nfft(3) ) Nfz = Nfz + 1

   energy = 0.d0
   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do k1 = Nas, Nfft(1), NProcsTemp

     Bspx = BsplineModuleX(k1)

     do k2 = 1, Nfft(2)

       Bspy = BsplineModuleY(k2)
       ini = 1
       if(k1==1.and.k2==1) ini = 2

       do k3 = ini, Nfft(3)

         ihx = k1 - 1
         if ( k1 > Nfx ) ihx = k1 - 1 - Nfft(1)

         ihy = k2 - 1
         if ( k2 > Nfy ) ihy = k2 - 1 - Nfft(2)

         ihz = k3 - 1
         if ( k3 > Nfz ) ihz = k3 - 1 - Nfft(3)

         knx = IHxx*ihx + IHyx*ihy + IHzx*ihz
         kny = IHxy*ihx + IHyy*ihy + IHzy*ihz
         knz = IHxz*ihx + IHyz*ihy + IHzz*ihz
         kn2 = knx*knx + kny*kny + knz*knz

         denom = pi * Volume * Bspx * Bspy * BsplineModuleZ(k3) * kn2
         eterm = exp( pref * kn2 ) / denom
         vterm = 2.d0 * ( -pref * kn2 + 1.d0 ) / kn2

         struc2 = gridQ(1,k1,k2,k3) * gridQ(1,k1,k2,k3) + gridQ(2,k1,k2,k3) * gridQ(2,k1,k2,k3)

         gridQ(1,k1,k2,k3) = eterm * gridQ(1,k1,k2,k3)
         gridQ(2,k1,k2,k3) = eterm * gridQ(2,k1,k2,k3)

         energy = energy + eterm * struc2

         est2 = eterm * struc2
         vtmx = vterm * knx
         vtmy = vterm * kny
         vtmz = vterm * knz

         Vxx = Vxx + est2 * (vtmx * knx - 1.d0)
         Vxy = Vxy + est2 * vtmx * kny
         Vxz = Vxz + est2 * vtmx * knz
         Vyy = Vyy + est2 * (vtmy * kny - 1.d0)
         Vyz = Vyz + est2 * vtmy * knz
         Vzz = Vzz + est2 * (vtmz * knz - 1.d0)

       end do
     end do
   end do

   Ene_Eksp = 0.5d0 * energy
   if(Qdebug) then
     print *, Ene_Eksp
   end if

   Vxx = - 0.50 * Vxx
   Vxy = - 0.50 * Vxy
   Vxz = - 0.50 * Vxz
   Vyy = - 0.50 * Vyy
   Vyz = - 0.50 * Vyz
   Vzz = - 0.50 * Vzz

   Vir_Eksp(1,1) = Vxx
   Vir_Eksp(1,2) = Vxy
   Vir_Eksp(1,3) = Vxz
   Vir_Eksp(2,1) = Vxy
   Vir_Eksp(2,2) = Vyy
   Vir_Eksp(2,3) = Vyz
   Vir_Eksp(3,1) = Vxz
   Vir_Eksp(3,2) = Vyz
   Vir_Eksp(3,3) = Vzz

end subroutine PME_energy


!######################################################################
!######################################################################


subroutine PME_force(gridQ)

use Numbers, only : N
use EwaldParam, only : Nel, Nelist, PCh, Frc_Eksp, Vir_Eksp
use CommonBlocks, only : QRigidBody, QPathInt
use PMEparam, only : Nfft, Nfftdim, ScRs, Bsp_order, &
&   BthetaX, BthetaY, BthetaZ, dBthetaX, dBthetaY, dBthetaZ
use RBparam, only : NumRB, QSingle, RBType, Rmolec, NumRBAtom
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess
use CellParam, only : InvH

implicit none

real(8), dimension(2,NfftDim(1),NfftDim(2),NfftDim(3)) :: gridQ
real(8), dimension(NfftDim(1),NfftDim(2),NfftDim(3)) :: Q

integer :: l,ll,ith1,ith2,ith3,ii,jj,kk,i,j,k, Nas
real(8) :: term, zz
real(8) :: dBx, dBy, dBz, Bx, By, Bz
real(8) :: fcx, fcy, fcz, Rx, Ry, Rz, Fx, Fy, Fz
real(8) :: IHxx, IHxy, IHxz, IHyx, IHyy, IHyz, IHzx, IHzy, IHzz
real(8) :: Saxx, Saxy, Saxz, Sayx, Sayy, Sayz, Sazx, Sazy, Sazz
real(8), dimension(3,N) :: RinMOL
integer :: Nc, MyType
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   IHxx = InvH(1,1)
   IHxy = InvH(1,2)
   IHxz = InvH(1,3)
   IHyx = InvH(2,1)
   IHyy = InvH(2,2)
   IHyz = InvH(2,3)
   IHzx = InvH(3,1)
   IHzy = InvH(3,2)
   IHzz = InvH(3,3)

   do k = Nas, NfftDim(3), NProcsTemp

     do j = 1, NfftDim(2)

       do i = 1, NfftDim(1)

         Q(i,j,k) = gridQ(1,i,j,k)

       end do

     end do

   end do

   if(NProcsTemp/=1) call DistChargeDens(Q)

   do l = Nas, Nel, NProcsTemp

     zz = PCh(l)
     ll = Nelist(l)

     fcx = 0.d0
     fcy = 0.d0
     fcz = 0.d0

     kk = int(ScRs(3,l)) - Bsp_order

     do ith3 = 1,Bsp_order

       kk = kk + 1
       k = kk + 1 + (Nfft(3) - isign(Nfft(3),kk))/2
       jj = int(ScRs(2,l)) - Bsp_order
       dBz = dBthetaZ(ith3,l)
       Bz  =  BthetaZ(ith3,l)

       do ith2 = 1,Bsp_order

         jj = jj + 1
         j = jj + 1 + (Nfft(2) - isign(Nfft(2),jj))/2
         ii = int(ScRs(1,l)) - Bsp_order
         dBy = dBthetaY(ith2,l)
         By  =  BthetaY(ith2,l)

         do ith1 = 1,Bsp_order

           ii = ii + 1
           i = ii + 1 + (Nfft(1) - isign(Nfft(1),ii))/2
           term = zz * Q(i,j,k)
           dBx = dBthetaX(ith1,l)
           Bx  =  BthetaX(ith1,l)

           fcx = fcx - Nfft(1) * term * dBx *  By *  Bz
           fcy = fcy - Nfft(2) * term *  Bx * dBy *  Bz
           fcz = fcz - Nfft(3) * term *  Bx *  By * dBz

         end do

       end do

     end do

     Frc_Eksp(1,ll) = Frc_Eksp(1,ll) + IHxx*fcx + IHyx*fcy + IHzx*fcz
     Frc_Eksp(2,ll) = Frc_Eksp(2,ll) + IHxy*fcx + IHyy*fcy + IHzy*fcz
     Frc_Eksp(3,ll) = Frc_Eksp(3,ll) + IHxz*fcx + IHyz*fcy + IHzz*fcz

   end do

! ## Virial correction for Rigid Bodies 

   if(QRigidBody) then

     k = 0

     do i = 1 , NumRB

       if(QSingle(i)) then

         k = k + 1
         RinMOL(:,k) = 0.d0

       else

         MyType = RBType(i)
         Nc     = NumRBAtom(MyType)
         do j = 1 , Nc
           ii = k + j
           RinMOL(:,ii) = Rmolec(:,j,i)
         end do
         k = k + Nc

       end if

     end do

     Saxx = 0.d0
     Saxy = 0.d0
     Saxz = 0.d0
     Sayx = 0.d0
     Sayy = 0.d0
     Sayz = 0.d0
     Sazx = 0.d0
     Sazy = 0.d0
     Sazz = 0.d0

     do l = Nas, N, NProcsTemp

       Fx = Frc_Eksp(1,l)
       Fy = Frc_Eksp(2,l)
       Fz = Frc_Eksp(3,l)
       Rx = RinMOL(1,l)
       Ry = RinMOL(2,l)
       Rz = RinMOL(3,l)
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


end subroutine PME_force


!######################################################################
!######################################################################


subroutine Get_Bspline_coeffs

use EwaldParam, only : Nel
use CommonBlocks, only : QPathInt
use PMEparam, only : ScRs, Bsp_order, BthetaX, BthetaY, BthetaZ, &
&                     dBthetaX, dBthetaY, dBthetaZ
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess

!---------------------------------------------------------------------
! INPUT:
!      N: number of atoms
!      ScRs: the scaled and shifted fractional coords
!      Bsp_order: the order of spline interpolation
! OUTPUT
!      BthetaX,BthetaY,BthetaZ: the spline coeff arrays
!      dBthetaX,dBthetaY,dBthetaZ: the 1st deriv of spline coeff arrays
!---------------------------------------------------------------------

implicit none

real(8) :: x, y, z
integer :: i, Nas
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   do i = Nas, Nel, NProcsTemp

     x = ScRs(1,i) - int(ScRs(1,i))
     y = ScRs(2,i) - int(ScRs(2,i))
     z = ScRs(3,i) - int(ScRs(3,i))

     call fill_bspline(x,Bsp_order,BthetaX(1,i),dBthetaX(1,i))
     call fill_bspline(y,Bsp_order,BthetaY(1,i),dBthetaY(1,i))
     call fill_bspline(z,Bsp_order,BthetaZ(1,i),dBthetaZ(1,i))

   end do


end subroutine Get_Bspline_coeffs


!######################################################################
!######################################################################


subroutine fill_bspline(w,order,array,darray)

!---------- use standard B-spline recursions: see doc file

implicit none
integer :: order
real(8) :: w
real(8), dimension(order) :: array,darray

integer :: i

! do linear case

   call init(array,w,order)

! compute standard b-spline recursion

   do i = 3,order-1

     call one_pass(array,w,i)

   end do

! perform standard b-spline differentiation

   call diff(array,darray,order)

! one more recursion

   call one_pass(array,w,order)


end subroutine fill_bspline


!######################################################################
!######################################################################


subroutine init(c,x,order)

implicit none

integer :: order
real(8), dimension(order) :: c
real(8) :: x

   c(order) = 0.d0
   c(2) = x
   c(1) = 1.d0 - x

end subroutine init


!######################################################################
!######################################################################


subroutine one_pass(c,x,i)

implicit none

real(8), dimension(*) :: c
real(8) :: x, div

integer :: i, j

   div  = 1.d0 / (i-1)
   c(i) = div * x * c(i-1)

   do j = 1, i-2
     c(i-j) = div * ( (x+j) * c(i-j-1) + (i-j-x) * c(i-j) )
   end do

   c(1) = div * ( 1.d0 - x ) * c(1)

end subroutine one_pass


!######################################################################
!######################################################################


subroutine diff(c,d,order)

implicit none

real(8), dimension(*) :: c,d

integer :: order, j

   d(1) = -c(1)

   do j = 2, order

     d(j) = c(j-1) - c(j)

   end do

end subroutine diff


!######################################################################
!######################################################################


subroutine Get_FFTdimension(SizeFFTtable, SizeFFTwork)

use PMEparam

implicit none

integer :: NfftMax
integer, dimension(3) :: NLS
integer :: SizeFFTtable, SizeFFTwork
#ifdef FFTW
integer :: n1, n2, n3
#endif

   NfftMax = max( Nfft(1) , Nfft(2) , Nfft(3) )

   Nfftdim(:) = Nfft(:)

   NLS(:) = Nfft(:) / 2

   if( Nfft(1) == 2 * NLS(1) ) Nfftdim(1) = Nfft(1) + 1
   if( Nfft(2) == 2 * NLS(2) ) Nfftdim(2) = Nfft(2) + 1
   if( Nfft(3) == 2 * NLS(3) ) Nfftdim(3) = Nfft(3) + 1

#ifdef SGIFFT
!   NffTable    = 2*(Nfftdim(1)+Nfftdim(2)+Nfftdim(3)+50)
   NffTable    = 2*(NfftMax+15)
   NffWork     = 0
   SizeFFTtable = NffTable
   SizeFFTwork  = NffWork
#endif
#ifdef PUBFFT
   NffTable     = 4 * NfftMax + 15
   NffWork      = NfftMax
   SizeFFTtable = 3 * NffTable
   SizeFFTwork  = 2 * NfftMax
#endif

#ifdef FFTW
   n1 = Nfft(1)
   n2 = Nfft(2)
   n3 = Nfft(3)
   allocate( Forward_in_x(n1) )
   allocate( Forward_out_x(n1) )
   allocate( Backward_in_x(n1) )
   allocate( Backward_out_x(n1) )
   allocate( Forward_in_y(n2) )
   allocate( Forward_out_y(n2) )
   allocate( Backward_in_y(n2) )
   allocate( Backward_out_y(n2) )
   allocate( Forward_in_z(n3) )
   allocate( Forward_out_z(n3) )
   allocate( Backward_in_z(n3) )
   allocate( Backward_out_z(n3) )
   SizeFFTtable = 0
   SizeFFTwork  = 0
#else
   allocate( FFTtable(SizeFFTtable) )
   if(SizeFFTwork/=0) then
     allocate( FFTwork(SizeFFTwork) )
   else
     allocate( FFTwork(1) )
   end if
#endif

end subroutine Get_FFTdimension


!######################################################################
!######################################################################


subroutine FFTsetup

use PMEparam, only : Nfft, FFTtable, NffTable

implicit none

integer :: n1,n2,n3

   n1 = Nfft(1)
   n2 = Nfft(2)
   n3 = Nfft(3)

#ifdef SGIFFT
   call sgifft3di(n1,n2,n3,FFTtable,NffTable)
!      call ZFFT3DI(n1,n2,n3,FFTtable)
#endif
#ifdef PUBFFT
   call pubz3di(n1,n2,n3,FFTtable,NffTable)
#endif
#ifdef FFTW
   call fftw_ini(n1,n2,n3)
#endif

end subroutine FFTsetup


!######################################################################
!######################################################################


subroutine FFT_forward(array)

use PMEparam, only : Nfft, Nfftdim, FFTtable, FFTwork, NffTable, NffWork

implicit none

real(8), dimension(*) :: array
integer :: nnsign, n1, n2, n3, d1, d2

   n1 = Nfft(1)
   n2 = Nfft(2)
   n3 = Nfft(3)
   d1 = Nfftdim(1)
   d2 = Nfftdim(2)

   nnsign = 1

#ifdef SGIFFT
!      call ZFFT3D(nnsign,n1,n2,n3,array,d1,d2,FFTtable)
   call sgifft3d(nnsign,n1,n2,n3,array,d1,d2,FFTtable,NffTable,FFTwork,NffWork)
#endif
#ifdef PUBFFT
   call pubz3d(nnsign,n1,n2,n3,array,d1,d2,FFTtable,NffTable,FFTwork,NffWork)
#endif
#ifdef FFTW
   call fftw3d(nnsign,n1,n2,n3,array)
#endif

end subroutine FFT_forward


!######################################################################
!######################################################################


subroutine FFT_back(array)

use PMEparam, only : Nfft, Nfftdim, FFTtable, FFTwork, NffTable, NffWork

implicit none

real(8), dimension(*) :: array

integer :: nnsign, n1, n2, n3, d1, d2

   n1 = Nfft(1)
   n2 = Nfft(2)
   n3 = Nfft(3)
   d1 = Nfftdim(1)
   d2 = Nfftdim(2)

   nnsign = -1

#ifdef SGIFFT
!      call ZFFT3D(nnsign,n1,n2,n3,array,d1,d2,FFTtable)
   call sgifft3d(nnsign,n1,n2,n3,array,d1,d2,FFTtable,NffTable,FFTwork,NffWork)
#endif
#ifdef PUBFFT
   call pubz3d(nnsign,n1,n2,n3,array,d1,d2,FFTtable,NffTable,FFTwork,NffWork)
#endif
#ifdef FFTW
   call fftw3d(nnsign,n1,n2,n3,array)
#endif

end subroutine FFT_back


#ifdef PUBFFT


!######################################################################
!######################################################################


!*****************************************************************************
!*
!*	3D (slow) Fourier Transform
!*   this 1d->3d code is brute force approach
!*   the 1d code is a double precision version of fftpack from netlib
!*   due to Paul N Swartztrauber at NCAR Boulder Coloraso
!*
!*****************************************************************************

subroutine pubz3di(n1,n2,n3,table,ntable)

implicit none

integer :: n1,n2,n3,ntable
real(8), dimension(ntable,3) :: table

! ntable should be 4*max(n1,n2,n3) +15

   call cffti(n1,table(1,1))
   call cffti(n2,table(1,2))
   call cffti(n3,table(1,3))

end subroutine pubz3di


#endif


#ifdef SGIFFT


!######################################################################
!######################################################################


!*****************************************************************************
!*
!*	3D (slow) Fourier Transform
!*   this 1d->3d code is brute force approach
!*   the 1d code is a double precision version of fftpack from netlib
!*   due to Paul N Swartztrauber at NCAR Boulder Coloraso
!*
!*****************************************************************************

subroutine sgifft3di(n1,n2,n3,table,ntable)

implicit none

integer :: n1,n2,n3,ntable
real(8), dimension(ntable,3) :: table

! ntable should be 4*max(n1,n2,n3) +15

   call ZFFT1DI(n1,table(1,1))
   call ZFFT1DI(n2,table(1,2))
   call ZFFT1DI(n3,table(1,3))

end subroutine sgifft3di


#endif


#ifdef FFTW


!######################################################################
!######################################################################


subroutine fftw_ini(n1,n2,n3)

use PMEparam, only : Nfft, planFx, planFy, planFz, planBx, planBy, planBz, &
                   &       Forward_in_x, Forward_out_x, Forward_in_y, &
                   &       Forward_out_y, Forward_in_z, Forward_out_z,&
                   &       Backward_in_x, Backward_out_x, Backward_in_y, &
                   &       Backward_out_y, Backward_in_z, Backward_out_z

implicit none

include 'fftw3.f'

integer :: n1, n2, n3

   n1 = Nfft(1)
   n2 = Nfft(2)
   n3 = Nfft(3)

   call dfftw_plan_dft_1d(planFx, n1, Forward_in_x, Forward_out_x, FFTW_FORWARD, FFTW_MEASURE)
   call dfftw_plan_dft_1d(planFy, n2, Forward_in_y, Forward_out_y, FFTW_FORWARD, FFTW_MEASURE)
   call dfftw_plan_dft_1d(planFz, n3, Forward_in_z, Forward_out_z, FFTW_FORWARD, FFTW_MEASURE)

   call dfftw_plan_dft_1d(planBx, n1, Backward_in_x, Backward_out_x, &
   &                      FFTW_BACKWARD, FFTW_MEASURE)
   call dfftw_plan_dft_1d(planBy, n2, Backward_in_y, Backward_out_y, &
   &                      FFTW_BACKWARD, FFTW_MEASURE)
   call dfftw_plan_dft_1d(planBz, n3, Backward_in_z, Backward_out_z, &
   &                      FFTW_BACKWARD, FFTW_MEASURE)

end subroutine fftw_ini


#endif


#ifdef PUBFFT


!######################################################################
!######################################################################


subroutine pubz3d(nnsign,n1,n2,n3,w,ld1,ld2,table,ntable,work,nwork)

use CommonBlocks, only : QPathInt
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess

implicit none

integer :: n1,n2,n3,ld1,ld2,nnsign,ntable,nwork
complex(8), dimension(ld1,ld2,n3) :: w
complex(8), dimension(nwork) :: work
real(8), dimension(ntable,3) :: table

integer :: i,j,k
integer :: Nas
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   if(nnsign == -1) then

! ntable should be 4*max(n1,n2,n3) +15
! nwork should be max(n1,n2,n3)
!
!   transform along X  first ...
!

   do k = Nas, n3, NProcsTemp

     do j = 1, n2

       do i = 1,n1
         work(i) = w(i,j,k)
       end do

       call cfftf(n1,work,table(1,1))

       do i = 1,n1
         w(i,j,k) = work(i)
       end do

     end do

   end do
!
!   transform along Y then ...
!
   do k = Nas, n3, NProcsTemp

     do i = 1,n1

       do j = 1,n2
         work(j) = w(i,j,k)
       end do

       call cfftf(n2,work,table(1,2))

       do j = 1,n2
         w(i,j,k) = work(j)
       end do

     end do

   end do

! ---------------------------------------------------------
   if(NProcsTemp/=1) call FFT_ChAxisF(w,n1,n2,n3,ld1,ld2)
! ---------------------------------------------------------
!
!   transform along Z finally ...
!
   do i = Nas, n1, NProcsTemp

     do j = 1, n2

       do k = 1,n3
         work(k) = w(i,j,k)
       end do

       call cfftf(n3,work,table(1,3))

       do k = 1,n3
         w(i,j,k) = work(k)
       end do

     end do

   end do

   else

!
!   transform along Z finally ...
!
   do i = Nas, n1, NProcsTemp

     do j = 1, n2

       do k = 1, n3
         work(k) = w(i,j,k)
       end do

       call cfftb(n3,work,table(1,3))

       do k = 1, n3
         w(i,j,k) = work(k)
       end do

     end do

   end do

! ---------------------------------------------------------
   if(NProcsTemp/=1) call FFT_ChAxisB(w,n1,n2,n3,ld1,ld2)
! ---------------------------------------------------------

! ntable should be 4*max(n1,n2,n3) +15
! nwork should be max(n1,n2,n3)
!
!   transform along X  first ...
!

   do k = Nas, n3, NProcsTemp

     do j = 1, n2

       do i = 1,n1
         work(i) = w(i,j,k)
       end do

       call cfftb(n1,work,table(1,1))

       do i = 1,n1
         w(i,j,k) = work(i)
       end do

     end do

   end do
!
!   transform along Y then ...
!
   do k = Nas, n3, NProcsTemp

     do i = 1,n1

       do j = 1,n2
         work(j) = w(i,j,k)
       end do

       call cfftb(n2,work,table(1,2))

       do j = 1,n2
         w(i,j,k) = work(j)
       end do

     end do

   end do

   end if

end subroutine pubz3d


#endif


#ifdef SGIFFT


!######################################################################
!######################################################################


subroutine sgifft3d(nnsign,n1,n2,n3,w,ld1,ld2,table,ntable,work,nwork)

use CommonBlocks, only : QPathInt
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess

implicit none

integer :: n1,n2,n3,ld1,ld2,nnsign,ntable,nwork
complex(8), dimension(ld1,ld2,n3) :: w
complex(8), dimension(nwork) :: work
real(8), dimension(ntable,3) :: table

integer :: i,j,k
integer :: NProcsTemp, MyRankTemp, Nas

! ntable should be 4*max(n1,n2,n3) +15
! nwork should be max(n1,n2,n3)
!
!   transform along X  first ...
!

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   if(nnsign == -1) then

   do k = Nas, n3, NProcsTemp

     do j = 1, n2

       do i = 1,n1
         work(i) = w(i,j,k)
       end do

       call cfft1d(nnsign,n1,work,1,table(1,1))

       do i = 1,n1
         w(i,j,k) = work(i)
       end do

     end do

   end do
!
!   transform along Y then ...
!
   do k = Nas, n3, NProcsTemp

     do i = 1,n1

       do j = 1,n2
         work(j) = w(i,j,k)
       end do

       call cfft1d(nnsign,n2,work,1,table(1,2))

       do j = 1,n2
         w(i,j,k) = work(j)
       end do

     end do

   end do

! --------------------------------------------------------
   if(NProcsTemp/=1) call FFT_ChAxisF(w,n1,n2,n3,ld1,ld2)
! --------------------------------------------------------
!
!   transform along Z finally ...
!
   do i = Nas, n1, NProcsTemp

     do j = 1, n2

       do k = 1,n3
         work(k) = w(i,j,k)
       end do

       call cfft1d(nnsign,n3,work,1,table(1,3))

       do k = 1,n3
         w(i,j,k) = work(k)
       end do

     end do

   end do

   else

!
!   transform along Z finally ...
!
   do i = Nas, n1, NProcsTemp

     do j = 1, n2

       do k = 1, n3
         work(k) = w(i,j,k)
       end do

       call cfft1d(nnsign,n3,work,1,table(1,3))

       do k = 1, n3
         w(i,j,k) = work(k)
       end do

     end do

   end do

! ---------------------------------------------------------
   if(NProcsTemp/=1) call FFT_ChAxisB(w,n1,n2,n3,ld1,ld2)
! ---------------------------------------------------------

! ntable should be 4*max(n1,n2,n3) +15
! nwork should be max(n1,n2,n3)
!
!   transform along X  first ...
!

   do k = Nas, n3, NProcsTemp

     do j = 1, n2

       do i = 1,n1
         work(i) = w(i,j,k)
       end do

       call cfft1d(nnsign,n1,work,1,table(1,1))

       do i = 1,n1
         w(i,j,k) = work(i)
       end do

     end do

   end do
!
!   transform along Y then ...
!
   do k = Nas, n3, NProcsTemp

     do i = 1,n1

       do j = 1,n2
         work(j) = w(i,j,k)
       end do

       call cfft1d(nnsign,n2,work,1,table(1,2))

       do j = 1,n2
         w(i,j,k) = work(j)
       end do

     end do

   end do

   end if

end subroutine sgifft3d


#endif


#ifdef FFTW


!######################################################################
!######################################################################


subroutine fftw3d(nnsign,n1,n2,n3,w)

use CommonBlocks, only : QPathInt
use CommonMPI
use CommonPI, only : MyRankPI, NumProcess
use PMEparam, only : planFx, planFy, planFz, planBx, planBy, planBz, &
                   & Forward_in_x, Forward_out_x, Forward_in_y, &
                   & Forward_out_y, Forward_in_z, Forward_out_z,&
                   & Backward_in_x, Backward_out_x, Backward_in_y, &
                   & Backward_out_y, Backward_in_z, Backward_out_z

implicit none

integer :: n1,n2,n3,ld1,ld2,nnsign
complex(8), dimension(n1,n2,n3) :: w

integer :: i, j, k
integer :: Nas
integer(8) :: plan
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

if(nnsign == -1) then ! Forward 

!   transform along X  first ...
!

   do k = Nas, n3, NProcsTemp

     do j = 1, n2

       do i = 1, n1
         Forward_in_x(i) = w(i,j,k)
       end do

       call dfftw_execute(planFx)

       do i = 1, n1
         w(i,j,k) = Forward_out_x(i)
       end do

     end do

   end do

!   transform along Y then ...
!

   do k = Nas, n3, NProcsTemp

     do i = 1, n1

       do j = 1, n2
         Forward_in_y(j) = w(i,j,k)
       end do

       call dfftw_execute(planFy)

       do j = 1, n2
         w(i,j,k) = Forward_out_y(j)
       end do

     end do

   end do

! --------------------------------------
   if(NProcsTemp/=1) call FFT_ChAxisF(w,n1,n2,n3,ld1,ld2)
! --------------------------------------
!
!   transform along Z finally ...
!

   do i = Nas, n1, NProcsTemp

     do j = 1, n2

       do k = 1, n3
         Forward_in_z(k) = w(i,j,k)
       end do

       call dfftw_execute(planFz)

       do k = 1, n3
         w(i,j,k) = Forward_out_z(k)
       end do

     end do

   end do

else ! Backward 

!   transform along X  first ...
!

   do k = Nas, n3, NProcsTemp

     do j = 1, n2

       do i = 1, n1
         Backward_in_x(i) = w(i,j,k)
       end do

       call dfftw_execute(planBx)

       do i = 1, n1
         w(i,j,k) = Backward_out_x(i)
       end do

     end do

   end do

! --------------------------------------
   if(NProcsTemp/=1) call FFT_ChAxisB(w,n1,n2,n3,ld1,ld2)
! --------------------------------------

!   transform along Y then ...
!

   do k = Nas, n3, NProcsTemp

     do i = 1, n1

       do j = 1, n2
         Backward_in_y(j) = w(i,j,k)
       end do

       call dfftw_execute(planBy)

       do j = 1, n2
         w(i,j,k) = Backward_out_y(j)
       end do

     end do

   end do
!
!   transform along Z finally ...
!

   do i = Nas, n1, NProcsTemp

     do j = 1, n2

       do k = 1, n3
         Backward_in_z(k) = w(i,j,k)
       end do

       call dfftw_execute(planBz)

       do k = 1, n3
         w(i,j,k) = Backward_out_z(k)
       end do

     end do

   end do

end if

end subroutine fftw3d


#endif
