! ############################
! ## SUBROUTINE LIST 
! ## -- Transform_Cell 
! ## -- Transform 
! ## -- ExtendChangeCell 
! ############################


program Transform_Cell

   call Setup

   call Transform
   call ExtendChangeCell

   call SaveParam

end program Transform_Cell


subroutine Transform

use Numbers, only : N
use Configuration
use BathParam, only : Rss, Vss, Vg
use CellParam, only : H, InvH

implicit none

integer :: i
real(8), dimension(3,N) :: ScR, ScV
real(8), dimension(3,3) :: Tfm, Ht, Ht1, Ht2, Ht3
real(8), dimension(3) :: Va
real(8) :: ct, st, cp, sp
real(8) :: ABx, ABy, ABz, Z, La

   Ht = H

   do i = 1 , N
     ScR(:,i) = matmul( InvH, R(:,i) )
     ScV(:,i) = matmul( InvH, Vel(:,i) )
   end do

   ABx = H(2,1)*H(3,2) - H(3,1)*H(2,2)
   ABy = H(3,1)*H(1,2) - H(1,1)*H(3,2)
   ABz = H(1,1)*H(2,2) - H(2,1)*H(1,2)
   Z = sqrt(ABx*ABx + ABy*ABy + ABz*ABz)
   ABx = ABx / Z
   ABy = ABy / Z
   ABz = ABz / Z

   ct = ABz
   st = sqrt(1.d0 - ABz*ABz)
   cp = ABx / sqrt(1.d0 - ABz*ABz)
   sp = ABy / sqrt(1.d0 - ABz*ABz)

! ##

   Tfm(1,1) = cp
   Tfm(1,2) = sp
   Tfm(1,3) = 0.d0
   Tfm(2,1) = -sp
   Tfm(2,2) = cp
   Tfm(2,3) = 0.d0
   Tfm(3,1) = 0.d0
   Tfm(3,2) = 0.d0
   Tfm(3,3) = 1.d0

   call matmat(Tfm,Ht,Ht1)

! ##

   Tfm(1,1) = ct
   Tfm(1,2) = 0.d0
   Tfm(1,3) = -st
   Tfm(2,1) = 0.d0
   Tfm(2,2) = 1.d0
   Tfm(2,3) = 0.d0
   Tfm(3,1) = st
   Tfm(3,2) = 0.d0
   Tfm(3,3) = ct

   call matmat(Tfm,Ht1,Ht2)

! ##

   Tfm(1,1) = cp
   Tfm(1,2) = -sp
   Tfm(1,3) = 0.d0
   Tfm(2,1) = sp
   Tfm(2,2) = cp
   Tfm(2,3) = 0.d0
   Tfm(3,1) = 0.d0
   Tfm(3,2) = 0.d0
   Tfm(3,3) = 1.d0

   call matmat(Tfm,Ht2,Ht3)

! ##

   Va = Ht3(:,1)

   La = sqrt( dot_product( Va, Va ) )
   cp = Va(1) / La
   sp = Va(2) / La

   Tfm(1,1) = cp
   Tfm(1,2) = sp
   Tfm(1,3) = 0.d0
   Tfm(2,1) = -sp
   Tfm(2,2) = cp
   Tfm(2,3) = 0.d0
   Tfm(3,1) = 0.d0
   Tfm(3,2) = 0.d0
   Tfm(3,3) = 1.d0

   call MatMat(Tfm,Ht3,H)

   do i = 1 , N
     R(:,i) = matmul( H, ScR(:,i) )
     Vel(:,i) = matmul( H, ScV(:,i) )
   end do

   Rss = 0.d0
   Vss = 0.d0
   Vg  = 0.d0

end subroutine Transform



subroutine ExtendChangeCell

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QRigidbody
use Configuration, only : R
use RBparam
use CellParam, only : H, InvH
use AtomParam, only : Mass

implicit none

integer :: i, j, Ng, k, moltype
real(8), dimension(3) :: Va, Vb, Vc
real(8) :: Lb2
integer :: ix, iy, ia, ib, nx, ny, nco, Num, MyType
real(8), dimension(3) :: Rg, RgT, Si
integer, dimension(3) :: Ic
real(8) :: SumM

   call InversMatrix(H,InvH)

   Va = H(:,1)
   Vb = H(:,2)
   Vc = H(:,3)

   write(*,'(3(3d16.8/))') H(1,:),H(2,:),H(3,:)

   if(Va(3) /= 0.) then
     write(*,*) 'ERROR va3 = ', Va(3)
     Va(3)     = 0.d0
     H(3,1)    = 0.d0
   end if

   ia = nint(Vc(1)/Va(1))

   if(ia /= 0) then
     Vc(1) = Vc(1) - ia * Va(1)
   end if

   Lb2 = dot_product( Vb, Vb )
   ib  = nint( dot_product( Vc, Vb ) / Lb2 )

   if(ib /= 0) then
     Vc = Vc - ib * Vb
   end if

   if((ia/=0).and.(ib/=0)) then
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     write(*,*) 'WARNING : transformed twice!'
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
   end if

   H(:,3) = Vc

   call InversMatrix(H,InvH)

   if(QRigidBody) then

     Ng = 0
     k  = 0

     do moltype = 1 , NumSpec

       do i = 1 , NumMol(moltype)

         Rg   = 0.d0
         SumM = 0.d0

         do j = 1 , NumRBinMol(moltype)

           Ng = Ng + 1

           if(QSingle(Ng)) then

             k  = k + 1
             Rg   = Rg   + Mass(k) * R_RB(:,Ng)
             SumM = SumM + Mass(k)

           else

             MyType = RBType(Ng)
             Rg     = Rg   + MassRB(MyType) * R_RB(:,Ng)
             SumM   = SumM + MassRB(Mytype)

             k = k + NumRBAtom(MyType)

           end if

         end do

         Rg = Rg / SumM

         nco = 0
         do ix = -2, 2
         do iy = -2, 2

           RgT = Rg + ix*Va + iy*Vb

           Si = matmul( InvH, RgT )
           Ic = nint(Si)
           if(Ic(1)==0.and.Ic(2)==0) then
             nco = nco + 1
             nx = ix
             ny = iy
           end if

         end do
         end do

         if(nco/=1) then
           write(*,*) 'ERROR : nco = ', nco
         end if

         do j = Ng - NumRBinMol(moltype) + 1, Ng

           R_RB(:,j) = R_RB(:,j) + nx*Va + ny*Vb

         end do

       end do

     end do

   else

     Num = 0

     do moltype = 1 , NumSpec

       do i = 1 , NumMol(moltype)

         Rg   = 0.d0
         SumM = 0.d0

         do j = 1 , NumAtm(moltype)

           Num  = Num  + 1
           Rg   = Rg   + Mass(Num) * R(:,Num)
           SumM = SumM + Mass(Num)

         end do

         Rg = Rg / SumM

         nco = 0
         do ix = -2, 2
         do iy = -2, 2

           RgT = Rg + ix*Va + iy*Vb

           Si = matmul( InvH, RgT )
           Ic = nint(Si)
           if(Ic(1)==0.and.Ic(2)==0) then
             nco = nco + 1
             nx = ix
             ny = iy
           end if
           if(Ic(3)/=0) then
             write(*,'(a,2i3)') 'ERROR : Ic(3)/=0 , ix,iy = ',ix,iy
             stop
           end if

         end do
         end do

         if(nco/=1) then
           write(*,'(a,i2)') 'ERROR : nco = ', nco
         end if

         do j = Num - NumAtm(moltype) + 1 , Num

           R(:,j) = R(:,j) + nx*Va + ny*Vb

         end do

       end do

     end do

   end if

end subroutine ExtendChangeCell
