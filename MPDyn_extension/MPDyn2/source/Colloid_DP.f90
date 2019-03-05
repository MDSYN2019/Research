! ############################
! ## SUBROUTINE LIST 
! ## -- RotationMatrix 
! ## -- SiteInfo_Colloid 
! ## -- SiteInfo_ColloidSC 
! ## -- SiteVeloc_Colloid 
! ## -- ForceTorque_COM_Colloid 
! ## -- ForceTorque_COM_Dissipative 
! ############################


!######################################################################
!######################################################################


subroutine RotationMatrix

use CommonDPD
use RBparam

implicit none

real(8), dimension(10) :: qt2
integer :: k
real(8) :: qn, qn2

   do k = 1 , NumColloid

     qn2 = dot_product(Quaternion(:,k),Quaternion(:,k))
     qn  = 1.d0 / sqrt(qn2)

     Quaternion(:,k) = Quaternion(:,k) * qn

   end do

   do k = 1 , NumColloid

     qt2( 1) = Quaternion(1,k) * Quaternion(1,k)
     qt2( 2) = Quaternion(2,k) * Quaternion(2,k)
     qt2( 3) = Quaternion(3,k) * Quaternion(3,k)
     qt2( 4) = Quaternion(4,k) * Quaternion(4,k)

     qt2( 5) = Quaternion(1,k) * Quaternion(2,k)
     qt2( 6) = Quaternion(1,k) * Quaternion(3,k)
     qt2( 7) = Quaternion(1,k) * Quaternion(4,k)
     qt2( 8) = Quaternion(2,k) * Quaternion(3,k)
     qt2( 9) = Quaternion(2,k) * Quaternion(4,k)
     qt2(10) = Quaternion(3,k) * Quaternion(4,k)

     Rotation(1,1,k) = qt2(1) + qt2(2) - qt2(3) - qt2(4)
     Rotation(2,2,k) = qt2(1) - qt2(2) + qt2(3) - qt2(4)
     Rotation(3,3,k) = qt2(1) - qt2(2) - qt2(3) + qt2(4)

     Rotation(1,2,k) = 2.d0 * (qt2( 8) + qt2( 7))
     Rotation(2,1,k) = 2.d0 * (qt2( 8) - qt2( 7))
     Rotation(1,3,k) = 2.d0 * (qt2( 9) - qt2( 6))
     Rotation(3,1,k) = 2.d0 * (qt2( 9) + qt2( 6))
     Rotation(2,3,k) = 2.d0 * (qt2(10) + qt2( 5))
     Rotation(3,2,k) = 2.d0 * (qt2(10) - qt2( 5))

   end do

end subroutine RotationMatrix


!######################################################################
!######################################################################


! ## INPUT  :: Rotation, R_RB, V_RBt, Omegat
! ## OUTPUT :: R, Rmolec, Velt


subroutine SiteInfo_Colloid

use Numbers, only : NumSpec, NumMol, NumAtm
use Configuration, only : R
use CommonDPD
use RBparam

implicit none

integer :: k, j, i, Num
real(8), dimension(3) :: RGi, VGi, Va, Om


! ## body fixed -> space fixed

   do k = 1 , NumColloid

     do j = 1 , NumCollAtm

       Rmolec(1,j,k) = dot_product( Rotation(:,1,k),R_onMol(:,j,1) )
       Rmolec(2,j,k) = dot_product( Rotation(:,2,k),R_onMol(:,j,1) )
       Rmolec(3,j,k) = dot_product( Rotation(:,3,k),R_onMol(:,j,1) )

     end do

   end do

   Num=0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do k = 1 , NumColloid

         RGi = R_RB(:,k)
         VGi = V_RBt(:,k)

! ## body fixed
         Om(:) = Omegat(:,k)

         do j = 1 , NumCollAtm

           Num = Num + 1

! ##### location of colloid components (Space-fixed)

           R(:,Num) = RGi + Rmolec(:,j,k)

! ##### velocity of colloid components (Space-fixed)
!       prepared for calculation of stochastic Force

           Va(1) = Om(2)*R_onMol(3,j,1) - Om(3)*R_onMol(2,j,1)
           Va(2) = Om(3)*R_onMol(1,j,1) - Om(1)*R_onMol(3,j,1)
           Va(3) = Om(1)*R_onMol(2,j,1) - Om(2)*R_onMol(1,j,1)

           Velt(:,Num) = VGi + Va

         end do

       end do

     else

       Num = Num + NumMol(i)*NumAtm(i)

     end if

   end do


end subroutine SiteInfo_Colloid


!######################################################################
!######################################################################


! ## INPUT  :: Rotation, R_RB, V_RB, Omega
! ## OUTPUT :: R, Rmolec, Vel


subroutine SiteInfo_ColloidSC

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonDPD
use RBparam
use Configuration

implicit none

integer :: k, j, i, Num
real(8), dimension(3) :: RGi, VGi, Va, Om


! ## body fixed -> space fixed

   do k = 1 , NumColloid

     do j = 1 , NumCollAtm

       Rmolec(1,j,k) = dot_product( Rotation(:,1,k),R_onMol(:,j,1) )
       Rmolec(2,j,k) = dot_product( Rotation(:,2,k),R_onMol(:,j,1) )
       Rmolec(3,j,k) = dot_product( Rotation(:,3,k),R_onMol(:,j,1) )

     end do

   end do

   Num=0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do k = 1 , NumColloid

         RGi = R_RB(:,k)
         VGi = V_RB(:,k)

! ## body fixed
         Om(:) = Omega(:,k)

         do j = 1 , NumCollAtm

           Num = Num + 1

! ##### location of colloid components (Space-fixed)

           R(:,Num) = RGi + Rmolec(:,j,k)

! ##### velocity of colloid components (Space-fixed)
!       prepared for calculation of stochastic Force

           Va(1) = Om(2)*R_onMol(3,j,1) - Om(3)*R_onMol(2,j,1)
           Va(2) = Om(3)*R_onMol(1,j,1) - Om(1)*R_onMol(3,j,1)
           Va(3) = Om(1)*R_onMol(2,j,1) - Om(2)*R_onMol(1,j,1)

           Vel(:,Num) = VGi + Va

         end do

       end do

     else

       Num = Num + NumMol(i)*NumAtm(i)

     end if

   end do


end subroutine SiteInfo_ColloidSC


!######################################################################
!######################################################################


! ## INPUT  :: V_RB, Omega
! ## OUTPUT :: Vel

subroutine SiteVeloc_Colloid

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonDPD
use RBparam
use Configuration, only : Vel

implicit none

integer :: k, j, i, Num
real(8), dimension(3) :: VGi, Va, Om

   Num=0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do k = 1 , NumColloid

         VGi = V_RB(:,k)

! ## body fixed
         Om(:) = Omega(:,k)

         do j = 1 , NumCollAtm

           Num = Num + 1

! ##### velocity of colloid components (Space-fixed)
!       prepared for calculation of stochastic Force

           Va(1) = Om(2)*R_onMol(3,j,1) - Om(3)*R_onMol(2,j,1)
           Va(2) = Om(3)*R_onMol(1,j,1) - Om(1)*R_onMol(3,j,1)
           Va(3) = Om(1)*R_onMol(2,j,1) - Om(2)*R_onMol(1,j,1)

           Vel(:,Num) = VGi + Va

         end do

       end do

     else

       Num = Num + NumMol(i)*NumAtm(i)

     end if

   end do


end subroutine SiteVeloc_Colloid


!######################################################################
!######################################################################


subroutine ForceTorque_COM_Colloid

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonDPD
use RBparam

implicit none

real(8), dimension(3,NumCollAtm,NumColloid) :: Fc
real(8), dimension(3) :: Trq
integer :: Num, i, j, k

   Num = 0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do k = 1 , NumColloid

         do j = 1 , NumCollAtm

           Num = Num + 1

           Fc(:,j,k) = FrcDPt(:,Num)

         end do

       end do

     else

       Num = Num + NumMol(i)*NumAtm(i)

     end if

   end do


   FrcCot = 0.d0

   do k = 1 , NumColloid

     Trq = 0.d0

     do j = 1 , NumCollAtm

       FrcCot(:,k) = FrcCot(:,k) + Fc(:,j,k)

       Trq(1) = Trq(1) + Rmolec(2,j,k) * Fc(3,j,k) - Rmolec(3,j,k) * Fc(2,j,k)
       Trq(2) = Trq(2) + Rmolec(3,j,k) * Fc(1,j,k) - Rmolec(1,j,k) * Fc(3,j,k)
       Trq(3) = Trq(3) + Rmolec(1,j,k) * Fc(2,j,k) - Rmolec(2,j,k) * Fc(1,j,k)

     end do

! space fixed ---> body fixed
     Torque(1,k) = dot_product(Rotation(1,:,k),Trq)
     Torque(2,k) = dot_product(Rotation(2,:,k),Trq)
     Torque(3,k) = dot_product(Rotation(3,:,k),Trq)

   end do

end subroutine ForceTorque_COM_Colloid


!######################################################################
!######################################################################


subroutine ForceTorque_COM_Dissipative

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonDPD
use RBparam

implicit none

real(8), dimension(3,NumCollAtm,NumColloid) :: Fc
real(8), dimension(3) :: Trq
integer :: Num, i, j, k

   Num = 0

   do i = 1 , NumSpec

     if(ColloidFlag(i)) then

       do k = 1 , NumColloid

         do j = 1 , NumCollAtm

           Num = Num + 1

           Fc(:,j,k) = FrcDPd(:,Num)

         end do

       end do

     else

       Num = Num + NumMol(i)*NumAtm(i)

     end if

   end do


   FrcCod = 0.d0

   do k = 1 , NumColloid

     Trq = 0.d0

     do j = 1 , NumCollAtm

       FrcCod(:,k) = FrcCod(:,k) + Fc(:,j,k)

       Trq(1) = Trq(1) + Rmolec(2,j,k) * Fc(3,j,k) - Rmolec(3,j,k) * Fc(2,j,k)
       Trq(2) = Trq(2) + Rmolec(3,j,k) * Fc(1,j,k) - Rmolec(1,j,k) * Fc(3,j,k)
       Trq(3) = Trq(3) + Rmolec(1,j,k) * Fc(2,j,k) - Rmolec(2,j,k) * Fc(1,j,k)

     end do

! space fixed ---> body fixed
     Torqued(1,k) = dot_product(Rotation(1,:,k),Trq)
     Torqued(2,k) = dot_product(Rotation(2,:,k),Trq)
     Torqued(3,k) = dot_product(Rotation(3,:,k),Trq)

   end do

end subroutine ForceTorque_COM_Dissipative
