

!######################################################################
!######################################################################


! *****************************************************
! **  subroutine for calculating intermolecular force **
! **  isolated system : No cut-off for interaction    **
! ******************************************************

subroutine Force_iso

use CommonBlocks, only : QPathInt, ForceField, QRigidBody, QEps_r
use Configuration, only : R
use CommonMPI
use CommonPI
use NonbondParam, only : Frc_Ersp, Frc_Elec, Ene_LJ, Ene_Ersp, &
& Ene_Elec

implicit none

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

! -----------
!  Zero clear
! -----------
   Frc_Ersp = 0.d0
   Ene_LJ   = 0.d0
   Ene_Ersp = 0.d0

if( ForceField(1:5) == 'CHARM' ) then

! ---------- Bond pair subtraction -------------

   call Force_Subt_14

! ----------------------------------------------

! ####  Calculation of LJ & Coulomb interaction

   if(QRigidBody) then

     if(QEps_r) then

       call Force_iso_ChRB_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_ChRB(Nas,NProcsTemp)

     end if

   else

     if(QEps_r) then

       call Force_iso_ChFX_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_ChFX(Nas,NProcsTemp)

     end if

   end if

else if(ForceField(1:4) == 'OPLS') then

! ---------- Bond pair subtraction -------------

   if(QEps_r) then

     call Force_Subt_14_Eps_OPLS

   else

     call Force_Subt_14_OPLS

   end if

! ----------------------------------------------

! ####  Calculation of LJ & Coulomb interaction

   if(QRigidBody) then

     if(QEps_r) then

       call Force_iso_OPRB_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_OPRB(Nas,NProcsTemp)

     end if

   else

     if(QEps_r) then

       call Force_iso_OPFX_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_OPFX(Nas,NProcsTemp)

     end if

   end if

end if

   Frc_Elec = Frc_Ersp
   Ene_Elec = Ene_Ersp

end subroutine Force_iso

!######################################################################

subroutine Force_iso_ChRB_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)
         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
inn0:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn0
             end if
           end do inn0

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
           R2  = dot_product( Rij, Rij )

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
inn1:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn1
             end if
           end do inn1

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_ChRB_EpsR

!######################################################################

subroutine Force_iso_ChRB(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R1, R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)
         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
inn2:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn2
             end if
           end do inn2

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
inn3:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn3
             end if
           end do inn3

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_ChRB

!######################################################################

subroutine Force_iso_ChFX_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           noa = 0
inn4:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn4
             end if
           end do inn4

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           noa = 0
inn5:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn5
             end if
           end do inn5

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_ChFX_EpsR

!######################################################################

subroutine  Force_iso_ChFX(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R1, R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           noa = 0
inn6:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn6
             end if
           end do inn6

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           noa = 0
inn7:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit inn7
             end if
           end do inn7

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + Eps * ( SR12 - 2.d0 * SR6 )

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine  Force_iso_ChFX

!######################################################################

subroutine Force_iso_OPRB_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)
         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
jnn0:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn0
             end if
           end do jnn0

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
jnn1:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn1
             end if
           end do jnn1

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_OPRB_EpsR

!######################################################################

subroutine Force_iso_OPRB(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R1, R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB

integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)
         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
jnn2:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn2
             end if
           end do jnn2

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           noa = 0
jnn3:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn3
             end if
           end do jnn3

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_OPRB

!######################################################################

subroutine Force_iso_OPFX_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use RBparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           noa = 0
jnn4:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn4
             end if
           end do jnn4

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           noa = 0
jnn5:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn5
             end if
           end do jnn5

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_OPFX_EpsR

!######################################################################

subroutine Force_iso_OPFX(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NoLJparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, non, noa
real(8) :: R1, R2, InvR2
real(8) :: fkLJ
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, cf
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         non = NumNoLJ(i)

         do j = i-2 , 1, -2

           noa = 0
jnn6:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn6
             end if
           end do jnn6

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

         do j = i+1 , N, 2

           noa = 0
jnn7:      do k = 1, non
             if(j == NoLJ(i,k)) then
               noa=1
               exit jnn7
             end if
           end do jnn7

           if(noa == 1) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif

           InvR2 = 1.d0 / R2

           Eps   = EpsLJ(i) * EpsLJ(j)

           if(Eps /= 0.) then

             Sgm2  = SgmLJ(i) * SgmLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

             Fij = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_LJ = Ene_LJ + 4.d0 * Eps * ( SR12 - SR6 )

           end if

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

end subroutine Force_iso_OPFX


!######################################################################
!######################################################################


! ************************************************************************
! **  subroutine for calculating intermolecular force : isolated system **
! **  Switching function is used in the region from "Ron" to "Roff"     **
! ************************************************************************

subroutine Force_iso_SW

use CommonBlocks, only : QPathInt, ForceField, QRigidBody, QEps_r
use Configuration, only : R
use CommonMPI
use CommonPI
use NonbondParam, only : Frc_Ersp, Frc_Elec, Ene_LJ, Ene_Ersp, &
& Ene_Elec

implicit none

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

! -----------
!  Zero clear
! -----------
   Frc_Ersp = 0.d0
   Ene_LJ   = 0.d0
   Ene_Ersp = 0.d0

if( ForceField(1:5) == 'CHARM' ) then

! ---------- Bond pair subtraction -------------

   if(QEps_r) then

     call Force_Subt_12_Eps
     call Force_Subt_13_Eps

   else

     call Force_Subt_12
     call Force_Subt_13

   end if

   call Force_Subt_14

! ----------------------------------------------

! ####  Calculation of LJ & Coulomb interaction

   if(QRigidBody) then

     if(QEps_r) then

       call Force_iso_SW_ChRB_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_SW_ChRB(Nas,NProcsTemp)

     end if

   else

     if(QEps_r) then

       call Force_iso_SW_ChFX_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_SW_ChFX(Nas,NProcsTemp)

     end if

   end if

else if(ForceField(1:4) == 'OPLS') then

! ---------- Bond pair subtraction -------------

   if(QEps_r) then

     call Force_Subt_12_Eps_OPLS
     call Force_Subt_13_Eps_OPLS
     call Force_Subt_14_Eps_OPLS

   else

     call Force_Subt_12_OPLS
     call Force_Subt_13_OPLS
     call Force_Subt_14_OPLS

   end if

! ----------------------------------------------

! ####  Calculation of LJ & Coulomb interaction

   if(QRigidBody) then

     if(QEps_r) then

       call Force_iso_SW_OPRB_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_SW_OPRB(Nas,NProcsTemp)

     end if

   else

     if(QEps_r) then

       call Force_iso_SW_OPFX_EpsR(Nas,NProcsTemp)

     else

       call Force_iso_SW_OPFX(Nas,NProcsTemp)

     end if

   end if

end if

   Frc_Elec = Frc_Ersp
   Ene_Elec = Ene_Ersp

end subroutine Force_iso_SW

!######################################################################

subroutine Force_iso_SW_ChRB_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use RBparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

           end if

         end do

       end do

end subroutine Force_iso_SW_ChRB_EpsR

!######################################################################

subroutine Force_iso_SW_ChRB(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use RBparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

           end if

         end do

       end do

end subroutine Force_iso_SW_ChRB

!######################################################################

subroutine Force_iso_SW_ChFX_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         do j = i-2 , 1, -2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

           end if

         end do

         do j = i+1 , N, 2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if

           end if

         end do

       end do

end subroutine Force_iso_SW_ChFX_EpsR

!######################################################################

subroutine Force_iso_SW_ChFX(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         do j = i-2 , 1, -2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

         do j = i+1 , N, 2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

       end do

end subroutine Force_iso_SW_ChFX

!######################################################################

subroutine Force_iso_SW_OPRB_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use RBparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm2, Eps

integer :: i, j
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

       end do

end subroutine Force_iso_SW_OPRB_EpsR

!######################################################################

subroutine Force_iso_SW_OPRB(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use RBparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm2, Eps

integer :: i, j
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: IRB
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         IRB = AtomUnitNum(i)

         do j = i-2 , 1, -2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

         do j = i+1 , N, 2

           if(IRB == AtomUnitNum(j)) cycle

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

       end do

end subroutine Force_iso_SW_OPRB

!######################################################################

subroutine Force_iso_SW_OPFX_EpsR(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm2, Eps

integer :: i, j
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         do j = i-2 , 1, -2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

         do j = i+1 , N, 2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             fk1 =  cf * InvR2 * 0.25d0
             fk  = fk1 * InvR2 * 2.d0

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

       end do

end subroutine Force_iso_SW_OPFX_EpsR

!######################################################################

subroutine Force_iso_SW_OPFX(Nas,NProcsTemp)

use Numbers, only : N
use Configuration, only : R
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

real(8) :: Sgm2, Eps

integer :: i, j
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fkLJ
real(8) :: fk1, fk, cf, ek
real(8), dimension(3) :: Rij
real(8), dimension(3) :: Fij

integer :: Nas
integer :: NProcsTemp

       do i = Nas, N, NProcsTemp

         do j = i-2 , 1, -2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

         do j = i+1 , N, 2

           Rij = R(:,i) - R(:,j)
#ifdef PCC
           R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
           R2  = dot_product( Rij, Rij )
#endif
           if(R2<=Rcutoff2) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
           call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Fij = fkLJ * Rij

           Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
           Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           if( cf /= 0.d0 ) then

             R1  = sqrt( R2 )

             fk1 = cf / R1
             fk  = fk1 * InvR2

             Fij = fk * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

             Ene_Ersp = Ene_Ersp + fk1

           end if
           end if

         end do

       end do

end subroutine Force_iso_SW_OPFX
