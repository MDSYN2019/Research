! ############################
! ## SUBROUTINE LIST 
! ## -- Force_Subt_12 
! ## -- Force_Subt_13 
! ## -- Force_Subt_14 
! ## -- Force_Subt_12_Eps 
! ## -- Force_Subt_13_Eps 
! ## -- Force_Subt_12_OPLS 
! ## -- Force_Subt_13_OPLS 
! ## -- Force_Subt_14_OPLS 
! ## -- Force_Subt_12_Eps_OPLS 
! ## -- Force_Subt_13_Eps_OPLS 
! ## -- Force_Subt_14_Eps_OPLS 
! ## -- Force_Subt_CG 
! ## -- Force_Subt_CG12 
! ## -- Force_Subt_CG13 
! ## -- COULOMB_SUB
! ## -- CorrectPotential 
! ## -- Viradd 
! ############################


!######################################################################
!######################################################################


subroutine Force_Subt_12

use CommonBlocks, only : QMaster
use Configuration, only : R
use RBparam
use SHAKEparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use BondedParam, only : NumBond, BondI, BondJ

implicit NONE

integer :: i, j, k, ii, jj, kk
real(8) :: Sgm, Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, FijLJ, FijEL

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm   = Rminh(i) + Rminh(j)
     Sgm2  = Sgm * Sgm
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2              !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2           !         ^6
     SR12 = SR6 * SR6                 !         ^12
     fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     FijLJ = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

     Ene_LJ = Ene_LJ - Eps * ( SR12 - 2.d0 * SR6 )

     cf = Charge(i) * Charge(j)

     if( cf /= 0. ) then

       fk1 = cf * sqrt( InvR2 )
       fk  = fk1 * InvR2

       FijEL = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

       Ene_Ersp = Ene_Ersp - fk1

     else

       FijEL = 0.d0

     end if

     call Viradd(FijLJ+FijEL,Rij,i,j)

   end do

   if(QMaster) then

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleBond(k)

       ii = CouplePair(k,kk,1)
       jj = CouplePair(k,kk,2)

       i  = CoupleAtom(k,ii)
       j  = CoupleAtom(k,jj)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       Sgm   = Rminh(i) + Rminh(j)
       Sgm2  = Sgm * Sgm
       Eps   = EpsLJ(i) * EpsLJ(j)
       InvR2 = 1.d0 / R2

       SR2  = Sgm2 * InvR2            !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2         !         ^6
       SR12 = SR6 * SR6               !         ^12
       fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

       FijLJ = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

       Ene_LJ = Ene_LJ - Eps * ( SR12 - 2.d0 * SR6 )

       cf = Charge(i) * Charge(j)

       if( cf /= 0. ) then

         fk1 = cf * sqrt( InvR2 )
         fk  = fk1 * InvR2

         FijEL = - fk * Rij

         Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
         Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

         Ene_Ersp = Ene_Ersp - fk1

       else

         FijEL = 0.d0

       end if

       call Viradd(FijLJ+FijEL,Rij,i,j)

     end do

   end do

   end if

end subroutine Force_Subt_12


!######################################################################
!######################################################################


subroutine Force_Subt_13

use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ
use BondedParam, only : NumAngle, AngleI, AngleK

implicit NONE

integer :: i, j, k
real(8) :: Sgm, Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, FijLJ, FijEL

   do k = 1 , NumAngle

     i = AngleI(k)
     j = AngleK(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm   = Rminh(i) + Rminh(j)
     Sgm2  = Sgm * Sgm
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2             !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2          !         ^6
     SR12 = SR6 * SR6                !         ^12
     fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     FijLJ = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

     Ene_LJ = Ene_LJ - Eps * ( SR12 - 2.d0 * SR6 )

     cf = Charge(i) * Charge(j)

     if(cf/=0.) then

       fk1 = cf * sqrt( InvR2 )
       fk  = fk1 * InvR2

       FijEL = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

       Ene_Ersp = Ene_Ersp - fk1

     else

       FijEL = 0.d0

     end if

     call Viradd(FijLJ+FijEL,Rij,i,j)

   end do

end subroutine Force_Subt_13


!######################################################################
!######################################################################


subroutine Force_Subt_14

use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_LJ, &
& Rminh14, Eps14
use BondedParam, only : NumDihedral, DihedI, DihedL, vdWSubtDih

implicit NONE

integer :: i, j, k
real(8) :: Sgm, Sgm2, Eps, Sgm14, Sgm14_2, Eps14i
real(8) :: SR2, SR6, SR12, fk
real(8) :: SC2, SC6, SC12, fk14
real(8) :: R2, InvR2
real(8), dimension(3) :: Rij, FijLJ

   do k = 1 , NumDihedral

     if(vdWSubtDih(k)) cycle

     i = DihedI(k)
     j = DihedL(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm   = Rminh(i) + Rminh(j)
     Sgm2  = Sgm * Sgm
     Eps   = EpsLJ(i) * EpsLJ(j)

     Sgm14   = Rminh14(i) + Rminh14(j)
     Sgm14_2 = Sgm14 * Sgm14
     Eps14i  = Eps14(i) * Eps14(j)

     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2          !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2       !         ^6
     SR12 = SR6 * SR6             !         ^12

     SC2  = Sgm14_2 * InvR2
     SC6  = SC2 * SC2 * SC2
     SC12 = SC6 * SC6

     fk   = Eps    * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
     fk14 = Eps14i * 12.d0 * ( SC12 - SC6 ) * InvR2

     FijLJ = ( fk14 - fk ) * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

     Ene_LJ = Ene_LJ + Eps14i * ( SC12 - 2.d0 * SC6 ) &
     &               - Eps    * ( SR12 - 2.d0 * SR6 )

     call Viradd(FijLJ,Rij,i,j)

   end do

end subroutine Force_Subt_14


!######################################################################
!######################################################################


subroutine Force_Subt_12_Eps

use CommonBlocks, only : QMaster
use Configuration, only : R
use SHAKEparam
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumBond, BondI, BondJ

implicit NONE

integer :: i, j, k, ii, jj, kk
real(8) :: Sgm, Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, Fij

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm   = Rminh(i) + Rminh(j)
     Sgm2  = Sgm * Sgm
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2              !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2           !         ^6
     SR12 = SR6 * SR6                 !         ^12
     fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     Fij = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

     Ene_LJ = Ene_LJ - Eps * ( SR12 - 2.d0 * SR6 )

     cf = Charge(i) * Charge(j)

     if( cf /= 0.d0 ) then

       fk1 =  cf * InvR2 * 0.25d0
       fk  = fk1 * InvR2 * 2.d0

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_Ersp = Ene_Ersp - fk1

     end if

   end do

   if(QMaster) then

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleBond(k)

       ii = CouplePair(k,kk,1)
       jj = CouplePair(k,kk,2)

       i  = CoupleAtom(k,ii)
       j  = CoupleAtom(k,jj)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       Sgm   = Rminh(i) + Rminh(j)
       Sgm2  = Sgm * Sgm
       Eps   = EpsLJ(i) * EpsLJ(j)
       InvR2 = 1.d0 / R2

       SR2  = Sgm2 * InvR2            !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2         !         ^6
       SR12 = SR6 * SR6               !         ^12
       fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_LJ = Ene_LJ - Eps * ( SR12 - 2.d0 * SR6 )

       cf = Charge(i) * Charge(j)

       if( cf /= 0.d0 ) then

         fk1 =  cf * InvR2 * 0.25d0
         fk  = fk1 * InvR2 * 2.d0

         Fij = - fk * Rij

         Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
         Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

         Ene_Ersp = Ene_Ersp - fk1

       end if

     end do

   end do

   end if

end subroutine Force_Subt_12_Eps


!######################################################################
!######################################################################


subroutine Force_Subt_13_Eps

use Configuration, only : R
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumAngle, AngleI, AngleK

implicit NONE

integer :: i, j, k
real(8) :: Sgm, Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, Fij

   do k = 1 , NumAngle

     i = AngleI(k)
     j = AngleK(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm   = Rminh(i) + Rminh(j)
     Sgm2  = Sgm * Sgm
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2             !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2          !         ^6
     SR12 = SR6 * SR6                !         ^12
     fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     Fij = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

     Ene_LJ = Ene_LJ - Eps * ( SR12 - 2.d0 * SR6 )

     cf = Charge(i) * Charge(j)

     if( cf /= 0.d0 ) then

       fk1 =  cf * InvR2 * 0.25
       fk  = fk1 * InvR2 * 2.d0

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_Ersp = Ene_Ersp - fk1

     end if

   end do

end subroutine Force_Subt_13_Eps


!######################################################################
!######################################################################


subroutine Force_Subt_12_OPLS

use CommonBlocks, only : QMaster
use Configuration, only : R
use RBparam
use SHAKEparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumBond, BondI, BondJ

implicit NONE

integer :: i, j, k, ii, jj, kk
real(8) :: Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, FijLJ, FijEL

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm2  = SgmLJ(i) * SgmLJ(j)
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2              !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2           !         ^6
     SR12 = SR6 * SR6                 !         ^12
     fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     FijLJ = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

     Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

     cf = Charge(i) * Charge(j)

     if( cf /= 0. ) then

       fk1 = cf * sqrt( InvR2 )
       fk  = fk1 * InvR2

       FijEL = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

       Ene_Ersp = Ene_Ersp - fk1

     else

       FijEL = 0.d0

     end if

     call Viradd(FijLJ+FijEL,Rij,i,j)

   end do

   if(QMaster) then

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleBond(k)

       ii = CouplePair(k,kk,1)
       jj = CouplePair(k,kk,2)

       i  = CoupleAtom(k,ii)
       j  = CoupleAtom(k,jj)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       Sgm2  = SgmLJ(i) * SgmLJ(j)
       Eps   = EpsLJ(i) * EpsLJ(j)
       InvR2 = 1.d0 / R2

       SR2  = Sgm2 * InvR2            !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2         !         ^6
       SR12 = SR6 * SR6               !         ^12
       fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

       FijLJ = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

       Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

       cf = Charge(i) * Charge(j)

       if( cf /= 0. ) then

         fk1 = cf * sqrt( InvR2 )
         fk  = fk1 * InvR2

         FijEL = - fk * Rij

         Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
         Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

         Ene_Ersp = Ene_Ersp - fk1

       else

         FijEL = 0.d0

       end if

       call Viradd(FijLJ+FijEL,Rij,i,j)

     end do

   end do

   end if

end subroutine Force_Subt_12_OPLS


!######################################################################
!######################################################################


subroutine Force_Subt_13_OPLS

use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumAngle, AngleI, AngleK

implicit NONE

integer :: i, j, k
real(8) :: Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, FijLJ, FijEL

   do k = 1 , NumAngle

     i = AngleI(k)
     j = AngleK(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm2  = SgmLJ(i) * SgmLJ(j)
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2             !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2          !         ^6
     SR12 = SR6 * SR6                !         ^12
     fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     FijLJ = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

     Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

     cf = Charge(i) * Charge(j)

     if(cf/=0.) then

       fk1 = cf * sqrt( InvR2 )
       fk  = fk1 * InvR2

       FijEL = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

       Ene_Ersp = Ene_Ersp - fk1

     else

       FijEL = 0.d0

     end if

     call Viradd(FijLJ+FijEL,Rij,i,j)

   end do

end subroutine Force_Subt_13_OPLS


!######################################################################
!######################################################################


subroutine Force_Subt_14_OPLS

use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumDihedral, DihedI, DihedL, vdWSubtDih

implicit NONE

integer :: i, j, k
real(8) :: Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk, fk1
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, FijLJ, FijEL
real(8), parameter :: Scale = 0.5d0

   do k = 1 , NumDihedral

     if(vdWSubtDih(k)) cycle

     i = DihedI(k)
     j = DihedL(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm2  = SgmLJ(i) * SgmLJ(j)
     Eps   = EpsLJ(i) * EpsLJ(j)

     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2          !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2       !         ^6
     SR12 = SR6 * SR6             !         ^12

     fk   = Eps    * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     FijLJ = - Scale * fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

     Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 ) * Scale

     cf = Charge(i) * Charge(j) * Scale

     if(cf /= 0.) then

       fk1 = cf * sqrt( InvR2 )
       fk  = fk1 * InvR2

       FijEL = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

       Ene_Ersp = Ene_Ersp - fk1

     else

       FijEL = 0.d0

     end if

     call Viradd(FijLJ+FijEL,Rij,i,j)

   end do

end subroutine Force_Subt_14_OPLS


!######################################################################
!######################################################################


subroutine Force_Subt_12_Eps_OPLS

use CommonBlocks, only : QMaster
use Configuration, only : R
use SHAKEparam
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumBond, BondI, BondJ

implicit NONE

integer :: i, j, k, ii, jj, kk
real(8) :: Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, Fij

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm2  = SgmLJ(i) * SgmLJ(j)
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2              !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2           !         ^6
     SR12 = SR6 * SR6                 !         ^12
     fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     Fij = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

     Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

     cf = Charge(i) * Charge(j)

     if( cf /= 0.d0 ) then

       fk1 =  cf * InvR2 * 0.25d0
       fk  = fk1 * InvR2 * 2.d0

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_Ersp = Ene_Ersp - fk1

     end if

   end do

   if(QMaster) then

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleBond(k)

       ii = CouplePair(k,kk,1)
       jj = CouplePair(k,kk,2)

       i  = CoupleAtom(k,ii)
       j  = CoupleAtom(k,jj)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       Sgm2  = SgmLJ(i) * SgmLJ(j)
       Eps   = EpsLJ(i) * EpsLJ(j)
       InvR2 = 1.d0 / R2

       SR2  = Sgm2 * InvR2            !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2         !         ^6
       SR12 = SR6 * SR6               !         ^12
       fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

       cf = Charge(i) * Charge(j)

       if( cf /= 0.d0 ) then

         fk1 =  cf * InvR2 * 0.25d0
         fk  = fk1 * InvR2 * 2.d0

         Fij = - fk * Rij

         Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
         Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

         Ene_Ersp = Ene_Ersp - fk1

       end if

     end do

   end do

   end if

end subroutine Force_Subt_12_Eps_OPLS


!######################################################################
!######################################################################


subroutine Force_Subt_13_Eps_OPLS

use Configuration, only : R
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumAngle, AngleI, AngleK

implicit NONE

integer :: i, j, k
real(8) :: Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, Fij

   do k = 1 , NumAngle

     i = AngleI(k)
     j = AngleK(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm2  = SgmLJ(i) * SgmLJ(j)
     Eps   = EpsLJ(i) * EpsLJ(j)
     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2             !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2          !         ^6
     SR12 = SR6 * SR6                !         ^12
     fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     Fij = - fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

     Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

     cf = Charge(i) * Charge(j)

     if( cf /= 0.d0 ) then

       fk1 =  cf * InvR2 * 0.25
       fk  = fk1 * InvR2 * 2.d0

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_Ersp = Ene_Ersp - fk1

     end if

   end do

end subroutine Force_Subt_13_Eps_OPLS


!######################################################################
!######################################################################


subroutine Force_Subt_14_Eps_OPLS

use Configuration, only : R
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_LJ, Ene_Ersp
use BondedParam, only : NumDihedral, DihedI, DihedL, vdWSubtDih

implicit NONE

integer :: i, j, k
real(8) :: Sgm2, Eps
real(8) :: SR2, SR6, SR12, fk1, fk
real(8) :: R2, InvR2, cf
real(8), dimension(3) :: Rij, Fij
real(8), parameter :: Scale = 0.5d0

   do k = 1 , NumDihedral

     if(vdWSubtDih(k)) cycle

     i = DihedI(k)
     j = DihedL(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     Sgm2  = SgmLJ(i) * SgmLJ(j)
     Eps   = EpsLJ(i) * EpsLJ(j)

     InvR2 = 1.d0 / R2

     SR2  = Sgm2 * InvR2             !(sigma/r)^2
     SR6  = SR2 * SR2 * SR2          !         ^6
     SR12 = SR6 * SR6                !         ^12
     fk   = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)

     Fij = - Scale * fk * Rij

     Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
     Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

     Ene_LJ = Ene_LJ - 4.d0 * Eps * ( SR12 - SR6 )

     cf = Charge(i) * Charge(j) * Scale

     if( cf /= 0.d0 ) then

       fk1 =  cf * InvR2 * 0.25
       fk  = fk1 * InvR2 * 2.d0

       Fij = - fk * Rij

       Frc_Ersp(:,i) = Frc_Ersp(:,i) + Fij
       Frc_Ersp(:,j) = Frc_Ersp(:,j) - Fij

       Ene_Ersp = Ene_Ersp - fk1

     end if

   end do

end subroutine Force_Subt_14_Eps_OPLS


!######################################################################
!######################################################################


subroutine Force_Subt_CG

use CommonBlocks, only : cCOULOMB, QPBC
use Configuration, only : R

implicit none

   if(QPBC.and.(cCOULOMB=='EWALD'.or.cCOULOMB=='PME')) then
     call Force_Subt_CG12
   end if
   call Force_Subt_CG13
!  call Force_Subt_14

end subroutine Force_Subt_CG


!######################################################################
!######################################################################


subroutine Force_Subt_CG12

use CommonBlocks, only : QMaster
use Configuration, only : R
use RBparam
use SHAKEparam
use CGdata
use NonbondParam, only : Charge, Frc_NBshrt, Ene_ELshrt, Vir_NBshrt
use BondedParam, only : NumBond, BondI, BondJ

implicit NONE

integer :: i, j, k, ii, jj, kk
real(8) :: R2, cf, Ueij, fke
real(8), dimension(3) :: Rij, Fij

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     Rij = R(:,i) - R(:,j)
     R2  = dot_product(Rij,Rij)

     cf = Charge(i) * Charge(j)
     if(cf/=0.) then
       call COULOMB_SUB(Ueij,fke,R2,cf)
       Fij = fke * Rij
       Ene_ELshrt = Ene_ELshrt + Ueij
       Frc_NBshrt(:,i) = Frc_NBshrt(:,i) + Fij
       Frc_NBshrt(:,j) = Frc_NBshrt(:,j) - Fij
       Vir_NBshrt(1,1) = Vir_NBshrt(1,1) + Fij(1) * Rij(1)
       Vir_NBshrt(1,2) = Vir_NBshrt(1,2) + Fij(1) * Rij(2)
       Vir_NBshrt(1,3) = Vir_NBshrt(1,3) + Fij(1) * Rij(3)
       Vir_NBshrt(2,2) = Vir_NBshrt(2,2) + Fij(2) * Rij(2)
       Vir_NBshrt(2,3) = Vir_NBshrt(2,3) + Fij(2) * Rij(3)
       Vir_NBshrt(3,3) = Vir_NBshrt(3,3) + Fij(3) * Rij(3)
     end if

   end do

   if(QMaster) then

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleBond(k)

       ii = CouplePair(k,kk,1)
       jj = CouplePair(k,kk,2)

       i  = CoupleAtom(k,ii)
       j  = CoupleAtom(k,jj)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       cf = Charge(i) * Charge(j)
       if(cf/=0.) then
         call COULOMB_SUB(Ueij,fke,R2,cf)
         Fij = fke * Rij
         Ene_ELshrt = Ene_ELshrt + Ueij
         Frc_NBshrt(:,i) = Frc_NBshrt(:,i) + Fij
         Frc_NBshrt(:,j) = Frc_NBshrt(:,j) - Fij
         Vir_NBshrt(1,1) = Vir_NBshrt(1,1) + Fij(1) * Rij(1)
         Vir_NBshrt(1,2) = Vir_NBshrt(1,2) + Fij(1) * Rij(2)
         Vir_NBshrt(1,3) = Vir_NBshrt(1,3) + Fij(1) * Rij(3)
         Vir_NBshrt(2,2) = Vir_NBshrt(2,2) + Fij(2) * Rij(2)
         Vir_NBshrt(2,3) = Vir_NBshrt(2,3) + Fij(2) * Rij(3)
         Vir_NBshrt(3,3) = Vir_NBshrt(3,3) + Fij(3) * Rij(3)
       end if

     end do

   end do

   end if

end subroutine Force_Subt_CG12


!######################################################################
!######################################################################


subroutine Force_Subt_CG13

use CommonBlocks, only : cCOULOMB, QPBC
use Configuration, only : R
use CGdata
use NonbondParam, only : Charge, Frc_NBshrt, Ene_ELshrt, Ene_NBshrt, Vir_NBshrt
use BondedParam, only : NumAngle, AngleI, AngleK, vdWSubtAng
#ifdef GEN
use TableFuncs, only : Invgs
#endif

implicit NONE

integer :: i, j, k, itype, jtype, ftype
real(8) :: R2, cf, Uij, Ueij, aij, bij, ek
real(8) :: InvR2, InvR1, InvR3, term1, term2
real(8) :: InvR4, fk, fke
#ifdef GEN
real(8) :: InvR6, xx, Rm, eps, R1
integer :: ii
real(8) :: InterPolate
#endif
real(8), dimension(3) :: Rij, Fij

   if(cCOULOMB=='EWALD'.or.cCOULOMB=='PME') then

     do k = 1 , NumAngle

       if(vdWSubtAng(k)) cycle

       i = AngleI(k)
       j = AngleK(k)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       itype = NBAtomType(i)
       jtype = NBAtomType(j)

       fk = 0.d0

       if(R2<=Rc13(itype,jtype)) then
!      if(R2>=100000000000000.) then

         ftype = NBFuncType(itype,jtype)

         if(ftype==1) then ! LJ9-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR1 = sqrt(InvR2)
           InvR3 = InvR2 * InvR1
           term1 = aij * InvR3 * InvR3 * InvR3
           term2 = bij * InvR3 * InvR3
           ek = term1 - term2
           fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==7) then ! LJ12-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR4
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

#ifdef GEN
         else if(ftype==8) then ! LJ12-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR6 = InvR2 * InvR2 * InvR2
           term1 = aij * InvR6 * InvR6
           term2 = bij * InvR6
           ek = term1 - term2
           fk = ( 12.d0 * term1 - 6.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==2) then ! LJ6-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR2
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==3) then ! tabulated

           ii  = IDtable(itype,jtype)
           Rm  = Rmin2(itype,jtype)
           xx  = Invgs(ii)
           Uij = InterPolate(R2,Rm,xx,1,ii)
           fk  = - InterPolate(R2,Rm,xx,2,ii)

         else if(ftype==4) then ! morse+switching func.

           eps = CoefAtype(itype,jtype)
           aij = CoefBtype(itype,jtype)
           R1  = sqrt(R2)
           term1 = aij/CoefCtype(itype,jtype)
           term2 = exp(- term1*R1 + aij)
           xx = 1.d0 - term2
           ek = eps * ( xx*xx - 1.d0 )
           fk = - 2.d0 * eps * term1 * xx * term2 / R1
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==5) then ! LJ8-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==6) then ! LJ10-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR2
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)
#endif
         end if

         Ene_NBshrt = Ene_NBshrt + Uij

       end if

       cf = Charge(i) * Charge(j)
       if(cf/=0.) then
         call COULOMB_SUB(Ueij,fke,R2,cf)
         fk = fk + fke
         Ene_ELshrt = Ene_ELshrt + Ueij
       end if
       Fij(:) = fk * Rij(:)
       Frc_NBshrt(:,i) = Frc_NBshrt(:,i) + Fij(:)
       Frc_NBshrt(:,j) = Frc_NBshrt(:,j) - Fij(:)
       Vir_NBshrt(1,1) = Vir_NBshrt(1,1) + Fij(1) * Rij(1)
       Vir_NBshrt(1,2) = Vir_NBshrt(1,2) + Fij(1) * Rij(2)
       Vir_NBshrt(1,3) = Vir_NBshrt(1,3) + Fij(1) * Rij(3)
       Vir_NBshrt(2,2) = Vir_NBshrt(2,2) + Fij(2) * Rij(2)
       Vir_NBshrt(2,3) = Vir_NBshrt(2,3) + Fij(2) * Rij(3)
       Vir_NBshrt(3,3) = Vir_NBshrt(3,3) + Fij(3) * Rij(3)

     end do

   else

     do k = 1 , NumAngle

       if(vdWSubtAng(k)) cycle

       i = AngleI(k)
       j = AngleK(k)

       Rij = R(:,i) - R(:,j)
       R2  = dot_product(Rij,Rij)

       itype = NBAtomType(i)
       jtype = NBAtomType(j)

       if(R2<=Rc13(itype,jtype)) then

         ftype = NBFuncType(itype,jtype)

         if(ftype==1) then ! LJ9-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR1 = sqrt(InvR2)
           InvR3 = InvR2 * InvR1
           term1 = aij * InvR3 * InvR3 * InvR3
           term2 = bij * InvR3 * InvR3
           ek = term1 - term2
           fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==7) then ! LJ12-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR4
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

#ifdef GEN
         else if(ftype==8) then ! LJ12-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR6 = InvR2 * InvR2 * InvR2
           term1 = aij * InvR6 * InvR6
           term2 = bij * InvR6
           ek = term1 - term2
           fk = ( 12.d0 * term1 - 6.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==2) then ! LJ6-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR2
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==3) then ! tabulated

           ii  = IDtable(itype,jtype)
           Rm  = Rmin2(itype,jtype)
           xx  = Invgs(ii)
           Uij = InterPolate(R2,Rm,xx,1,ii)
           fk  = - InterPolate(R2,Rm,xx,2,ii)

         else if(ftype==4) then ! morse+switching func.

           eps = CoefAtype(itype,jtype)
           aij = CoefBtype(itype,jtype)
           R1  = sqrt(R2)
           term1 = aij/CoefCtype(itype,jtype)
           term2 = exp(- term1*R1 + aij)
           xx = 1.d0 - term2
           ek = eps * ( xx*xx - 1.d0 )
           fk = - 2.d0 * eps * term1 * xx * term2 / R1
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==5) then ! LJ8-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)

         else if(ftype==6) then ! LJ10-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR2
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
           Uij = ek + Eps13(itype,jtype)
#endif
         end if

         Ene_NBshrt = Ene_NBshrt + Uij

         Fij = fk * Rij
         Frc_NBshrt(:,i) = Frc_NBshrt(:,i) + Fij
         Frc_NBshrt(:,j) = Frc_NBshrt(:,j) - Fij
         Vir_NBshrt(1,1) = Vir_NBshrt(1,1) + Fij(1) * Rij(1)
         Vir_NBshrt(1,2) = Vir_NBshrt(1,2) + Fij(1) * Rij(2)
         Vir_NBshrt(1,3) = Vir_NBshrt(1,3) + Fij(1) * Rij(3)
         Vir_NBshrt(2,2) = Vir_NBshrt(2,2) + Fij(2) * Rij(2)
         Vir_NBshrt(2,3) = Vir_NBshrt(2,3) + Fij(2) * Rij(3)
         Vir_NBshrt(3,3) = Vir_NBshrt(3,3) + Fij(3) * Rij(3)

       end if

     end do

   end if

end subroutine Force_Subt_CG13


!######################################################################
!######################################################################


subroutine COULOMB_SUB(Uij,fk,R2,cf)

use EwaldParam, only : Alpha, ar2

implicit none

real(8) :: Uij, R2, cf, R1, x, xtm, fk1, fk2, fk
real(8) :: ErrorFunc, Error_Function
external Error_Function

   R1  = sqrt( R2 )
   x   = Alpha * R1

   ErrorFunc = Error_Function(x)

   xtm = -x * x
   fk1 = cf * (ErrorFunc-1.d0) / R1
   fk2 = cf * ar2 * exp(xtm)
   fk  = ( fk1 + fk2 ) / R2

   Uij = fk1

end subroutine COULOMB_SUB


!#####################################################################
!#####################################################################


! ****************************************
! **  Calculation of the correction of  **
! **  Potential Cutoff                  **
! ****************************************

subroutine CorrectPotential

use Numbers, only : N, NumSpec, NumMol
use CommonBlocks, only : ForceField, Qdebug
use CGdata
use UnitExParam, only : pi
use NonbondParam, only : EpsLJ, Rminh, SgmLJ, &
& PairBack, CoeA, CoeB, CoeC
use TailCorrect, only: CorrectE, CorrectV
use CutoffParam, only : Rcutoff2

implicit NONE

integer :: i, j, ij
integer :: itype, jtype, ftype
real(8) :: pref, Inv, Inv3, Inv6, Inv9
real(8) :: prE1, prE2, prV1, prV2
real(8) :: eps, sgm, sg3, sg6, sg12
real(8) :: Aij, Bij, Cij
real(8) :: R1, R2, R3, BR, InvB, InvB2, InvB3
real(8) :: E_co
real(8) :: V_co

   if( ForceField(1:5) == 'CHARM' ) then

   pref = 2.d0 * pi
   Inv  = 1.d0 / sqrt(Rcutoff2)
   Inv3 = Inv  * Inv  * Inv
   Inv9 = Inv3 * Inv3 * Inv3

   prE1 =  1.d0 / 9.d0 * pref * Inv9
   prE2 = -2.d0 / 3.d0 * pref * Inv3
   prV1 = -4.d0 / 3.d0 * pref * Inv9
   prV2 =  4.d0        * pref * Inv3

   CorrectE = 0.d0
   CorrectV = 0.d0

   do i = 1 , N

     do j = 1 , N

       eps = EpsLJ(i) * EpsLJ(j)
       sgm = Rminh(i) + Rminh(j)
       sg3 = sgm * sgm * sgm 
       sg6 = sg3 * sg3
       sg12= sg6 * sg6
       CorrectE = CorrectE + eps * ( prE1 * sg12 + prE2 * sg6 )
       CorrectV = CorrectV + eps * ( prV1 * sg12 + prV2 * sg6 )

     end do

   end do

   else if( ForceField(1:4) == 'OPLS' ) then

   pref = 8.d0 * pi
   Inv  = 1.d0 / sqrt(Rcutoff2)
   Inv3 = Inv  * Inv  * Inv
   Inv9 = Inv3 * Inv3 * Inv3

   prE1 =  1.d0 / 9.d0 * pref * Inv9
   prE2 = -1.d0 / 3.d0 * pref * Inv3
   prV1 = -4.d0 / 3.d0 * pref * Inv9
   prV2 =  2.d0        * pref * Inv3

   CorrectE = 0.d0
   CorrectV = 0.d0

   do i = 1 , N

     do j = 1 , N

       eps = EpsLJ(i) * EpsLJ(j)
       sgm = sqrt( SgmLJ(i) * SgmLJ(j) )
       sg3 = sgm * sgm * sgm
       sg6 = sg3 * sg3
       sg12= sg6 * sg6
       CorrectE = CorrectE + eps * ( prE1 * sg12 + prE2 * sg6 )
       CorrectV = CorrectV + eps * ( prV1 * sg12 + prV2 * sg6 )

     end do

   end do

   else if( ForceField(1:3) == 'BKS' ) then

   R1 = sqrt(Rcutoff2)
   R2 = R1 * R1
   R3 = R1 * R2

   CorrectE = 0.d0
   CorrectV = 0.d0

   do i = 1, NumSpec

     do j = 1, NumSpec

       ij = PairBack(i,j)

       Aij = CoeA(ij)
       Bij = CoeB(ij)
       Cij = CoeC(ij)

       if(Bij==0.) cycle

       InvB  = 1.d0 / Bij
       InvB2 = InvB * InvB
       InvB3 = InvB * InvB2
       BR = Bij * R1

       pref = 2.d0 * pi * NumMol(i) * NumMol(j)

       E_co = pref * ( Aij*InvB* exp(-BR)* (R2 + 2.d0*InvB*R1 + 2.d0*InvB2) &
       &             - Cij / R3 / 3.d0 )

       V_co = pref * ( -Aij * exp(-BR) * (R3 + 3.d0*InvB*R2 + 6.d0*InvB2*R1 + 6.d0*InvB3) &
       &             + 2.d0 * Cij / R3 )

       CorrectE = CorrectE + E_co
       CorrectV = CorrectV + V_co

     end do

   end do

   if(Qdebug) then
   print *, 'Correct=',CorrectE,CorrectV
   end if

   else if( ForceField(1:2) == 'CG' ) then

   pref = 2.d0 * pi

   CorrectE = 0.d0
   CorrectV = 0.d0

   do i = 1 , N

     itype = NBAtomType(i)

     do j = 1 , N

       jtype = NBAtomType(j)
       ftype = NBFuncType(itype,jtype)

       if(ftype==1) then ! LJ9-6

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)

         Inv  = 1.d0 / sqrt(Rcut2(itype,jtype))
         Inv3 = Inv * Inv * Inv
         Inv6 = Inv3 * Inv3

         prE1 =  1.d0 / 6.d0 * pref * Inv6 * aij
         prE2 = -1.d0 / 3.d0 * pref * Inv3 * bij
         prV1 = -1.5d0       * pref * Inv6 * aij
         prV2 =  2.d0        * pref * Inv3 * bij

       else if(ftype==2) then ! LJ6-4

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)

         Inv  = 1.d0 / sqrt(Rcut2(itype,jtype))
         Inv3 = Inv * Inv * Inv

         prE1 =  1.d0 / 3.d0 * pref * Inv3 * aij
         prE2 = -              pref * Inv  * bij
         prV1 = -2.d0        * pref * Inv3 * aij
         prV2 =  4.d0        * pref * Inv  * bij

       else if(ftype==3) then ! Tabular function 

         prE1 = 0.d0
         prE2 = 0.d0
         prV1 = 0.d0
         prV2 = 0.d0

       end if

       CorrectE = CorrectE + (prE1 + prE2)
       CorrectV = CorrectV + (prV1 + prV2)

     end do

   end do

   end if


end subroutine CorrectPotential


!######################################################################
!######################################################################


subroutine Viradd(Fij,Rij,i,j)

use CommonBlocks, only : QRigidBody, QPINPT
use RBparam
use SHAKEparam
use CommonPI, only: Rcentroid
use NonbondParam, only : Vir_Ersp

implicit none

real(8), dimension(3) :: Fij, Rij, RGij
integer :: IRB, JRB, i, j

   if(QRigidBody) then

     IRB = AtomUnitNum(i)
     JRB = AtomUnitNum(j)

     RGij = R_RB(:,IRB) - R_RB(:,JRB)

     Vir_Ersp(1,1) = Vir_Ersp(1,1) + Fij(1) * RGij(1)
     Vir_Ersp(1,2) = Vir_Ersp(1,2) + Fij(1) * RGij(2)
     Vir_Ersp(1,3) = Vir_Ersp(1,3) + Fij(1) * RGij(3)
     Vir_Ersp(2,1) = Vir_Ersp(2,1) + Fij(2) * RGij(1)
     Vir_Ersp(2,2) = Vir_Ersp(2,2) + Fij(2) * RGij(2)
     Vir_Ersp(2,3) = Vir_Ersp(2,3) + Fij(2) * RGij(3)
     Vir_Ersp(3,1) = Vir_Ersp(3,1) + Fij(3) * RGij(1)
     Vir_Ersp(3,2) = Vir_Ersp(3,2) + Fij(3) * RGij(2)
     Vir_Ersp(3,3) = Vir_Ersp(3,3) + Fij(3) * RGij(3)

   else if(QPINPT) then

     RGij = Rcentroid(:,i) - Rcentroid(:,j)

     Vir_Ersp(1,1) = Vir_Ersp(1,1) + Fij(1) * RGij(1)
     Vir_Ersp(1,2) = Vir_Ersp(1,2) + Fij(1) * RGij(2)
     Vir_Ersp(1,3) = Vir_Ersp(1,3) + Fij(1) * RGij(3)
     Vir_Ersp(2,1) = Vir_Ersp(2,1) + Fij(2) * RGij(1)
     Vir_Ersp(2,2) = Vir_Ersp(2,2) + Fij(2) * RGij(2)
     Vir_Ersp(2,3) = Vir_Ersp(2,3) + Fij(2) * RGij(3)
     Vir_Ersp(3,1) = Vir_Ersp(3,1) + Fij(3) * RGij(1)
     Vir_Ersp(3,2) = Vir_Ersp(3,2) + Fij(3) * RGij(2)
     Vir_Ersp(3,3) = Vir_Ersp(3,3) + Fij(3) * RGij(3)

   else

     Vir_Ersp(1,1) = Vir_Ersp(1,1) + Fij(1) * Rij(1)
     Vir_Ersp(1,2) = Vir_Ersp(1,2) + Fij(1) * Rij(2)
     Vir_Ersp(1,3) = Vir_Ersp(1,3) + Fij(1) * Rij(3)
     Vir_Ersp(2,2) = Vir_Ersp(2,2) + Fij(2) * Rij(2)
     Vir_Ersp(2,3) = Vir_Ersp(2,3) + Fij(2) * Rij(3)
     Vir_Ersp(3,3) = Vir_Ersp(3,3) + Fij(3) * Rij(3)

   end if

end subroutine Viradd
