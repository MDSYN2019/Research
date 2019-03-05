! ############################
! ## SUBROUTINE LIST 
! ## -- Force_Bond 
! ## -- Force_Angle 
! ## -- Force_UB 
! ## -- Force_Dihedral 
! ## -- Force_Improper 
! ############################


!######################################################################
!######################################################################


! ***************************
! ** Bond Stretching Force **
! ***************************

subroutine Force_Bond

use Numbers, only : N
use CommonBlocks, only : QRigidBody, QPINPT
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use BondedParam, only : NumBond, BondI, BondJ, kBond, rBond, &
&   Frc_Bond, Vir_Bond, Ene_Bond

implicit NONE

integer :: i, j, k
real(8) :: R2, R1, dR, Fc
real(8), dimension(3) :: Rij, Fij

   Ene_Bond = 0.d0
   Frc_Bond = 0.d0
   Vir_Bond = 0.d0

   if( NumBond == 0 ) Return

   do k = 1 , NumBond

     i= BondI(k)
     j= BondJ(k)

     Rij = R(:,i) - R(:,j)

#ifdef PCC
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
     R2  = dot_product(Rij,Rij)
#endif

     R1 = sqrt(R2)
     dR = R1 - rBond(k)
     Fc = kBond(k) * dR

     Ene_Bond = Ene_Bond + Fc * dR

     Fij = - 2.d0 * Fc * Rij / R1

     Frc_Bond(:,i) = Frc_Bond(:,i) + Fij
     Frc_Bond(:,j) = Frc_Bond(:,j) - Fij

   end do

   if(QRigidBody) then

     do i = 1  , N

       j = AtomUnitNum(i)

       Vir_Bond(1,1) = Vir_Bond(1,1) + Frc_Bond(1,i) * R_RB(1,j)
       Vir_Bond(1,2) = Vir_Bond(1,2) + Frc_Bond(1,i) * R_RB(2,j)
       Vir_Bond(1,3) = Vir_Bond(1,3) + Frc_Bond(1,i) * R_RB(3,j)
       Vir_Bond(2,1) = Vir_Bond(2,1) + Frc_Bond(2,i) * R_RB(1,j)
       Vir_Bond(2,2) = Vir_Bond(2,2) + Frc_Bond(2,i) * R_RB(2,j)
       Vir_Bond(2,3) = Vir_Bond(2,3) + Frc_Bond(2,i) * R_RB(3,j)
       Vir_Bond(3,1) = Vir_Bond(3,1) + Frc_Bond(3,i) * R_RB(1,j)
       Vir_Bond(3,2) = Vir_Bond(3,2) + Frc_Bond(3,i) * R_RB(2,j)
       Vir_Bond(3,3) = Vir_Bond(3,3) + Frc_Bond(3,i) * R_RB(3,j)

     end do

   else if(QPINPT) then

     do i = 1  , N

       Vir_Bond(1,1) = Vir_Bond(1,1) + Frc_Bond(1,i) * Rcentroid(1,i)
       Vir_Bond(1,2) = Vir_Bond(1,2) + Frc_Bond(1,i) * Rcentroid(2,i)
       Vir_Bond(1,3) = Vir_Bond(1,3) + Frc_Bond(1,i) * Rcentroid(3,i)
       Vir_Bond(2,1) = Vir_Bond(2,1) + Frc_Bond(2,i) * Rcentroid(1,i)
       Vir_Bond(2,2) = Vir_Bond(2,2) + Frc_Bond(2,i) * Rcentroid(2,i)
       Vir_Bond(2,3) = Vir_Bond(2,3) + Frc_Bond(2,i) * Rcentroid(3,i)
       Vir_Bond(3,1) = Vir_Bond(3,1) + Frc_Bond(3,i) * Rcentroid(1,i)
       Vir_Bond(3,2) = Vir_Bond(3,2) + Frc_Bond(3,i) * Rcentroid(2,i)
       Vir_Bond(3,3) = Vir_Bond(3,3) + Frc_Bond(3,i) * Rcentroid(3,i)

     end do

   else

     do i = 1 , N

       Vir_Bond(1,1) = Vir_Bond(1,1) + Frc_Bond(1,i) * R(1,i)
       Vir_Bond(1,2) = Vir_Bond(1,2) + Frc_Bond(1,i) * R(2,i)
       Vir_Bond(1,3) = Vir_Bond(1,3) + Frc_Bond(1,i) * R(3,i)
       Vir_Bond(2,2) = Vir_Bond(2,2) + Frc_Bond(2,i) * R(2,i)
       Vir_Bond(2,3) = Vir_Bond(2,3) + Frc_Bond(2,i) * R(3,i)
       Vir_Bond(3,3) = Vir_Bond(3,3) + Frc_Bond(3,i) * R(3,i)

     end do

     Vir_Bond(2,1) = Vir_Bond(1,2)
     Vir_Bond(3,1) = Vir_Bond(1,3)
     Vir_Bond(3,2) = Vir_Bond(2,3)

   end if

end subroutine Force_Bond


!######################################################################
!######################################################################


! ********************
! **  Bending Force **
! ********************

subroutine Force_Angle

use Numbers, only : N
use CommonBlocks, only : QRigidBody, QPINPT
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use BondedParam, only : NumAngle, AngleI, AngleJ, AngleK, kTheta, Theta0, &
&   Frc_Angle, Vir_Angle, Ene_Angle

implicit NONE

integer :: i, j, k, l
real(8) :: Ra2, Rb2
real(8) :: rRa2, rRb2, rRab
real(8) :: Cst, Theta
real(8) :: dT
real(8) :: cf
real(8), dimension(3) :: Rij, Rkj
real(8), dimension(3) :: Fi, Fj, Fk

   Ene_Angle = 0.d0
   Frc_Angle = 0.d0
   Vir_Angle = 0.d0

   if( NumAngle == 0 ) Return

   do l = 1 , NumAngle

   if( kTheta(l) == 0. ) cycle

     i = AngleI(l)
     j = AngleJ(l)
     k = AngleK(l)

     Rij = R(:,i) - R(:,j)
     Rkj = R(:,k) - R(:,j)

#ifdef PCC
     Ra2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
     Rb2  = Rkj(1)*Rkj(1) + Rkj(2)*Rkj(2) + Rkj(3)*Rkj(3)
#else
     Ra2 = dot_product(Rij,Rij)
     Rb2 = dot_product(Rkj,Rkj)
#endif

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt( rRa2 * rRb2 )

#ifdef PCC
     Cst = ( Rij(1)*Rkj(1) + Rij(2)*Rkj(2) + Rij(3)*Rkj(3) ) * rRab
#else
     Cst = dot_product(Rij,Rkj) * rRab
#endif

     Theta = acos(Cst)

     dT = Theta - Theta0(l)

     cf = 2.d0 * kTheta(l) * dT / sin(Theta)

     Fi = cf * ( Rkj * rRab - Cst * Rij * rRa2 )
     Fk = cf * ( Rij * rRab - Cst * Rkj * rRb2 )
     Fj = - ( Fi + Fk )

     Frc_Angle(:,i) = Frc_Angle(:,i) + Fi
     Frc_Angle(:,j) = Frc_Angle(:,j) + Fj
     Frc_Angle(:,k) = Frc_Angle(:,k) + Fk

     Ene_Angle = Ene_Angle + kTheta(l) * dT * dT

   end do

   if(QRigidBody) then

     do i = 1, N

       j = AtomUnitNum(i)

       Vir_Angle(1,1) = Vir_Angle(1,1) + Frc_Angle(1,i) * R_RB(1,j)
       Vir_Angle(1,2) = Vir_Angle(1,2) + Frc_Angle(1,i) * R_RB(2,j)
       Vir_Angle(1,3) = Vir_Angle(1,3) + Frc_Angle(1,i) * R_RB(3,j)
       Vir_Angle(2,1) = Vir_Angle(2,1) + Frc_Angle(2,i) * R_RB(1,j)
       Vir_Angle(2,2) = Vir_Angle(2,2) + Frc_Angle(2,i) * R_RB(2,j)
       Vir_Angle(2,3) = Vir_Angle(2,3) + Frc_Angle(2,i) * R_RB(3,j)
       Vir_Angle(3,1) = Vir_Angle(3,1) + Frc_Angle(3,i) * R_RB(1,j)
       Vir_Angle(3,2) = Vir_Angle(3,2) + Frc_Angle(3,i) * R_RB(2,j)
       Vir_Angle(3,3) = Vir_Angle(3,3) + Frc_Angle(3,i) * R_RB(3,j)

     end do

   else if(QPINPT) then

     do i = 1, N

       Vir_Angle(1,1) = Vir_Angle(1,1) + Frc_Angle(1,i) * Rcentroid(1,i)
       Vir_Angle(1,2) = Vir_Angle(1,2) + Frc_Angle(1,i) * Rcentroid(2,i)
       Vir_Angle(1,3) = Vir_Angle(1,3) + Frc_Angle(1,i) * Rcentroid(3,i)
       Vir_Angle(2,1) = Vir_Angle(2,1) + Frc_Angle(2,i) * Rcentroid(1,i)
       Vir_Angle(2,2) = Vir_Angle(2,2) + Frc_Angle(2,i) * Rcentroid(2,i)
       Vir_Angle(2,3) = Vir_Angle(2,3) + Frc_Angle(2,i) * Rcentroid(3,i)
       Vir_Angle(3,1) = Vir_Angle(3,1) + Frc_Angle(3,i) * Rcentroid(1,i)
       Vir_Angle(3,2) = Vir_Angle(3,2) + Frc_Angle(3,i) * Rcentroid(2,i)
       Vir_Angle(3,3) = Vir_Angle(3,3) + Frc_Angle(3,i) * Rcentroid(3,i)

     end do

   else

     do i = 1, N

       Vir_Angle(1,1) = Vir_Angle(1,1) + Frc_Angle(1,i) * R(1,i)
       Vir_Angle(1,2) = Vir_Angle(1,2) + Frc_Angle(1,i) * R(2,i)
       Vir_Angle(1,3) = Vir_Angle(1,3) + Frc_Angle(1,i) * R(3,i)
       Vir_Angle(2,2) = Vir_Angle(2,2) + Frc_Angle(2,i) * R(2,i)
       Vir_Angle(2,3) = Vir_Angle(2,3) + Frc_Angle(2,i) * R(3,i)
       Vir_Angle(3,3) = Vir_Angle(3,3) + Frc_Angle(3,i) * R(3,i)

     end do

     Vir_Angle(2,1) = Vir_Angle(1,2)
     Vir_Angle(3,1) = Vir_Angle(1,3)
     Vir_Angle(3,2) = Vir_Angle(2,3)

   end if

end subroutine Force_Angle


!######################################################################
!######################################################################


! ******************************
! ** Urey-Bradley Angle Force **
! ******************************

subroutine Force_UB

use Numbers, only : N
use CommonBlocks, only : QRigidBody, ForceField, QPINPT
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use BondedParam, only : NumUB, UB_I, UB_J, Kub, S0, &
&   Frc_UB, Vir_UB, Ene_UB

implicit NONE

integer :: i, j, k
real(8) :: R2, R1, dR, Fc
real(8), dimension(3) :: Rij, Fij

   Ene_UB = 0.d0
   Frc_UB = 0.d0
   Vir_UB = 0.d0

   if( ForceField(1:4) == 'OPLS' ) Return

   if( NumUB == 0 ) Return

   do k = 1 , NumUB

     i= UB_I(k)
     j= UB_J(k)

     Rij = R(:,i) - R(:,j)
#ifdef PCC
     R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
     R2  = dot_product(Rij,Rij)
#endif

     R1 = sqrt(R2)
     dR = R1 - S0(k)
     Fc = Kub(k) * dR

     Ene_UB = Ene_UB + Fc * dR

     Fij = - 2.d0 * Fc * Rij / R1

     Frc_UB(:,i) = Frc_UB(:,i) + Fij
     Frc_UB(:,j) = Frc_UB(:,j) - Fij

   end do

   if(QRigidBody) then

     do i = 1, N

       j = AtomUnitNum(i)

       Vir_UB(1,1) = Vir_UB(1,1) + Frc_UB(1,i) * R_RB(1,j)
       Vir_UB(1,2) = Vir_UB(1,2) + Frc_UB(1,i) * R_RB(2,j)
       Vir_UB(1,3) = Vir_UB(1,3) + Frc_UB(1,i) * R_RB(3,j)
       Vir_UB(2,1) = Vir_UB(2,1) + Frc_UB(2,i) * R_RB(1,j)
       Vir_UB(2,2) = Vir_UB(2,2) + Frc_UB(2,i) * R_RB(2,j)
       Vir_UB(2,3) = Vir_UB(2,3) + Frc_UB(2,i) * R_RB(3,j)
       Vir_UB(3,1) = Vir_UB(3,1) + Frc_UB(3,i) * R_RB(1,j)
       Vir_UB(3,2) = Vir_UB(3,2) + Frc_UB(3,i) * R_RB(2,j)
       Vir_UB(3,3) = Vir_UB(3,3) + Frc_UB(3,i) * R_RB(3,j)

     end do

   else if(QPINPT) then

     do i = 1, N

       Vir_UB(1,1) = Vir_UB(1,1) + Frc_UB(1,i) * Rcentroid(1,i)
       Vir_UB(1,2) = Vir_UB(1,2) + Frc_UB(1,i) * Rcentroid(2,i)
       Vir_UB(1,3) = Vir_UB(1,3) + Frc_UB(1,i) * Rcentroid(3,i)
       Vir_UB(2,1) = Vir_UB(2,1) + Frc_UB(2,i) * Rcentroid(1,i)
       Vir_UB(2,2) = Vir_UB(2,2) + Frc_UB(2,i) * Rcentroid(2,i)
       Vir_UB(2,3) = Vir_UB(2,3) + Frc_UB(2,i) * Rcentroid(3,i)
       Vir_UB(3,1) = Vir_UB(3,1) + Frc_UB(3,i) * Rcentroid(1,i)
       Vir_UB(3,2) = Vir_UB(3,2) + Frc_UB(3,i) * Rcentroid(2,i)
       Vir_UB(3,3) = Vir_UB(3,3) + Frc_UB(3,i) * Rcentroid(3,i)

     end do

   else

     do i = 1 , N

       Vir_UB(1,1) = Vir_UB(1,1) + Frc_UB(1,i) * R(1,i)
       Vir_UB(1,2) = Vir_UB(1,2) + Frc_UB(1,i) * R(2,i)
       Vir_UB(1,3) = Vir_UB(1,3) + Frc_UB(1,i) * R(3,i)
       Vir_UB(2,2) = Vir_UB(2,2) + Frc_UB(2,i) * R(2,i)
       Vir_UB(2,3) = Vir_UB(2,3) + Frc_UB(2,i) * R(3,i)
       Vir_UB(3,3) = Vir_UB(3,3) + Frc_UB(3,i) * R(3,i)

     end do

     Vir_UB(2,1) = Vir_UB(1,2)
     Vir_UB(3,1) = Vir_UB(1,3)
     Vir_UB(3,2) = Vir_UB(2,3)

   end if

end subroutine Force_UB


!######################################################################
!######################################################################


! ********************
! **  Torsion Force **
! ********************

subroutine Force_Dihedral

use Numbers, only : N
use CommonBlocks, only : QMaster, QRigidBody, QPINPT, ForceField
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use UnitExParam, only : InvPi
use BondedParam, only : NumDihedral, DihedI, DihedJ, DihedK, DihedL,    &
&   Ts, kChi, DeltaDih, NDih, DupFlag, DupkChi, DupDeltaDih, Frc_Dihed, &
&   Vir_Dihed, Ene_Dihed, CsDelDih, DupCsDelDih, DupNDih

implicit NONE

integer :: i, j, k
integer :: l, m, mm
real(8) :: Ra2, Rb2
real(8) :: rRa2, rRb2, rRab
real(8) :: CsPhi, Phi
real(8) :: Ds, Dir
real(8) :: Deg
real(8) :: cf
real(8), dimension(3) :: Rij, Rkj, Rlj
real(8), dimension(3) :: Fi, Fj, Fk, Fl
real(8), dimension(3) :: Pla, Plb, Ori

   Ene_Dihed = 0.d0
   Frc_Dihed = 0.d0
   Vir_Dihed = 0.d0

   if( NumDihedral == 0 ) Return

   if( ForceField(1:5) == 'CHARM' ) then

   do m = 1 , NumDihedral

     i = DihedI(m)
     j = DihedJ(m)
     k = DihedK(m)
     l = DihedL(m)

     Rij = R(:,i) - R(:,j)
     Rkj = R(:,k) - R(:,j)
     Rlj = R(:,l) - R(:,j)

     Pla = VecProd(Rij,Rkj)
     Plb = VecProd(Rlj,Rkj)

#ifdef PCC
     Ra2 = Pla(1)*Pla(1) + Pla(2)*Pla(2) + Pla(3)*Pla(3)
     Rb2 = Plb(1)*Plb(1) + Plb(2)*Plb(2) + Plb(3)*Plb(3)
#else
     Ra2 = dot_product(Pla,Pla)
     Rb2 = dot_product(Plb,Plb)
#endif

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt(rRa2*rRb2)

#ifdef PCC
     CsPhi = (Pla(1)*Plb(1) + Pla(2)*Plb(2) + Pla(3)*Plb(3)) * rRab
#else
     CsPhi = dot_product(Pla,Plb) * rRab
#endif

     if(CsPhi >  1.d0) CsPhi =  1.d0
     if(CsPhi < -1.d0) CsPhi = -1.d0

     Ori = VecProd(Pla,Plb)
#ifdef PCC
     Ds  = Ori(1)*Rkj(1) + Ori(2)*Rkj(2) + Ori(3)*Rkj(3)
#else
     Ds  = dot_product(Ori,Rkj)
#endif

     Dir = 1.d0
     if(Ds < 0.d0) Dir = -1.d0

     Phi = acos(CsPhi)
     Deg = Dir * Phi * 180.d0 * InvPi
     Ts(m) = Deg

     cf = - VPh(kChi(m),CsPhi,DeltaDih(m),NDih(m))

     Fi = cf * ( VecProd(Rkj,Plb) * rRab + CsPhi * VecProd(Pla,Rkj) * rRa2 )

     Fk = cf * ( ( VecProd(Plb,Rij) + VecProd(Pla,Rlj) ) * rRab          &
     &  + CsPhi * ( VecProd(Rij,Pla) * rRa2 + VecProd(Rlj,Plb) * rRb2 )  )

     Fl = cf * ( VecProd(Rkj,Pla) * rRab + CsPhi * VecProd(Plb,Rkj) * rRb2 )

     Fj = - ( Fi + Fk + Fl )

     Frc_Dihed(:,i) = Frc_Dihed(:,i) + Fi
     Frc_Dihed(:,j) = Frc_Dihed(:,j) + Fj
     Frc_Dihed(:,k) = Frc_Dihed(:,k) + Fk
     Frc_Dihed(:,l) = Frc_Dihed(:,l) + Fl

     Ene_Dihed = Ene_Dihed + &
     &           kChi(m) * ( 1.d0 + cos( NDih(m) * Phi - DeltaDih(m) ) )

     do mm = 1, DupFlag(m)

       cf = -VPh(DupkChi(mm,m),CsPhi,DupDeltaDih(mm,m),DupNDih(mm,m))

       Fi = cf * ( VecProd(Rkj,Plb) * rRab + CsPhi * VecProd(Pla,Rkj) * rRa2 )

       Fk = cf * ( ( VecProd(Plb,Rij) + VecProd(Pla,Rlj) ) * rRab          &
       &  + CsPhi * ( VecProd(Rij,Pla) * rRa2 + VecProd(Rlj,Plb) * rRb2 )  )

       Fl = cf * ( VecProd(Rkj,Pla) * rRab + CsPhi * VecProd(Plb,Rkj) * rRb2 )

       Fj = - ( Fi + Fk + Fl )

       Frc_Dihed(:,i) = Frc_Dihed(:,i) + Fi
       Frc_Dihed(:,j) = Frc_Dihed(:,j) + Fj
       Frc_Dihed(:,k) = Frc_Dihed(:,k) + Fk
       Frc_Dihed(:,l) = Frc_Dihed(:,l) + Fl

       Ene_Dihed = Ene_Dihed + DupkChi(mm,m) &
       &           * ( 1.d0 + cos( DupNDih(mm,m) * Phi - DupDeltaDih(mm,m) ) )

     end do

   end do

   else if( ForceField(1:4) == 'OPLS' ) then

   do m = 1 , NumDihedral

     i = DihedI(m)
     j = DihedJ(m)
     k = DihedK(m)
     l = DihedL(m)

     Rij = R(:,i) - R(:,j)
     Rkj = R(:,k) - R(:,j)
     Rlj = R(:,l) - R(:,j)

     Pla = VecProd(Rij,Rkj)
     Plb = VecProd(Rlj,Rkj)

#ifdef PCC
     Ra2 = Pla(1)*Pla(1) + Pla(2)*Pla(2) + Pla(3)*Pla(3)
     Rb2 = Plb(1)*Plb(1) + Plb(2)*Plb(2) + Plb(3)*Plb(3)
#else
     Ra2 = dot_product(Pla,Pla)
     Rb2 = dot_product(Plb,Plb)
#endif

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt(rRa2*rRb2)

#ifdef PCC
     CsPhi = (Pla(1)*Plb(1) + Pla(2)*Plb(2) + Pla(3)*Plb(3)) * rRab
#else
     CsPhi = dot_product(Pla,Plb) * rRab
#endif

     if(CsPhi >  1.d0) CsPhi =  1.d0
     if(CsPhi < -1.d0) CsPhi = -1.d0

     Ori = VecProd(Pla,Plb)
#ifdef PCC
     Ds  = Ori(1)*Rkj(1) + Ori(2)*Rkj(2) + Ori(3)*Rkj(3)
#else
     Ds  = dot_product(Ori,Rkj)
#endif

     Dir = 1.d0
     if(Ds < 0.d0) Dir = -1.d0

     Phi = acos(CsPhi)
     Deg = Dir * Phi * 180.d0 * InvPi
     Ts(m) = Deg

     cf = - VPo(kChi(m),CsPhi,CsDelDih(m),NDih(m))

     Fi = cf * ( VecProd(Rkj,Plb) * rRab + CsPhi * VecProd(Pla,Rkj) * rRa2 )

     Fk = cf * ( ( VecProd(Plb,Rij) + VecProd(Pla,Rlj) ) * rRab          &
     &  + CsPhi * ( VecProd(Rij,Pla) * rRa2 + VecProd(Rlj,Plb) * rRb2 )  )

     Fl = cf * ( VecProd(Rkj,Pla) * rRab + CsPhi * VecProd(Plb,Rkj) * rRb2 )

     Fj = - ( Fi + Fk + Fl )

     Frc_Dihed(:,i) = Frc_Dihed(:,i) + Fi
     Frc_Dihed(:,j) = Frc_Dihed(:,j) + Fj
     Frc_Dihed(:,k) = Frc_Dihed(:,k) + Fk
     Frc_Dihed(:,l) = Frc_Dihed(:,l) + Fl

     Ene_Dihed = Ene_Dihed + UPo(kChi(m),CsPhi,CsDelDih(m),NDih(m))

     do mm = 1, DupFlag(m)

       cf = -VPo(DupkChi(mm,m),CsPhi,DupCsDelDih(mm,m),DupNDih(mm,m))

       Fi = cf * ( VecProd(Rkj,Plb) * rRab + CsPhi * VecProd(Pla,Rkj) * rRa2 )

       Fk = cf * ( ( VecProd(Plb,Rij) + VecProd(Pla,Rlj) ) * rRab          &
       &  + CsPhi * ( VecProd(Rij,Pla) * rRa2 + VecProd(Rlj,Plb) * rRb2 )  )

       Fl = cf * ( VecProd(Rkj,Pla) * rRab + CsPhi * VecProd(Plb,Rkj) * rRb2 )

       Fj = - ( Fi + Fk + Fl )

       Frc_Dihed(:,i) = Frc_Dihed(:,i) + Fi
       Frc_Dihed(:,j) = Frc_Dihed(:,j) + Fj
       Frc_Dihed(:,k) = Frc_Dihed(:,k) + Fk
       Frc_Dihed(:,l) = Frc_Dihed(:,l) + Fl

       Ene_Dihed = Ene_Dihed + &
       &           UPo(DupkChi(mm,m),CsPhi,DupCsDelDih(mm,m),DupNDih(mm,m))

     end do

   end do

   end if

   if(QRigidBody) then

     do i = 1, N

       j = AtomUnitNum(i)

       Vir_Dihed(1,1) = Vir_Dihed(1,1) + Frc_Dihed(1,i) * R_RB(1,j)
       Vir_Dihed(1,2) = Vir_Dihed(1,2) + Frc_Dihed(1,i) * R_RB(2,j)
       Vir_Dihed(1,3) = Vir_Dihed(1,3) + Frc_Dihed(1,i) * R_RB(3,j)
       Vir_Dihed(2,1) = Vir_Dihed(2,1) + Frc_Dihed(2,i) * R_RB(1,j)
       Vir_Dihed(2,2) = Vir_Dihed(2,2) + Frc_Dihed(2,i) * R_RB(2,j)
       Vir_Dihed(2,3) = Vir_Dihed(2,3) + Frc_Dihed(2,i) * R_RB(3,j)
       Vir_Dihed(3,1) = Vir_Dihed(3,1) + Frc_Dihed(3,i) * R_RB(1,j)
       Vir_Dihed(3,2) = Vir_Dihed(3,2) + Frc_Dihed(3,i) * R_RB(2,j)
       Vir_Dihed(3,3) = Vir_Dihed(3,3) + Frc_Dihed(3,i) * R_RB(3,j)

     end do

   else if(QPINPT) then

     do i = 1, N

       Vir_Dihed(1,1) = Vir_Dihed(1,1) + Frc_Dihed(1,i) * Rcentroid(1,i)
       Vir_Dihed(1,2) = Vir_Dihed(1,2) + Frc_Dihed(1,i) * Rcentroid(2,i)
       Vir_Dihed(1,3) = Vir_Dihed(1,3) + Frc_Dihed(1,i) * Rcentroid(3,i)
       Vir_Dihed(2,1) = Vir_Dihed(2,1) + Frc_Dihed(2,i) * Rcentroid(1,i)
       Vir_Dihed(2,2) = Vir_Dihed(2,2) + Frc_Dihed(2,i) * Rcentroid(2,i)
       Vir_Dihed(2,3) = Vir_Dihed(2,3) + Frc_Dihed(2,i) * Rcentroid(3,i)
       Vir_Dihed(3,1) = Vir_Dihed(3,1) + Frc_Dihed(3,i) * Rcentroid(1,i)
       Vir_Dihed(3,2) = Vir_Dihed(3,2) + Frc_Dihed(3,i) * Rcentroid(2,i)
       Vir_Dihed(3,3) = Vir_Dihed(3,3) + Frc_Dihed(3,i) * Rcentroid(3,i)

     end do

   else

     do i = 1, N

       Vir_Dihed(1,1) = Vir_Dihed(1,1) + Frc_Dihed(1,i) * R(1,i)
       Vir_Dihed(1,2) = Vir_Dihed(1,2) + Frc_Dihed(1,i) * R(2,i)
       Vir_Dihed(1,3) = Vir_Dihed(1,3) + Frc_Dihed(1,i) * R(3,i)
       Vir_Dihed(2,2) = Vir_Dihed(2,2) + Frc_Dihed(2,i) * R(2,i)
       Vir_Dihed(2,3) = Vir_Dihed(2,3) + Frc_Dihed(2,i) * R(3,i)
       Vir_Dihed(3,3) = Vir_Dihed(3,3) + Frc_Dihed(3,i) * R(3,i)

     end do

     Vir_Dihed(2,1) = Vir_Dihed(1,2)
     Vir_Dihed(3,1) = Vir_Dihed(1,3)
     Vir_Dihed(3,2) = Vir_Dihed(2,3)

   end if

! #############

Contains

! ###############################################

   function VecProd(x,y) Result(z)

   real(8), dimension(3) :: x, y, z

     z(1) = x(2) * y(3) - y(2) * x(3)
     z(2) = x(3) * y(1) - y(3) * x(1)
     z(3) = x(1) * y(2) - y(1) * x(2)

   end function VecProd

! ###############################################

   function VPh(k,csp,del,nn) Result(x)

   integer :: nn
   real(8) :: x, csp, del, k

     if(nn==1) then

       x = k * cos(del)

     else if(nn==2) then

       x = 4.d0 * k * csp * cos(del)

     else if(nn==3) then

       x = 3.d0 * k * (4.d0 * csp * csp - 1.d0 ) * cos(del)

     else if(nn==4) then

       x = 16.d0 * k * csp * ( 2.d0 * csp * csp - 1.d0 ) * cos(del)

     else if(nn==5) then

       x = 5.d0 * k * ( 8.d0 * csp * csp * ( csp * csp + csp - 1.d0 ) &
       &   - 4.d0 * csp + 1.d0 ) * cos(del)

     else if(nn==6) then

       x = 12.d0 * k * (4.d0 * csp * csp - 1.d0 ) * csp * &
       &   ( 4.d0 * csp * csp - 3.d0 ) * cos(del)

     else

       if(QMaster) write(*,*) 'ERROR : Dihedral Calc.'
       call Finalize

     end if

   end function VPh

! ###############################################

   function VPo(k,csp,csdel,nn) Result(x)

   integer :: nn
   real(8) :: x, csp, csdel, k, csp2

     if(nn==1) then

       x = 0.5d0 * k * csdel

     else if(nn==2) then

       x = 2.d0 * k * csp * csdel

     else if(nn==3) then

       x = 3.d0 * k * (2.d0 * csp * csp - 0.5d0 ) * csdel

     else if(nn==4) then

       x = 8.d0 * k * (2.d0 * csp * csp - 1.d0 ) * csp * csdel

     else if(nn==6) then

       csp2 = csp * csp
       x = k * ( 96.d0 * ( csp2 - 1.d0 ) * csp2 + 18.d0 ) * csp * csdel

     else

       if(QMaster) write(*,*) 'ERROR : Dihedral Calc.'
       call Finalize

     end if

   end function VPo

! ###############################################

   function UPo(k,csp,csdel,nn) Result(x)

   integer :: nn
   real(8) :: x, csp, csdel, k, csp2

     if(nn==1) then

       x = 0.5d0 * k * ( 1.d0 + csp * csdel )

     else if(nn==2) then

       x = 0.5d0 * k * ( 1.d0 + (2.d0 * csp * csp - 1.d0 ) * csdel )

     else if(nn==3) then

       x = 0.5d0 * k * ( 1.d0 + (4.d0 * csp * csp - 3.d0 ) * csp * csdel )

     else if(nn==4) then

       csp2 = csp * csp
       x = 0.5d0 * k * ( 1.d0 + ( (8.d0 * csp2 - 8.d0 ) * csp2 + 1.d0 ) * csdel )

     else if(nn==6) then

       csp2 = csp * csp
       x = 0.5d0 * k * ( 1.d0 + ( ( (32.d0 * csp2 - 48.d0 ) * csp2 + 18.d0 ) &
       &             * csp2 - 1.d0 ) * csdel )

     else

       if(QMaster) write(*,*) 'ERROR : Dihedral Calc.'
       call Finalize

     end if

   end function UPo

! ###############################################

end subroutine Force_Dihedral


!######################################################################
!######################################################################


! *********************
! **  Improper Force **
! *********************

subroutine Force_Improper

use Numbers, only : N
use CommonBlocks, only : QMaster, QRigidBody, QPINPT, ForceField
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use UnitExParam, only : InvPi
use BondedParam, only : NumImproper, ImproI, ImproJ, ImproK, ImproL, &
&   kPsi, PsiImp, kImp, CsDelImp, NImp, Frc_Impro, Vir_Impro, Ene_Impro

implicit NONE

real(8), parameter :: k3 = - 1.d0 /      6.d0
real(8), parameter :: k5 =   1.d0 /    120.d0
real(8), parameter :: k7 = - 1.d0 /   5040.d0
real(8), parameter :: k9 =   1.d0 / 362880.d0

integer :: i, j, k, l, m
real(8) :: Ra2, Rb2
real(8) :: rRa2, rRb2, rRab
real(8) :: CsPsi, Ds
real(8) :: Dir
real(8) :: Deg
real(8) :: Psi, dPsi
real(8) :: cf
real(8) :: Psi2, Psi3, Psi5
real(8) :: Psi7, Psi9
! real(8) :: Psi4, Psi6, Psi8
real(8), dimension(3) :: Rij, Rkj, Rlj
real(8), dimension(3) :: Rik, Rjk, Rlk
real(8), dimension(3) :: Fi, Fj, Fk, Fl
real(8), dimension(3) :: Ori
real(8), dimension(3) :: Pla, Plb

   Ene_Impro = 0.d0
   Frc_Impro = 0.d0
   Vir_Impro = 0.d0

   if( NumImproper == 0 ) Return

   if( ForceField(1:5) == 'CHARM' ) then

   do m = 1 , NumImproper
                            !       (l)
     i = ImproI(m)          !        |
     j = ImproJ(m)          !       (i)
     k = ImproK(m)          !       / \
     l = ImproL(m)          !    (j)   (k)

     Rij = R(:,i) - R(:,j)
     Rkj = R(:,k) - R(:,j)
     Rlj = R(:,l) - R(:,j)

     Pla = VecProd(Rij,Rkj)  ! Plane (ijk)
     Plb = VecProd(Rlj,Rkj)  ! Plane (jkl)

#ifdef PCC
     Ra2 = Pla(1)*Pla(1) + Pla(2)*Pla(2) + Pla(3)*Pla(3)
     Rb2 = Plb(1)*Plb(1) + Plb(2)*Plb(2) + Plb(3)*Plb(3)
#else
     Ra2 = dot_product(Pla,Pla)
     Rb2 = dot_product(Plb,Plb)
#endif

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt(rRa2*rRb2)

#ifdef PCC
     CsPsi = (Pla(1)*Plb(1) + Pla(2)*Plb(2) + Pla(3)*Plb(3)) * rRab
#else
     CsPsi = dot_product(Pla,Plb) * rRab
#endif
     if((CsPsi >= 1.d0).or.(CsPsi <= -1.d0)) cycle

     Ori = VecProd(Pla,Plb)
#ifdef PCC
     Ds  = Ori(1)*Rkj(1) + Ori(2)*Rkj(2) + Ori(3)*Rkj(3)
#else
     Ds  = dot_product(Ori,Rkj)
#endif

     Dir = 1.d0
     if(Ds < 0.d0) Dir = -1.d0

     Psi = acos(CsPsi)
     Deg = Dir * Psi * 180.d0 * InvPi

     dPsi = Psi - PsiImp(m)
     Psi2 = Psi*Psi
     Psi3 = Psi2*Psi
     Psi5 = Psi3*Psi2
     Psi7 = Psi5*Psi2
     Psi9 = Psi7*Psi2

     cf = 2.d0 * kPsi(m) * dPsi / &
     &    ( Psi + k3*Psi3 + k5*Psi5 + k7*Psi7 + k9*Psi9 )

!     Psi4 = Psi2*Psi2
!     Psi6 = Psi4*Psi2
!     Psi8 = Psi4*Psi4
!     cf = 2.d0 * kPsi(m) / ( 1.d0 + k3*Psi2 + k5*Psi4 + k7*Psi6 + k9*Psi8 )

     Fi = cf * ( VecProd(Rkj,Plb) * rRab + CsPsi * VecProd(Pla,Rkj) * rRa2 )

     Fk = cf * ( ( VecProd(Plb,Rij) + VecProd(Pla,Rlj) ) * rRab          &
     &  + CsPsi * ( VecProd(Rij,Pla) * rRa2 + VecProd(Rlj,Plb) * rRb2 )  )

     Fl = cf * ( VecProd(Rkj,Pla) * rRab + CsPsi * VecProd(Plb,Rkj) * rRb2 )

     Fj = - ( Fi + Fk + Fl )

     Frc_Impro(:,i) = Frc_Impro(:,i) + Fi
     Frc_Impro(:,j) = Frc_Impro(:,j) + Fj
     Frc_Impro(:,k) = Frc_Impro(:,k) + Fk
     Frc_Impro(:,l) = Frc_Impro(:,l) + Fl

     Ene_Impro = Ene_Impro + kPsi(m) * dPsi * dPsi

   end do

   else if(ForceField(1:4) == 'OPLS') then

   do m = 1 , NumImproper
                            !       (i)
     i = ImproI(m)          !        |
     j = ImproJ(m)          !       (k)
     k = ImproK(m)          !       / \
     l = ImproL(m)          !    (j)   (l)

     Rik = R(:,i) - R(:,k)  ! a
     Rjk = R(:,j) - R(:,k)  ! b
     Rlk = R(:,l) - R(:,k)  ! c

     Pla = VecProd(Rik,Rlk)  ! Plane (ikl) a*c
     Plb = VecProd(Rjk,Rlk)  ! Plane (jkl) b*c

#ifdef PCC
     Ra2 = Pla(1)*Pla(1) + Pla(2)*Pla(2) + Pla(3)*Pla(3)
     Rb2 = Plb(1)*Plb(1) + Plb(2)*Plb(2) + Plb(3)*Plb(3)
#else
     Ra2 = dot_product(Pla,Pla)
     Rb2 = dot_product(Plb,Plb)
#endif

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt(rRa2*rRb2)

#ifdef PCC
     CsPsi = (Pla(1)*Plb(1) + Pla(2)*Plb(2) + Pla(3)*Plb(3)) * rRab
#else
     CsPsi = dot_product(Pla,Plb) * rRab
#endif

     if(CsPsi >  1.d0) CsPsi =  1.d0
     if(CsPsi < -1.d0) CsPsi = -1.d0

!     Ori = VecProd(Pla,Plb)
!     Ds  = dot_product(Ori,Rlk)

!     Dir = 1.d0
!     if(Ds < 0.d0) Dir = -1.d0

!     Psi = acos(CsPsi)
!     Deg = Dir * Psi * 180.d0 * InvPi

     cf = - VPo(kImp(m),CsPsi,CsDelImp(m),NImp(m))

     Fi = cf * ( VecProd(Rlk,Plb) * rRab + CsPsi * VecProd(Pla,Rlk) * rRa2 )

     Fl = cf * ( ( VecProd(Plb,Rik) + VecProd(Pla,Rjk) ) * rRab          &
     &  + CsPsi * ( VecProd(Rik,Pla) * rRa2 + VecProd(Rjk,Plb) * rRb2 )  )

     Fj = cf * ( VecProd(Rlk,Pla) * rRab + CsPsi * VecProd(Plb,Rlk) * rRb2 )

     Fk = - ( Fi + Fj + Fl )

     Frc_Impro(:,i) = Frc_Impro(:,i) + Fi
     Frc_Impro(:,j) = Frc_Impro(:,j) + Fj
     Frc_Impro(:,k) = Frc_Impro(:,k) + Fk
     Frc_Impro(:,l) = Frc_Impro(:,l) + Fl

     Ene_Impro = Ene_Impro + Upo(kImp(m),CsPsi,CsDelImp(m),NImp(m))

   end do

   end if

   if(QRigidBody) then

     do i = 1, N

       j = AtomUnitNum(i)

       Vir_Impro(1,1) = Vir_Impro(1,1) + Frc_Impro(1,i) * R_RB(1,j)
       Vir_Impro(1,2) = Vir_Impro(1,2) + Frc_Impro(1,i) * R_RB(2,j)
       Vir_Impro(1,3) = Vir_Impro(1,3) + Frc_Impro(1,i) * R_RB(3,j)
       Vir_Impro(2,1) = Vir_Impro(2,1) + Frc_Impro(2,i) * R_RB(1,j)
       Vir_Impro(2,2) = Vir_Impro(2,2) + Frc_Impro(2,i) * R_RB(2,j)
       Vir_Impro(2,3) = Vir_Impro(2,3) + Frc_Impro(2,i) * R_RB(3,j)
       Vir_Impro(3,1) = Vir_Impro(3,1) + Frc_Impro(3,i) * R_RB(1,j)
       Vir_Impro(3,2) = Vir_Impro(3,2) + Frc_Impro(3,i) * R_RB(2,j)
       Vir_Impro(3,3) = Vir_Impro(3,3) + Frc_Impro(3,i) * R_RB(3,j)

     end do

   else if(QPINPT) then

     do i = 1, N

       Vir_Impro(1,1) = Vir_Impro(1,1) + Frc_Impro(1,i) * Rcentroid(1,i)
       Vir_Impro(1,2) = Vir_Impro(1,2) + Frc_Impro(1,i) * Rcentroid(2,i)
       Vir_Impro(1,3) = Vir_Impro(1,3) + Frc_Impro(1,i) * Rcentroid(3,i)
       Vir_Impro(2,1) = Vir_Impro(2,1) + Frc_Impro(2,i) * Rcentroid(1,i)
       Vir_Impro(2,2) = Vir_Impro(2,2) + Frc_Impro(2,i) * Rcentroid(2,i)
       Vir_Impro(2,3) = Vir_Impro(2,3) + Frc_Impro(2,i) * Rcentroid(3,i)
       Vir_Impro(3,1) = Vir_Impro(3,1) + Frc_Impro(3,i) * Rcentroid(1,i)
       Vir_Impro(3,2) = Vir_Impro(3,2) + Frc_Impro(3,i) * Rcentroid(2,i)
       Vir_Impro(3,3) = Vir_Impro(3,3) + Frc_Impro(3,i) * Rcentroid(3,i)

     end do

   else

     do i = 1, N

       Vir_Impro(1,1) = Vir_Impro(1,1) + Frc_Impro(1,i) * R(1,i)
       Vir_Impro(1,2) = Vir_Impro(1,2) + Frc_Impro(1,i) * R(2,i)
       Vir_Impro(1,3) = Vir_Impro(1,3) + Frc_Impro(1,i) * R(3,i)
       Vir_Impro(2,2) = Vir_Impro(2,2) + Frc_Impro(2,i) * R(2,i)
       Vir_Impro(2,3) = Vir_Impro(2,3) + Frc_Impro(2,i) * R(3,i)
       Vir_Impro(3,3) = Vir_Impro(3,3) + Frc_Impro(3,i) * R(3,i)

     end do

     Vir_Impro(2,1) = Vir_Impro(1,2)
     Vir_Impro(3,1) = Vir_Impro(1,3)
     Vir_Impro(3,2) = Vir_Impro(2,3)

   end if

! #############

Contains

   function VecProd(x,y) Result(z)

   real(8), dimension(3) :: x, y, z

     z(1) = x(2) * y(3) - y(2) * x(3)
     z(2) = x(3) * y(1) - y(3) * x(1)
     z(3) = x(1) * y(2) - y(1) * x(2)

   end function VecProd

! ###############################################

   function VPo(k,csp,csdel,nn) Result(x)

   integer :: nn
   real(8) :: x, csp, csdel, k

     if(nn==1) then

       x = 0.5d0 * k * csdel

     else if(nn==2) then

       x = 2.d0 * k * csp * csdel

     else if(nn==3) then

       x = 3.d0 * k * (2.d0 * csp * csp - 0.5d0 ) * csdel

     else

       if(QMaster) write(*,*) 'ERROR : Dihedral Calc.'
       call Finalize

     end if

   end function VPo

! ###############################################

   function UPo(k,csp,csdel,nn) Result(x)

   integer :: nn
   real(8) :: x, csp, csdel, k

     if(nn==1) then

       x = 0.5d0 * k * ( 1.d0 + csp * csdel )

     else if(nn==2) then

       x = 0.5d0 * k * ( 1.d0 + (2.d0 * csp * csp - 1.d0 ) * csdel )

     else if(nn==3) then

       x = 0.5d0 * k * ( 1.d0 + (4.d0 * csp * csp - 3.d0 ) * csp * csdel )

     else

       if(QMaster) write(*,*) 'ERROR : Dihedral Calc.'
       call Finalize

     end if

   end function UPo

! ###############################################

end subroutine Force_Improper
