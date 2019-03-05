! ############################
! ## SUBROUTINE LIST 
! ## -- Force_Bond_CG 
! ## -- Force_Angle_CG 
! ## -- Force_Dihed_CG
! ## -- Force_Improper_CG
! ############################


!######################################################################
!######################################################################


! ***************************
! ** Bond Stretching Force **
! ***************************

subroutine Force_Bond_CG

use Numbers, only : N
use CommonBlocks, only : QRigidBody
use Configuration, only : R
use RBparam
use BondedParam, only : NumBond, BondI, BondJ, kBond, rBond, FTypeBond, &
&   Frc_Bond, Vir_Bond, Ene_Bond

implicit NONE

integer :: i, j, k
real(8) :: R2, R1, dR, Fc, pref
real(8) :: Rx, Ry, Rz, Fx, Fy, Fz

   Ene_Bond = 0.d0
   Frc_Bond = 0.d0
   Vir_Bond = 0.d0

   if( NumBond == 0 ) Return

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     Rx = R(1,i) - R(1,j)
     Ry = R(2,i) - R(2,j)
     Rz = R(3,i) - R(3,j)

     R2  = Rx*Rx + Ry*Ry + Rz*Rz

     if (FTypeBond(k) == 1) then

       R1 = sqrt(R2)
       dR = R1 - rBond(k)
       Fc = kBond(k) * dR

       Ene_Bond = Ene_Bond + Fc * dR

       pref = - 2.d0 * Fc / R1

       Fx = pref * Rx
       Fy = pref * Ry
       Fz = pref * Rz

     else

       write(*,*) 'error: bond type in (Force_Bond_CG)'
       call Finalize

     end if

     Frc_Bond(1,i) = Frc_Bond(1,i) + Fx
     Frc_Bond(2,i) = Frc_Bond(2,i) + Fy
     Frc_Bond(3,i) = Frc_Bond(3,i) + Fz
     Frc_Bond(1,j) = Frc_Bond(1,j) - Fx
     Frc_Bond(2,j) = Frc_Bond(2,j) - Fy
     Frc_Bond(3,j) = Frc_Bond(3,j) - Fz

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

end subroutine Force_Bond_CG


!######################################################################
!######################################################################


! ********************
! **  Bending Force **
! ********************

subroutine Force_Angle_CG

use Numbers, only : N
use CommonBlocks, only : QRigidBody
use Configuration, only : R
use RBparam
use UnitExParam, only : pi
use BondedParam, only : NumAngle, AngleI, AngleJ, AngleK, &
&   kTheta, Theta0, FTypeAngle, Frc_Angle, Vir_Angle, Ene_Angle

implicit NONE

integer :: i, j, k, l
real(8) :: Ra2, Rb2
real(8) :: rRa2, rRb2, rRab
real(8) :: Cst, Theta, Thpi, Thpi2
real(8) :: cf, ek, dT, Snt, t1
real(8) :: Rijx, Rkjx
real(8) :: Rijy, Rkjy
real(8) :: Rijz, Rkjz
real(8) :: Fix, Fjx, Fkx
real(8) :: Fiy, Fjy, Fky
real(8) :: Fiz, Fjz, Fkz
real(8) :: pref, pref1, pref2
real(8) :: poly_quart
external poly_quart

   Ene_Angle = 0.d0
   Frc_Angle = 0.d0
   Vir_Angle = 0.d0

   if( NumAngle == 0 ) Return

   do l = 1 , NumAngle

     i = AngleI(l)
     j = AngleJ(l)
     k = AngleK(l)

     Rijx = R(1,i) - R(1,j)
     Rijy = R(2,i) - R(2,j)
     Rijz = R(3,i) - R(3,j)
     Rkjx = R(1,k) - R(1,j)
     Rkjy = R(2,k) - R(2,j)
     Rkjz = R(3,k) - R(3,j)

     Ra2  = Rijx * Rijx + Rijy * Rijy + Rijz * Rijz
     Rb2  = Rkjx * Rkjx + Rkjy * Rkjy + Rkjz * Rkjz

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt( rRa2 * rRb2 )

     Cst = ( Rijx * Rkjx + Rijy * Rkjy + Rijz * Rkjz ) * rRab

     if(FTypeAngle(l) == 1) then ! cosine function

       cf = - kTheta(l)
       ek = kTheta(l) * ( 1.d0 + Cst )

     else if(FTypeAngle(l) == 2) then ! quartic function 

       if(Cst < -1.) Cst  = -1.d0

       Theta = acos(Cst)
       Thpi  = Theta - pi
       Snt   = sin(Theta)
       if( Snt > 1.d0 ) then
         Snt = 1.d0
       else if( Snt < -1.d0 ) then
         Snt = -1.d0
       end if

       if(abs(Thpi) < 0.01) then
         t1 = -1.d0 * poly_quart(Snt, 20)
       else
         t1 = -Thpi/Snt
       end if

       Thpi2 = Theta0(l) - Thpi * Thpi

       cf = kTheta(l) * Thpi2 * t1
       ek = 0.25d0 * kTheta(l) * Thpi2 * Thpi2

     else if(FTypeAngle(l) == 3) then

       if(Cst <= -1.d0) then

         dT = pi - Theta0(l)
         ek = kTheta(l) * dT * dT
         cf = 0.d0

       else

         Theta = acos(Cst)

         dT = Theta - Theta0(l)
         cf = 2.d0 * kTheta(l) * dT / sin(Theta)
         ek = kTheta(l) * dT * dT

       end if

     end if

     pref  = cf * rRab
     pref1 = cf * Cst * rRa2
     pref2 = cf * Cst * rRb2

     Fix = pref * Rkjx - pref1 * Rijx
     Fiy = pref * Rkjy - pref1 * Rijy
     Fiz = pref * Rkjz - pref1 * Rijz
     Fkx = pref * Rijx - pref2 * Rkjx
     Fky = pref * Rijy - pref2 * Rkjy
     Fkz = pref * Rijz - pref2 * Rkjz
     Fjx = - ( Fix + Fkx )
     Fjy = - ( Fiy + Fky )
     Fjz = - ( Fiz + Fkz )

     Frc_Angle(1,i) = Frc_Angle(1,i) + Fix
     Frc_Angle(2,i) = Frc_Angle(2,i) + Fiy
     Frc_Angle(3,i) = Frc_Angle(3,i) + Fiz
     Frc_Angle(1,j) = Frc_Angle(1,j) + Fjx
     Frc_Angle(2,j) = Frc_Angle(2,j) + Fjy
     Frc_Angle(3,j) = Frc_Angle(3,j) + Fjz
     Frc_Angle(1,k) = Frc_Angle(1,k) + Fkx
     Frc_Angle(2,k) = Frc_Angle(2,k) + Fky
     Frc_Angle(3,k) = Frc_Angle(3,k) + Fkz

     Ene_Angle = Ene_Angle + ek

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

end subroutine Force_Angle_CG

function fact_quart(m)

integer :: i, m
real(8) :: fact_quart

   fact_quart = 1.d0
   i = 2
   do while(i <= m)
     fact_quart = fact_quart * i
     i = i + 1
   end do

end function fact_quart

function poly_quart( s_theta, order )

integer :: order, i
real(8) :: s_theta
real(8) :: poly_quart, fqi, fact_quart

   poly_quart = 1.d0
   do i = 1, order
     fqi = fact_quart(i)
     poly_quart = poly_quart + fact_quart(2*i) / (2**(2*i)*fqi*fqi*(2*i+1)) &
     &          * s_theta**(2*i)
   end do

end function poly_quart


!######################################################################
!######################################################################


! *********************
! **  Torsion Force **
! *********************

subroutine Force_Dihed_CG

use Numbers, only : N
use CommonBlocks, only : QMaster, QRigidBody
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use UnitExParam, only : InvPi
use BondedParam, only : NumDihedral, DihedI, DihedJ, DihedK, DihedL, &
&   FTypeDihed,   &
&   kChi, DeltaDih, NDih, DupFlag, DupkChi, DupDeltaDih, Frc_Dihed, &
&   Vir_Dihed, Ene_Dihed, CsDelDih, DupCsDelDih, DupNDih

implicit NONE

integer :: i, j, k
integer :: l, m, mm
real(8) :: Ra2, Rb2
real(8) :: rRa2, rRb2, rRab
real(8) :: CsPhi, Phi
real(8) :: Ds, Dir
real(8) :: Deg
real(8) :: cf, SnPhi, xx
real(8), dimension(3) :: Rij, Rkj, Rlj
real(8), dimension(3) :: Fi, Fj, Fk, Fl
real(8), dimension(3) :: Pla, Plb, Ori

   Ene_Dihed = 0.d0
   Frc_Dihed = 0.d0
   Vir_Dihed = 0.d0

   if( NumDihedral == 0 ) Return

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
     Phi = Dir * Phi
     SnPhi = sin(Phi)

     if(FTypeDihed(m)==1) then
       cf = - VPh(kChi(m),CsPhi,SnPhi,DeltaDih(m),NDih(m))
       Ene_Dihed = Ene_Dihed + kChi(m) * ( 1.d0 + cos( NDih(m) * Phi - DeltaDih(m) ) )
     else if(FTypeDihed(m)==2) then
       xx = Phi - DeltaDih(m)
       cf = 2.d0 * kChi(m) * xx / SnPhi
       Ene_Dihed = Ene_Dihed + kChi(m) * xx * xx
     end if

     Fi = cf * ( VecProd(Rkj,Plb) * rRab + CsPhi * VecProd(Pla,Rkj) * rRa2 )

     Fk = cf * ( ( VecProd(Plb,Rij) + VecProd(Pla,Rlj) ) * rRab          &
     &  + CsPhi * ( VecProd(Rij,Pla) * rRa2 + VecProd(Rlj,Plb) * rRb2 )  )

     Fl = cf * ( VecProd(Rkj,Pla) * rRab + CsPhi * VecProd(Plb,Rkj) * rRb2 )

     Fj = - ( Fi + Fk + Fl )

     Frc_Dihed(:,i) = Frc_Dihed(:,i) + Fi
     Frc_Dihed(:,j) = Frc_Dihed(:,j) + Fj
     Frc_Dihed(:,k) = Frc_Dihed(:,k) + Fk
     Frc_Dihed(:,l) = Frc_Dihed(:,l) + Fl

     do mm = 1, DupFlag(m)

       cf = -VPh(DupkChi(mm,m),CsPhi,SnPhi,DupDeltaDih(mm,m),DupNDih(mm,m))

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

   function VPh(k,csp,snp,del,nn) Result(x)

   integer :: nn
   real(8) :: x, csp, snp, del, k

     if(nn==1) then

       x = k * (cos(del) - csp/snp*sin(del))

     else if(nn==2) then

       x = k * (4.d0* csp * cos(del) + 2.d0*sin(del)*(snp-csp*csp/snp))

     else if(nn==3) then

       x = k*((4.d0*csp*csp-1.d0)*3.d0*cos(del)-3.d0*sin(del)*(csp/snp-snp*csp))

     else if(nn==4) then

       x = k*(16.d0*csp*(2.d0*csp*csp-1.d0)*cos(del)-&
&          2.d0*sin(del)*(csp*csp/snp*(2.d0*csp*csp-1)-snp*(6.d0*csp*csp-1.d0)))

     else

       if(QMaster) write(*,*) 'ERROR : Dihedral Calc.'
       call Finalize

     end if

   end function VPh


end subroutine Force_Dihed_CG

!######################################################################
!######################################################################


! *********************
! **  Improper Force **
! *********************

subroutine Force_Improper_CG

use Numbers, only : N
use CommonBlocks, only : QMaster, QRigidBody, QPINPT, ForceField
use Configuration, only : R
use RBparam
use CommonPI, only: Rcentroid
use UnitExParam, only : InvPi
use BondedParam, only : NumImproper, ImproI, ImproJ, ImproK, ImproL, &
&   FTypeImpro, kImp, NImp, DeltaImp, Frc_Impro, Vir_Impro, Ene_Impro

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
real(8) :: Rix, Riy, Riz
real(8) :: Rjx, Rjy, Rjz
real(8) :: Rkx, Rky, Rkz
real(8) :: Rlx, Rly, Rlz
real(8) :: Rijx, Rijy, Rijz
real(8) :: Rkjx, Rkjy, Rkjz
real(8) :: Rljx, Rljy, Rljz
real(8) :: Fix, Fjx, Fkx, Flx
real(8) :: Fiy, Fjy, Fky, Fly
real(8) :: Fiz, Fjz, Fkz, Flz
real(8) :: Orix, Oriy, Oriz
real(8) :: Plax, Play, Plaz
real(8) :: Plbx, Plby, Plbz

   Ene_Impro = 0.d0
   Frc_Impro = 0.d0
   Vir_Impro = 0.d0

   if( NumImproper == 0 ) Return

   do m = 1 , NumImproper
                            !       (l)
     i = ImproI(m)          !        |
     j = ImproJ(m)          !       (i)
     k = ImproK(m)          !       / \
     l = ImproL(m)          !    (j)   (k)

     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)
     Rjx = R(1,j)
     Rjy = R(2,j)
     Rjz = R(3,j)
     Rkx = R(1,k)
     Rky = R(2,k)
     Rkz = R(3,k)
     Rlx = R(1,l)
     Rly = R(2,l)
     Rlz = R(3,l)

     Rijx = Rix - Rjx
     Rijy = Riy - Rjy
     Rijz = Riz - Rjz
     Rkjx = Rkx - Rjx
     Rkjy = Rky - Rjy
     Rkjz = Rkz - Rjz
     Rljx = Rlx - Rjx
     Rljy = Rly - Rjy
     Rljz = Rlz - Rjz

     Plax = Rijy * Rkjz - Rkjy * Rijz
     Play = Rijz * Rkjx - Rkjz * Rijx
     Plaz = Rijx * Rkjy - Rkjx * Rijy
     Plbx = Rljy * Rkjz - Rkjy * Rljz
     Plby = Rljz * Rkjx - Rkjz * Rljx
     Plbz = Rljx * Rkjy - Rkjx * Rljy

     Ra2 = Plax*Plax + Play*Play + Plaz*Plaz
     Rb2 = Plbx*Plbx + Plby*Plby + Plbz*Plbz

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt(rRa2*rRb2)

     CsPsi = (Plax*Plbx + Play*Plby + Plaz*Plbz) * rRab
     if((CsPsi >= 1.d0).or.(CsPsi <= -1.d0)) cycle

     Orix = Play*Plbz - Plby*Plaz
     Oriy = Plaz*Plbx - Plbz*Plax
     Oriz = Plax*Plby - Plbx*Play

     Ds  = Orix*Rkjx + Oriy*Rkjy + Oriz*Rkjz

     Dir = 1.d0
     if(Ds < 0.d0) Dir = -1.d0

     Psi = acos(CsPsi)
     Deg = Dir * Psi * 180.d0 * InvPi

     dPsi = Psi - DeltaImp(m)
     Psi2 = Psi*Psi
     Psi3 = Psi2*Psi
     Psi5 = Psi3*Psi2
     Psi7 = Psi5*Psi2
     Psi9 = Psi7*Psi2

     cf = 2.d0 * kImp(m) * dPsi / &
     &    ( Psi + k3*Psi3 + k5*Psi5 + k7*Psi7 + k9*Psi9 )

!     Psi4 = Psi2*Psi2
!     Psi6 = Psi4*Psi2
!     Psi8 = Psi4*Psi4
!     cf = 2.d0 * kPsi(m) / ( 1.d0 + k3*Psi2 + k5*Psi4 + k7*Psi6 + k9*Psi8 )

     Fix = cf * ( (Rkjy*Plbz-Plby*Rkjz) * rRab + CsPsi * (Play*Rkjz-Rkjy*Plaz) * rRa2 )
     Fiy = cf * ( (Rkjz*Plbx-Plbz*Rkjx) * rRab + CsPsi * (Plaz*Rkjx-Rkjz*Plax) * rRa2 )
     Fiz = cf * ( (Rkjx*Plby-Plbx*Rkjy) * rRab + CsPsi * (Plax*Rkjy-Rkjx*Play) * rRa2 )

     Fkx = cf * ( ( (Plby*Rijz-Rijy*Plbz) +        (Play*Rljz-Rljy*Plaz) ) * rRab    &
     &  + CsPsi * ( (Rijy*Plaz-Play*Rijz) * rRa2 + (Rljy*Plbz-Plby*Rljz) * rRb2 )  )
     Fky = cf * ( ( (Plbz*Rijx-Rijz*Plbx) +        (Plaz*Rljx-Rljz*Plax) ) * rRab    &
     &  + CsPsi * ( (Rijz*Plax-Plaz*Rijx) * rRa2 + (Rljz*Plbx-Plbz*Rljx) * rRb2 )  )
     Fkz = cf * ( ( (Plbx*Rijy-Rijx*Plby) +        (Plax*Rljy-Rljx*Play) ) * rRab    &
     &  + CsPsi * ( (Rijx*Play-Plax*Rijy) * rRa2 + (Rljx*Plby-Plbx*Rljy) * rRb2 )  )

     Flx = cf * ( (Rkjy*Plaz-Play*Rkjz) * rRab + CsPsi * (Plby*Rkjz-Rkjy*Plbz) * rRb2 )
     Fly = cf * ( (Rkjz*Plax-Plaz*Rkjx) * rRab + CsPsi * (Plbz*Rkjx-Rkjz*Plbx) * rRb2 )
     Flz = cf * ( (Rkjx*Play-Plax*Rkjy) * rRab + CsPsi * (Plbx*Rkjy-Rkjx*Plby) * rRb2 )

     Fjx = - ( Fix + Fkx + Flx )
     Fjy = - ( Fiy + Fky + Fly )
     Fjz = - ( Fiz + Fkz + Flz )

     Frc_Impro(1,i) = Frc_Impro(1,i) + Fix
     Frc_Impro(2,i) = Frc_Impro(2,i) + Fiy
     Frc_Impro(3,i) = Frc_Impro(3,i) + Fiz
     Frc_Impro(1,j) = Frc_Impro(1,j) + Fjx
     Frc_Impro(2,j) = Frc_Impro(2,j) + Fjy
     Frc_Impro(3,j) = Frc_Impro(3,j) + Fjz
     Frc_Impro(1,k) = Frc_Impro(1,k) + Fkx
     Frc_Impro(2,k) = Frc_Impro(2,k) + Fky
     Frc_Impro(3,k) = Frc_Impro(3,k) + Fkz
     Frc_Impro(1,l) = Frc_Impro(1,l) + Flx
     Frc_Impro(2,l) = Frc_Impro(2,l) + Fly
     Frc_Impro(3,l) = Frc_Impro(3,l) + Flz

     Ene_Impro = Ene_Impro + kImp(m) * dPsi * dPsi

   end do


   if(QRigidBody) then

     do i = 1  , N

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

end subroutine Force_Improper_CG
