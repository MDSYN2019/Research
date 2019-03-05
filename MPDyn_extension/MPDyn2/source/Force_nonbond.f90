! ############################
! ## SUBROUTINE LIST 
! ## -- Force_rspace 
! ## -- Force_rspace_ChRB_Fixed 
! ## -- Force_rspace_ChRB
! ## -- Force_rspace_ChRB_INTRASUBT
! ## -- Force_rspace_ChFX_Fixed
! ## -- Force_rspace_ChFX
! ## -- Force_rspace_OPRB_Fixed
! ## -- Force_rspace_OPRB
! ## -- Force_rspace_OPRB_INTRASUBT
! ## -- Force_rspace_OPFX_Fixed
! ## -- Force_rspace_OPFX
! ## -- Force_kspace 
! ## -- Force_rcut 
! ## -- Force_rcut_ChRB 
! ## -- Force_rcut_ChFX 
! ## -- Force_rcut_OPRB 
! ## -- Force_rcut_OPFX 
! ############################


!######################################################################
!######################################################################


! ******************************************************
! **  subroutine for calculating intermolecular force **
! **  real-space part of Ewald sum                    **
! **  interaction based cut-off method (LJ part )     **
! ******************************************************

subroutine Force_rspace

use Numbers, only : N
use CommonBlocks, only : QRigidBody, ForceField, QMacro
use NonbondParam, only : Frc_Ersp, Ene_LJ, Ene_Ersp, Vir_Ersp

implicit none

! -----------
!  Zero clear
! -----------
   Vir_Ersp = 0.d0
   Frc_Ersp = 0.d0
   Ene_LJ   = 0.d0
   Ene_Ersp = 0.d0

if(ForceField(1:3) == 'EAM') Return

if(ForceField(1:5) == 'CHARM') then

! ---------- Bond pair subtraction -------------

   call Force_Subt_12
   call Force_Subt_13
   call Force_Subt_14

! ----------------------------------------------
! ####  Calculation of LJ & Coulomb interaction
!

! ## for the system composed by the Rigid-Body molecules

   if(QRigidBody) then

     call Force_rspace_ChRB
!------------(intramolecule-subtraction)-------------------------
     call Force_rspace_ChRB_INTRASUBT

! ## flexible molecules

   else

     call Force_rspace_ChFX

     Vir_Ersp(2,1) = Vir_Ersp(1,2)
     Vir_Ersp(3,1) = Vir_Ersp(1,3)
     Vir_Ersp(3,2) = Vir_Ersp(2,3)

   end if

else if(ForceField(1:4) == 'OPLS') then

! ---------- Bond pair subtraction -------------
   call Force_Subt_12_OPLS
   call Force_Subt_13_OPLS
   call Force_Subt_14_OPLS
! ----------------------------------------------

! ####  Calculation of LJ & Coulomb interaction
! ## for the system composed by the Rigid-Body molecules
   if(QRigidBody) then

     call Force_rspace_OPRB
     call Force_rspace_OPRB_INTRASUBT ! intramolecule subtraction

! ## flexible molecules
   else

     call Force_rspace_OPFX

     Vir_Ersp(2,1) = Vir_Ersp(1,2)
     Vir_Ersp(3,1) = Vir_Ersp(1,3)
     Vir_Ersp(3,2) = Vir_Ersp(2,3)

   end if

else if(ForceField(1:3) == 'BKS') then

   call ForceR_BKS

end if

end subroutine Force_rspace

!######################################################################

subroutine Force_rspace_ChRB

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QSwitch
use BookParam, only : Npair, ListIJ
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, &
&  Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, l
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: ErrorFunc
real(8) :: fk1, fk2, fk, fkLJ, fkt
real(8) :: ek
real(8) :: cf
real(8) :: xtm, x
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk
real(8) :: Error_Function
external Error_Function

       Forcer = 0.d0
       Gk = 0.d0

       do l = 1 , Npair

         i = ListIJ(1,l)
         j = ListIJ(2,l)
         k = ListIJ(3,l)

         Rx = R(1,i) - R(1,j) + CellShft(1,k)
         Ry = R(2,i) - R(2,j) + CellShft(2,k)
         Rz = R(3,i) - R(3,j) + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if(R2 <= Rcutoff2) then

           InvR2 = 1.d0 / R2

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fk   = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           R1  = sqrt( R2 )
           x   = Alpha * R1

           ErrorFunc = Error_Function(x)

           xtm = -x * x
           fk1 = cf * ErrorFunc * R1 * InvR2
           fk2 = cf * ar2 * exp(xtm)
           fk  = ( fk1 + fk2 ) * InvR2

           Ene_Ersp = Ene_Ersp + fk1

           fkt = fkLJ + fk

           Fx = fkt * Rx
           Fy = fkt * Ry
           Fz = fkt * Rz

           Forcer(1,i) = Forcer(1,i) + Fx
           Forcer(2,i) = Forcer(2,i) + Fy
           Forcer(3,i) = Forcer(3,i) + Fz
           Forcer(1,j) = Forcer(1,j) - Fx
           Forcer(2,j) = Forcer(2,j) - Fy
           Forcer(3,j) = Forcer(3,j) - Fz
           Gk(1,k)     = Gk(1,k)     + Fx
           Gk(2,k)     = Gk(2,k)     + Fy
           Gk(3,k)     = Gk(3,k)     + Fz

         end if

       end do

       call VirialBekkerRB(Forcer,Vir_Ersp,Gk)

       Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rspace_ChRB

!######################################################################

subroutine Force_rspace_ChRB_INTRASUBT

use Numbers, only : N
use Configuration, only : R
use CommonMPI, only : NProcs, MyRank
use RBparam, only : NumRB, QSingle, InitAtomRB, RBType, NumRBAtom
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

integer :: i, j
real(8) :: R1, R2, InvR2
real(8) :: ErrorFunc
real(8) :: fk1, fk2, fk
real(8) :: cf
real(8) :: xtm, x
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
integer :: MyType, Nc
integer :: i1, i2, k1, k2
integer :: Nas
real(8) :: Error_Function
external Error_Function

     Nas = NProcs - MyRank

     do i = Nas , NumRB, NProcs

       if(QSingle(i)) cycle

       j = InitAtomRB(i)
       MyType = RBType(i)
       Nc = NumRBAtom(MyType)

       do k1 = 1 , Nc-1

         do k2 = k1 + 1 , Nc

           i1 = j + k1
           i2 = j + k2

           cf = Charge(i1) * Charge(i2)

           if( cf /= 0. ) then

             Rx= R(1,i1) - R(1,i2)
             Ry= R(2,i1) - R(2,i2)
             Rz= R(3,i1) - R(3,i2)
             R2 = Rx*Rx + Ry*Ry + Rz*Rz

             InvR2 = 1.d0 / R2

             R1  = sqrt( R2 )
             x   = alpha * R1

             ErrorFunc = Error_Function(x)

             xtm = -x * x
             fk1 = cf * (ErrorFunc-1.d0) * R1 * InvR2
             fk2 = cf * ar2 * exp(xtm)
             fk  = ( fk1 + fk2 ) * InvR2

             Fx = fk * Rx
             Fy = fk * Ry
             Fz = fk * Rz

             Frc_Ersp(1,i1) = Frc_Ersp(1,i1) + Fx
             Frc_Ersp(2,i1) = Frc_Ersp(2,i1) + Fy
             Frc_Ersp(3,i1) = Frc_Ersp(3,i1) + Fz
             Frc_Ersp(1,i2) = Frc_Ersp(1,i2) - Fx
             Frc_Ersp(2,i2) = Frc_Ersp(2,i2) - Fy
             Frc_Ersp(3,i2) = Frc_Ersp(3,i2) - Fz

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

     end do

end subroutine Force_rspace_ChRB_INTRASUBT

!######################################################################

subroutine Force_rspace_ChFX

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QPINPT, QSwitch
use BookParam, only : Npair, ListIJ
use EwaldParam, only : Alpha, ar2, msh, EFlist
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, &
&  Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, l, ii
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: ErrorFunc
real(8) :: fk1, fk2, fk, fkLJ, fkt
real(8) :: ek
real(8) :: cf
real(8) :: xtm, x
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk
real(8) :: yy, dx, y1, y2, y3, z1, z2
real(8) :: Xoff, Xon
real(8) :: Switch, Dswitch

       Gk = 0.d0
       Forcer = 0.d0

       do l = 1 , Npair

         i = ListIJ(1,l)
         j = ListIJ(2,l)
         k = ListIJ(3,l)

         Rx = R(1,i) - R(1,j) + CellShft(1,k)
         Ry = R(2,i) - R(2,j) + CellShft(2,k)
         Rz = R(3,i) - R(3,j) + CellShft(3,k)
         R2 = Rx * Rx + Ry * Ry + Rz * Rz

         if(R2 <= Rcutoff2) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)
           InvR2 = 1.d0 / R2

           SR2  = Sgm2 * InvR2
           SR6  = SR2 * SR2 * SR2
           SR12 = SR6 * SR6
           fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2
           ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
           if(QSwitch.and.R2>Ron2) then
             Xoff = R2 - Rcutoff2      ! x-xoff
             Xon  = R2 - Ron2          ! x-xon
             Switch  = ( 3.d0 * Xon - Xoff ) * Xoff * Xoff * swf1
             Dswitch = 12.d0 * swf1 * Xoff * Xon
             fk = fk * Switch + Dswitch * ek
             ek = ek * Switch
           end if
! --------------------------------------------------------

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           R1  = sqrt( R2 )
           x   = Alpha * R1

           yy = x * msh
           ii = nint( yy )

           dx  = yy - dble(ii)

           y1  = EFList( ii-1 )
           y2  = EFList( ii   )
           y3  = EFList( ii+1 )

           z2  = y1 - 2.d0 * y2 + y3
           z1  = - y1 + y3
           ErrorFunc  = 0.5d0 * (z1 + z2 * dx ) * dx + y2

           xtm = -x * x
           fk1 = cf * ErrorFunc * R1 * InvR2
           fk2 = cf * ar2 * exp(xtm)
           fk  = ( fk1 + fk2 ) * InvR2

           Ene_Ersp = Ene_Ersp + fk1

           fkt = fkLJ + fk

           Fx = fkt * Rx
           Fy = fkt * Ry
           Fz = fkt * Rz

           Forcer(1,i) = Forcer(1,i) + Fx
           Forcer(2,i) = Forcer(2,i) + Fy
           Forcer(3,i) = Forcer(3,i) + Fz
           Forcer(1,j) = Forcer(1,j) - Fx
           Forcer(2,j) = Forcer(2,j) - Fy
           Forcer(3,j) = Forcer(3,j) - Fz
           Gk(1,k)     = Gk(1,k)     + Fx
           Gk(2,k)     = Gk(2,k)     + Fy
           Gk(3,k)     = Gk(3,k)     + Fz

         end if

       end do

       if(QPINPT) then
         call VirialBekkerPI(Forcer,Vir_Ersp,Gk)
       else
         call VirialBekker(Forcer,Vir_Ersp,Gk)
       end if

       Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rspace_ChFX

!######################################################################

subroutine Force_rspace_OPRB

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QSwitch
use BookParam, only : Npair, ListIJ
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, &
&  Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, l
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: ErrorFunc
real(8) :: fk1, fk2, fk, fkLJ, fkt
real(8) :: ek
real(8) :: cf
real(8) :: xtm, x
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk
real(8) :: Error_Function
external Error_Function

       Forcer = 0.d0
       Gk = 0.d0

       do l = 1 , Npair

         i = ListIJ(1,l)
         j = ListIJ(2,l)
         k = ListIJ(3,l)

         Rx = R(1,i) - R(1,j) + CellShft(1,k)
         Ry = R(2,i) - R(2,j) + CellShft(2,k)
         Rz = R(3,i) - R(3,j) + CellShft(3,k)
         R2 = Rx*Rx + Ry*Ry + Rz*Rz

         if(R2 <= Rcutoff2) then

           InvR2 = 1.d0 / R2
           Eps   = EpsLJ(i) * EpsLJ(j)

           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           R1  = sqrt( R2 )
           x   = Alpha * R1

           ErrorFunc = Error_Function(x)

           xtm = -x * x
           fk1 = cf * ErrorFunc * R1 * InvR2
           fk2 = cf * ar2 * exp(xtm)
           fk  = ( fk1 + fk2 ) * InvR2

           Ene_Ersp = Ene_Ersp + fk1

           fkt = fkLJ + fk

           Fx = fkt * Rx
           Fy = fkt * Ry
           Fz = fkt * Rz

           Forcer(1,i) = Forcer(1,i) + Fx
           Forcer(2,i) = Forcer(2,i) + Fy
           Forcer(3,i) = Forcer(3,i) + Fz
           Forcer(1,j) = Forcer(1,j) - Fx
           Forcer(2,j) = Forcer(2,j) - Fy
           Forcer(3,j) = Forcer(3,j) - Fz
           Gk(1,k)     = Gk(1,k)     + Fx
           Gk(2,k)     = Gk(2,k)     + Fy
           Gk(3,k)     = Gk(3,k)     + Fz

         end if

       end do

       call VirialBekkerRB(Forcer,Vir_Ersp,Gk)

       Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rspace_OPRB

!######################################################################

subroutine Force_rspace_OPRB_INTRASUBT

use Configuration, only : R
use CommonMPI, only : NProcs, MyRank
use RBparam, only : NumRB, QSingle, InitAtomRB, RBType, NumRBAtom
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, Ene_LJ

implicit none

integer :: i, j
real(8) :: R1, R2, InvR2
real(8) :: ErrorFunc
real(8) :: fk1, fk2, fk
real(8) :: cf
real(8) :: xtm, x
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
integer :: MyType, Nc
integer :: i1, i2, k1, k2
integer :: Nas
real(8) :: Error_Function
external Error_Function

     Nas = NProcs - MyRank

     do i = Nas , NumRB, NProcs

       if(QSingle(i)) cycle

       j = InitAtomRB(i)
       MyType = RBType(i)
       Nc = NumRBAtom(MyType)

       do k1 = 1 , Nc-1

         do k2 = k1 + 1 , Nc

           i1 = j + k1
           i2 = j + k2

           cf = Charge(i1) * Charge(i2)

           if( cf /= 0. ) then

             Rx= R(1,i1) - R(1,i2)
             Ry= R(2,i1) - R(2,i2)
             Rz= R(3,i1) - R(3,i2)
             R2 = Rx*Rx + Ry*Ry + Rz*Rz

             InvR2 = 1.d0 / R2

             R1  = sqrt( R2 )
             x   = alpha * R1

             ErrorFunc = Error_Function(x)

             xtm = -x * x
             fk1 = cf * (ErrorFunc-1.d0) * R1 * InvR2
             fk2 = cf * ar2 * exp(xtm)
             fk  = ( fk1 + fk2 ) * InvR2

             Fx = fk * Rx
             Fy = fk * Ry
             Fz = fk * Rz

             Frc_Ersp(1,i1) = Frc_Ersp(1,i1) + Fx
             Frc_Ersp(2,i1) = Frc_Ersp(2,i1) + Fy
             Frc_Ersp(3,i1) = Frc_Ersp(3,i1) + Fz
             Frc_Ersp(1,i2) = Frc_Ersp(1,i2) - Fx
             Frc_Ersp(2,i2) = Frc_Ersp(2,i2) - Fy
             Frc_Ersp(3,i2) = Frc_Ersp(3,i2) - Fz

             Ene_Ersp = Ene_Ersp + fk1

           end if

         end do

       end do

     end do



end subroutine Force_rspace_OPRB_INTRASUBT

!######################################################################

subroutine Force_rspace_OPFX

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QPINPT, QSwitch
use BookParam, only : Npair, ListIJ
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, &
&  Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, l
real(8) :: R1, R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: ErrorFunc
real(8) :: fk1, fk2, fk, fkLJ, fkt
real(8) :: ek
real(8) :: cf
real(8) :: xtm, x
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk
real(8) :: Error_Function
external Error_Function

       Forcer = 0.d0
       Gk = 0.d0

       do l = 1 , Npair

         i = ListIJ(1,l)
         j = ListIJ(2,l)
         k = ListIJ(3,l)

         Rx = R(1,i) - R(1,j) + CellShft(1,k)
         Ry = R(2,i) - R(2,j) + CellShft(2,k)
         Rz = R(3,i) - R(3,j) + CellShft(3,k)
         R2 = Rx * Rx + Ry * Ry + Rz * Rz

         if(R2 <= Rcutoff2) then

           InvR2 = 1.d0 / R2
           Eps   = EpsLJ(i) * EpsLJ(j)
           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           ek   = Eps * 4.d0 * ( SR12 - SR6 )

! --------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)

           R1  = sqrt( R2 )
           x   = Alpha * R1

           ErrorFunc = Error_Function(x)

           xtm = -x * x
           fk1 = cf * ErrorFunc * R1 * InvR2
           fk2 = cf * ar2 * exp(xtm)
           fk  = ( fk1 + fk2 ) * InvR2

           Ene_Ersp = Ene_Ersp + fk1

           fkt = fkLJ + fk

           Fx = fkt * Rx
           Fy = fkt * Ry
           Fz = fkt * Rz

           Forcer(1,i) = Forcer(1,i) + Fx
           Forcer(2,i) = Forcer(2,i) + Fy
           Forcer(3,i) = Forcer(3,i) + Fz
           Forcer(1,j) = Forcer(1,j) - Fx
           Forcer(2,j) = Forcer(2,j) - Fy
           Forcer(3,j) = Forcer(3,j) - Fz
           Gk(1,k)     = Gk(1,k)     + Fx
           Gk(2,k)     = Gk(2,k)     + Fy
           Gk(3,k)     = Gk(3,k)     + Fz

         end if

       end do

       if(QPINPT) then
         call VirialBekkerPI(Forcer,Vir_Ersp,Gk)
       else
         call VirialBekker(Forcer,Vir_Ersp,Gk)
       end if

       Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rspace_OPFX


!#####################################################################
!#####################################################################


! ******************************************************
! **  subroutine for calculating intermolecular force **
! **  interaction in reciprocal space (Ewald sum)     **
! ******************************************************

subroutine Force_kspace

use CommonBlocks, only : ForceField, cCOULOMB
use EwaldParam, only : Frc_Eksp, Vir_Eksp, Ene_Eksp

implicit NONE

!  Zero clear

   Ene_Eksp = 0.d0
   Vir_Eksp = 0.d0

   Frc_Eksp = 0.d0

   if(ForceField(1:3) == 'EAM') Return

   if(cCOULOMB == 'PME') then

     call PME_Reciprocal

   else if(cCOULOMB == 'EWALD') then

     call Ewald_Reciprocal

   end if

end subroutine Force_kspace


!######################################################################
!######################################################################


! ******************************************************
! **  subroutine for calculating intermolecular force **
! **  interaction based cut-off method                **
! ******************************************************

subroutine Force_rcut

use Numbers, only : N
use CommonBlocks, only : QRigidBody, ForceField, QMacro
use Configuration, only : R
use NonbondParam, only : Frc_Ersp, Frc_Elec, Ene_Ersp, Ene_LJ, &
& Ene_Elec, Vir_Ersp

implicit none

! -----------
!  Zero clear
! -----------
   Vir_Ersp = 0.d0
   Frc_Ersp = 0.d0
   Ene_LJ   = 0.d0
   Ene_Ersp = 0.d0

   if( ForceField(1:5) == 'CHARM' ) then

! ---------- Bond pair subtraction -------------

     call Force_Subt_12
     call Force_Subt_13
     call Force_Subt_14

! ----------------------------------------------
! ####  Calculation of LJ & Coulomb interaction
!
     if(QRigidBody) then
       call Force_rcut_ChRB
     else
       call Force_rcut_ChFX
     end if

   else if(ForceField(1:4) == 'OPLS') then

! ---------- Bond pair subtraction -------------

     call Force_Subt_12_OPLS
     call Force_Subt_13_OPLS
     call Force_Subt_14_OPLS

! ----------------------------------------------
! ####  Calculation of LJ & Coulomb interaction
!
     if(QRigidBody) then
       call Force_rcut_OPRB
     else
       print *,"EE#"
       call Force_rcut_OPFX
     end if

   end if

   Frc_Elec = Frc_Ersp
   Ene_Elec = Ene_Ersp


end subroutine Force_rcut


!######################################################################
!######################################################################


subroutine Force_rcut_ChRB

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QSwitch
use BookParam, only : Npair, ListIJ
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, &
&  Vir_Ersp, Ene_LJ
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, l
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, fkLJ, fkt
real(8) :: ek, cf
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk

   Forcer = 0.d0
   Gk = 0.d0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)
     k = ListIJ(3,l)

     Rx = R(1,i) - R(1,j) + CellShft(1,k)
     Ry = R(2,i) - R(2,j) + CellShft(2,k)
     Rz = R(3,i) - R(3,j) + CellShft(3,k)
     R2 = Rx * Rx + Ry * Ry + Rz * Rz

     if(R2 <= Rcutoff2) then

       InvR2 = 1.d0 / R2

       Sgm   = Rminh(i) + Rminh(j)
       Sgm2  = Sgm * Sgm
       Eps   = EpsLJ(i) * EpsLJ(j)

       SR2  = Sgm2 * InvR2                    !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2                 !         ^6
       SR12 = SR6 * SR6                       !         ^12
       fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
       ek   = Eps * ( SR12 - 2.d0 * SR6 )

! --------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

       Ene_LJ = Ene_LJ + ek

       cf = Charge(i) * Charge(j)

       fk1 = cf * sqrt( InvR2 )
       fk  = fk1 * InvR2

       fkt = fkLJ + fk

       Fx = fkt * Rx
       Fy = fkt * Ry
       Fz = fkt * Rz

       Forcer(1,i) = Forcer(1,i) + Fx
       Forcer(2,i) = Forcer(2,i) + Fy
       Forcer(3,i) = Forcer(3,i) + Fz
       Forcer(1,j) = Forcer(1,j) - Fx
       Forcer(2,j) = Forcer(2,j) - Fy
       Forcer(3,j) = Forcer(3,j) - Fz
       Gk(1,k)     = Gk(1,k)     + Fx
       Gk(2,k)     = Gk(2,k)     + Fy
       Gk(3,k)     = Gk(3,k)     + Fz

       Ene_Ersp = Ene_Ersp + fk1

     end if

   end do

   call VirialBekkerRB(Forcer,Vir_Ersp,Gk)

   Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rcut_ChRB

! ##################################################

subroutine Force_rcut_ChFX

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QPINPT, QSwitch
use BookParam, only : Npair, ListIJ
use NonbondParam, only : Charge, Rminh, EpsLJ, Frc_Ersp, Ene_Ersp, &
& Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm, Sgm2, Eps

integer :: i, j, k, l
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, fkLJ, fkt
real(8) :: ek, cf
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk

   Forcer = 0.d0
   Gk = 0.d0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)
     k = ListIJ(3,l)

     Rx = R(1,i) - R(1,j) + CellShft(1,k)
     Ry = R(2,i) - R(2,j) + CellShft(2,k)
     Rz = R(3,i) - R(3,j) + CellShft(3,k)
     R2 = Rx * Rx + Ry * Ry + Rz * Rz

     if(R2 <= Rcutoff2) then

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
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

       Ene_LJ = Ene_LJ + ek

       cf = Charge(i) * Charge(j)

       fk1 = cf * sqrt(InvR2)
       fk  = fk1 * InvR2

       Ene_Ersp = Ene_Ersp + fk1

       fkt = fkLJ + fk

       Fx = fkt * Rx
       Fy = fkt * Ry
       Fz = fkt * Rz

       Forcer(1,i) = Forcer(1,i) + Fx
       Forcer(2,i) = Forcer(2,i) + Fy
       Forcer(3,i) = Forcer(3,i) + Fz
       Forcer(1,j) = Forcer(1,j) - Fx
       Forcer(2,j) = Forcer(2,j) - Fy
       Forcer(3,j) = Forcer(3,j) - Fz
       Gk(1,k)     = Gk(1,k)     + Fx
       Gk(2,k)     = Gk(2,k)     + Fy
       Gk(3,k)     = Gk(3,k)     + Fz

     end if

   end do

   if(QPINPT) then
     call VirialBekkerPI(Forcer,Vir_Ersp,Gk)
   else
     call VirialBekker(Forcer,Vir_Ersp,Gk)
   end if

   Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rcut_ChFX

! ##################################################

subroutine Force_rcut_OPRB

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QSwitch
use BookParam, only : Npair, ListIJ
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, &
& Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, l
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, fkLJ, fkt
real(8) :: ek, cf
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk

   Forcer = 0.d0
   Gk = 0.d0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)
     k = ListIJ(3,l)

     Rx = R(1,i) - R(1,j) + CellShft(1,k)
     Ry = R(2,i) - R(2,j) + CellShft(2,k)
     Rz = R(3,i) - R(3,j) + CellShft(3,k)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     if(R2 <= Rcutoff2) then

       InvR2 = 1.d0 / R2

       Eps   = EpsLJ(i) * EpsLJ(j)
       Sgm2  = SgmLJ(i) * SgmLJ(j)

       SR2  = Sgm2 * InvR2                    !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2                 !         ^6
       SR12 = SR6 * SR6                       !         ^12
       fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
       ek   = Eps * 4.d0  * ( SR12 - SR6 )

! --------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

       Ene_LJ = Ene_LJ + ek

       cf = Charge(i) * Charge(j)

       fk1 = cf * sqrt(InvR2)
       fk  = fk1 * InvR2

       Ene_Ersp = Ene_Ersp + fk1

       fkt = fkLJ + fk

       Fx = fkt * Rx
       Fy = fkt * Ry
       Fz = fkt * Rz

       Forcer(1,i) = Forcer(1,i) + Fx
       Forcer(2,i) = Forcer(2,i) + Fy
       Forcer(3,i) = Forcer(3,i) + Fz
       Forcer(1,j) = Forcer(1,j) - Fx
       Forcer(2,j) = Forcer(2,j) - Fy
       Forcer(3,j) = Forcer(3,j) - Fz
       Gk(1,k)     = Gk(1,k)     + Fx
       Gk(2,k)     = Gk(2,k)     + Fy
       Gk(3,k)     = Gk(3,k)     + Fz

     end if

   end do

   call VirialBekkerRB(Forcer,Vir_Ersp,Gk)

   Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rcut_OPRB

! ##################################################

subroutine Force_rcut_OPFX

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QPINPT, QSwitch
use BookParam, only : Npair, ListIJ
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Frc_Ersp, Ene_Ersp, &
& Ene_LJ, Vir_Ersp
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CellParam, only : CellShft

implicit none

real(8) :: Sgm2, Eps

integer :: i, j, k, l
real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: fk1, fk, fkLJ, fkt
real(8) :: ek, cf
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8), dimension(3,N) :: Forcer
real(8), dimension(3,-13:13) :: Gk

   Gk = 0.d0
   Forcer = 0.d0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)
     k = ListIJ(3,l)

     Rx = R(1,i) - R(1,j) + CellShft(1,k)
     Ry = R(2,i) - R(2,j) + CellShft(2,k)
     Rz = R(3,i) - R(3,j) + CellShft(3,k)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     if(R2 <= Rcutoff2) then

       InvR2 = 1.d0 / R2

       Eps   = EpsLJ(i) * EpsLJ(j)
       Sgm2  = SgmLJ(i) * SgmLJ(j)

       SR2  = Sgm2 * InvR2                    !(sigma/r)^2
       SR6  = SR2 * SR2 * SR2                 !         ^6
       SR12 = SR6 * SR6                       !         ^12
       fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
       ek   = Eps * 4.d0  * ( SR12 - SR6 )

! --------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------

       Ene_LJ = Ene_LJ + ek

       cf = Charge(i) * Charge(j)

       fk1 = cf * sqrt(InvR2)
       fk  = fk1 * InvR2

       Ene_Ersp = Ene_Ersp + fk1

       fkt = fkLJ + fk

       Fx = fkt * Rx
       Fy = fkt * Ry
       Fz = fkt * Rz

       Forcer(1,i) = Forcer(1,i) + Fx
       Forcer(2,i) = Forcer(2,i) + Fy
       Forcer(3,i) = Forcer(3,i) + Fz
       Forcer(1,j) = Forcer(1,j) - Fx
       Forcer(2,j) = Forcer(2,j) - Fy
       Forcer(3,j) = Forcer(3,j) - Fz
       Gk(1,k)     = Gk(1,k)     + Fx
       Gk(2,k)     = Gk(2,k)     + Fy
       Gk(3,k)     = Gk(3,k)     + Fz

     end if

   end do

   if(QPINPT) then
     call VirialBekkerPI(Forcer,Vir_Ersp,Gk)
   else
     call VirialBekker(Forcer,Vir_Ersp,Gk)
   end if

   Frc_Ersp = Frc_Ersp + Forcer

end subroutine Force_rcut_OPFX
