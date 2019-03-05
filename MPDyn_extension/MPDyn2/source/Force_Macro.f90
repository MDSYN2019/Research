! ############################
! ## SUBROUTINE LIST 
! ## -- Force_MacroSphAA
! ## -- Force_MacroSphCG
! ## -- Force_MacroMacro
! ## -- InterPolateMacroSph
! ## -- InterPolateMacroMacro
! ## -- Set_Sphere
! ############################


!######################################################################
!######################################################################


subroutine Force_MacroSphAA

use CGball
use Numbers, only : N
use Configuration, only : R
use OptConstraintParam, only : Frc_OptC, Ene_OptC, Vir_OptC
use CellParam, only : CellShft

implicit none

integer :: i, j, k, ii, jj
real(8) :: Rc2, Rm, xx
real(8) :: Rx, Ry, Rz, FF, R2, Uij
real(8) :: InterPolateMacroSph, Fc, Fr
real(8) :: Fx, Fy, Fz
real(8) :: Rix, Riy, Riz, Fix, Fiy, Fiz
real(8), dimension(:,:), allocatable :: Forcer
real(8), dimension(:,:), allocatable :: Gk
external InterPolateMacroSph

   allocate(Forcer(3,N))
   allocate(Gk(3,-13:13))

   Rc2 = Rcutmacro2(1,1)
   Rm  = Rminmacro2(1,1)
   xx = InvddSph(1,1)

   Forcer = 0.d0
   Gk = 0.d0

   do ii = 1, NumSphere

     i = IdSph(ii)
     FF = 0.d0
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)
     Fix = 0.d0
     Fiy = 0.d0
     Fiz = 0.d0

     do jj = 1, NmacList(ii)

       j = Listmac(1,jj,ii)
       k = Listmac(2,jj,ii)

       Rx = Rix - R(1,j) + CellShft(1,k)
       Ry = Riy - R(2,j) + CellShft(2,k)
       Rz = Riz - R(3,j) + CellShft(3,k)
       R2 = Rx * Rx + Ry * Ry + Rz * Rz

       if(R2 <= Rc2) then

         Uij =   InterPolateMacroSph(R2,Rm,xx,1,1,1)
         Fc  = - InterPolateMacroSph(R2,Rm,xx,2,1,1)
         Fr  =   InterPolateMacroSph(R2,Rm,xx,3,1,1)

         Ene_OptC = Ene_OptC + Uij

         Fx = Fc * Rx
         Fy = Fc * Ry
         Fz = Fc * Rz

         Fix = Fix + Fx
         Fiy = Fiy + Fy
         Fiz = Fiz + Fz
         Forcer(1,j) = Forcer(1,j) - Fx
         Forcer(2,j) = Forcer(2,j) - Fy
         Forcer(3,j) = Forcer(3,j) - Fz
         Gk(1,k) = Gk(1,k) + Fx
         Gk(2,k) = Gk(2,k) + Fy
         Gk(3,k) = Gk(3,k) + Fz

         FF = FF + Fr

       end if

     end do

     FSphRs(ii) = FSphRs(ii) + FF
     Forcer(1,i) = Fix
     Forcer(2,i) = Fiy
     Forcer(3,i) = Fiz

   end do

   call VirialBekker(Forcer,Vir_OptC,Gk)

   Frc_OptC = Frc_OptC + Forcer

   deallocate(Forcer,Gk)

end subroutine Force_MacroSphAA


!######################################################################
!######################################################################


subroutine Force_MacroSphCG

use Numbers, only : N
use Configuration, only : R
use OptConstraintParam, only : Frc_OptC, Ene_OptC, Vir_OptC
use CellParam, only : CellShft
use CGdata, only : NBAtomType
use CGball

implicit none

integer :: i, j, k, ii, jj, itype
integer :: jtype
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: Rix, Riy, Riz
real(8) :: Fix, Fiy, Fiz
real(8) :: R2, Uij, FF, Fc, Fr
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
real(8), dimension(:,:), allocatable :: Forcer
real(8), dimension(:,:), allocatable :: Gk
real(8) :: Rm, InterPolateMacroSph, xx
external InterPolateMacroSph

   allocate(Forcer(3,N))
   allocate(Gk(3,-13:13))

   Forcer = 0.d0
   Gk = 0.d0

if(NumSphere > 1) then

   do ii = 1, NumSphere

     i = IdSph(ii)
     itype = ItypeSph(ii)
     FF = 0.d0
     Fix = 0.d0
     Fiy = 0.d0
     Fiz = 0.d0
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     do jj = 1 , NmacList(ii)

       j = Listmac(1,jj,ii)
       k = Listmac(2,jj,ii)

       Rx = Rix - R(1,j) + CellShft(1,k)
       Ry = Riy - R(2,j) + CellShft(2,k)
       Rz = Riz - R(3,j) + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       jtype = NBAtomType(j)

       if(R2 <= Rcutmacro2(jtype,itype)) then

         Rm  = Rminmacro2(jtype,itype)
         xx  = InvddSph(jtype,itype)
         Uij =   InterPolateMacroSph(R2,Rm,xx,1,jtype,itype)
         Fc  = - InterPolateMacroSph(R2,Rm,xx,2,jtype,itype)

         Ene_OptC = Ene_OptC + Uij

         Fx  = Fc * Rx
         Fy  = Fc * Ry
         Fz  = Fc * Rz

         Fix = Fix + Fx
         Fiy = Fiy + Fy
         Fiz = Fiz + Fz
         Forcer(1,j) = Forcer(1,j) - Fx
         Forcer(2,j) = Forcer(2,j) - Fy
         Forcer(3,j) = Forcer(3,j) - Fz
         Gk(1,k) = Gk(1,k) + Fx
         Gk(2,k) = Gk(2,k) + Fy
         Gk(3,k) = Gk(3,k) + Fz

       end if

     end do

     Forcer(1,i) = Forcer(1,i) + Fix
     Forcer(2,i) = Forcer(2,i) + Fiy
     Forcer(3,i) = Forcer(3,i) + Fiz

   end do

else

   i = IdSph(1)
   itype = ItypeSph(1)
   FF = 0.d0
   Fix = 0.d0
   Fiy = 0.d0
   Fiz = 0.d0
   Rix = R(1,i)
   Riy = R(2,i)
   Riz = R(3,i)

   do jj = 1 , NmacList(1)

     j = Listmac(1,jj,1)
     k = Listmac(2,jj,1)

     Rx = Rix - R(1,j) + CellShft(1,k)
     Ry = Riy - R(2,j) + CellShft(2,k)
     Rz = Riz - R(3,j) + CellShft(3,k)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     jtype = NBAtomType(j)

     if(R2 <= Rcutmacro2(jtype,itype)) then

       Rm  = Rminmacro2(jtype,itype)
       xx  = InvddSph(jtype,itype)
       Uij =   InterPolateMacroSph(R2,Rm,xx,1,jtype,itype)
       Fc  = - InterPolateMacroSph(R2,Rm,xx,2,jtype,itype)
       Fr  =   InterPolateMacroSph(R2,Rm,xx,3,jtype,itype)

       Ene_OptC = Ene_OptC + Uij

       Fx  = Fc * Rx
       Fy  = Fc * Ry
       Fz  = Fc * Rz

       Fix = Fix + Fx
       Fiy = Fiy + Fy
       Fiz = Fiz + Fz
       Forcer(1,j) = Forcer(1,j) - Fx
       Forcer(2,j) = Forcer(2,j) - Fy
       Forcer(3,j) = Forcer(3,j) - Fz
       Gk(1,k) = Gk(1,k) + Fx
       Gk(2,k) = Gk(2,k) + Fy
       Gk(3,k) = Gk(3,k) + Fz

       FF = FF + Fr

     end if

   end do

   Forcer(1,i) = Forcer(1,i) + Fix
   Forcer(2,i) = Forcer(2,i) + Fiy
   Forcer(3,i) = Forcer(3,i) + Fiz
   FSphRs(1) = FSphRs(1) + FF

end if

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do i = 1, N
     Fx = Forcer(1,i)
     Fy = Forcer(2,i)
     Fz = Forcer(3,i)
     Rx = R(1,i)
     Ry = R(2,i)
     Rz = R(3,i)
     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz
   end do

   do k = -13, 13
     if(k==0) cycle
     Fx = Gk(1,k)
     Fy = Gk(2,k)
     Fz = Gk(3,k)
     Rx = CellShft(1,k)
     Ry = CellShft(2,k)
     Rz = CellShft(3,k)
     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz
   end do

   Frc_OptC = Frc_OptC + Forcer

   Vir_OptC(1,1) = Vir_OptC(1,1) + Vxx
   Vir_OptC(1,2) = Vir_OptC(1,2) + Vxy
   Vir_OptC(1,3) = Vir_OptC(1,3) + Vxz
   Vir_OptC(2,2) = Vir_OptC(2,2) + Vyy
   Vir_OptC(2,3) = Vir_OptC(2,3) + Vyz
   Vir_OptC(3,3) = Vir_OptC(3,3) + Vzz

   Vir_OptC(2,1) = Vir_OptC(1,2)
   Vir_OptC(3,1) = Vir_OptC(1,3)
   Vir_OptC(3,2) = Vir_OptC(2,3)

   deallocate(Forcer,Gk)

   if(NumSphere>1) call Force_MacroMacro

end subroutine Force_MacroSphCG


!######################################################################
!######################################################################


subroutine Force_MacroMacro

use Numbers, only : N
use Configuration, only : R
use OptConstraintParam, only : Frc_OptC, Ene_OptC, Vir_OptC
use CellParam, only : CellShft
use CGball

implicit none

integer :: i, j, k, ii, itype, jtype
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: R2, Uij, Fc
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
real(8), dimension(:,:), allocatable :: Rmc
real(8), dimension(:,:), allocatable :: Forcer
real(8), dimension(:,:), allocatable :: Gk
real(8) :: InterPolateMacroMacro
external InterPolateMacroMacro

   allocate(Rmc(3,NumSphere))
   allocate(Forcer(3,NumSphere))
   allocate(Gk(3,-13:13))

   do ii = 1, NumSphere
     i = IdSph(ii)
     Rmc(:,ii) = R(:,i)
   end do

   Forcer = 0.d0
   Gk = 0.d0

   do ii = 1, NumPairMM

     i = ListMM(1,ii)
     j = ListMM(2,ii)
     k = ListMM(3,ii)

     itype = ItypeSph(i)
     jtype = ItypeSph(j)

     Rx = Rmc(1,i) - Rmc(1,j) + CellShft(1,k)
     Ry = Rmc(2,i) - Rmc(2,j) + CellShft(2,k)
     Rz = Rmc(3,i) - Rmc(3,j) + CellShft(3,k)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     if(R2 <= Rcutmm2(itype,jtype)) then

       Uij =   InterPolateMacroMacro(R2,1,itype,jtype)
       Fc  = - InterPolateMacroMacro(R2,2,itype,jtype)

       Ene_OptC = Ene_OptC + Uij

       Fx  = Fc * Rx
       Fy  = Fc * Ry
       Fz  = Fc * Rz

       Forcer(1,i) = Forcer(1,i) + Fx
       Forcer(2,i) = Forcer(2,i) + Fy
       Forcer(3,i) = Forcer(3,i) + Fz
       Forcer(1,j) = Forcer(1,j) - Fx
       Forcer(2,j) = Forcer(2,j) - Fy
       Forcer(3,j) = Forcer(3,j) - Fz
       Gk(1,k) = Gk(1,k) + Fx
       Gk(2,k) = Gk(2,k) + Fy
       Gk(3,k) = Gk(3,k) + Fz

     end if

   end do

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do i = 1, NumSphere
     Fx = Forcer(1,i)
     Fy = Forcer(2,i)
     Fz = Forcer(3,i)
     Rx = Rmc(1,i)
     Ry = Rmc(2,i)
     Rz = Rmc(3,i)
     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz
   end do

   do k = -13, 13
     if(k==0) cycle
     Fx = Gk(1,k)
     Fy = Gk(2,k)
     Fz = Gk(3,k)
     Rx = CellShft(1,k)
     Ry = CellShft(2,k)
     Rz = CellShft(3,k)
     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz
   end do

   do ii = 1, NumSphere
     i = IdSph(ii)
     Frc_OptC(:,i) = Frc_OptC(:,i) + Forcer(:,ii)
   end do

   Vir_OptC(1,1) = Vir_OptC(1,1) + Vxx
   Vir_OptC(1,2) = Vir_OptC(1,2) + Vxy
   Vir_OptC(1,3) = Vir_OptC(1,3) + Vxz
   Vir_OptC(2,2) = Vir_OptC(2,2) + Vyy
   Vir_OptC(2,3) = Vir_OptC(2,3) + Vyz
   Vir_OptC(3,3) = Vir_OptC(3,3) + Vzz

   Vir_OptC(2,1) = Vir_OptC(1,2)
   Vir_OptC(3,1) = Vir_OptC(1,3)
   Vir_OptC(3,2) = Vir_OptC(2,3)

   deallocate(Rmc,Forcer,Gk)

end subroutine Force_MacroMacro


!######################################################################
!######################################################################


Function InterPolateMacroSph(x,drmin,invdr,j,k,i) Result(rs)

use CGball

implicit none

integer :: ii, j, k, i
real(8) :: x, xx, dx, y1, y2, y3
real(8) :: z1, z2
real(8) :: rs, drmin, invdr

#ifdef CUBIC
real(8) :: y4, z3
real(8), parameter :: sb = 0.3333333333333333d0
real(8), parameter :: sb2 = 0.6666666666666667d0

   xx  = (x - drmin) * invdr

   ii  = int( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = SphFunc( ii  , j, k, i )
   y2  = SphFunc( ii+1, j, k, i )
   y3  = SphFunc( ii+2, j, k, i )
   y4  = SphFunc( ii+3, j, k, i )

   z3  = y2 - y3 + sb * (y4 - y1)
   z2  = y1 - 2.d0 * y2 + y3
   z1  = -sb2 * y1 - y2 + 2.d0 * y3 - sb * y4

   rs = 0.5d0 * dx * ( z1 + dx * ( z2 + dx * z3 ) ) + y2

   else

   rs = SphFunc( 1, j, k, i )
   print *, 'Macrosphere Bad contact'

   end if

#else

   xx  = (x - drmin) * invdr

   ii  = nint( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = SphFunc( ii  , j, k, i )
   y2  = SphFunc( ii+1, j, k, i )
   y3  = SphFunc( ii+2, j, k, i )

   z2  = y1 - 2.d0 * y2 + y3
   z1  = - y1 + y3
   rs  = 0.5d0 * (z1 + z2 * dx ) * dx + y2

   else

   rs = SphFunc( 1, j, k, i )
   print *, 'Macrosphere Bad contact'

   end if

#endif

end Function InterPolateMacroSph


!######################################################################
!######################################################################


Function InterPolateMacroMacro(x,j,itype,jtype) Result(rs)

use CGball

implicit none

integer :: ii, j, itype, jtype
real(8) :: x, xx, dx, y1, y2, y3
real(8) :: z1, z2
real(8) :: rs

#ifdef CUBIC
real(8) :: y4, z3
real(8), parameter :: sb = 0.3333333333333333d0
real(8), parameter :: sb2 = 0.6666666666666667d0

   xx  = (x - Rminmm2(itype,jtype)) * InvddMM(itype,jtype)

   ii  = int( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = UFtableMM( ii  , j, itype, jtype )
   y2  = UFtableMM( ii+1, j, itype, jtype )
   y3  = UFtableMM( ii+2, j, itype, jtype )
   y4  = UFtableMM( ii+3, j, itype, jtype )

   z3  = y2 - y3 + sb * (y4 - y1)
   z2  = y1 - 2.d0 * y2 + y3
   z1  = -sb2 * y1 - y2 + 2.d0 * y3 - sb * y4

   rs = 0.5d0 * dx * ( z1 + dx * ( z2 + dx * z3 ) ) + y2

   else

   rs = UFtableMM( 1, j, itype, jtype )
   print *, 'Macrosphere Bad contact'

   end if

#else

   xx  = (x - Rminmm2(itype,jtype)) * InvddMM(itype,jtype)

   ii  = nint( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = UFtableMM( ii  , j, itype, jtype )
   y2  = UFtableMM( ii+1, j, itype, jtype )
   y3  = UFtableMM( ii+2, j, itype, jtype )

   z2  = y1 - 2.d0 * y2 + y3
   z1  = - y1 + y3
   rs  = 0.5d0 * (z1 + z2 * dx ) * dx + y2

   else

   rs = UFtableMM( 1, j, itype, jtype )
   print *, 'Macrosphere Bad contact'

   end if

#endif

end Function InterPolateMacroMacro


!######################################################################
!######################################################################


subroutine Set_Sphere(iflag)

use CGball
use UnitExParam, only : pi, kb
use CommonBlocks, only : ForceField, QMaster, Qstdout
use Numbers, only : N
use CommonMPI, only : NProcs, MyRank
use CGdata, only : NumAtype, AtomTypeList
use AtomParam, only : Mass, ResidName, AtomName
use IOparam, only : Sphere_file

implicit none

integer :: i, ii, j, k, iflag, Mac, itype
integer :: eofile, Nas
integer, parameter :: Ngrid = 1000
! ## assuming carbon
real(8), parameter :: MassCarbon = 12.011d0
real(8) :: Rls, deltaD, dd, dr, guesn
real(8) :: U, F, Fr
real(8), dimension(:), allocatable :: eps_ws
real(8), dimension(Ngrid) :: distance
character(len=100) :: String, String1
character(len=80), dimension(:,:), allocatable :: filename_sph
character(len=80) :: temp_filename
#ifdef DEBUG
character(len=80) :: filename_type
#endif
character(len=6) :: temp_atomtype
logical, dimension(:,:), allocatable :: flag_assigned
real(8) :: Rcutmm, Rbkmm

if(iflag == 1) then

   allocate( IDsphere(N) )
   IDsphere = 0
   allocate( IdSph(NumSphere) )
   allocate( ItypeSph(NumSphere) )

   ii = 0
   do i = 1, N
     if(ResidName(i)=='CGSP') then
       ii = ii + 1
       read(AtomName(i),*) ItypeSph(ii)
       Mass(i) = MassCarbon * RhoSphere * 4.d0 * pi * RadiusSphere(ItypeSph(ii))**3 / 3.d0
       IDsphere(i) = ii
       IdSph(ii) = i
     end if
   end do

   if(ii/=NumSphere) then
     if(QMaster) then
       write(*,*) "Number of colloid (CGSP) particles is not consistent with the number"
       write(*,*) "given in 'input.data'"
     end if
     call Finalize
   end if

   MaxPairMac = 0
   do i = 1, NumTypeSphere
     guesn = 4.d0 * pi * ((RadiusSphere(i)+20.d0)**3 - RadiusSphere(i)**3) / 3.d0 * 0.6022045 / 18.d0
     Mac = int(guesn) + 1000
     if(Mac > MaxPairMac) MaxPairMac = Mac
   end do
   allocate( NmacList(NumSphere) )
   allocate( Listmac(2,MaxPairMac,NumSphere) )
   allocate( FSphRs(NumSphere) )
   allocate( PMFball(NumSphere) )

   PMFball(:) = 0.d0

   if(ForceField(1:4)=='OPLS') then

     NumSphTab = 1
     allocate( InvddSph(NumSphTab,NumTypeSphere) )
     allocate( SphFunc(Ngrid,3,NumSphTab,NumTypeSphere) )
     allocate( Rcutmacro(NumSphTab,NumTypeSphere) )
     allocate( Rcutmacro2(NumSphTab,NumTypeSphere) )
     allocate( Rlstmacro2(NumSphTab,NumTypeSphere) )
     allocate( Rminmacro(NumSphTab,NumTypeSphere) )
     allocate( Rminmacro2(NumSphTab,NumTypeSphere) )

     Rcutmacro(1,1) = RadiusSphere(1) + 20.d0
     Rcutmacro2(1,1) = Rcutmacro(1,1) * Rcutmacro(1,1)
     Rls = Rcutmacro(1,1) + 4.d0
     Rlstmacro2(1,1) = Rls * Rls
     Rminmacro(1,1) = RadiusSphere(1)
     Rminmacro2(1,1) = Rminmacro(1,1) * Rminmacro(1,1)

     deltaD = (Rcutmacro2(1,1) - Rminmacro2(1,1)) / (Ngrid-3.d0)
     InvddSph(1,1) = 1.d0 / deltaD

     dd = Rminmacro2(1,1) - deltaD
     do i = 1, Ngrid
       dd = dd + deltaD
       dr = sqrt(dd)
! ## just to prevent from going to infinity
       if((dr-RadiusSphere(1))<1.5) then
         dr = RadiusSphere(1) + 1.5
       end if
       call AAwall_SPCE(U,F,Fr,dr,eps_CO)
       SphFunc(i,1,1,1) = U
       SphFunc(i,2,1,1) = F / dr
       SphFunc(i,3,1,1) = Fr
     end do

     Nas = NProcs - MyRank

     allocate( IdOx(N) )
     Nox = 0
     do i = Nas, N, NProcs
       if(AtomName(i)=='OH2') then
         Nox = Nox + 1
         IdOx(Nox) = i
       end if
     end do

   end if

else if(iflag == 2) then
! ## for CG only

     NumSphTab = NumAtype
     allocate( InvddSph(NumSphTab,NumTypeSphere) )
     allocate( SphFunc(Ngrid,3,NumSphTab,NumTypeSphere) )
     allocate( Rcutmacro(NumSphTab,NumTypeSphere) )
     allocate( Rcutmacro2(NumSphTab,NumTypeSphere) )
     allocate( Rlstmacro2(NumSphTab,NumTypeSphere) )
     allocate( Rminmacro(NumSphTab,NumTypeSphere) )
     allocate( Rminmacro2(NumSphTab,NumTypeSphere) )

     allocate( filename_sph(NumSphTab,NumTypeSphere) )
     allocate( flag_assigned(NumSphTab,NumTypeSphere) )

     allocate( eps_ws(NumTypeSphere) )

     flag_assigned(:,:) = .False.

     open(71,file=trim(Sphere_file),status='old')

     itype = 0

     do

       read(71,'(a100)',iostat=eofile) String1
       if(eofile == -1) exit
       String = trim(adjustl(String1))

       if(String(1:5)=='#TYPE') then
         read(String(6:),*) itype
       end if
       if(String(1:1)=='#'.or.String(1:1)=='!'.or.String(1:1)==' ') cycle

       if(itype/=0) then

         read(String,*) temp_atomtype
         if(temp_atomtype=='W') then
           read(String,*) temp_atomtype, eps_ws(itype)
         else
           read(String,*) temp_atomtype, temp_filename
         end if

innsel:  do i = 1, NumAtype
           if(AtomTypeList(i) == temp_atomtype) then
             write(filename_sph(i,itype),'(a)') trim(temp_filename)
             flag_assigned(i,itype) = .True.
             exit innsel
           end if
         end do innsel

       end if

     end do

     close(71)

     if(Qstdout.and.QMaster) print *, 'eps_ws=', eps_ws(:), '[K]'

     do i = 1, NumAtype

       do j = 1, NumTypeSphere

         if(.not.flag_assigned(i,j)) then
           if(QMaster) then
             write(*,*) 'error: no table potential file is assigned for '
             write(*,*) 'macroparticle type ',j,'.vs.', AtomTypeList(i)
           end if
           call Finalize
         end if

         if(AtomTypeList(i)=='W') then ! for water

           Rcutmacro(i,j) = RadiusSphere(j) + 20.d0
           Rcutmacro2(i,j) = Rcutmacro(i,j) * Rcutmacro(i,j)
           Rls = Rcutmacro(i,j) + 2.d0
           Rlstmacro2(i,j) = Rls * Rls
           Rminmacro(i,j) = RadiusSphere(j)
           Rminmacro2(i,j) = Rminmacro(i,j) * Rminmacro(i,j)

           deltaD = (Rcutmacro2(i,j) - Rminmacro2(i,j)) / (Ngrid-3.d0)
           InvddSph(i,j) = 1.d0 / deltaD

           eps_ws(j) = eps_ws(j) * kb
           dd = Rminmacro2(i,j) - deltaD

#ifdef DEBUG
           write(filename_type,'(a,i1,a)') 'interaction_Macro',i,'Water.dat'
           open(121,file=trim(filename_type))
#endif
           do k = 1, Ngrid
             dd = dd + deltaD
             dr = sqrt(dd)
! ##   just to prevent from going to infinity
             if((dr-RadiusSphere(j))<1.5) then
               dr = RadiusSphere(j) + 1.5
             end if
             call CGwall_water(U,F,Fr,dr,eps_ws(j),RadiusSphere(j))

#ifdef DEBUG
             write(121,*) dr, U/kb, F,kb
#endif

             SphFunc(k,1,i,j) = U
             SphFunc(k,2,i,j) = F / dr
             SphFunc(k,3,i,j) = Fr
           end do

#ifdef DEBUG
           close(121)
#endif

         else ! not water

           open(71,file=trim(filename_sph(i,j)),status='old')
           k = 0
           do while(k<Ngrid)
             read(71,'(a100)') String
             if(String(1:1)=='!'.or.String(1:1)=='#') cycle
             k = k + 1
! ## assuming those energies are given in Kelvin
             read(String,*) dr, U, F, Fr
             distance(k) = dr
             SphFunc(k,1,i,j) = U * kb
             SphFunc(k,2,i,j) = F / dr  * kb
             SphFunc(k,3,i,j) = Fr  * kb
           end do
           close(71)

           Rcutmacro(i,j) = distance(Ngrid-2)
           Rminmacro(i,j) = distance(1)
           Rminmacro2(i,j) = Rminmacro(i,j) * Rminmacro(i,j)
           Rcutmacro2(i,j) = Rcutmacro(i,j) * Rcutmacro(i,j)

           Rls = Rcutmacro(i,j) + 2.d0
           Rlstmacro2(i,j) = Rls * Rls

           deltaD = (Rcutmacro2(i,j) - Rminmacro2(i,j)) / (Ngrid-3.d0)
           InvddSph(i,j) = 1.d0 / deltaD

         end if

         if(QMaster.and.Qstdout) then
           print *, 'Type', j, AtomTypeList(i), Rminmacro(i,j), Rcutmacro(i,j)
         end if

       end do

     end do

     deallocate(flag_assigned, filename_sph)

! ## Macro-Macro

     MaxPairMM = 200 * NumSphere
     allocate( UFtableMM(Ngrid,2,NumTypeSphere,NumTypeSphere) )
     allocate( ListMM(3,MaxPairMM) )
     allocate( Rminmm2(NumTypeSphere,NumTypeSphere) )
     allocate( Rcutmm2(NumTypeSphere,NumTypeSphere) )
     allocate( Rbkmm2(NumTypeSphere,NumTypeSphere) )
     allocate( InvddMM(NumTypeSphere,NumTypeSphere) )

     do i = 1, NumTypeSphere
       do j = 1, NumTypeSphere

         Rcutmm = RadiusSphere(i) + RadiusSphere(j) + 20.d0
         Rcutmm2(i,j) = Rcutmm * Rcutmm
         Rbkmm = Rcutmm + 4.d0
         Rbkmm2(i,j) = Rbkmm * Rbkmm
         Rminmm2(i,j) = (RadiusSphere(i) + RadiusSphere(j) + 1.0d0)**2

         deltaD = (Rcutmm2(i,j) - Rminmm2(i,j)) / (Ngrid-3)
         InvddMM(i,j) = 1.d0 / deltaD
         dd = Rminmm2(i,j) - deltaD

#ifdef DEBUG
         write(filename_type,'(a,i1,a,i1,a)') 'interaction_Macro',i,'Macro',j,'.dat'
         open(121,file=trim(filename_type))
#endif

         do k = 1, Ngrid
           dd = dd + deltaD
           dr = sqrt(dd)
           call Set_MacroMacroID(U,F,RadiusSphere(i),RadiusSphere(j),dr)

#ifdef DEBUG
           write(121,*) dr, U/kb, F/kb
#endif

           UFtableMM(k,1,i,j) = U
           UFtableMM(k,2,i,j) = F / dr
         end do

#ifdef DEBUG
         close(121)
#endif
       end do
     end do

end if

Contains

   subroutine CGwall_water(U,F,Fr,d0,ep,rs1)

   implicit none

   real(8), parameter :: a = 0.0338d0! [-]
   real(8), parameter :: b = 0.118d0! [-]
   real(8), parameter :: sg = 4.d0 ! [A]
   real(8) :: ep

   real(8) :: d0, d, U, F, Fr
   real(8) :: pU1, pU2, sg3, sg6, sg9
   real(8) :: rs1, rs2, rs3, rs4, rs5, rs6
   real(8) :: d2, d3, d4, d5, d6, d7
   real(8) :: dd2, dd4, dd6, pT1, pT2
   real(8) :: ds, ds3, ds6, ds7, bU1, bU2
   real(8) :: U1, bF1, bF2, F1, bR1, Fr1, bR2, Fr2

      rs2 = rs1 * rs1
      rs3 = rs2 * rs1
      rs4 = rs2 * rs2
      rs5 = rs4 * rs1
      rs6 = rs4 * rs2

      d  = d0 - rs1
      d2 = d * d
      d3 = d2 * d
      d4 = d2 * d2
      d5 = d3 * d2
      d6 = d3 * d3
      d7 = d6 * d

      dd2 = d0 * d0
      dd4 = dd2 * dd2
      dd6 = dd4 * dd2

      sg3 = sg * sg * sg
      sg6 = sg3 * sg3
      sg9 = sg6 * sg3

      pT1 = 4.d0 * a * ep * sg9 * rs2
      pT2 = 8.d0 * b * ep * sg6 * rs2
      pU1 = pT1 * rs1
      pU2 = pT2 * rs1

      ds  = d0 + rs1
      ds3 = ds * ds * ds
      ds6 = ds3 * ds3
      ds7 = ds6 * ds
      bU1 = 5.d0 * d6 * d0 * ds6
      bU2 = d3 * ds3

      U1 = 35.d0 * d4 + 140.d0 * rs1 * d3 + 252.d0 * rs2 * d2 &
      &  + 224.d0 * rs3 * d + 80.d0 * rs4
      U = pU1 * U1 / bU1 - pU2 / bU2

      bF1 = bU1 * d * d0 * ds
      bF2 = 1.d0 / (d4 * ds3 * ds)
      F1 = 105.d0 * d6 + 630.d0 * d5 * rs1 + 1764.d0 * d4 * rs2 + 2856.d0 * d3 * rs3 &
      &  + 2736.d0 * d2 * rs4 + 1440.d0 * d * rs5 + 320.d0 * rs6
      F  = - 3.d0 * pU1 * F1 / bF1 + 6.d0 * pU2 * d0 * bF2

      bR1 = d0 * d7 * ds7
      Fr1 = 7.d0 * dd6 + 35.d0 * dd4 * rs2 + 21.d0 * dd2 * rs4 + rs6
      bR2 = bF2
      Fr2 = 3.d0 * pT2 * (dd2 + rs2)
      Fr  = 3.d0 * pT1 * Fr1 / bR1 - Fr2 * bR2

   end subroutine CGwall_water

   subroutine AAwall_SPCE(U,F,Fr,d0,ep)

   use UnitExParam, only : pi
   use CGball, only : RadiusSphere, RhoSphere

   implicit none

   real(8) :: ep
   real(8), parameter :: sg = 3.190d0 ! [A]

   real(8) :: d0, d, U, F, Fr
   real(8) :: pU1, pU2, sg3, sg6
   real(8) :: rs1, rs2, rs3, rs4, rs5, rs6, rs7, rs8
   real(8) :: d2, d3, d4, d5, d6, d7, d9, d10
   real(8) :: dd2, dd4, dd6, dd8, pr1, pr2
   real(8) :: ds, ds3, ds4, ds6, ds7, ds9, ds10, bU1, bU2
   real(8) :: U1, bF1, bF2, F1, bR1, Fr1, bR2, Fr2

      rs1 = RadiusSphere(1)
      rs2 = rs1 * rs1
      rs3 = rs2 * rs1
      rs4 = rs2 * rs2
      rs5 = rs4 * rs1
      rs6 = rs4 * rs2
      rs7 = rs4 * rs3
      rs8 = rs4 * rs4

      d  = d0 - rs1
      d2 = d * d
      d3 = d2 * d
      d4 = d2 * d2
      d5 = d3 * d2
      d6 = d3 * d3
      d7 = d6 * d
      d9 = d6 * d3
      d10 = d6 * d4

      dd2 = d0 * d0
      dd4 = dd2 * dd2
      dd6 = dd4 * dd2
      dd8 = dd4 * dd4

      sg3 = sg * sg * sg
      sg6 = sg3 * sg3

      pU2 = 16.d0 * pi * ep * RhoSphere * sg6  * rs3
      pU1 = pU2 * sg6

      ds  = d0 + rs1
      ds3 = ds * ds * ds
      ds4 = ds3 * ds
      ds6 = ds3 * ds3
      ds7 = ds6 * ds
      ds9 = ds6 * ds3
      ds10 = ds9 * ds
      bU1 = 45.d0 * d9 * ds9
      bU2 = 3.d0 * d3 * ds3

      U1 = 15.d0 * d6 + 90.d0 * rs1 * d5 + 288.d0 * rs2 * d4 &
      &  + 552.d0 * rs3 * d3 + 648.d0 * rs4 * d2 + 432.d0 * rs5 * d + 128.d0 * rs6
      U = pU1 * U1 / bU1 - pU2 / bU2

      bF1 = 5.d0 * d10 * ds10
      bF2 = d4 * ds4
      F1 = 5.d0 * d7 + 35.d0 * d6 * rs1 + 132.d0 * d5 * rs2 + 310.d0 * d4 * rs3 &
      &  + 472.d0 * d3 * rs4 + 456.d0 * d2 * rs5 + 256.d0 * d * rs6 + 64.d0 * rs7
      F  = - 4.d0 * pU1 * F1 / bF1 + 2.d0 * pU2 * d0 / bF2

      bR1 = bF1
      Fr1 = 5.d0 * dd8 + 60.d0 * dd6 * rs2 + 126.d0 * dd4 * rs4 &
      &   + 60.d0 * dd2 * rs6 + 5.d0 * rs8
      bR2 = bF2
      pr2 = 16.d0 * pi * ep * RhoSphere * sg6 * rs2
      pr1 = pr2 * sg6
      Fr2 = pr2 * (dd2 + rs2)
      Fr  = pr1 * Fr1 / bR1 - Fr2 / bR2

   end subroutine AAwall_SPCE

   subroutine Set_MacroMacroID(U,F,Rad1,Rad2,dist)

   use UnitExParam, only : pi
   use CGball, only : epsMM, sigMM, RhoSphere

   implicit none

   real(8) :: U, F, Rad1, Rad2
   real(8) :: dist, Rs
   real(8) :: Rs2, Rs3, Rs4, Rs5, Rs6, Rs7, Rs8, Rs9
   real(8) :: Rs10, Rs11, Rs12, Rs13, Rs14, Rs15, Rs16
   real(8) :: Ra, Ra2, Ra3, Ra4, Ra5, Ra6, Ra7, Ra8, Ra9, Ra10, Ra11
   real(8) :: Rb, Rb2, Rb3, Rb4, Rb5, Rb6, Rb7, Rb8, Rb9, Rb10, Rb11
   real(8), dimension(-2:2) :: d, d2, d3, d4, d5, d6, d7, d8, d9
   real(8), dimension(-2:2) :: d10, d11, d12, d13, d14, d15, d16
   real(8), dimension(-2:2) :: p2, p2a, p2b, p3, t1, t2, t3, t4, t5
   real(8), dimension(-2:2) :: a1, a2, a
   real(8), parameter :: delta = 0.01
   real(8), parameter :: delta2 = 2.d0 * delta
   real(8) :: p1, s2, s6, s12
   real(8) :: hh

   if(Rad1==Rad2) then

      Rs  = Rad1
      Rs2 = Rs  * Rs
      Rs3 = Rs2 * Rs
      Rs4 = Rs2 * Rs2
      Rs5 = Rs2 * Rs3
      Rs6 = Rs3 * Rs3
      Rs7 = Rs3 * Rs4
      Rs8 = Rs4 * Rs4
      Rs9 = Rs4 * Rs5
      Rs10= Rs5 * Rs5
      Rs11= Rs5 * Rs6
      Rs12= Rs6 * Rs6
      Rs13= Rs6 * Rs7
      Rs14= Rs7 * Rs7
      Rs15= Rs7 * Rs8
      Rs16= Rs8 * Rs8

      p1  = pi * pi * RhoSphere * RhoSphere * epsMM
      s2  = sigMM * sigMM
      s6  = s2 * s2 * s2
      s12 = s6 * s6

      hh = dist - 2.d0*Rs
      d(-2) = hh - delta2
      d(-1) = hh - delta
      d( 0) = hh
      d( 1) = hh + delta
      d( 2) = hh + delta2

      p2(:)  = (2.d0 * Rs + d(:)) * (2.d0 * Rs + d(:))
      p3(:)  = 4.d0 * Rs + d(:)
      d2(:)  = d(:)  * d(:)
      d3(:)  = d2(:) * d(:)
      d4(:)  = d2(:) * d2(:)
      d5(:)  = d3(:) * d2(:)
      d6(:)  = d3(:) * d3(:)
      d7(:)  = d3(:) * d4(:)
      d8(:)  = d4(:) * d4(:)
      d9(:)  = d5(:) * d4(:)
      d10(:) = d5(:) * d5(:)
      d11(:) = d6(:) * d5(:)
      d12(:) = d6(:) * d6(:)
      d13(:) = d7(:) * d6(:)
      d14(:) = d7(:) * d7(:)
      d15(:) = d8(:) * d7(:)
      d16(:) = d8(:) * d8(:)

      t1(:) = 2.d0 * p1 * s6 / (3.d0 * d(:) * p2(:) * p3(:))
      t2(:) = - 4.d0 * Rs2 * (d2(:) + 4.d0 * Rs * d(:) + 2.d0 * Rs2)
      t3(:) = d(:) * ( d3(:) + 8.d0 * Rs * d2(:) + 20.d0 * Rs2 * d(:) &
      &     + 16.d0 * Rs3 ) * log( p2(:) / (d(:) * p3(:)) )

      a1(:) = t1(:) * (t2(:) + t3(:))

      t4(:) = 2.0 * p1 * s12 * 32.0 * Rs6 / &
      &       ( 4725.0 * d7(:) * p2(:)**7 * p3(:)**7 )

      t5(:) = 525.0 * d16(:) + 16800.0 * Rs * d15(:) + 251160.0 * Rs2 * d14(:) &
      &     + 2328480.0 * Rs3 * d13(:) + 14986776.0 * Rs4 * d12(:)             &
      &     + 71045184.0 * Rs5 * d11(:) + 256802336.0 * Rs6 * d10(:)           &
      &     + 722727040.0 * Rs7 * d9(:) + 1602330496.0 * Rs8 * d8(:)           &
      &     + 2811136000.0 * Rs9 * d7(:) + 3893728768.0 * Rs10 * d6(:)         &
      &     + 4216440832.0 * Rs11 * d5(:) + 3500892160.0 * Rs12 * d4(:)        &
      &     + 2154790912.0 * Rs13 * d3(:) + 927072256.0 * Rs14 * d2(:)         &
      &     + 249036800.0 * Rs15 * d(:) + 31457280.0 * Rs16

      a2(:) = t4(:)*t5(:)

      a(:) = a1(:) + a2(:)

      U = a(0)
      F = (0.25d0 * a(-2) - 2.d0 * a(-1) + 2.d0 * a(1) - 0.25d0 * a(2)) &
      & / (3.d0 * delta)

   else

      Ra  = Rad1
      Ra2 = Ra  * Ra
      Ra3 = Ra2 * Ra
      Ra4 = Ra2 * Ra2
      Ra5 = Ra2 * Ra3
      Ra6 = Ra3 * Ra3
      Ra7 = Ra3 * Ra4
      Ra8 = Ra4 * Ra4
      Ra9 = Ra4 * Ra5
      Ra10= Ra5 * Ra5
      Ra11= Ra5 * Ra6
      Rb  = Rad2
      Rb2 = Rb  * Rb
      Rb3 = Rb2 * Rb
      Rb4 = Rb2 * Rb2
      Rb5 = Rb2 * Rb3
      Rb6 = Rb3 * Rb3
      Rb7 = Rb3 * Rb4
      Rb8 = Rb4 * Rb4
      Rb9 = Rb4 * Rb5
      Rb10= Rb5 * Rb5
      Rb11= Rb5 * Rb6

      p1  = pi * pi * RhoSphere * RhoSphere * epsMM
      s2  = sigMM * sigMM
      s6  = s2 * s2 * s2
      s12 = s6 * s6

      hh = dist - Rad1 - Rad2
      d(-2) = hh - delta2
      d(-1) = hh - delta
      d( 0) = hh
      d( 1) = hh + delta
      d( 2) = hh + delta2

      d2(:)  = d(:)  * d(:)
      d3(:)  = d2(:) * d(:)
      d4(:)  = d2(:) * d2(:)
      d5(:)  = d3(:) * d2(:)
      d6(:)  = d3(:) * d3(:)
      d7(:)  = d3(:) * d4(:)
      d8(:)  = d4(:) * d4(:)
      d9(:)  = d5(:) * d4(:)
      d10(:) = d5(:) * d5(:)
      d11(:) = d6(:) * d5(:)
      d12(:) = d6(:) * d6(:)
      d13(:) = d7(:) * d6(:)
      d14(:) = d7(:) * d7(:)
      d15(:) = d8(:) * d7(:)
      d16(:) = d8(:) * d8(:)

      p2a(:)  = 2.d0 * Ra + d(:)
      p2b(:)  = 2.d0 * Rb + d(:)
      p3(:)   = 2.d0 * Ra + 2.d0 * Rb + d(:)

      t2(:) = 4.d0 * Ra * Rb * ( 2.d0 * Ra * Rb + 2.d0 * Ra * d(:) &
      &     + 2.d0 * Rb * d(:) + d2(:) ) / ( d(:) * p2a(:) * p2b(:) * p3(:) )

      t3(:) = log( d(:) * p3(:) / ( p2a(:) * p2b(:) ) )

      a1(:) = (-2.d0 * p1 * s6 / 3.d0) * (t2(:) + t3(:))

      t4(:) = p1 * s12 * 64.0 * Ra3 * Rb3 / &
      &       ( 4725.0 * d7(:) * p2a(:)**7 * p2b(:)**7 * p3(:)**7 )

      t5(:) = 491520.0 * Ra11 * Rb5 + 2949120.0 * Ra10 * Rb6                      &
      &     + 7372800.0 * Ra9 * Rb7 + 9830400.0 * Ra8 * Rb8                       &
      &     + 7372800.0 * Ra7 * Rb9 + 2949120.0 * Ra6 * Rb10                      &
      &     + 491520.0 * Ra5 * Rb11 + 1638400.0 * Ra11 * Rb4 * d(:)               &
      &     + 12697600.0 * Ra10 * Rb5 * d(:) + 40550400.0 * Ra9 * Rb6 * d(:)      &
      &     + 69632000.0 * Ra8 * Rb7 * d(:) + 69632000.0 * Ra7 * Rb8 * d(:)       &
      &     + 40550400.0 * Ra6 * Rb9 * d(:) + 12697600.0 * Ra5 * Rb10 * d(:)      &
      &     + 1638400.0 * Ra4 * Rb11 * d(:) + 2293760.0 * Ra11 * Rb3 * d2(:)      &
      &     + 23322624.0 * Ra10 * Rb4 * d2(:) + 95412224.0 * Ra9 * Rb5 * d2(:)    &
      &     + 208445440.0 * Ra8 * Rb6 * d2(:) + 268124160.0 * Ra7 * Rb7 * d2(:)   &
      &     + 208445440.0 * Ra6 * Rb8 * d2(:) + 95412224.0 * Ra5 * Rb9 * d2(:)    &
      &     + 23322624.0 * Ra4 * Rb10 * d2(:) + 2293760.0 * Ra3 * Rb11 * d2(:)    &
      &     + 1720320.0 * Ra11 * Rb2 * d3(:) + 23711744.0 * Ra10 * Rb3 * d3(:)    &
      &     + 126230528.0 * Ra9 * Rb4 * d3(:) + 351596544.0 * Ra8 * Rb5 * d3(:)   &
      &     + 574136320.0 * Ra7 * Rb6 * d3(:) + 574136320.0 * Ra6 * Rb7 * d3(:)   &
      &     + 351596544.0 * Ra5 * Rb8 * d3(:) + 126230528.0 * Ra4 * Rb9 * d3(:)   &
      &     + 23711744.0 * Ra3 * Rb10 * d3(:) + 1720320.0 * Ra2 * Rb11 * d3(:)    &
      &     + 716800.0 * Ra11 * Rb * d4(:) + 14350336.0 * Ra10 * Rb2 * d4(:)      &
      &     + 102932480.0 * Ra9 * Rb3 * d4(:) + 371501056.0 * Ra8 * Rb4 * d4(:)   &
      &     + 771573760.0 * Ra7 * Rb5 * d4(:) + 978743296.0 * Ra6 * Rb6 * d4(:)   &
      &     + 771573760.0 * Ra5 * Rb7 * d4(:) + 371501056.0 * Ra4 * Rb8 * d4(:)   &
      &     + 102932480.0 * Ra3 * Rb9 * d4(:) + 14350336.0 * Ra2 * Rb10 * d4(:)   &
      &     + 716800.0 * Ra * Rb11 * d4(:) + 143360.0 * Ra11 * d5(:)              &
      &     + 5053440.0 * Ra10 * Rb * d5(:) + 52699136.0 * Ra9 * Rb2 * d5(:)      &
      &     + 255761408.0 * Ra8 * Rb3 * d5(:) + 687083520.0 * Ra7 * Rb4 * d5(:)   &
      &     + 1107479552.0 * Ra6 * Rb5 * d5(:) + 1107479552.0 * Ra5 * Rb6 * d5(:) &
      &     + 687083520.0 * Ra4 * Rb7 * d5(:) + 255761408.0 * Ra3 * Rb8 * d5(:)   &
      &     + 52699136.0 * Ra2 * Rb9 * d5(:) + 5053440.0 * Ra * Rb10 * d5(:)      &
      &     + 143360.0 * Rb11 * d5(:) + 842240.0 * Ra10 * d6(:)                   &
      &     + 16056320.0 * Ra9 * Rb * d6(:) + 113960448.0 * Ra8 * Rb2 * d6(:)     &
      &     + 411944960.0 * Ra7 * Rb3 * d6(:) + 858627584.0 * Ra6 * Rb4 * d6(:)   &
      &     + 1090865664.0 * Ra5 * Rb5 * d6(:) + 858627584.0 * Ra4 * Rb6 * d6(:)  &
      &     + 411944960.0 * Ra3 * Rb7 * d6(:) + 113960448.0 * Ra2 * Rb8 * d6(:)   &
      &     + 16056320.0 * Ra * Rb9 * d6(:) + 842240.0 * Rb10 * d6(:)             &
      &     + 2293760.0 * Ra9 * d7(:) + 30571520.0 * Ra8 * Rb * d7(:)             &
      &     + 162821120.0 * Ra7 * Rb2 * d7(:) + 457630208.0 * Ra6 * Rb3 * d7(:)   &
      &     + 752251392.0 * Ra5 * Rb4 * d7(:) + 752251392.0 * Ra4 * Rb5 * d7(:)   &
      &     + 457630208.0 * Ra3 * Rb6 * d7(:) + 162821120.0 * Ra2 * Rb7 * d7(:)   &
      &     + 30571520.0 * Ra * Rb8 * d7(:) + 2293760.0 * Rb9 * d7(:)             &
      &     + 3821440.0 * Ra8 * d8(:) + 38989440.0 * Ra7 * Rb * d8(:)             &
      &     + 162674176.0 * Ra6 * Rb2 * d8(:) + 361589760.0 * Ra5 * Rb3 * d8(:)   &
      &     + 468180864.0 * Ra4 * Rb4 * d8(:) + 361589760.0 * Ra3 * Rb5 * d8(:)   &
      &     + 162674176.0 * Ra2 * Rb6 * d8(:) + 38989440.0 * Ra * Rb7 * d8(:)     &
      &     + 3821440.0 * Rb8 * d8(:) + 4332160.0 * Ra7 * d9(:)                   &
      &     + 35156800.0 * Ra6 * Rb * d9(:) + 116807040.0 * Ra5 * Rb2 * d9(:)     &
      &     + 205067520.0 * Ra4 * Rb3 * d9(:) + 205067520.0 * Ra3 * Rb4 * d9(:)   &
      &     + 116807040.0 * Ra2 * Rb5 * d9(:) + 35156800.0 * Ra * Rb6 * d9(:)     &
      &     + 4332160.0 * Rb7 * d9(:) + 3515680.0 * Ra6 * d10(:)                  &
      &     + 22989120.0 * Ra5 * Rb * d10(:) + 60682272.0 * Ra4 * Rb2 * d10(:)    &
      &     + 82428192.0 * Ra3 * Rb3 * d10(:) + 60682272.0 * Ra2 * Rb4 * d10(:)   &
      &     + 22989120.0 * Ra * Rb5 * d10(:) + 3515680.0 * Rb6 * d10(:)           &
      &     + 2089920.0 * Ra5 * d11(:) + 10956960.0 * Ra4 * Rb * d11(:)           &
      &     + 22475712.0 * Ra3 * Rb2 * d11(:) + 22475712.0 * Ra2 * Rb3 * d11(:)   &
      &     + 10956960.0 * Ra * Rb4 * d11(:) + 2089920.0 * Rb5 * d11(:)           &
      &     + 913080.0 * Ra4 * d12(:) + 3745560.0 * Ra3 * Rb * d12(:)             &
      &     + 5669496.0 * Ra2 * Rb2 * d12(:) + 3745560.0 * Ra * Rb3 * d12(:)      &
      &     + 913080.0 * Rb4 * d12(:) + 288120.0 * Ra3 * d13(:)                   &
      &     + 876120.0 * Ra2 * Rb * d13(:) + 876120.0 * Ra * Rb2 * d13(:)         &
      &     + 288120.0 * Rb3 * d13(:) + 62580.0 * Ra2 * d14(:)                    &
      &     + 126000.0 * Ra * Rb * d14(:) + 62580.0 * Rb2 * d14(:)                &
      &     + 8400.0 * Ra * d15(:) + 8400.0 * Rb * d15(:) + 525.0 * d16(:)

      a2(:) = t4(:)*t5(:)

      a(:) = a1(:) + a2(:)

      U = a(0)
      F = (0.25d0 * a(-2) - 2.d0 * a(-1) + 2.d0 * a(1) - 0.25d0 * a(2)) &
      & / (3.d0 * delta)

   end if

   end subroutine Set_MacroMacroID

end subroutine Set_Sphere


!######################################################################
!######################################################################


subroutine MacroPMFsample(istep)

use CGball
use TimeParam, only : Nstep, lk
use UnitExParam, only : cvol
use CommonBlocks, only : Job_name, QMaster

implicit none

integer :: istep
real(8), dimension(NumSphere) :: PMFtemp
character(len=58) :: Filename

   if(istep==lk.and.QMaster) then
     NMacSample = 0
     write(Filename,'(a,a)') trim(adjustl(Job_name)),'_PMF.dat'
     open(111,file=trim(Filename))
   end if

   NMacSample = NMacSample + 1
   PMFball(:) = PMFball(:) + FSphRs(:)

   if(mod(NMacSample,100)==0) then
     PMFtemp(:) = PMFball(:) / NMacSample
     call SumPMF_CGball(PMFtemp,NumSphere)
     if(QMaster) write(111,'(i10,5e16.8)') NMacSample, PMFtemp(:) * cvol
   else if(istep==Nstep) then
     PMFtemp(:) = PMFball(:) / NMacSample
     call SumPMF_CGball(PMFtemp,NumSphere)
     if(QMaster) write(111,'(i10,5e16.8)') NMacSample, PMFtemp(:) * cvol
   end if

   if(istep==Nstep.and.QMaster) close(111)

end subroutine MacroPMFsample
