! ForCG 
! ############################
! ## SUBROUTINE LIST 
! ## -- Force_nonbond_CG_long  
! ## -- Force_nonbond_CG_short 
! ## -- Force_nonbond_CG_PBC 
! ## -- Force_nonbond_Ortho_CELL 
! ## -- Force_nonbond_short 
! ## -- PreventContact 
! ## -- NONBONDF 
! ## -- COULOMBFR 
! ## -- COULOMBFiso 
! ## -- ScaledCoordinate 
! ## -- Force_nonbond_iso_long 
! ## -- Force_nonbond_iso_short 
! ############################


!######################################################################
!######################################################################


subroutine Force_nonbond_CG_long(istep)

use Numbers, only : N
use CommonBlocks, only : QPBC, QCellList, QMacro
use NonbondParam, only : Frc_NBlong, Vir_NBlong, Ene_NBlong, Ene_ELlong
use CellParam, only : H, CellShape
use TimeParam, only : BookFreq
use CGball, only : NumSphere
use CommonMPI

implicit none

integer :: istep, Nas
#ifdef GEN
real(8), dimension(3,N) :: ScR
#endif
real(8) :: clhx, clhy, clhz

   Frc_NBlong = 0.d0
   Ene_NBlong = 0.d0
   Ene_ELlong = 0.d0
   Vir_NBlong = 0.d0

   if(QPBC) then
     if(QCellList) then
       call PBC
       if(QMacro) then
         Nas = NProcs - MyRank
#ifdef GEN
         call ScaledCoordinate(ScR)
         call PairListMacrovsCGParticle(ScR,Nas)
         if(NumSphere>1) call PairListMacroMacro(ScR,Nas)
#else
         clhx = H(1,1)*0.5d0
         clhy = H(2,2)*0.5d0
         clhz = H(3,3)*0.5d0
         call PairListMacrovsCGParticle(clhx,clhy,clhz,Nas)
         if(NumSphere>1) call PairListMacroMacro(clhx,clhy,clhz,Nas)
#endif
       end if
       call Force_nonbond_Ortho_CELL
     else
       if( mod(istep,BookFreq) == 0 ) then
         call PBC
         call PairListCG
       end if
       call Force_nonbond_CG_PBC
     end if
   else
     call Force_nonbond_iso_long
   end if

end subroutine Force_nonbond_CG_long


!######################################################################
!######################################################################


subroutine Force_nonbond_CG_short

use CommonBlocks, only : QPBC, QMacro
use NonbondParam, only : Frc_NBshrt, Vir_NBshrt, Ene_NBshrt, Ene_ELshrt
use CellParam, only : CellShape

implicit none

   Frc_NBshrt = 0.d0
   Ene_NBshrt = 0.d0
   Ene_ELshrt = 0.d0
   Vir_NBshrt = 0.d0

   if(QPBC) then
!     call Force_Subt_CG
     call Force_nonbond_short
   else
!     call Force_Subt_CG
     call Force_nonbond_iso_short
   end if

end subroutine Force_nonbond_CG_short


!######################################################################
!######################################################################


subroutine Force_nonbond_CG_PBC

use Numbers, only : N
use Configuration, only : R
use BookParam, only : Npair, ListIJ, Npair_short, List_shortIJ
use NonbondParam, only : Frc_NBlong, Vir_NBlong, Charge, Ene_ELlong, Ene_NBlong
use CellParam, only : CellShft
use CutoffParam, only : Rcutoff2
#ifdef GEN
use CommonBlocks, only : Qcoulomb, QSwitch
use CGdata, only : NBAtomType, Rcut2, Rcut_short2, Rbk_short2, Rheal2, &
&   NBFuncType, Fsw1, CoefAtype, CoefBtype, CoefCtype, Rsw2, Swch, &
&   IDtable, Rmin2
use TableFuncs
#else
use CommonBlocks, only : Qcoulomb
use CGdata, only : NBAtomType, Rcut2, Rcut_short2, Rbk_short2, Rheal2, &
&   NBFuncType, Fsw1, CoefAtype, CoefBtype
#endif

implicit none

integer :: i,j,k,l
integer :: itype, jtype, ftype
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: InvR2, InvR1, InvR3, InvR4, term1, term2, ek, fk, fke
real(8) :: aij, bij
real(8) :: R2, Uij, Ueij, cf, ScF, yr
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
real(8), dimension(:,:),   allocatable :: Gk
#ifdef GEN
real(8) :: Rm, InvR6, InterPolate, xx
real(8) :: eps, R1
integer :: ii
external InterPolate
#endif

   allocate(Gk(3,-13:13))
   Gk = 0.d0

   Npair_short = 0

   do l = 1 , Npair

     i = ListIJ(1,l)
     j = ListIJ(2,l)
     k = ListIJ(3,l)

     Rx = R(1,i) - R(1,j) + CellShft(1,k)
     Ry = R(2,i) - R(2,j) + CellShft(2,k)
     Rz = R(3,i) - R(3,j) + CellShft(3,k)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     itype = NBAtomType(i)
     jtype = NBAtomType(j)

     if(R2 <= Rcut2(itype,jtype)) then

       if(R2 <= Rbk_short2) then
         Npair_short = Npair_short + 1
         List_shortIJ(1,Npair_short) = i
         List_shortIJ(2,Npair_short) = j
         List_shortIJ(3,Npair_short) = k
       end if

       if(R2 > Rcut_short2) then

         ftype = NBFuncType(itype,jtype)

         if(ftype==7) then ! LJ12-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR4
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
#ifdef GEN
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
#endif
           Uij = ek

         else if(ftype==1) then ! LJ9-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR1 = sqrt(InvR2)
           InvR3 = InvR2 * InvR1
           term1 = aij * InvR3 * InvR3 * InvR3
           term2 = bij * InvR3 * InvR3
           ek = term1 - term2
           fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
#ifdef GEN
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
#endif
           Uij = ek

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
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
           Uij = ek

         else if(ftype==2) then ! LJ6-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR2
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
           Uij = ek

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
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
           Uij = ek

         else if(ftype==5) then ! LJ8-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
           Uij = ek

         else if(ftype==6) then ! LJ10-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR2
           term2 = bij * InvR4
           ek = term1 - term2
           fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
           if(QSwitch) then
             call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
           end if
           Uij = ek
#endif
         end if

         if(Qcoulomb) then
           if(R2 <= Rcutoff2) then
             cf = Charge(i) * Charge(j)
             if(cf/=0.) then
               call COULOMBFR(Ueij,fke,R2,cf)
               fk = fk + fke
               Ene_ELlong = Ene_ELlong + Ueij
             end if
           end if
         end if

         Ene_NBlong = Ene_NBlong + Uij

         if(R2 <= Rheal2) then
           yr  = (R2 - Rcut_short2) * Fsw1
           ScF = - yr * yr * (2.d0 * yr - 3.d0)
           fk  = fk * ScF
         end if

         Fx  = fk * Rx
         Fy  = fk * Ry
         Fz  = fk * Rz

         Frc_NBlong(1,i) = Frc_NBlong(1,i) + Fx
         Frc_NBlong(2,i) = Frc_NBlong(2,i) + Fy
         Frc_NBlong(3,i) = Frc_NBlong(3,i) + Fz
         Frc_NBlong(1,j) = Frc_NBlong(1,j) - Fx
         Frc_NBlong(2,j) = Frc_NBlong(2,j) - Fy
         Frc_NBlong(3,j) = Frc_NBlong(3,j) - Fz
         Gk(1,k) = Gk(1,k) + Fx
         Gk(2,k) = Gk(2,k) + Fy
         Gk(3,k) = Gk(3,k) + Fz

       end if

     end if

   end do

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do i = 1, N
     Fx = Frc_NBlong(1,i)
     Fy = Frc_NBlong(2,i)
     Fz = Frc_NBlong(3,i)
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

   Vir_NBlong(1,1) = Vxx
   Vir_NBlong(1,2) = Vxy
   Vir_NBlong(1,3) = Vxz
   Vir_NBlong(2,2) = Vyy
   Vir_NBlong(2,3) = Vyz
   Vir_NBlong(3,3) = Vzz

   Vir_NBlong(2,1) = Vxy
   Vir_NBlong(3,1) = Vxz
   Vir_NBlong(3,2) = Vyz

   deallocate(Gk)

end subroutine Force_nonbond_CG_PBC


!######################################################################
!######################################################################


subroutine Force_nonbond_Ortho_CELL

use Numbers, only : N
use Configuration, only : R
use CommonMPI, only : MyRank, NProcs
use CellListMethod
use NoLJparam
use BookParam, only : MaxPair, Npair_short, List_shortIJ
use NonbondParam, only : Charge, Ene_ELlong, Ene_NBlong, Frc_NBlong, Vir_NBlong
use CellParam, only : CellShft, CellL, InvCL, H, InvH
use CutoffParam, only : Rcutoff2
#ifdef GEN
use CommonBlocks, only : Qcoulomb, QSwitch, QMaster
use CGdata, only : NBAtomType, Rcut2, Rcut_short2, Rbk_short2, Rheal2, &
&   NBFuncType, Fsw1, CoefAtype, CoefBtype, CoefCtype, Rsw2, Swch, &
&   IDtable, Rmin2
use TableFuncs
#else
use CommonBlocks, only : Qcoulomb, QMaster
use CGdata, only : NBAtomType, Rcut2, Rcut_short2, Rbk_short2, Rheal2, &
&   NBFuncType, Fsw1, CoefAtype, CoefBtype
#endif

implicit none

integer :: i, j, k, itype, jtype, Nas, jj
integer :: icell, jcell, jcell0, nabor, ftype
integer :: non, check_bonded
real(8) :: Rix, Riy, Riz
real(8) :: Rx, Ry, Rz
real(8) :: Risx, Risy, Risz
real(8) :: Fx, Fy, Fz
real(8) :: R2, Uij, Ueij, cf, ScF, yr
real(8) :: aij, bij
real(8) :: InvR2, InvR1, InvR3, InvR4, term1, term2, ek, fk, fke
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
real(8) :: Gvx, Gvy, Gvz
real(8), dimension(:,:), allocatable :: Gk
real(8), dimension(:,:), allocatable :: TmpR
#ifdef GEN
integer :: ii
real(8) :: Rm, InvR6, InterPolate, xx
real(8) :: eps, R1
external InterPolate
#endif
external check_bonded

   allocate(Gk(3,-13:13))
   allocate(TmpR(3,N))

   Npair_short = 0
   Gk  = 0.d0

   Nas = NProcs - MyRank

   TmpR = R

! ---------------------
   do i = 1, 3
     CellL(i) = H(i,i)
     InvCL(i) = InvH(i,i)
   end do
! ---------------------
   call LinkCell
! ---------------------

   do icell = Nas , Ncell, NProcs

     i = Head(icell)

     do while( i > 0 )

       non   = NumNoLJ(i)

       Rix   = R(1,i)
       Riy   = R(2,i)
       Riz   = R(3,i)
       itype = NBAtomType(i)

       j = NextP(i)

       do while( j > 0 )

         if(check_bonded(i,j,non) == 1) then
           j = NextP(j)
           cycle
         end if

         jtype = NBAtomType(j)

         Rx = Rix - R(1,j)
         Ry = Riy - R(2,j)
         Rz = Riz - R(3,j)
         R2  = Rx*Rx + Ry*Ry + Rz*Rz

         if(R2 <= Rcut2(itype,jtype)) then

           if(R2 <= Rbk_short2) then
             Npair_short = Npair_short + 1
             List_shortIJ(1,Npair_short) = i
             List_shortIJ(2,Npair_short) = j
             List_shortIJ(3,Npair_short) = SelfShft(i) - SelfShft(j)
           end if

           if(R2 > Rcut_short2) then

             ftype = NBFuncType(itype,jtype)

             if(ftype==7) then ! LJ12-4

               aij = CoefAtype(itype,jtype)
               bij = CoefBtype(itype,jtype)
               InvR2 = 1.d0 / R2
               InvR4 = InvR2 * InvR2
               term1 = aij * InvR4 * InvR4 * InvR4
               term2 = bij * InvR4
               ek = term1 - term2
               fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
#ifdef GEN
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
#endif
               Uij = ek

             else if(ftype==1) then ! LJ9-6

               aij = CoefAtype(itype,jtype)
               bij = CoefBtype(itype,jtype)
               InvR2 = 1.d0 / R2
               InvR1 = sqrt(InvR2)
               InvR3 = InvR2 * InvR1
               term1 = aij * InvR3 * InvR3 * InvR3
               term2 = bij * InvR3 * InvR3
               ek = term1 - term2
               fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
#ifdef GEN
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
#endif
               Uij = ek

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
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
               Uij = ek

             else if(ftype==2) then ! LJ6-4

               aij = CoefAtype(itype,jtype)
               bij = CoefBtype(itype,jtype)
               InvR2 = 1.d0 / R2
               InvR4 = InvR2 * InvR2
               term1 = aij * InvR4 * InvR2
               term2 = bij * InvR4
               ek = term1 - term2
               fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
               Uij = ek

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
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
               Uij = ek

             else if(ftype==5) then ! LJ8-4

               aij = CoefAtype(itype,jtype)
               bij = CoefBtype(itype,jtype)
               InvR2 = 1.d0 / R2
               InvR4 = InvR2 * InvR2
               term1 = aij * InvR4 * InvR4
               term2 = bij * InvR4
               ek = term1 - term2
               fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
               Uij = ek

             else if(ftype==6) then ! LJ10-4

               aij = CoefAtype(itype,jtype)
               bij = CoefBtype(itype,jtype)
               InvR2 = 1.d0 / R2
               InvR4 = InvR2 * InvR2
               term1 = aij * InvR4 * InvR4 * InvR2
               term2 = bij * InvR4
               ek = term1 - term2
               fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
               if(QSwitch) then
               call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
               end if
               Uij = ek
#endif
             end if

             if(Qcoulomb) then
               if(R2 <= Rcutoff2) then
                 cf = Charge(i) * Charge(j)
                 if(cf/=0.) then
                   call COULOMBFR(Ueij,fke,R2,cf)
                   fk = fk + fke
                   Ene_ELlong = Ene_ELlong + Ueij
                 end if
               end if
             end if

             Ene_NBlong = Ene_NBlong + Uij

             if(R2 <= Rheal2) then
               yr  = (R2 - Rcut_short2) * Fsw1
               ScF = - yr * yr * (2.d0 * yr - 3.d0)
               fk  = fk * ScF
             end if

             Fx  = fk * Rx
             Fy  = fk * Ry
             Fz  = fk * Rz

             Frc_NBlong(1,i) = Frc_NBlong(1,i) + Fx
             Frc_NBlong(2,i) = Frc_NBlong(2,i) + Fy
             Frc_NBlong(3,i) = Frc_NBlong(3,i) + Fz
             Frc_NBlong(1,j) = Frc_NBlong(1,j) - Fx
             Frc_NBlong(2,j) = Frc_NBlong(2,j) - Fy
             Frc_NBlong(3,j) = Frc_NBlong(3,j) - Fz

           end if

         end if

         j = NextP(j)

       end do

       jcell0 = 13 * ( icell - 1 )

       do nabor = 1 , 13

         jj = jcell0 + nabor
         jcell = Map( jj )

         if(jcell > 0) then

           k    = MapShft( jj )
           Risx = Rix + CellShft(1,k)
           Risy = Riy + CellShft(2,k)
           Risz = Riz + CellShft(3,k)

           Gvx  = 0.d0
           Gvy  = 0.d0
           Gvz  = 0.d0

           j = Head(jcell)

           do while( j /= 0 )

             if(check_bonded(i,j,non) == 1) then
               j = NextP(j)
               cycle
             end if

             jtype = NBAtomType(j)

             Rx = Risx - R(1,j)
             Ry = Risy - R(2,j)
             Rz = Risz - R(3,j)
             R2  = Rx*Rx + Ry*Ry + Rz*Rz

             if(R2 <= Rcut2(itype,jtype)) then

               if(R2 <= Rbk_short2) then
                 Npair_short = Npair_short + 1
                 List_shortIJ(1,Npair_short) = i
                 List_shortIJ(2,Npair_short) = j
                 List_shortIJ(3,Npair_short) = k + SelfShft(i) - SelfShft(j)
               end if

               if(R2 > Rcut_short2) then

                 ftype = NBFuncType(itype,jtype)

                 if(ftype==7) then ! LJ12-4

                   aij = CoefAtype(itype,jtype)
                   bij = CoefBtype(itype,jtype)
                   InvR2 = 1.d0 / R2
                   InvR4 = InvR2 * InvR2
                   term1 = aij * InvR4 * InvR4 * InvR4
                   term2 = bij * InvR4
                   ek = term1 - term2
                   fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
#ifdef GEN
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
#endif
                   Uij = ek

                 else if(ftype==1) then ! LJ9-6

                   aij = CoefAtype(itype,jtype)
                   bij = CoefBtype(itype,jtype)
                   InvR2 = 1.d0 / R2
                   InvR1 = sqrt(InvR2)
                   InvR3 = InvR2 * InvR1
                   term1 = aij * InvR3 * InvR3 * InvR3
                   term2 = bij * InvR3 * InvR3
                   ek = term1 - term2
                   fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
#ifdef GEN
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
#endif
                   Uij = ek

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
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
                   Uij = ek

                 else if(ftype==2) then ! LJ6-4

                   aij = CoefAtype(itype,jtype)
                   bij = CoefBtype(itype,jtype)
                   InvR2 = 1.d0 / R2
                   InvR4 = InvR2 * InvR2
                   term1 = aij * InvR4 * InvR2
                   term2 = bij * InvR4
                   ek = term1 - term2
                   fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
                   Uij = ek

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
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
                   Uij = ek

                 else if(ftype==5) then ! LJ8-4

                   aij = CoefAtype(itype,jtype)
                   bij = CoefBtype(itype,jtype)
                   InvR2 = 1.d0 / R2
                   InvR4 = InvR2 * InvR2
                   term1 = aij * InvR4 * InvR4
                   term2 = bij * InvR4
                   ek = term1 - term2
                   fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
                   Uij = ek

                 else if(ftype==6) then ! LJ10-4

                   aij = CoefAtype(itype,jtype)
                   bij = CoefBtype(itype,jtype)
                   InvR2 = 1.d0 / R2
                   InvR4 = InvR2 * InvR2
                   term1 = aij * InvR4 * InvR4 * InvR2
                   term2 = bij * InvR4
                   ek = term1 - term2
                   fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
                   if(QSwitch) then
                   call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),&
                   &    Swch(itype,jtype))
                   end if
                   Uij = ek
#endif
                 end if

                 if(Qcoulomb) then
                   if(R2 <= Rcutoff2) then
                     cf = Charge(i) * Charge(j)
                     if(cf/=0.) then
                       call COULOMBFR(Ueij,fke,R2,cf)
                       fk = fk + fke
                       Ene_ELlong = Ene_ELlong + Ueij
                     end if
                   end if
                 end if

                 Ene_NBlong = Ene_NBlong + Uij

                 if(R2 <= Rheal2) then
                   yr  = (R2 - Rcut_short2) * Fsw1
                   ScF = - yr * yr * (2.d0 * yr - 3.d0)
                   fk  = fk * ScF
                 end if

                 Fx = fk * Rx
                 Fy = fk * Ry
                 Fz = fk * Rz

                 Frc_NBlong(1,i) = Frc_NBlong(1,i) + Fx
                 Frc_NBlong(2,i) = Frc_NBlong(2,i) + Fy
                 Frc_NBlong(3,i) = Frc_NBlong(3,i) + Fz
                 Frc_NBlong(1,j) = Frc_NBlong(1,j) - Fx
                 Frc_NBlong(2,j) = Frc_NBlong(2,j) - Fy
                 Frc_NBlong(3,j) = Frc_NBlong(3,j) - Fz
                 Gvx = Gvx + Fx
                 Gvy = Gvy + Fy
                 Gvz = Gvz + Fz

               end if

             end if

             j = NextP(j)

           end do

           Gk(1,k) = Gk(1,k) + Gvx
           Gk(2,k) = Gk(2,k) + Gvy
           Gk(3,k) = Gk(3,k) + Gvz

         end if

       end do

       i = NextP(i)

     end do

   end do

   if(Npair_short > MaxPair) then
     if(QMaster) write(*,*) 'ERROR : increase MAXPAIR in PBC'
     call Finalize
   end if

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do i = 1, N
     Fx = Frc_NBlong(1,i)
     Fy = Frc_NBlong(2,i)
     Fz = Frc_NBlong(3,i)
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

   Vir_NBlong(1,1) = Vxx
   Vir_NBlong(1,2) = Vxy
   Vir_NBlong(1,3) = Vxz
   Vir_NBlong(2,2) = Vyy
   Vir_NBlong(2,3) = Vyz
   Vir_NBlong(3,3) = Vzz

   Vir_NBlong(2,1) = Vxy
   Vir_NBlong(3,1) = Vxz
   Vir_NBlong(3,2) = Vyz

   R = TmpR

   deallocate(Gk,TmpR)

end subroutine Force_nonbond_Ortho_CELL


!######################################################################
!######################################################################


subroutine Force_nonbond_short

use Numbers, only : N
use Configuration, only : R
use BookParam, only : Npair_short, List_shortIJ
use NonbondParam, only : Frc_NBshrt, Vir_NBshrt
use CellParam, only : CellShft
use Numbers, only : N
use NonbondParam, only : Charge, Ene_ELshrt, Ene_NBshrt
#ifdef GEN
use CommonBlocks, only : Qcoulomb, QSwitch
use CGdata, only : NBAtomType, Rcut_short2, Rheal2, &
&   NBFuncType, Fsw1, CoefAtype, CoefBtype, Rsw2, Swch, &
&   IDtable, Rmin2, Rcut2, CoefCtype
use TableFuncs
#else
use CommonBlocks, only : Qcoulomb
use CGdata, only : NBAtomType, NBFuncType, Rcut_short2, Rmin2, Rheal2, Fsw1, &
&   CoefAtype, CoefBtype
#endif

implicit none

real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: cf, ScF, yr, y3
integer :: i,j,k,l,itype,jtype, ftype
real(8) :: Uij, Ueij, R2, aij, bij, InvR2, InvR1
real(8) :: term1, term2, fk, InvR4, InvR3, ek, fke
real(8), dimension(:,:), allocatable :: Forcer
real(8), dimension(:,:), allocatable :: Gk
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
#ifdef GEN
integer :: ii
real(8) :: Rm, InvR6, InterPolate, xx, eps, R1
external InterPolate
#endif

   allocate(Forcer(3,N))
   allocate(Gk(3,-13:13))

   Forcer = 0.d0
   Gk = 0.d0

   do l = 1 , Npair_short

     i = List_shortIJ(1,l)
     j = List_shortIJ(2,l)
     k = List_shortIJ(3,l)

     Rx = R(1,i) - R(1,j) + CellShft(1,k)
     Ry = R(2,i) - R(2,j) + CellShft(2,k)
     Rz = R(3,i) - R(3,j) + CellShft(3,k)
     R2 = Rx * Rx + Ry * Ry + Rz * Rz

     if(R2 <= Rheal2) then

       itype = NBAtomType(i)
       jtype = NBAtomType(j)

#ifdef EQUILIBRATION
       if(R2 <= Rmin2(itype,jtype)) call PreventContact(Rx,Ry,Rz,R2,Rmin2(itype,jtype))
#endif

       ftype = NBFuncType(itype,jtype)

       if(ftype==7) then ! LJ12-4

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)
         InvR2 = 1.d0 / R2
         InvR4 = InvR2 * InvR2
         term1 = aij * InvR4 * InvR4 * InvR4
         term2 = bij * InvR4
         ek = term1 - term2
         fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
#ifdef GEN
         if(QSwitch) then
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
#endif
         Uij = ek

       else if(ftype==1) then ! LJ9-6

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)
         InvR2 = 1.d0 / R2
         InvR1 = sqrt(InvR2)
         InvR3 = InvR2 * InvR1
         term1 = aij * InvR3 * InvR3 * InvR3
         term2 = bij * InvR3 * InvR3
         ek = term1 - term2
         fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
#ifdef GEN
         if(QSwitch) then
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
#endif
         Uij = ek

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
         if(QSwitch) then
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
         Uij = ek

       else if(ftype==2) then ! LJ6-4

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)
         InvR2 = 1.d0 / R2
         InvR4 = InvR2 * InvR2
         term1 = aij * InvR4 * InvR2
         term2 = bij * InvR4
         ek = term1 - term2
         fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
         if(QSwitch) then
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
         Uij = ek

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
         if(QSwitch) then 
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
         Uij = ek

       else if(ftype==5) then ! LJ8-4

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)
         InvR2 = 1.d0 / R2
         InvR4 = InvR2 * InvR2
         term1 = aij * InvR4 * InvR4
         term2 = bij * InvR4
         ek = term1 - term2
         fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
         if(QSwitch) then
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
         Uij = ek

       else if(ftype==6) then ! LJ10-4

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)
         InvR2 = 1.d0 / R2
         InvR4 = InvR2 * InvR2
         term1 = aij * InvR4 * InvR4 * InvR2
         term2 = bij * InvR4
         ek = term1 - term2
         fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
         if(QSwitch) then
         call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
         end if
         Uij = ek
#endif
       end if

       if(Qcoulomb) then
         cf = Charge(i) * Charge(j)
         if(cf/=0.) then
           call COULOMBFR(Ueij,fke,R2,cf)
           fk = fk + fke
           if(R2 <= Rcut_short2) Ene_ELshrt = Ene_ELshrt + Ueij
         end if
       end if

       if(R2 > Rcut_short2) then
         yr  = (R2 - Rcut_short2) * Fsw1
         y3  = yr * yr * (2.d0 * yr - 3.d0)
         ScF = 1.d0 + y3
         fk  = fk * ScF
       else
         Ene_NBshrt = Ene_NBshrt + Uij
       end if

       Fx = fk * Rx
       Fy = fk * Ry
       Fz = fk * Rz

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

   Frc_NBshrt = Frc_NBshrt + Forcer

   Vir_NBshrt(1,1) = Vir_NBshrt(1,1) + Vxx
   Vir_NBshrt(1,2) = Vir_NBshrt(1,2) + Vxy
   Vir_NBshrt(1,3) = Vir_NBshrt(1,3) + Vxz
   Vir_NBshrt(2,2) = Vir_NBshrt(2,2) + Vyy
   Vir_NBshrt(2,3) = Vir_NBshrt(2,3) + Vyz
   Vir_NBshrt(3,3) = Vir_NBshrt(3,3) + Vzz

   Vir_NBshrt(2,1) = Vir_NBshrt(1,2)
   Vir_NBshrt(3,1) = Vir_NBshrt(1,3)
   Vir_NBshrt(3,2) = Vir_NBshrt(2,3)

   deallocate(Forcer,Gk)

end subroutine Force_nonbond_short


!######################################################################
!######################################################################


#ifdef EQUILIBRATION

subroutine PreventContact(Rx,Ry,Rz,R2,Rm)

implicit none

real(8) :: R2, R1, R0, Rm
real(8) :: Rx, Ry, Rz

   R0 = sqrt(R2)
   R2 = Rm       ! R2 = Rmin2 
   R1 = sqrt(Rm) / R0
   Rx = Rx * R1
   Ry = Ry * R1
   Rz = Rz * R1

end subroutine PreventContact

#endif

!######################################################################
!######################################################################


subroutine NONBONDF(Uij,fk,R2,ftype,itype,jtype)

use CGdata, only : CoefAtype, CoefBtype, CoefCtype, Rcut2, Rsw2, Swch, &
&   IDtable, Rmin2
use CommonBlocks, only : QSwitch
use TableFuncs

implicit none

real(8) :: Uij, R2, aij, bij
real(8) :: InvR2, InvR1, InvR3, InvR4, term1, term2, ek, fk
integer :: ftype, itype, jtype
#ifdef GEN
real(8) :: xx, eps, R1, InvR6
real(8) :: Rm, yy, InterPolate
integer :: ii
external InterPolate
#endif

   if(ftype==7) then ! LJ12-4

     aij = CoefAtype(itype,jtype)
     bij = CoefBtype(itype,jtype)
     InvR2 = 1.d0 / R2
     InvR4 = InvR2 * InvR2
     term1 = aij * InvR4 * InvR4 * InvR4
     term2 = bij * InvR4
     ek = term1 - term2
     fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2
#ifdef GEN
     if(QSwitch) then
       call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
#endif
     Uij = ek

   else if(ftype==1) then ! LJ9-6

     aij = CoefAtype(itype,jtype)
     bij = CoefBtype(itype,jtype)
     InvR2 = 1.d0 / R2
     InvR1 = sqrt(InvR2)
     InvR3 = InvR2 * InvR1
     term1 = aij * InvR3 * InvR3 * InvR3
     term2 = bij * InvR3 * InvR3
     ek = term1 - term2
     fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2
#ifdef GEN
     if(QSwitch) then
       call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
#endif
     Uij = ek

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
     if(QSwitch) then
       call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
     Uij = ek

   else if(ftype==2) then ! LJ6-4

     aij = CoefAtype(itype,jtype)
     bij = CoefBtype(itype,jtype)
     InvR2 = 1.d0 / R2
     InvR4 = InvR2 * InvR2
     term1 = aij * InvR4 * InvR2
     term2 = bij * InvR4
     ek = term1 - term2
     fk = ( 6.d0 * term1 - 4.d0 * term2 ) * InvR2
     if(QSwitch) then
       call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
     Uij = ek

   else if(ftype==3) then ! tabulated

     ii  = IDtable(itype,jtype)
     Rm  = Rmin2(itype,jtype)
     xx  = Invgs(ii)
     Uij = InterPolate(R2,Rm,xx,1,ii)
     fk  = InterPolate(R2,Rm,xx,2,ii)

   else if(ftype==4) then ! morse+switching func.

     eps = CoefAtype(itype,jtype)
     aij = CoefBtype(itype,jtype)
     R1  = sqrt(R2)
     term1 = aij/CoefCtype(itype,jtype)
     term2 = exp(- term1*R1 + aij)
     xx = 1.d0 - term2
     ek = eps * ( xx*xx - 1.d0 )
     fk = - 2.d0 * eps * term1 * xx * term2 / R1
     if(QSwitch) then
     call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
     Uij = ek

   else if(ftype==5) then ! LJ8-4

     aij = CoefAtype(itype,jtype)
     bij = CoefBtype(itype,jtype)
     InvR2 = 1.d0 / R2
     InvR4 = InvR2 * InvR2
     term1 = aij * InvR4 * InvR4
     term2 = bij * InvR4
     ek = term1 - term2
     fk = ( 8.d0 * term1 - 4.d0 * term2 ) * InvR2
     if(QSwitch) then
       call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
     Uij = ek

   else if(ftype==6) then ! LJ10-4

     aij = CoefAtype(itype,jtype)
     bij = CoefBtype(itype,jtype)
     InvR2 = 1.d0 / R2
     InvR4 = InvR2 * InvR2
     term1 = aij * InvR4 * InvR4 * InvR2
     term2 = bij * InvR4
     ek = term1 - term2
     fk = ( 10.d0 * term1 - 4.d0 * term2 ) * InvR2
     if(QSwitch) then
       call SwitchFunc(R2,fk,ek,Rsw2(itype,jtype),Rcut2(itype,jtype),Swch(itype,jtype))
     end if
     Uij = ek
#endif

   end if

end subroutine NONBONDF


!######################################################################
!######################################################################


subroutine COULOMBFR(Uij,fk,R2,cf)

use CommonBlocks, only : cCOULOMB
use EwaldParam, only : Alpha, ar2, msh, EFlist
use CGdata, only : Kappa, Ush0, Ush3, Ush4, Fsh2, Fsh3

implicit none

integer :: ii
real(8) :: Uij, R2, cf, R1, x, xtm, fk1, fk2, fk
real(8) :: InvR2, R3, R4
real(8) :: ErrorFunc
real(8) :: yy, yz, zz, dx, y1, y2, y3, z1, z2

   InvR2 = 1.d0 / R2
   R1    = sqrt( R2 )

if(cCOULOMB=='SCREEN') then

   x   = Kappa * R1

   yy  = exp(-x)
   yz  = cf * yy
   zz  = R1 * InvR2
   fk1 = yz * zz
   fk2 = Kappa * yz

   fk  = ( fk1 + fk2 ) * InvR2

else if(cCOULOMB=='SHIFT') then

   x   = InvR2 * R1
   R3  = R2 * R1
   R4  = R2 * R2
   yy  = x + Ush0 + Ush3 * R3 + Ush4 * R4
   zz  = x + Fsh2 * R3 + Fsh3 * R4
   fk1 = cf * yy
   fk  = cf * zz * InvR2

else if(cCOULOMB=='PME'.or.cCOULOMB=='EWALD') then

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

end if

   Uij = fk1

end subroutine COULOMBFR


!######################################################################
!######################################################################


subroutine COULOMBFiso(Uij,fk,R2,cf)

implicit none

real(8) :: Uij, R2, cf, InvR1, fk1, fk, InvR2

   InvR2 = 1.d0 / R2
   InvR1 = sqrt( InvR2 )

   fk1 = cf * InvR1
   fk  = fk1* InvR2

   Uij = fk1

end subroutine COULOMBFiso


!######################################################################
!######################################################################


subroutine ScaledCoordinate(ScR)

use Numbers, only : N
use Configuration, only : R
use CellParam, only : InvH

implicit none

integer :: i
real(8), dimension(3,N) :: ScR
real(8) :: IHxx, IHxy, IHxz
real(8) :: IHyx, IHyy, IHyz
real(8) :: IHzx, IHzy, IHzz
real(8) :: Rx, Ry, Rz

   IHxx = InvH(1,1)
   IHyx = InvH(2,1)
   IHzx = InvH(3,1)
   IHxy = InvH(1,2)
   IHyy = InvH(2,2)
   IHzy = InvH(3,2)
   IHxz = InvH(1,3)
   IHyz = InvH(2,3)
   IHzz = InvH(3,3)

   do i = 1 , N

     Rx = R(1,i)
     Ry = R(2,i)
     Rz = R(3,i)
     ScR(1,i) = IHxx * Rx + IHxy * Ry + IHxz * Rz
     ScR(2,i) = IHyx * Rx + IHyy * Ry + IHyz * Rz
     ScR(3,i) = IHzx * Rx + IHzy * Ry + IHzz * Rz

   end do

end subroutine ScaledCoordinate


!######################################################################
!######################################################################


subroutine Force_nonbond_iso_long

use Numbers, only : N
use CommonBlocks, only : QRigidBody, QCoulomb
use Configuration, only : R
use NoLJparam, only : NumNoLJ
use CommonMPI, only : NProcs, MyRank
use CGdata, only : NBAtomType, Rcut2, Rcut_short2, NBFuncType, &
&   Rbk_short2, Rheal2, Fsw1
use RBparam, only : AtomUnitNum
use BookParam, only : Npair_short, List_shortIJ
use NonbondParam, only : Charge, Frc_NBlong, Ene_ELlong, Ene_NBlong

implicit none
 
integer :: Nas, IRB, i, j, check_bonded, non, itype, jtype
external check_bonded

   Nas = NProcs - MyRank

   Npair_short = 0

   if(QRigidBody) then

     do i = Nas, N, NProcs

       IRB   = AtomUnitNum(i)
       non   = NumNoLJ(i)
       itype = NBAtomType(i)

       do j = i-2 , 1, -2

         if(IRB == AtomUnitNum(j)) cycle

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         call calc_long_interactions

       end do

       do j = i+1 , N, 2

         if(IRB == AtomUnitNum(j)) cycle

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         call calc_long_interactions

       end do

     end do

   else

     do i = Nas, N, NProcs

       non   = NumNoLJ(i)
       itype = NBAtomType(i)

       do j = i-2 , 1, -2

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         call calc_long_interactions

       end do

       do j = i+1 , N, 2

         if(check_bonded(i,j,non) == 1) cycle

         jtype = NBAtomType(j)

         call calc_long_interactions

       end do

     end do

   end if

Contains 


   subroutine calc_long_interactions

   implicit none

   real(8) :: Rx, Ry, Rz, Fx, Fy, Fz, fk, fke
   real(8) :: R2, Uij, Ueij, cf, ScF, yr

      Rx = R(1,i) - R(1,j)
      Ry = R(2,i) - R(2,j)
      Rz = R(3,i) - R(3,j)
      R2  = Rx*Rx + Ry*Ry + Rz*Rz

      if((R2 <= Rcut2(itype,jtype)).and.(R2 > Rcut_short2)) then

        call NONBONDF(Uij,fk,R2,NBFuncType(itype,jtype),itype,jtype)

        if(QCoulomb) then
          cf = Charge(i) * Charge(j)
          if(cf/=0.) then
            call COULOMBFiso(Ueij,fke,R2,cf)
            fk = fk + fke
            Ene_ELlong = Ene_ELlong + Ueij
          end if
        end if

        Ene_NBlong = Ene_NBlong + Uij

        if(R2 <= Rheal2) then
          yr  = (R2 - Rcut_short2) * Fsw1
          ScF = - yr * yr * (2.d0 * yr - 3.d0)
          fk = fk * ScF
        end if

        Fx = fk * Rx
        Fy = fk * Ry
        Fz = fk * Rz

        Frc_NBlong(1,i) = Frc_NBlong(1,i) + Fx
        Frc_NBlong(2,i) = Frc_NBlong(2,i) + Fy
        Frc_NBlong(3,i) = Frc_NBlong(3,i) + Fz
        Frc_NBlong(1,j) = Frc_NBlong(1,j) - Fx
        Frc_NBlong(2,j) = Frc_NBlong(2,j) - Fy
        Frc_NBlong(3,j) = Frc_NBlong(3,j) - Fz

      end if

     if(R2 <= Rbk_short2) then
       Npair_short = Npair_short + 1
       List_shortIJ(1,Npair_short) = i
       List_shortIJ(2,Npair_short) = j
     end if

   end subroutine calc_long_interactions


end subroutine Force_nonbond_iso_long


!######################################################################
!######################################################################


subroutine Force_nonbond_iso_short

use CommonBlocks, only : QCoulomb
use Configuration, only : R
use CGdata, only : NBAtomType, NBFuncType, Rheal2, Rcut_short2, Fsw1
use BookParam, only : Npair_short, List_shortIJ
use NonbondParam, only : Charge, Frc_NBshrt, Ene_ELshrt, Ene_NBshrt

implicit none

integer :: i, j, l
real(8) :: Rx, Ry, Rz, Fx, Fy, Fz, fk, fke
real(8) :: R2, Uij, Ueij, cf, ScF, yr, y3
integer :: itype, jtype

   do l = 1, Npair_short
     i = List_shortIJ(1,l)
     j = List_shortIJ(2,l)

     Rx = R(1,i) - R(1,j)
     Ry = R(2,i) - R(2,j)
     Rz = R(3,i) - R(3,j)
     R2  = Rx*Rx + Ry*Ry + Rz*Rz

     if(R2 <= Rheal2) then

       itype = NBAtomType(i)
       jtype = NBAtomType(j)

       call NONBONDF(Uij,fk,R2,NBFuncType(itype,jtype),itype,jtype)

       if(QCoulomb) then
         cf = Charge(i) * Charge(j)
         if(cf/=0.) then
           call COULOMBFiso(Ueij,fke,R2,cf)
           fk = fk + fke
           if(R2 <= Rcut_short2) Ene_ELshrt = Ene_ELshrt + Ueij
         end if
       end if

       if(R2 > Rcut_short2) then
         yr  = (R2 - Rcut_short2) * Fsw1
         y3  = yr * yr * (2.d0 * yr - 3.d0)
         ScF = 1.d0 + y3
         fk  = fk * ScF
       else
         Ene_NBshrt = Ene_NBshrt + Uij
       end if

       Fx = fk * Rx
       Fy = fk * Ry
       Fz = fk * Rz

       Frc_NBshrt(1,i) = Frc_NBshrt(1,i) + Fx
       Frc_NBshrt(2,i) = Frc_NBshrt(2,i) + Fy
       Frc_NBshrt(3,i) = Frc_NBshrt(3,i) + Fz
       Frc_NBshrt(1,j) = Frc_NBshrt(1,j) - Fx
       Frc_NBshrt(2,j) = Frc_NBshrt(2,j) - Fy
       Frc_NBshrt(3,j) = Frc_NBshrt(3,j) - Fz

     end if

   end do

end subroutine Force_nonbond_iso_short
