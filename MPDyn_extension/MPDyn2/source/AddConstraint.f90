! ############################
! ## SUBROUTINE LIST 
! ## -- ConstPrepare 
! ## -- AddConstR 
! ## -- AddConstV 
! ## -- ConstSingleR 
! ## -- ConstPairR 
! ## -- ConstSingleV 
! ## -- ConstPairV 
! ############################


!######################################################################
!######################################################################


subroutine ConstPrepare

use Configuration
use CommonBlocks, only : Job_name, QRigidBody, QPBC
use RBparam, only : MassRB, AtomUnitNum, V_RB, R_RB
use SMDparam, only : NumOpFixS, NumFixF, NumFixL, FixDir, MassRG, NumOpFixP, &
&   NumFixFi, NumFixLi, NumFixFj, NumFixLj, MassRGi, MassRGj, Rcomij
use CellParam, only : H, InvH
use AtomParam , only : Mass
use IOparam, only : DirectoryName

implicit none

integer :: i, j, k
real(8) :: VcomZ
real(8), dimension(3) :: Rcomi, Rcomj, Rij, Sij

   if(NumOpFixS/=0) then

     do i = 1, NumOpFixS

       MassRG(i) = 0.d0

       if(QRigidBody) then
         j = NumFixF(i)
         NumFixF(i) = AtomUnitNum(j)
         j = NumFixL(i)
         NumFixL(i) = AtomUnitNum(j)

         do j = NumFixF(i), NumFixL(i)
           MassRG(i) = MassRG(i) + MassRB(j)
         end do
       else
         do j = NumFixF(i), NumFixL(i)
           MassRG(i) = MassRG(i) + Mass(j)
         end do
       end if

     end do

     do i = 1, NumOpFixS

       VcomZ = 0.d0
       k = FixDir(i)

       if(QRigidBody) then
         do j = NumFixF(i), NumFixL(i)
           VcomZ = VcomZ + MassRB(j) * V_RB(k,j)
         end do
       else
         do j = NumFixF(i), NumFixL(i)
           VcomZ = VcomZ + Mass(j) * Vel(k,j)
         end do
       end if

       VcomZ = VcomZ / MassRG(i)

       if(QRigidBody) then
         do j = NumFixF(i), NumFixL(i)
           V_RB(k,j) = V_RB(k,j) - VcomZ
         end do
       else
         do j = NumFixF(i), NumFixL(i)
           Vel(k,j) = Vel(k,j) - VcomZ
         end do
       end if

     end do

   end if

   if(NumOpFixP/=0) then

     open(26,file=trim(DirectoryName)//trim(adjustl(Job_name))//'F_PAIR.dat',status='unknown')

     do i = 1, NumOpFixP

       MassRGi(i) = 0.d0
       MassRGj(i) = 0.d0

       if(QRigidBody) then
         j = NumFixFi(i)
         k = NumFixLi(i)
         NumFixFi(i) = AtomUnitNum(j)
         NumFixLi(i) = AtomUnitNum(k)
         j = NumFixFj(i)
         k = NumFixLj(i)
         NumFixFj(i) = AtomUnitNum(j)
         NumFixLj(i) = AtomUnitNum(k)

         Rcomi(:) = 0.d0
         do j = NumFixFi(i), NumFixLi(i)
           MassRGi(i) = MassRGi(i) + MassRB(j)
           Rcomi(:) = Rcomi(:) + MassRB(j) * R_RB(:,j)
         end do
         Rcomi(:) = Rcomi(:) / MassRGi(i)

         Rcomj(:) = 0.d0
         do j = NumFixFj(i), NumFixLj(i)
           MassRGj(i) = MassRGj(i) + MassRB(j)
           Rcomj(:) = Rcomj(:) + MassRB(j) * R_RB(:,j)
         end do
         Rcomj(:) = Rcomj(:) / MassRGj(i)

       else

         Rcomi(:) = 0.d0
         do j = NumFixFi(i), NumFixLi(i)
           MassRGi(i) = MassRGi(i) + Mass(j)
           Rcomi(:) = Rcomi(:) + Mass(j) * R(:,j)
         end do
         Rcomi(:) = Rcomi(:) / MassRGi(i)

         Rcomj(:) = 0.d0
         do j = NumFixFj(i), NumFixLj(i)
           MassRGj(i) = MassRGj(i) + Mass(j)
           Rcomj(:) = Rcomj(:) + Mass(j) * R(:,j)
         end do
         Rcomj(:) = Rcomj(:) / MassRGj(i)

       end if

       Rij = Rcomi - Rcomj

       if(QPBC) then
         Sij = matmul( InvH, Rij )
         Sij = Sij - nint( Sij )
         Rij = matmul( H, Sij )
       end if

       Rcomij(:,i) = Rij(:)

     end do

   end if

end subroutine ConstPrepare


!######################################################################
!######################################################################


subroutine AddConstR

use SMDparam, only : NumOpFixS, NumOpFixP

   if(NumOpFixS/=0) then
     call ConstSingleR
   end if
   if(NumOpFixP/=0) then
     call ConstPairR
   end if

end subroutine AddConstR


!######################################################################
!######################################################################


subroutine AddConstV(i)

use SMDparam, only : NumOpFixS, NumOpFixP

integer :: i

   if(NumOpFixS/=0) then
     call ConstSingleV(i)
   end if
   if(NumOpFixP/=0) then
     call ConstPairV
   end if

end subroutine AddConstV


!######################################################################
!######################################################################


subroutine ConstSingleR

use Configuration
use CommonBlocks, only : QRigidBody
use RBparam, only : MassRB, R_RB
use SMDparam, only : NumOpFixS, NumFixF, NumFixL, FixDir, RGConst, MassRG
use AtomParam , only : Mass

implicit none

integer :: i, j, k
real(8) :: Rcom

   do i = 1, NumOpFixS

     Rcom = 0.d0
     k = FixDir(i)

     if(QRigidBody) then
       do j = NumFixF(i), NumFixL(i)
         Rcom = Rcom + MassRB(i) * R_RB(k,j)
       end do
     else
       do j = NumFixF(i), NumFixL(i)
         Rcom = Rcom + Mass(i) * R(k,j)
       end do
     end if

     Rcom = Rcom / MassRG(i)
     Rcom = Rcom - RGConst(i)

     if(QRigidBody) then
       do j = NumFixF(i), NumFixL(i)
         R_RB(k,j) = R_RB(k,j) - Rcom
       end do
     else
       do j = NumFixF(i), NumFixL(i)
         R(k,j) = R(k,j) - Rcom
       end do
     end if

   end do

end subroutine ConstSingleR


!######################################################################
!######################################################################


subroutine ConstPairR

use Configuration
use CommonBlocks, only : QPBC, QRigidBody
use RBparam, only : MassRB, R_RB, V_RB
use SMDparam, only : NumOpFixP, NumFixFi, NumFixLi, NumFixFj, NumFixLj, &
&   MassRGi, MassRGj, ConstDis2, Rcomij, Vir_ConstR
use CellParam, only : H, InvH
use AtomParam, only : Mass
use TimeParam, only : deltat, dt2

implicit none

integer :: i, j
real(8), dimension(3,NumOpFixP) :: Rcomi, Rcomj
real(8), dimension(3) :: Rij, Sij, drs, dr
real(8), dimension(3,3) :: VircA
real(8), dimension(3) :: dRi, dRj
real(8) :: R2, d2, dq, drf2, gab, InvMi, InvMj, InvDt
real(8), parameter :: RpTol=1.0d-28

   Rcomi = 0.d0
   Rcomj = 0.d0

   InvDt = 1.d0 / deltat
   if(QPBC) VircA = 0.d0

   do i = 1, NumOpFixP

     InvMi = 1.d0/MassRGi(i)
     InvMj = 1.d0/MassRGj(i)

     if(QRigidBody) then

     do j = 1, NumFixFi(i), NumFixLi(i)
       Rcomi(:,i) = Rcomi(:,i) + MassRB(j) * R_RB(:,j)
     end do
     Rcomi(:,i) = Rcomi(:,i) * InvMi

     do j = 1, NumFixFj(i), NumFixLj(i)
       Rcomj(:,i) = Rcomj(:,i) + MassRB(j) * R_RB(:,j)
     end do
     Rcomj(:,i) = Rcomj(:,i) * InvMj

     else

     do j = 1, NumFixFi(i), NumFixLi(i)
       Rcomi(:,i) = Rcomi(:,i) + Mass(j) * R(:,j)
     end do
     Rcomi(:,i) = Rcomi(:,i) * InvMi

     do j = 1, NumFixFj(i), NumFixLj(i)
       Rcomj(:,i) = Rcomj(:,i) + Mass(j) * R(:,j)
     end do
     Rcomj(:,i) = Rcomj(:,i) * InvMj

     end if

     Rij = Rcomi(:,i) - Rcomj(:,i)
     if(QPBC) then
       Sij = matmul( InvH, Rij )
       Sij = Sij - nint(Sij)
       Rij = matmul( H, Sij )
     end if
     R2 = dot_product( Rij, Rij )
     d2 = ConstDis2(i)
     dq = d2 - R2

     dr(:) = Rcomij(:,i)
     drf2 = dot_product(dr,Rij)

     if( drf2 < (d2*RpTol) ) then
       write(*,'(a)') ' constraint failure '
       write(*,*) ' in subroutine "ComstPairR" ', i
       call Finalize
     end if

     gab = dq / (2.d0 * (InvMi + InvMj) * drf2)
     drs = dr * gab

     if(QPBC) then
       VircA(1,1) = VircA(1,1) + gab * dr(1) * dr(1)
       VircA(1,2) = VircA(1,2) + gab * dr(1) * dr(2)
       VircA(1,3) = VircA(1,3) + gab * dr(1) * dr(3)
       VircA(2,2) = VircA(2,2) + gab * dr(2) * dr(2)
       VircA(2,3) = VircA(2,3) + gab * dr(2) * dr(3)
       VircA(3,3) = VircA(3,3) + gab * dr(3) * dr(3)
     end if

     dRi =   InvMi * drs
     dRj = - InvMj * drs

     Rcomi(:,i) = Rcomi(:,i) + dRi
     Rcomj(:,i) = Rcomj(:,i) + dRj

     if(QRigidBody) then

     do j = 1, NumFixFi(i), NumFixLi(i)
       R_RB(:,j) = R_RB(:,j) + dRi
     end do

     do j = 1, NumFixFj(i), NumFixLj(i)
       R_RB(:,j) = R_RB(:,j) + dRj
     end do

     else

     do j = 1, NumFixFi(i), NumFixLi(i)
       R(:,j) = R(:,j) + dRi
     end do

     do j = 1, NumFixFj(i), NumFixLj(i)
       R(:,j) = R(:,j) + dRj
     end do

     end if

     dRi = dRi * InvDt
     dRj = dRj * InvDt

     if(QRigidBody) then

     do j = 1, NumFixFi(i), NumFixLi(i)
       V_RB(:,j) = V_RB(:,j) + dRi
     end do

     do j = 1, NumFixFj(i), NumFixLj(i)
       V_RB(:,j) = V_RB(:,j) + dRj
     end do

     else

     do j = 1, NumFixFi(i), NumFixLi(i)
       Vel(:,j) = Vel(:,j) + dRi
     end do

     do j = 1, NumFixFj(i), NumFixLj(i)
       Vel(:,j) = Vel(:,j) + dRj
     end do

     end if

     Rij = Rcomi(:,i) - Rcomj(:,i)
     if(QPBC) then
       Sij = matmul( InvH, Rij )
       Sij = Sij - nint(Sij)
       Rij = matmul( H, Sij )
     end if
     Rcomij(:,i) = Rij(:)

   end do

   if(QPBC) then

     VircA(2,1) = VircA(1,2)
     VircA(3,1) = VircA(1,3)
     VircA(3,2) = VircA(2,3)

     Vir_ConstR = VircA / dt2

   end if

end subroutine ConstPairR


!######################################################################
!######################################################################


subroutine ConstSingleV(step)

use Configuration
use CommonBlocks, only : Job_name, QMaster, QRigidBody
use RBparam, only : NumRBAtom, RBType, QSingle, InitAtomRB, MassRB, &
&   V_RB, R_RB
use EAM_param, only : Frc_EAM
use SMDparam, only : NumFreqConst, NumOpFixS, NumFixF, NumFixL, FixDir, &
& MassRG, FFsm
use UnitExParam, only : cvol
use EwaldParam, only : Frc_Eksp
use OptConstraintParam, only : Frc_OptC
use NonbondParam, only : Frc_Ersp, Frc_NBlong, Frc_NBshrt
use BondedParam, only : Frc_Bond, Frc_Angle, Frc_UB, Frc_Dihed, Frc_Impro
use AtomParam, only : Mass
use TimeParam, only : Timeps
use IOparam, only : DirectoryName

implicit none

integer :: i, j, k, step, ii, MyType, Nc, jj, kk
real(8) :: Vcom
real(8), dimension(NumOpFixS) :: Fconst, Fav
real(8), dimension(3,NumOpFixS) :: Rcom

   if(QMaster) then

     do i = 1, NumOpFixS

       Vcom = 0.d0
       k = FixDir(i)

       if(QRigidBody) then
         do j = NumFixF(i), NumFixL(i)
           Vcom = Vcom + MassRB(i) * V_RB(k,j)
         end do
       else
         do j = NumFixF(i), NumFixL(i)
           Vcom = Vcom + Mass(i) * Vel(k,j)
         end do
       end if

       Vcom = Vcom / MassRG(i)

       if(QRigidBody) then
         do j = NumFixF(i), NumFixL(i)
           V_RB(k,j) = V_RB(k,j) - Vcom
         end do
       else
         do j = NumFixF(i), NumFixL(i)
           Vel(k,j) = Vel(k,j) - Vcom
         end do
       end if

     end do

     if(step == NumFreqConst) then
       open(26,file=trim(DirectoryName)//trim(adjustl(Job_name))//'F_SNGL.dat',status='unknown')
       FFsm = 0.d0
     end if

   end if

   Fconst = 0.d0
   do i = 1, NumOpFixS

     k = FixDir(i)

     if(QRigidBody) then

       ii = InitAtomRB(NumFixF(i))

       do j = NumFixF(i), NumFixL(i)
         if(QSingle(j)) then
           ii = ii + 1
           Fconst(i) = Fconst(i) + Frc_Bond  (k,ii) &
           &                     + Frc_Angle (k,ii) &
           &                     + Frc_UB    (k,ii) &
           &                     + Frc_Dihed (k,ii) &
           &                     + Frc_Impro (k,ii) &
           &                     + Frc_OptC  (k,ii) &
           &                     + Frc_EAM   (k,ii) &
           &                     + Frc_Ersp  (k,ii) &
           &                     + Frc_NBshrt(k,ii) &
           &                     + Frc_NBlong(k,ii) &
           &                     + Frc_Eksp  (k,ii)
         else
           MyType = RBType(j)
           Nc = NumRBATom(MyType)
           do jj = 1, Nc
             kk = ii + jj
             Fconst(i) = Fconst(i) + Frc_Bond  (k,kk) &
             &                     + Frc_Angle (k,kk) &
             &                     + Frc_UB    (k,kk) &
             &                     + Frc_Dihed (k,kk) &
             &                     + Frc_Impro (k,kk) &
             &                     + Frc_OptC  (k,kk) &
             &                     + Frc_EAM   (k,kk) &
             &                     + Frc_Ersp  (k,kk) &
             &                     + Frc_NBshrt(k,kk) &
             &                     + Frc_NBlong(k,kk) &
             &                     + Frc_Eksp  (k,kk)
           end do
           ii = ii + Nc
         end if
       end do

     else

       do j = NumFixF(i), NumFixL(i)
         Fconst(i) = Fconst(i) + Frc_Bond  (k,j) &
         &                     + Frc_Angle (k,j) &
         &                     + Frc_UB    (k,j) &
         &                     + Frc_Dihed (k,j) &
         &                     + Frc_Impro (k,j) &
         &                     + Frc_OptC  (k,j) &
         &                     + Frc_EAM   (k,j) &
         &                     + Frc_Ersp  (k,j) &
         &                     + Frc_NBshrt(k,j) &
         &                     + Frc_NBlong(k,j) &
         &                     + Frc_Eksp  (k,j)
       end do

     end if

   end do

   Fconst(:) = Fconst(:) * cvol

   FFsm(:) = FFsm(:) + Fconst(:)

   if(mod(step,NumFreqConst)==0) then

     call SumConstF(NumOpFixS,FFsm,Fconst)

     if(QMaster) then

       do i = 1, NumOpFixS

         Rcom(:,i) = 0.d0

         if(QRigidBody) then
           do j = NumFixF(i), NumFixL(i)
             Rcom(:,i) = Rcom(:,i) + MassRB(j) * R_RB(:,j)
           end do
         else
           do j = NumFixF(i), NumFixL(i)
             Rcom(:,i) = Rcom(:,i) + Mass(j) * R(:,j)
           end do
         end if

         Rcom(:,i) = Rcom(:,i) / MassRG(i)

       end do

       Fav(:)  = FFsm(:) / dble(step)

       do i = 1, NumOpFixS
         if(i==1) then
           write(26,'(f12.4,5e16.8)') Timeps, Fconst(i), Fav(i), Rcom(:,i)
         else
           write(26,'(12x,5e16.8)') Fconst(i), Fav(i), Rcom(:,i)
         end if
       end do

     end if

   end if

end subroutine ConstSingleV


!######################################################################
!######################################################################


subroutine ConstPairV

use Configuration
use CommonBlocks, only : QPBC, QRigidBody
use RBparam, only : MassRB, V_RB
use SMDparam, only : NumOpFixP, NumFixFi, NumFixLi, NumFixFj, NumFixLj, &
& ConstDis2, MassRGi, MassRGj, Rcomij
use SMDparam, only : Vir_ConstV
use AtomParam , only : Mass
use TimeParam, only : dt2

implicit none

integer :: i, j
real(8), dimension(3,NumOpFixP) :: Vcomi, Vcomj
real(8), dimension(3) :: Rij, Vij, drs
real(8), dimension(3) :: dRi, dRj
real(8) :: gab, RV, InvMi, InvMj
real(8), dimension(3,3) :: VircB

   if(QPBC) VircB = 0.d0

   do i = 1, NumOpFixP

     InvMi = 1.d0/MassRGi(i)
     InvMj = 1.d0/MassRGj(i)

     if(QRigidBody) then

       Vcomi(:,i) = 0.d0
       do j = 1, NumFixFi(i), NumFixLi(i)
         Vcomi(:,i) = Vcomi(:,i) + MassRB(j) * V_RB(:,j)
       end do
       Vcomi(:,i) = Vcomi(:,i) * InvMi

       Vcomj(:,i) = 0.d0
       do j = 1, NumFixFj(i), NumFixLj(i)
         Vcomj(:,i) = Vcomj(:,i) + MassRB(j) * V_RB(:,j)
       end do
       Vcomj(:,i) = Vcomj(:,i) * InvMj

     else

       Vcomi(:,i) = 0.d0
       do j = 1, NumFixFi(i), NumFixLi(i)
         Vcomi(:,i) = Vcomi(:,i) + Mass(j) * Vel(:,j)
       end do
       Vcomi(:,i) = Vcomi(:,i) * InvMi

       Vcomj(:,i) = 0.d0
       do j = 1, NumFixFj(i), NumFixLj(i)
         Vcomj(:,i) = Vcomj(:,i) + Mass(j) * Vel(:,j)
       end do
       Vcomj(:,i) = Vcomj(:,i) * InvMj

     end if

     Rij(:) = Rcomij(:,i)
     Vij(:) = Vcomi(:,i) - Vcomj(:,i)
     RV = dot_product(Rij,Vij)
     gab = -RV / ((InvMi + InvMj) * ConstDis2(i))
     drs = Rij * gab

     if(QPBC) then
       VircB(1,1) = VircB(1,1) + gab * Rij(1) * Rij(1)
       VircB(1,2) = VircB(1,2) + gab * Rij(1) * Rij(2)
       VircB(1,3) = VircB(1,3) + gab * Rij(1) * Rij(3)
       VircB(2,2) = VircB(2,2) + gab * Rij(2) * Rij(2)
       VircB(2,3) = VircB(2,3) + gab * Rij(2) * Rij(3)
       VircB(3,3) = VircB(3,3) + gab * Rij(3) * Rij(3)
     end if

     dRi =   InvMi * drs
     dRj = - InvMj * drs

     if(QRigidBody) then
       do j = 1, NumFixFi(i), NumFixLi(i)
         V_RB(:,j) = V_RB(:,j) + dRi
       end do
       do j = 1, NumFixFj(i), NumFixLj(i)
         V_RB(:,j) = V_RB(:,j) + dRj
       end do
     else
       do j = 1, NumFixFi(i), NumFixLi(i)
         Vel(:,j) = Vel(:,j) + dRi
       end do
       do j = 1, NumFixFj(i), NumFixLj(i)
         Vel(:,j) = Vel(:,j) + dRj
       end do
     end if

   end do

   if(QPBC) then
     VircB(2,1) = VircB(1,2)
     VircB(3,1) = VircB(1,3)
     VircB(3,2) = VircB(2,3)
     Vir_ConstV = VircB / dt2
   end if

end subroutine ConstPairV
