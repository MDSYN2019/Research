! ############################
! ## SUBROUTINE LIST 
! ## -- SMD_reference
! ## -- SMD_pre
! ## -- SMD_sample
! ## -- Force_Jarzynski
! ############################


!######################################################################
!######################################################################


subroutine SMD_reference

use SMDparam, only : R0_steering, Vshift_steering

implicit none

   R0_steering = R0_steering + Vshift_steering

end subroutine SMD_reference


!######################################################################
!######################################################################


subroutine SMD_pre

use CommonBlocks, only : Job_name, QRigidBody, QOpJarz, QClust, QFixCOM
use AtomParam, only : Mass
use SMDparam, only : NumOpFixS, NumOpFixP, MassRG, MassRGi, MassRGj, &
&   NumFixF, NumFixL, NumFixFi, NumFixLi, NumFixFj, NumFixLj,        &
&   Nsample_steering, Workforce, FFsm, Vshift_steering, NSelCone
use RBparam, only : MassRB, AtomUnitNum
use IOparam, only : DirectoryName
use OptConstraintParam, only : fscom

implicit none

integer :: i, j, k
character(len=72) :: File1

   Nsample_steering = 0
   Workforce = 0.d0

   if(QOpJarz.or.(QClust.and.(Vshift_steering==0.))) then
     write(File1,'(a,a)') trim(adjustl(Job_name)),'_ConstraintF.dat'
     open(26,file=trim(DirectoryName)//trim(File1),form='formatted',status='unknown')
     write(26,'(a)') '# time[ps]  Force[kcal/mol/A] AveF[kcal/mol/A]  Rcom [A]'
     FFsm = 0.d0
   else
     write(File1,'(a,a)') trim(adjustl(Job_name)),'_SMD_WORK.dat'
     open(26,file=trim(DirectoryName)//trim(File1),form='formatted',status='unknown')
     write(26,'(a)') '# time[ps]  work[kcal/mol] R0(target)[A] R0(actual)[A]'
   end if

   if(NSelCone/=0) then
     write(File1,'(a,a)') trim(adjustl(Job_name)),'_FixedCOM_PMF.dat'
     open(91,file=trim(DirectoryName)//trim(File1),form='formatted',status='unknown')
     write(91,'(a)') '# time[ps]  AveF(x,y,z)[kcal/mol/A]  Rcom(x,y,z) [A]'
     fscom(:,:) = 0.d0
   end if

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

   end if

   if(NumOpFixP/=0) then

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

         do j = NumFixFi(i), NumFixLi(i)
           MassRGi(i) = MassRGi(i) + MassRB(j)
         end do

         do j = NumFixFj(i), NumFixLj(i)
           MassRGj(i) = MassRGj(i) + MassRB(j)
         end do

       else

         do j = NumFixFi(i), NumFixLi(i)
           MassRGi(i) = MassRGi(i) + Mass(j)
         end do

         do j = NumFixFj(i), NumFixLj(i)
           MassRGj(i) = MassRGj(i) + Mass(j)
         end do

       end if

     end do

   end if

end subroutine SMD_pre


!######################################################################
!######################################################################


subroutine SMD_sample(step)

use CommonBlocks, only : QRigidBody, QPBC, QMaster, QOpJarz, QClust, QFixCOM, &
&   QCorrCone
use TimeParam, only : deltat, lk, Timeps
use AtomParam, only : Mass
use CellParam, only : H, InvH, CellShft
use OptConstraintParam, only : NumClust, NAtomClus, Ipartner, Nccom, fcom, fscom, &
&   Rcon, ShiftCone
use SMDparam, only : NumOpFixS, NumOpFixP, NumFixF, NumFixL, NumFixFi, &
& NumFixLi, NumFixFj, NumFixLj, MassRG, MassRGi, MassRGj, NumFreqConst, FFsm, &
& R0_initial, Nsample_steering, Workforce, k_steering,R0_steering, Uvec_steering, &
& Vshift_steering, NSelCyl, NSelCylI, NSelCone, FFcyl, FFcone
use Configuration, only : R
use RBparam, only : MassRB, R_RB
use UnitExParam, only : cvol

implicit none

integer :: i, j, l, step
real(8), dimension(3) :: Rcomi, Rcomj
real(8), dimension(3) :: Rij, Sij, ScRi, ScRj
real(8) :: R2, Rpos, ff, Work
real(8) :: Fconst, Fav
real(8) :: clhx, clhy, clhz
real(8) :: Rix, Riy, Riz, Rmin2
real(8) :: Rx, Ry, Rz
integer :: Nx, Ny, Nz
real(8) :: Favx, Favy, Favz, xx

   if(.not.QMaster) return

   if(QClust) then

     if(QPBC) then
       clhx = H(1,1)*0.5d0
       clhy = H(2,2)*0.5d0
       clhz = H(3,3)*0.5d0
     end if

     Rix = R(1,Ipartner)
     Riy = R(2,Ipartner)
     Riz = R(3,Ipartner)

     Rmin2 = H(1,1)*H(1,1)

     do i = 1, NumClust

       j = NAtomClus(i)

       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
       if(QPBC) then
         if(Rx>clhx) then
           Nx = -9
         else if(Rx<-clhx) then
           Nx =  9
         else
           Nx =  0
         end if
         if(Ry>clhy) then
           Ny = -3
         else if(Ry<-clhy) then
           Ny =  3
         else
           Ny =  0
         end if
         if(Rz>clhz) then
           Nz = -1
         else if(Rz<-clhz) then
           Nz =  1
         else
           Nz =  0
         end if
         l = Nx + Ny + Nz
         Rx = Rx + CellShft(1,l)
         Ry = Ry + CellShft(2,l)
         Rz = Rz + CellShft(3,l)
       end if
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if(R2<Rmin2) then
         Rmin2 = R2
       end if

     end do

     Rpos = sqrt(Rmin2)
     ff = - k_steering * (Rpos - R0_steering)

   else if(NumOpFixS==1) then

     Rcomi(:) = 0.d0
     if(QRigidBody) then
       do j = NumFixF(1), NumFixL(1)
         Rcomi(:) = Rcomi(:) + MassRB(j)*R_RB(:,j)
       end do
     else
       do j = NumFixF(1), NumFixL(1)
         Rcomi(:) = Rcomi(:) + Mass(j)*R(:,j)
       end do
     end if
     Rcomi(:) = Rcomi(:) / MassRG(1)

     Rpos = dot_product(Rcomi(:),Uvec_steering(:))

     ff = - k_steering * (Rpos - R0_steering)

   else if(NumOpFixP==1) then

     Rcomi(:) = 0.d0
     Rcomj(:) = 0.d0
     if(QRigidBody) then
       do j = NumFixFi(1), NumFixLi(1)
         Rcomi(:) = Rcomi(:) + MassRB(j)*R_RB(:,j)
       end do
       do j = NumFixFj(1), NumFixLj(1)
         Rcomj(:) = Rcomj(:) + MassRB(j)*R_RB(:,j)
       end do
     else
       do j = NumFixFi(1), NumFixLi(1)
         Rcomi(:) = Rcomi(:) + Mass(j)*R(:,j)
       end do
       do j = NumFixFj(1), NumFixLj(1)
         Rcomj(:) = Rcomj(:) + Mass(j)*R(:,j)
       end do
     end if
     Rcomi(:) = Rcomi(:) / MassRGi(1)
     Rcomj(:) = Rcomj(:) / MassRGj(1)

     if(QPBC) then
       ScRi(:) = matmul(InvH,Rcomi(:))
       ScRj(:) = matmul(InvH,Rcomj(:))
       Sij(:) = ScRi(:) - ScRj(:)
       Sij(:) = Sij(:) - nint(Sij(:))
       Rij(:) = matmul(H,Sij(:))
     else
       Rij(:) = Rcomi(:) - Rcomj(:)
     end if
     R2 = dot_product(Rij,Rij)
     Rpos = sqrt(R2)

     ff = - k_steering * (Rpos - R0_steering)

   else if((NSelCyl>0).or.(NSelCylI>0)) then

     ff = FFcyl

   else if(NSelCone>0) then

     ff = FFcone

   else

     write(*,*) 'error : no constraint for Jarzynski (steered MD)'
     call Finalize

   end if

   if(QOpJarz.or.(QClust.and.(Vshift_steering==0.))) then ! ## Constraint

     Fconst = ff * cvol
     FFsm(1) = FFsm(1) + Fconst

     Fav = FFsm(1) / dble(step/lk)

     if(mod(step,NumFreqConst)==0) then
       if(NSelCyl>0.or.NSelCone>0.or.NSelCylI>0) then
         if(QCorrCone) then
           write(26,'(f12.4,2e16.8,f10.4)') Timeps, Fconst, Fav, -ShiftCone
         else
           write(26,'(f12.4,2e16.8)') Timeps, Fconst, Fav
         end if
       else
         write(26,'(f12.4,3e16.8)') Timeps, Fconst, Fav, Rpos
       end if
     end if

   else ! ## normal Jarzynski

     Nsample_steering = Nsample_steering + 1

     Workforce = Workforce + ff

     Work = Workforce / dble(Nsample_steering) * (R0_steering - R0_initial)

     if(mod(step,NumFreqConst)==0) then
       write(26,'(f9.4,e16.8,2f12.5)') deltat*lk*Nsample_steering, Work*cvol, R0_steering, Rpos
     end if

   end if

   if(NSelCone>0) then

     xx = cvol * dble(lk) / dble(step)

     fscom(:,1) = fscom(:,1) + fcom(:,1)
     Favx = fscom(1,1) * xx
     Favy = fscom(2,1) * xx
     Favz = fscom(3,1) * xx
     if(mod(step,NumFreqConst)==0) then
       write(91,'(f9.4,3e16.8,3f10.4)') deltat*lk*Nsample_steering, Favx, Favy, Favz, &
       &                               Rcon(:,1)
     end if

   end if

end subroutine SMD_sample


!######################################################################
!######################################################################


subroutine Force_Jarzynski

use CommonBlocks, only : QRigidBody, QPBC, QMaster, QClust, QCorrCone
use AtomParam, only : Mass
use CellParam, only : H, InvH, CellShft
use OptConstraintParam, only : Frc_OptC, Vir_OptC, Ene_OptC, NumClust, &
&   NAtomClus, Ipartner, Thick_Lipo, Vol_Lipo, facta, factb, factc,    &
&   Invfacta, kccom, Xcon, Rcon, fcom, InvMscom, ShiftCone
use SMDparam, only : NumOpFixS, NumOpFixP, NumFixF, NumFixL, NumFixFi, &
&   NumFixLi, NumFixFj, NumFixLj, MassRG, MassRGi, MassRGj, &
&   k_steering, R0_steering, Uvec_steering, Nsample_steering, &
&   NSelCyl, NSelCylI, NSelCone, FFcyl, FFcone, icaxis, icx1, icx2, IdSelC
use Configuration, only : R
use RBparam, only : MassRB, R_RB, InitAtomRB, RBType, NumRBAtom
use UnitExParam, only : pi

implicit none

integer :: i, j, k, i1, MyType, l, ii
real(8), dimension(3) :: Rcomi, Rcomj
real(8), dimension(3) :: Fi, Fj, ScRi, ScRj
real(8), dimension(3) :: Rij, Sij
real(8) :: InvMGi, InvMGj
real(8) :: Rpos, dR, ff, R2, R1, fc
real(8) :: clhx, clhy, clhz
real(8) :: Rx, Ry, Rz, Rix, Riy, Riz
real(8) :: Fx, Fy, Fz, Rmx, Rmy, Rmz, Rmin2
integer :: Nx, Ny, Nz
real(8) :: dR2, zz, theta, cst, snt, cst1, csti
real(8) :: fa, fb, fd, rg, pref1, x3, x4, x5, x6, x7

   if(.not.QMaster) return

   if(QClust) then

     clhx = H(1,1)*0.5d0
     clhy = H(2,2)*0.5d0
     clhz = H(3,3)*0.5d0

     Rix = R(1,Ipartner)
     Riy = R(2,Ipartner)
     Riz = R(3,Ipartner)

     Rmin2 = H(1,1)*H(1,1)

     do i = 1, NumClust

       j = NAtomClus(i)

       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
       if(Rx>clhx) then
         Nx = -9
       else if(Rx<-clhx) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Ry>clhy) then
         Ny = -3
       else if(Ry<-clhy) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Rz>clhz) then
         Nz = -1
       else if(Rz<-clhz) then
         Nz =  1
       else
         Nz =  0
       end if
       l = Nx + Ny + Nz
       Rx = Rx + CellShft(1,l)
       Ry = Ry + CellShft(2,l)
       Rz = Rz + CellShft(3,l)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if(R2<Rmin2) then
         Rmin2 = R2
         Rmx = Rx
         Rmy = Ry
         Rmz = Rz
         ii  = j
       end if

     end do

     R1 = sqrt(Rmin2)
     dR = R1 - R0_steering
     ff = - k_steering * dR
     fc = ff / R1

     Fx = fc * Rmx
     Fy = fc * Rmy
     Fz = fc * Rmz

     Frc_OptC(1,Ipartner) = Frc_OptC(1,Ipartner) + Fx
     Frc_OptC(2,Ipartner) = Frc_OptC(2,Ipartner) + Fy
     Frc_OptC(3,Ipartner) = Frc_OptC(3,Ipartner) + Fz

     Frc_OptC(1,ii) = Frc_OptC(1,ii) - Fx
     Frc_OptC(2,ii) = Frc_OptC(2,ii) - Fy
     Frc_OptC(3,ii) = Frc_OptC(3,ii) - Fz

     Ene_OptC = Ene_OptC - 0.5d0 * ff * dR

     Vir_OptC(1,1) = Vir_OptC(1,1) + Fx * Rmx
     Vir_OptC(1,2) = Vir_OptC(1,2) + Fx * Rmy
     Vir_OptC(1,3) = Vir_OptC(1,3) + Fx * Rmz
     Vir_OptC(2,1) = Vir_OptC(2,1) + Fy * Rmx
     Vir_OptC(2,2) = Vir_OptC(2,2) + Fy * Rmy
     Vir_OptC(2,3) = Vir_OptC(2,3) + Fy * Rmz
     Vir_OptC(3,1) = Vir_OptC(3,1) + Fz * Rmx
     Vir_OptC(3,2) = Vir_OptC(3,2) + Fz * Rmy
     Vir_OptC(3,3) = Vir_OptC(3,3) + Fz * Rmz

   else if(NumOpFixS==1) then

     InvMGi = 1.d0/MassRG(1)

     Rcomi(:) = 0.d0
     if(QRigidBody) then
       do j = NumFixF(1), NumFixL(1)
         Rcomi(:) = Rcomi(:) + MassRB(j)*R_RB(:,j)
       end do
     else
       do j = NumFixF(1), NumFixL(1)
         Rcomi(:) = Rcomi(:) + Mass(j)*R(:,j)
       end do
     end if
     Rcomi(:) = Rcomi(:) * InvMGi

     Rpos = dot_product(Rcomi(:),Uvec_steering(:))

     dR = Rpos - R0_steering
     ff = - k_steering * dR
     Fi(:) = ff * Uvec_steering(:)

     if(QRigidBody) then
       do i = NumFixF(1), NumFixL(1)
         j = InitAtomRB(i)
         MyType = RBType(i)
         do k = 1, NumRBAtom(MyType)
           i1 = j + k
           Frc_OptC(:,i1) = Frc_OptC(:,i1) + Fi(:) * Mass(i1) * InvMGi
         end do
       end do
     else
       do i = NumFixF(1), NumFixL(1)
         Frc_OptC(:,i) = Frc_OptC(:,i) + Fi(:) * Mass(i) * InvMGi
       end do
     end if

     Ene_OptC = Ene_OptC - 0.5d0 * ff * dR

     if(QPBC) then
       Vir_OptC(1,1) = Vir_OptC(1,1) + Fi(1) * Rcomi(1)
       Vir_OptC(1,2) = Vir_OptC(1,2) + Fi(1) * Rcomi(2)
       Vir_OptC(1,3) = Vir_OptC(1,3) + Fi(1) * Rcomi(3)
       Vir_OptC(2,1) = Vir_OptC(2,1) + Fi(2) * Rcomi(1)
       Vir_OptC(2,2) = Vir_OptC(2,2) + Fi(2) * Rcomi(2)
       Vir_OptC(2,3) = Vir_OptC(2,3) + Fi(2) * Rcomi(3)
       Vir_OptC(3,1) = Vir_OptC(3,1) + Fi(3) * Rcomi(1)
       Vir_OptC(3,2) = Vir_OptC(3,2) + Fi(3) * Rcomi(2)
       Vir_OptC(3,3) = Vir_OptC(3,3) + Fi(3) * Rcomi(3)
     end if

   else if(NumOpFixP==1) then

     InvMGi = 1.d0/MassRGi(1)
     InvMGj = 1.d0/MassRGj(1)

     Rcomi(:) = 0.d0
     Rcomj(:) = 0.d0
     if(QRigidBody) then
       do j = NumFixFi(1), NumFixLi(1)
         Rcomi(:) = Rcomi(:) + MassRB(j)*R_RB(:,j)
       end do
       do j = NumFixFj(1), NumFixLj(1)
         Rcomj(:) = Rcomj(:) + MassRB(j)*R_RB(:,j)
       end do
     else
       do j = NumFixFi(1), NumFixLi(1)
         Rcomi(:) = Rcomi(:) + Mass(j)*R(:,j)
       end do
       do j = NumFixFj(1), NumFixLj(1)
         Rcomj(:) = Rcomj(:) + Mass(j)*R(:,j)
       end do
     end if
     Rcomi(:) = Rcomi(:) * InvMGi
     Rcomj(:) = Rcomj(:) * InvMGj

     if(QPBC) then
       ScRi(:) = matmul(InvH,Rcomi(:))
       ScRj(:) = matmul(InvH,Rcomj(:))
       Sij(:) = ScRi(:) - ScRj(:)
       Sij(:) = Sij(:) - nint(Sij(:))
       Rij(:) = matmul(H,Sij(:))
     else
       Rij(:) = Rcomi(:) - Rcomj(:)
     end if

     R2 = dot_product(Rij,Rij)
     Rpos = sqrt(R2)

     dR = Rpos - R0_steering
     ff = - k_steering * dR
     Fi(:) = ff * Rij(:) / Rpos
     Fj(:) = -Fi(:)

     if(QRigidBody) then
       do i = NumFixFi(1), NumFixLi(1)
         j = InitAtomRB(i)
         MyType = RBType(i)
         do k = 1, NumRBAtom(MyType)
           i1 = j + k
           Frc_OptC(:,i1) = Frc_OptC(:,i1) + Fi(:) * Mass(i1) * InvMGi
         end do
       end do
       do i = NumFixFj(1), NumFixLj(1)
         j = InitAtomRB(i)
         MyType = RBType(i)
         do k = 1, NumRBAtom(MyType)
           i1 = j + k
           Frc_OptC(:,i1) = Frc_OptC(:,i1) + Fj(:) * Mass(i1) * InvMGj
         end do
       end do
     else
       do i = NumFixFi(1), NumFixLi(1)
         Frc_OptC(:,i) = Frc_OptC(:,i) + Fi(:) * Mass(i) * InvMGi
       end do
       do i = NumFixFj(1), NumFixLj(1)
         Frc_OptC(:,i) = Frc_OptC(:,i) + Fj(:) * Mass(i) * InvMGj
       end do
     end if

     Ene_OptC = Ene_OptC - 0.5d0 * ff * dR

     if(QPBC) then
       Vir_OptC(1,1) = Vir_OptC(1,1) + Fi(1) * Rij(1)
       Vir_OptC(1,2) = Vir_OptC(1,2) + Fi(1) * Rij(2)
       Vir_OptC(1,3) = Vir_OptC(1,3) + Fi(1) * Rij(3)
       Vir_OptC(2,1) = Vir_OptC(2,1) + Fi(2) * Rij(1)
       Vir_OptC(2,2) = Vir_OptC(2,2) + Fi(2) * Rij(2)
       Vir_OptC(2,3) = Vir_OptC(2,3) + Fi(2) * Rij(3)
       Vir_OptC(3,1) = Vir_OptC(3,1) + Fi(3) * Rij(1)
       Vir_OptC(3,2) = Vir_OptC(3,2) + Fi(3) * Rij(2)
       Vir_OptC(3,3) = Vir_OptC(3,3) + Fi(3) * Rij(3)
     end if

   else if(NSelCyl/=0) then

     FFcyl = 0.d0

     do i = 1, NSelCyl

       j = IdSelC(i)
       R2 = R(icx1,j)*R(icx1,j) + R(icx2,j)*R(icx2,j)
       Rpos = sqrt(R2)

       dR = Rpos - R0_steering
       if(dR < 0.d0) then
         dR2 = dR * dR
         ff = k_steering * dR2
         fc = ff / Rpos
         Ene_OptC = Ene_OptC - ff * dR * 0.3333333333333333d0
         Fi(icx1) = fc * R(icx1,j)
         Fi(icx2) = fc * R(icx2,j)
         Frc_OptC(icx1,j) = Frc_OptC(icx1,j) + Fi(icx1)
         Frc_OptC(icx2,j) = Frc_OptC(icx2,j) + Fi(icx2)
        if(QPBC) then
           Vir_OptC(icx1,1) = Vir_OptC(icx1,1) + Fi(icx1) * R(1,j)
           Vir_OptC(icx1,2) = Vir_OptC(icx1,2) + Fi(icx1) * R(2,j)
           Vir_OptC(icx1,3) = Vir_OptC(icx1,3) + Fi(icx1) * R(3,j)
           Vir_OptC(icx2,1) = Vir_OptC(icx2,1) + Fi(icx2) * R(1,j)
           Vir_OptC(icx2,2) = Vir_OptC(icx2,2) + Fi(icx2) * R(2,j)
           Vir_OptC(icx2,3) = Vir_OptC(icx2,3) + Fi(icx2) * R(3,j)
         end if
         FFcyl = FFcyl + ff
       end if

     end do

   else if(NSelCylI/=0) then

     FFcyl = 0.d0

     do i = 1, NSelCylI

       j = IdSelC(i)
       R2 = R(icx1,j)*R(icx1,j) + R(icx2,j)*R(icx2,j)
       Rpos = sqrt(R2)

       dR = Rpos - R0_steering
       if(dR > 0.d0) then
         dR2 = dR * dR
         ff = - k_steering * dR2
         fc = ff / Rpos
         Ene_OptC = Ene_OptC - ff * dR * 0.3333333333333333d0
         Fi(icx1) = fc * R(icx1,j)
         Fi(icx2) = fc * R(icx2,j)
         Frc_OptC(icx1,j) = Frc_OptC(icx1,j) + Fi(icx1)
         Frc_OptC(icx2,j) = Frc_OptC(icx2,j) + Fi(icx2)
        if(QPBC) then
           Vir_OptC(icx1,1) = Vir_OptC(icx1,1) + Fi(icx1) * R(1,j)
           Vir_OptC(icx1,2) = Vir_OptC(icx1,2) + Fi(icx1) * R(2,j)
           Vir_OptC(icx1,3) = Vir_OptC(icx1,3) + Fi(icx1) * R(3,j)
           Vir_OptC(icx2,1) = Vir_OptC(icx2,1) + Fi(icx2) * R(1,j)
           Vir_OptC(icx2,2) = Vir_OptC(icx2,2) + Fi(icx2) * R(2,j)
           Vir_OptC(icx2,3) = Vir_OptC(icx2,3) + Fi(icx2) * R(3,j)
         end if
         FFcyl = FFcyl + ff
       end if

     end do

   else if(NSelCone/=0) then

! ## constraint for COM
     Rcon(:,1) = 0.d0
     do i = 1, NSelCone
       j = IdSelC(i)
       Rcon(:,1) = Rcon(:,1) + Mass(j) * R(:,j)
     end do
     Rcon(:,1) = Rcon(:,1) * InvMscom(1)
     Rij(:) = Rcon(:,1) - Xcon(:,1)
     fcom(:,1) = - kccom * Rij(:)
     if(QPBC) then
       Vir_OptC(1,1) = Vir_OptC(1,1) + fcom(1,1) * Rcon(1,1)
       Vir_OptC(1,2) = Vir_OptC(1,2) + fcom(1,1) * Rcon(2,1)
       Vir_OptC(1,3) = Vir_OptC(1,3) + fcom(1,1) * Rcon(3,1)
       Vir_OptC(2,1) = Vir_OptC(2,1) + fcom(2,1) * Rcon(1,1)
       Vir_OptC(2,2) = Vir_OptC(2,2) + fcom(2,1) * Rcon(2,1)
       Vir_OptC(2,3) = Vir_OptC(2,3) + fcom(2,1) * Rcon(3,1)
       Vir_OptC(3,1) = Vir_OptC(3,1) + fcom(3,1) * Rcon(1,1)
       Vir_OptC(3,2) = Vir_OptC(3,2) + fcom(3,1) * Rcon(2,1)
       Vir_OptC(3,3) = Vir_OptC(3,3) + fcom(3,1) * Rcon(3,1)
     end if
     Fi(:) = fcom(:,1) * InvMscom(1)
     do i = 1, NSelCone
       j = IdSelC(i)
       Frc_OptC(:,j) = Frc_OptC(:,j) + Fi(:) * Mass(j)
     end do
     Ene_OptC = Ene_OptC + 0.5d0 * kccom * dot_product( Rij, Rij )

! ## cone potential
     FFcone = 0.d0

     theta = R0_steering/180.d0*pi
     cst = cos(theta)
     snt = sin(theta)

     cst1 = (1+cst)
     csti = 1.d0/cst1

     pref1 = 0.25 * pi / Vol_Lipo * (cst * cst - 1.d0)

     do i = 1, NSelCone

       j = IdSelC(i)
       zz  = R(icaxis,j)

       if(zz < 0.d0) cycle

       if(QCorrCone) then

         fa = facta * cst1
         fb = factb * cst1
         fc = factc * cst1 - Vol_Lipo
         fd = fb * fb - 4.d0 * fa * fc
         rg = 0.5d0*(-fb + sqrt(fd))*Invfacta*csti

         x3 = rg*rg
         x4 = x3*x3
         x5 = rg-Thick_Lipo
         x6 = x5*x5
         x7 = x6*x6
         ShiftCone = pref1 * (x4-x7)

         zz = zz + ShiftCone

       end if

       R2 = R(icx1,j)*R(icx1,j) + R(icx2,j)*R(icx2,j)
       Rpos = sqrt(R2)

       dR  = Rpos*cst - zz*snt

       if(dR < 0.d0) then
         dR2 = dR * dR
         ff = k_steering * dR2
         fc = ff / Rpos * cst
         Fi(icx1) = fc * R(icx1,j)
         Fi(icx2) = fc * R(icx2,j)
         Fi(icaxis) = - ff * snt
         Ene_OptC = Ene_OptC - ff * dR * 0.3333333333333333d0
         Frc_OptC(icx1,j)   = Frc_OptC(icx1,j)   + Fi(icx1)
         Frc_OptC(icx2,j)   = Frc_OptC(icx2,j)   + Fi(icx2)
         Frc_OptC(icaxis,j) = Frc_OptC(icaxis,j) + Fi(icaxis)
         if(QPBC) then
           Vir_OptC(icx1,1)   = Vir_OptC(icx1,1)   + Fi(icx1)   * R(1,j)
           Vir_OptC(icx1,2)   = Vir_OptC(icx1,2)   + Fi(icx1)   * R(2,j)
           Vir_OptC(icx1,3)   = Vir_OptC(icx1,3)   + Fi(icx1)   * R(3,j)
           Vir_OptC(icx2,1)   = Vir_OptC(icx2,1)   + Fi(icx2)   * R(1,j)
           Vir_OptC(icx2,2)   = Vir_OptC(icx2,2)   + Fi(icx2)   * R(2,j)
           Vir_OptC(icx2,3)   = Vir_OptC(icx2,3)   + Fi(icx2)   * R(3,j)
           Vir_OptC(icaxis,1) = Vir_OptC(icaxis,1) + Fi(icaxis) * R(1,j)
           Vir_OptC(icaxis,2) = Vir_OptC(icaxis,2) + Fi(icaxis) * R(2,j)
           Vir_OptC(icaxis,3) = Vir_OptC(icaxis,3) + Fi(icaxis) * R(3,j)
         end if
         FFcone = FFcone + ff * (Rpos*snt + zz*cst)
       end if

     end do

   end if

end subroutine Force_Jarzynski
