! ############################
! ## SUBROUTINE LIST 
! ## -- Force_OptC 
! ## -- Force_Harmonic 
! ## -- Force_Cluster
! ## -- Force_repulsivewall 
! ## -- Force_wconstraint 
! ## -- Force_steelewall
! ## -- Set_Steele
! ## -- Force_CGwall
! ## -- Force_cylinder
! ## -- Force_CGcylinder
! ## -- Force_AAcylinder
! ## -- CylPMFsample
! ## -- Force_FScylinder
! ## -- Force_FixonPlane
! ## -- Force_HamC 
! ## -- Force_Constraint 
! ## -- Force_Div_Component 
! ## -- SPCF_Intra
! ## -- Set_MembrPot
! ## -- Force_MembrPot
! ############################


!######################################################################
!######################################################################


subroutine Force_OptC

use Numbers, only : N
use CommonBlocks, only : QREPWALL, QSTWALL, QMP, QSPCF, QOption, &
&   QRigidBody, QJarzynski, QCGWALL, Qwcformol, QMaster, QMacro, &
&   ForceField, QPLC, QClust, QInert, QCyl, QFSCyl, QFixCOM, QEfield
use Configuration, only : R
use RBparam, only : R_RB, AtomUnitNum
use OptConstraintParam, only : Frc_OptC, Vir_OptC, Ene_OptC
use CGball, only : FSphRs

implicit NONE

integer :: i, j
real(8) :: Rx, Ry, Rz
real(8) :: Fx, Fy, Fz
real(8) :: Vxx, Vxy, Vxz, Vyx, Vyy, Vyz, Vzx, Vzy, Vzz

   Ene_OptC = 0.d0
   Frc_OptC = 0.d0
   Vir_OptC = 0.d0
   if(QMacro) then
     FSphRs(:) = 0.d0
   end if

   if(QREPWALL) call Force_repulsivewall
   if(QSTWALL)  call Force_steelewall
   if(QCGWALL)  call Force_CGwall
   if(QPLC)     call Force_FixonPlane
   if(QMP)      call Force_MembrPot
   if(QSPCF)    call SPCF_Intra
   if(QInert)   call ContInertia
   if(QCyl)     call Force_cylinder
   if(QFSCyl)   call Force_FScylinder
   if(QFixCOM)  call Force_ConstCOM
   if(QEfield)  call Force_Efield

   if(QREPWALL.or.QSTWALL.or.QCGWALL.or.QSPCF.or.QOption &
   & .or.QInert.or.QCyl.or.QFSCyl.or.QFixCOM) then

     Vxx = 0.d0
     Vxy = 0.d0
     Vxz = 0.d0
     Vyy = 0.d0
     Vyz = 0.d0
     Vzz = 0.d0

     if(QRigidBody) then

       Vyx = 0.d0
       Vzx = 0.d0
       Vzy = 0.d0

       do i = 1 , N

         j = AtomUnitNum(i)

         Fx = Frc_OptC(1,i)
         Fy = Frc_OptC(2,i)
         Fz = Frc_OptC(3,i)
         Rx = R_RB(1,j)
         Ry = R_RB(2,j)
         Rz = R_RB(3,j)

         Vxx = Vxx + Fx * Rx
         Vxy = Vxy + Fx * Ry
         Vxz = Vxz + Fx * Rz
         Vyx = Vyx + Fy * Rx
         Vyy = Vyy + Fy * Ry
         Vyz = Vyz + Fy * Rz
         Vzx = Vzx + Fz * Rx
         Vzy = Vzy + Fz * Ry
         Vzz = Vzz + Fz * Rz

       end do

     else

       do i = 1 , N

         Fx = Frc_OptC(1,i)
         Fy = Frc_OptC(2,i)
         Fz = Frc_OptC(3,i)
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

       Vyx = Vxy
       Vzx = Vxz
       Vzy = Vyz

     end if

     Vir_OptC(1,1) = Vxx
     Vir_OptC(1,2) = Vxy
     Vir_OptC(1,3) = Vxz
     Vir_OptC(2,1) = Vyx
     Vir_OptC(2,2) = Vyy
     Vir_OptC(2,3) = Vyz
     Vir_OptC(3,1) = Vzx
     Vir_OptC(3,2) = Vzy
     Vir_OptC(3,3) = Vzz

   end if

   if(QOption)  call Force_Harmonic
   if(QClust.and.QMaster)  call Force_Cluster
   if(Qwcformol) call Force_wconstraint
   if(QJarzynski) call Force_Jarzynski
   if(QMacro) then
     if(ForceField=='OPLS') then
       call Force_MacroSphAA
     else if(ForceField=='CG') then
       call Force_MacroSphCG
     end if
   end if

end subroutine Force_OptC


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! *******************************

subroutine Force_Harmonic

#ifdef GEN
use Numbers, only : N
#endif
use CommonBlocks, only : QMacro
use Configuration, only : R
use OptConstraintParam, only : NumOptC, OptCI, OptCJ, kOptC, rOptC, &
&   Frc_OptC, Ene_OptC, Vir_OptC
use CGball, only : IDsphere, FSphRs
use CellParam, only : H, CellShft

implicit NONE

integer :: i, j, k, l
real(8) :: R2, R1, dR, Fc, FF
real(8) :: Rx, Ry, Rz, Fx, Fy, Fz
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
integer :: Nx, Ny, Nz
integer :: id, jd
#ifdef GEN
real(8), dimension(3,N) :: ScR
real(8) :: Sx, Sy, Sz
#else
real(8) :: clhx, clhy, clhz
#endif

#ifdef GEN
   call ScaledCoordinate(ScR)
#else
   clhx = H(1,1)*0.5d0
   clhy = H(2,2)*0.5d0
   clhz = H(3,3)*0.5d0
#endif

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

if(QMacro) then

   do k = 1 , NumOptC

     i= OptCI(k)
     j= OptCJ(k)

     Rx = R(1,i) - R(1,j)
     Ry = R(2,i) - R(2,j)
     Rz = R(3,i) - R(3,j)
#ifdef GEN
     Sx = ScR(1,i) - ScR(1,j)
     Sy = ScR(2,i) - ScR(2,j)
     Sz = ScR(3,i) - ScR(3,j)
     if(Sx>0.5) then
       Nx = -9
     else if(Sx<-0.5) then
       Nx =  9
     else
       Nx =  0
     end if
     if(Sy>0.5) then
       Ny = -3
     else if(Sy<-0.5) then
       Ny =  3
     else
       Ny =  0
     end if
     if(Sz>0.5) then
       Nz = -1
     else if(Sz<-0.5) then
       Nz =  1
     else
       Nz =  0
     end if
#else
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
#endif
     l = Nx + Ny + Nz
     Rx = Rx + CellShft(1,l)
     Ry = Ry + CellShft(2,l)
     Rz = Rz + CellShft(3,l)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     R1 = sqrt(R2)
     dR = R1 - rOptC(k)
     Fc = kOptC(k) * dR

     Ene_OptC = Ene_OptC + Fc * dR

     FF = - 2.d0 * Fc / R1
     Fx = FF * Rx
     Fy = FF * Ry
     Fz = FF * Rz

     Frc_OptC(1,i) = Frc_OptC(1,i) + Fx
     Frc_OptC(2,i) = Frc_OptC(2,i) + Fy
     Frc_OptC(3,i) = Frc_OptC(3,i) + Fz
     Frc_OptC(1,j) = Frc_OptC(1,j) - Fx
     Frc_OptC(2,j) = Frc_OptC(2,j) - Fy
     Frc_OptC(3,j) = Frc_OptC(3,j) - Fz

     id = IDsphere(i)
     jd = IDsphere(j)
     if(id/=0) FSphRs(id) = FSphRs(id) - 2.d0 * Fc
     if(jd/=0) FSphRs(jd) = FSphRs(jd) - 2.d0 * Fc

     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz

   end do

else

   do k = 1 , NumOptC

     i= OptCI(k)
     j= OptCJ(k)

     Rx = R(1,i) - R(1,j)
     Ry = R(2,i) - R(2,j)
     Rz = R(3,i) - R(3,j)
#ifdef GEN
     Sx = ScR(1,i) - ScR(1,j)
     Sy = ScR(2,i) - ScR(2,j)
     Sz = ScR(3,i) - ScR(3,j)
     if(Sx>0.5) then
       Nx = -9
     else if(Sx<-0.5) then
       Nx =  9
     else
       Nx =  0
     end if
     if(Sy>0.5) then
       Ny = -3
     else if(Sy<-0.5) then
       Ny =  3
     else
       Ny =  0
     end if
     if(Sz>0.5) then
       Nz = -1
     else if(Sz<-0.5) then
       Nz =  1
     else
       Nz =  0
     end if
#else
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
#endif
     l = Nx + Ny + Nz
     Rx = Rx + CellShft(1,l)
     Ry = Ry + CellShft(2,l)
     Rz = Rz + CellShft(3,l)
     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     R1 = sqrt(R2)
     dR = R1 - rOptC(k)
     Fc = kOptC(k) * dR

     Ene_OptC = Ene_OptC + Fc * dR

     FF = - 2.d0 * Fc / R1
     Fx = FF * Rx
     Fy = FF * Ry
     Fz = FF * Rz

     Frc_OptC(1,i) = Frc_OptC(1,i) + Fx
     Frc_OptC(2,i) = Frc_OptC(2,i) + Fy
     Frc_OptC(3,i) = Frc_OptC(3,i) + Fz
     Frc_OptC(1,j) = Frc_OptC(1,j) - Fx
     Frc_OptC(2,j) = Frc_OptC(2,j) - Fy
     Frc_OptC(3,j) = Frc_OptC(3,j) - Fz

     Vxx = Vxx + Fx * Rx
     Vxy = Vxy + Fx * Ry
     Vxz = Vxz + Fx * Rz
     Vyy = Vyy + Fy * Ry
     Vyz = Vyz + Fy * Rz
     Vzz = Vzz + Fz * Rz

   end do

end if

   Vir_OptC(1,1) = Vir_OptC(1,1) + Vxx
   Vir_OptC(1,2) = Vir_OptC(1,2) + Vxy
   Vir_OptC(1,3) = Vir_OptC(1,3) + Vxz
   Vir_OptC(2,1) = Vir_OptC(2,1) + Vxy
   Vir_OptC(2,2) = Vir_OptC(2,2) + Vyy
   Vir_OptC(2,3) = Vir_OptC(2,3) + Vyz
   Vir_OptC(3,1) = Vir_OptC(3,1) + Vxz
   Vir_OptC(3,2) = Vir_OptC(3,2) + Vyz
   Vir_OptC(3,3) = Vir_OptC(3,3) + Vzz

end subroutine Force_Harmonic


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! *******************************

subroutine Force_Cluster

use CommonBlocks, only : QMacro
use Configuration, only : R
use OptConstraintParam, only : NumClust, NAtomClus, Rsh_clust, k_clust, &
&   Frc_OptC, Ene_OptC, Vir_OptC
#ifdef GEN
use CellParam, only : H, InvH, CellShft
#else
use CellParam, only : H, CellShft
#endif

implicit NONE

integer :: i, j, k, l, ii, jj
integer :: Nx, Ny, Nz
real(8) :: Rx, Ry, Rz, Fx, Fy, Fz
real(8) :: Rix, Riy, Riz
real(8) :: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
#ifdef GEN
real(8) :: Six, Siy, Siz
real(8), dimension(3,NumClust) :: ScR
real(8) :: Sx, Sy, Sz
#else
real(8) :: clhx, clhy, clhz
#endif
real(8) :: dmin, dr, aa, fc
real(8), dimension(3,NumClust) :: Rclust, Fclust
real(8), dimension(NumClust,NumClust) :: R2list
real(8), dimension(3,NumClust,NumClust) :: Rlist
logical, dimension(NumClust) :: Qdone

   do i = 1, NumClust
     j = NAtomClus(i)
     Rclust(:,i) = R(:,j)
   end do

#ifdef GEN
   do i = 1, NumClust
     ScR(:,i) = matmul(InvH, Rclust(:,i))
   end do
#else
   clhx = H(1,1)*0.5d0
   clhy = H(2,2)*0.5d0
   clhz = H(3,3)*0.5d0
#endif

   Fclust = 0.d0

   Vxx = 0.d0
   Vxy = 0.d0
   Vxz = 0.d0
   Vyy = 0.d0
   Vyz = 0.d0
   Vzz = 0.d0

   do i = 1, NumClust-1
     Rix = Rclust(1,i)
     Riy = Rclust(2,i)
     Riz = Rclust(3,i)
#ifdef GEN
     Six = ScR(1,i)
     Siy = ScR(2,i)
     Siz = ScR(3,i)
#endif
     do j = i+1, NumClust
       Rx = Rix - Rclust(1,j)
       Ry = Riy - Rclust(2,j)
       Rz = Riz - Rclust(3,j)
#ifdef GEN
       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       if(Sx>0.5) then
         Nx = -9
       else if(Sx<-0.5) then
         Nx =  9
       else
         Nx =  0
       end if
       if(Sy>0.5) then
         Ny = -3
       else if(Sy<-0.5) then
         Ny =  3
       else
         Ny =  0
       end if
       if(Sz>0.5) then
         Nz = -1
       else if(Sz<-0.5) then
         Nz =  1
       else
         Nz =  0
       end if
#else
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
#endif
       l = Nx + Ny + Nz
       Rx = Rx + CellShft(1,l)
       Ry = Ry + CellShft(2,l)
       Rz = Rz + CellShft(3,l)
       Rlist(1,i,j) = Rx
       Rlist(2,i,j) = Ry
       Rlist(3,i,j) = Rz
       R2list(i,j) = Rx*Rx + Ry*Ry + Rz*Rz
       Rlist(1,j,i) = -Rx
       Rlist(2,j,i) = -Ry
       Rlist(3,j,i) = -Rz
       R2list(j,i) = Rx*Rx + Ry*Ry + Rz*Rz
     end do
   end do

   Qdone(:) = .False.
   Qdone(1) = .True.

   do k = 1, NumClust - 1

     dmin = 1.d+9

     do i = 1, NumClust

       if(Qdone(i)) then

         do j = 1, NumClust

           if(Qdone(j)) cycle

           if(R2list(i,j)<dmin) then
             dmin=R2list(i,j)
             ii = i
             jj = j
           end if

         end do

       end if

     end do

     dmin = sqrt(dmin)
     if(dmin>Rsh_clust) then
       dr = dmin - Rsh_clust
       aa = k_clust * dr * dr
       fc = - 3.d0 * aa / dmin
       Ene_OptC = Ene_OptC + aa * dr
       Rx = Rlist(1,ii,jj)
       Ry = Rlist(2,ii,jj)
       Rz = Rlist(3,ii,jj)
       fx = fc * Rx
       fy = fc * Ry
       fz = fc * Rz
       Fclust(1,ii) = Fclust(1,ii) + fx
       Fclust(2,ii) = Fclust(2,ii) + fy
       Fclust(3,ii) = Fclust(3,ii) + fz
       Fclust(1,jj) = Fclust(1,jj) - fx
       Fclust(2,jj) = Fclust(2,jj) - fy
       Fclust(3,jj) = Fclust(3,jj) - fz
       Vxx = Vxx + Fx * Rx
       Vxy = Vxy + Fx * Ry
       Vxz = Vxz + Fx * Rz
       Vyy = Vyy + Fy * Ry
       Vyz = Vyz + Fy * Rz
       Vzz = Vzz + Fz * Rz
     end if
     Qdone(jj) = .True.

   end do

#ifdef DEBUG
!# for check
   do i = 1, NumClust
     if(.not.Qdone(i)) then
       write(*,*) "wrong calc in cluster"
       call Finalize
     end if
   end do
#endif

   do i = 1, NumClust
     j = NAtomClus(i)
     Frc_OptC(:,j) = Frc_OptC(:,j) + Fclust(:,i)
   end do

   Vir_OptC(1,1) = Vir_OptC(1,1) + Vxx
   Vir_OptC(1,2) = Vir_OptC(1,2) + Vxy
   Vir_OptC(1,3) = Vir_OptC(1,3) + Vxz
   Vir_OptC(2,1) = Vir_OptC(2,1) + Vxy
   Vir_OptC(2,2) = Vir_OptC(2,2) + Vyy
   Vir_OptC(2,3) = Vir_OptC(2,3) + Vyz
   Vir_OptC(3,1) = Vir_OptC(3,1) + Vxz
   Vir_OptC(3,2) = Vir_OptC(3,2) + Vyz
   Vir_OptC(3,3) = Vir_OptC(3,3) + Vzz

end subroutine Force_Cluster


!######################################################################
!######################################################################


subroutine Force_repulsivewall

use Numbers, only : N
use CommonBlocks, only : QPathInt
use Configuration, only : R
use CommonMPI
use CommonPI
use WallParam
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, j, Nas, ix, iy, k
real(8) :: dR, R2, Rcap2, R1, Fc, xx, Rpore2, dx
real(8), dimension(3) :: Rij
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   if(Ccap=='PLANER') then

     if(Rpore/=0.) then

       do j = 1, nwall

         xx = Rpcap(j)

         if(Ifcdir(j)==-1) then

           do k = Nas, NSelR, NProcsTemp

             i = IdSelR(k)
             dx = R(Lcap,i) - xx
             if((dx >0.).and.(dx<3.d0)) then

               dR = R(Lcap,i) - xx
               Frc_OptC(Lcap,i) = Frc_OptC(Lcap,i) - 2.d0 * kcap * dR
               Ene_OptC = Ene_OptC + kcap * dR * dR

             end if

           end do

         else if(Ifcdir(j)==1) then

           do k = Nas, NSelR, NProcsTemp

             i = IdSelR(k)
             dx = R(Lcap,i) - xx
             if((dx<0.).and.(dx>-3.)) then

               dR = R(Lcap,i) - xx
               Frc_OptC(Lcap,i) = Frc_OptC(Lcap,i) - 2.d0 * kcap * dR
               Ene_OptC = Ene_OptC + kcap * dR * dR

             end if

           end do

         end if

       end do

     else

       if(Lcap==1) then
         ix = 2
         iy = 3
       else if(Lcap==2) then
         ix = 1
         iy = 3
       else if(Lcap==3) then
         ix = 1
         iy = 2
       end if

       Rpore2 = Rpore*Rpore

       do j = 1, nwall

         xx = Rpcap(j)

         if(Ifcdir(j)==-1) then

           do k = Nas, NSelR, NProcsTemp

             i = IdSelR(k)
             R2 = R(ix,i)*R(ix,i)+R(iy,i)*R(iy,i)

             dx = R(Lcap,i) - xx
             if((dx >0.).and.(dx<3.d0).and.(R2>Rpore2)) then

               dR = R(Lcap,i) - xx
               Frc_OptC(Lcap,i) = Frc_OptC(Lcap,i) - 2.d0 * kcap * dR
               Ene_OptC = Ene_OptC + kcap * dR * dR

             end if

           end do

         else if(Ifcdir(j)==1) then

           do k = Nas, NSelR, NProcsTemp

             i = IdSelR(k)
             R2 = R(ix,i)*R(ix,i)+R(iy,i)*R(iy,i)

             dx = R(Lcap,i) - xx
             if((dx<0.).and.(dx>-3.).and.(R2>Rpore2)) then

               dR = R(Lcap,i) - xx
               Frc_OptC(Lcap,i) = Frc_OptC(Lcap,i) - 2.d0 * kcap * dR
               Ene_OptC = Ene_OptC + kcap * dR * dR

             end if

           end do

         end if

       end do

     end if

   else if(Ccap=='SPHERE') then

     Rcap2 = Rcap * Rcap

     do k = Nas, NSelR, NProcsTemp

       i = IdSelR(k)
       Rij = R(:,i) - Ocap
       R2  = dot_product(Rij,Rij)

       if(R2 > Rcap2) then

         R1 = sqrt(R2)
         dR = R1 - Rcap
         Fc = -kcap * dR / R1
         Frc_OptC(:,i) = Frc_OptC(:,i) + Fc * Rij
         Ene_OptC = Ene_OptC + kcap * dR * dR

       end if

     end do

   end if

end subroutine Force_repulsivewall


!######################################################################
!######################################################################


subroutine Force_wconstraint

use CommonBlocks, only : QMaster, QPathInt
use Configuration, only : R
use wcparam
use OptConstraintParam, only : Frc_OptC, Ene_OptC, Vir_OptC
use AtomParam, only : Mass
use CellParam, only : H
use CommonMPI

implicit NONE

integer :: i, j, Nas
real(8) :: dR, Fc, xx, CellLength, Rtmp

   Nas = NProcs - MyRank

   CellLength = H(Iaxis_const,Iaxis_const)

if(Nmolwc==0) then

   if(QMaster) then

   Ronc = 0.d0
   do i = atom_init, atom_fin
     Ronc = Ronc + R(Iaxis_const,i) * Mass(i)
   end do
   Ronc = Ronc * InvMconst_atom

   if(Qwc_edge) then
     dw_mid = CellLength * 0.5d0
     rc_low = dw_low - dw_mid
     rc_upp = - rc_low
   end if

   Rtmp = Ronc - dw_mid
   Rtmp = Rtmp - CellLength * nint(Rtmp/CellLength)

   Fc = 0.d0

   if(Rtmp > rc_upp) then

     dR = Rtmp - rc_upp
     xx = kc_wall * dR ** (corder-1)
     Fc = - dble(corder) * xx
     Ene_OptC = Ene_OptC + xx * dR

   else if(Rtmp < rc_low) then

     dR = rc_low - Rtmp
     xx = kc_wall * dR ** (corder-1)
     Fc = dble(corder) * xx
     Ene_OptC = Ene_OptC + xx * dR

   end if

   Vir_OptC(Iaxis_const,Iaxis_const) = &
   &   Vir_OptC(Iaxis_const,Iaxis_const) + Fc * Ronc

   Fc = Fc * InvMconst_atom
   do i = atom_init, atom_fin
     Frc_OptC(Iaxis_const,i) = Frc_OptC(Iaxis_const,i) + Fc * Mass(i)
   end do

   end if

else

   if(Qwc_edge) then
     dw_mid = CellLength * 0.5d0
     rc_low = dw_low - dw_mid
     rc_upp = - rc_low
   end if

   do i = Nas, Nmolwc, NProcs

     Ronc = 0.d0
     do j = idf(i), ide(i)
       Ronc = Ronc + R(Iaxis_const,j) * Mass(j)
     end do
     Ronc = Ronc * invMwc(i)

     Rtmp = Ronc - dw_mid
     Rtmp = Rtmp - CellLength * nint(Rtmp/CellLength)

     Fc = 0.d0

     if(Rtmp > rc_upp) then

       dR = Rtmp - rc_upp
       xx = kc_wall * dR ** (corder-1)
       Fc = - dble(corder) * xx
       Ene_OptC = Ene_OptC + xx * dR

     else if(Rtmp < rc_low) then

       dR = rc_low - Rtmp
       xx = kc_wall * dR ** (corder-1)
       Fc = dble(corder) * xx
       Ene_OptC = Ene_OptC + xx * dR

     end if

     Vir_OptC(Iaxis_const,Iaxis_const) = &
     &   Vir_OptC(Iaxis_const,Iaxis_const) + Fc * Ronc

     Fc = Fc * InvMconst_atom
     do j = idf(i), ide(i)
       Frc_OptC(Iaxis_const,j) = Frc_OptC(Iaxis_const,j) + Fc * Mass(j)
     end do

   end do

end if


end subroutine Force_wconstraint


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! *******************************

subroutine Force_steelewall

use Numbers, only : N
use CommonBlocks, only : QPathInt
use Configuration, only : R
use CommonMPI
use CommonPI
use WallParam
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, j, k, Nas
real(8) :: dR, xx
real(8) :: zz, rzz, rz2, rz4, rz5, rz10, rz11
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   do j = 1, nwall

     xx = Rstwall(j)

     if(FdirST(j)==-1) then

       do i = Nas, N, NProcsTemp

         dR = R(STdir,i) - xx

         do k = 1, nlayer

           zz  = dR - (k-1) * Delt_z
           rzz = 1.d0 / zz
           rz2 = rzz * rzz
           rz4 = rz2 * rz2
           rz5 = rz4 * rzz
           rz10 = rz5 * rz5
           rz11 = rz10 * rzz

           Frc_OptC(STdir,i) = Frc_OptC(STdir,i) + FAijST(i) * rz11 + FBijST(i) * rz5
           Ene_OptC = Ene_OptC + AijST(i) * rz10 + BijST(i) * rz4

         end do

       end do

     else if(FdirST(j)==1) then

       do i = Nas, N, NProcsTemp

         dR = R(STdir,i) - xx

         do k = 1, nlayer

           zz  = dR + (k-1) * Delt_z
           rzz = 1.d0 / zz
           rz2 = rzz * rzz
           rz4 = rz2 * rz2
           rz5 = rz4 * rzz
           rz10 = rz5 * rz5
           rz11 = rz10 * rzz

           Frc_OptC(STdir,i) = Frc_OptC(STdir,i) + FAijST(i) * rz11 + FBijST(i) * rz5
           Ene_OptC = Ene_OptC + AijST(i) * rz10 + BijST(i) * rz4

         end do

       end do

     end if

   end do

end subroutine Force_steelewall


!######################################################################
!######################################################################


subroutine Set_Steele

use CommonBlocks, only : ForceField, QMaster
use Numbers, only : N
use NonbondParam, only : EpsLJ, SgmLJ, Rminh
use UnitExParam, only : pi
use WallParam

implicit NONE

integer :: i
real(8) :: x, Inv_a1, Inv_a2, ras, fc
real(8), dimension(N) :: SigST, EpsST

   Inv_a1 = 1.d0 / dlattice_wall
   Inv_a2 = Inv_a1 * Inv_a1
   ras    = 1.154700538379252d0

   allocate( AijST(N) )
   allocate( BijST(N) )
   allocate( FAijST(N) )
   allocate( FBijST(N) )

   ep_wall = sqrt( ep_wall )
   sg_wall = sg_wall * 0.5d0

   do i = 1, N
     EpsST(i) = ep_wall * EpsLJ(i)
   end do

   if( ForceField(1:5) == 'CHARM' ) then
     fc = 1.d0 / ( 2.d0 ** (1.d0/6.d0) )
     do i = 1, N
       x = Rminh(i) * fc + sg_wall
       SigST(i) = x ** 6
     end do
   else if( ForceField(1:4) == 'OPLS' ) then
     do i = 1, N
       x = 0.5d0*SgmLJ(i) + sg_wall
       SigST(i) = x ** 6
     end do
   else
     if(QMaster) then
     write(*,*) 'error : whenever using the Steele wall, one should choose'
     write(*,*) '        (CHARMm or OPLS) force fileds! '
     end if
     call Finalize
   end if

   do i = 1, N
     BijST(i) = 4.d0 * pi * EpsST(i) * SigST(i) * Inv_a2 * ras
     AijST(i) = BijST(i) * 0.4d0 * SigST(i)
   end do

   BijST = - BijST

   FAijST = 10.d0 * AijST
   FBijST =  4.d0 * BijST

end subroutine Set_Steele


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! *******************************

subroutine Force_CGwall

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use WallParam, only : nstwall, Rstwall, Rwmin, Rwmax, invdz_size, STdir, &
&   FdirST, TabWall
use CGdata, only : NBAtomType
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, j, ii, Nas, itype
real(8) :: dR, xx, y1, y2, y3, z1, z2
real(8) :: zz, a1, a2, a3, b1, b2, dx, rsp, rsf
real(8) :: invdr

   Nas = NProcs - MyRank

   do j = 1, nstwall

     zz = Rstwall(j)

     if(FdirST(j)==-1) then

       do i = Nas, N, NProcs

         itype = NBAtomType(i)

         dR = zz - R(STdir,i) - Rwmin(itype)

         if(dR > Rwmax(itype)) cycle

         invdr = invdz_size(itype)

         xx = dR * invdr

         ii  = nint( xx )

         if(ii >= 1) then

         dx  = xx - dble(ii)

         y1  = TabWall( 1, ii  , itype )
         a1  = TabWall( 2, ii  , itype )
         y2  = TabWall( 1, ii+1, itype )
         a2  = TabWall( 2, ii+1, itype )
         y3  = TabWall( 1, ii+2, itype )
         a3  = TabWall( 2, ii+2, itype )

         z2  = y1 - 2.d0 * y2 + y3
         b2  = a1 - 2.d0 * a2 + a3
         z1  = - y1 + y3
         b1  = - a1 + a3
         rsp = 0.5d0 * (z1 + z2 * dx ) * dx + y2
         rsf = 0.5d0 * (b1 + b2 * dx ) * dx + a2

         else

         rsp = TabWall( 1, 1, itype )
         rsf = TabWall( 2, 1, itype )
         print *, 'Bad contact'

         end if

         Frc_OptC(STdir,i) = Frc_OptC(STdir,i) + rsf
         Ene_OptC = Ene_OptC + rsp

       end do

     else if(FdirST(j)==1) then

       do i = Nas, N, NProcs

         itype = NBAtomType(i)

         dR = R(STdir,i) - zz - Rwmin(itype)

         if(dR > Rwmax(itype)) cycle

         invdr = invdz_size(itype)

         xx = dR * invdr

         ii  = nint( xx )

         if(ii >= 1) then

         dx  = xx - dble(ii)

         y1  = TabWall( 1, ii  , itype )
         a1  = TabWall( 2, ii  , itype )
         y2  = TabWall( 1, ii+1, itype )
         a2  = TabWall( 2, ii+1, itype )
         y3  = TabWall( 1, ii+2, itype )
         a3  = TabWall( 2, ii+2, itype )

         z2  = y1 - 2.d0 * y2 + y3
         b2  = a1 - 2.d0 * a2 + a3
         z1  = - y1 + y3
         b1  = - a1 + a3
         rsp = 0.5d0 * (z1 + z2 * dx ) * dx + y2
         rsf = 0.5d0 * (b1 + b2 * dx ) * dx + a2

         else

         rsp = TabWall( 1, 1, itype )
         rsf = TabWall( 2, 1, itype )
         print *, 'Bad contact'

         end if

         Frc_OptC(STdir,i) = Frc_OptC(STdir,i) - rsf
         Ene_OptC = Ene_OptC + rsp

       end do

     end if

   end do

end subroutine Force_CGwall


!######################################################################
!######################################################################


subroutine Force_cylinder

use CommonBlocks, only : ForceField

   if(ForceField=='CG') then
     call Force_CGcylinder
   else
     call Force_AAcylinder
   end if

end subroutine Force_cylinder


!######################################################################
!######################################################################


subroutine Force_CGcylinder

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use CylParam, only : icx1, icx2, invdr_size, Rcylmin2, Rcylmax2, TabCyl, &
&   FCylRs
use CGdata, only : NBAtomType
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, Nas, itype
real(8) :: U, Fc, Fr, FF, F1, F2
real(8) :: dx1, dx2, dR2, Rm
real(8) :: InterPolateCylinder
external InterPolateCylinder

   Nas = NProcs - MyRank

   FF = 0.d0

   do i = Nas, N, NProcs

     itype = NBAtomType(i)

     dx1 = R(icx1,i)
     dx2 = R(icx2,i)

     dR2 = dx1*dx1 + dx2*dx2

     if(dR2 > Rcylmax2(itype)) cycle

     Rm = Rcylmin2(itype)

     U  =   InterPolateCylinder(dR2,Rm,invdr_size(itype),1,itype)
     Fc = - InterPolateCylinder(dR2,Rm,invdr_size(itype),2,itype)
     Fr =   InterPolateCylinder(dR2,Rm,invdr_size(itype),3,itype)

     F1 = Fc * dx1
     F2 = Fc * dx2

     Frc_OptC(icx1,i) = Frc_OptC(icx1,i) + F1
     Frc_OptC(icx2,i) = Frc_OptC(icx2,i) + F2

     Ene_OptC = Ene_OptC + U

     FF = FF + Fr

   end do

   FCylRs = FF

end subroutine Force_CGcylinder


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! *******************************

subroutine Force_AAcylinder

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use CylParam, only : icx1, icx2, invdr_size, Rcylmin2, Rcylmax2, &
&   TabCyl, FCylRs
use AtomParam, only : NBAAType
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, Nas, itype
real(8) :: U, Fc, Fr, FF, F1, F2
real(8) :: dx1, dx2, dR2, Rm
real(8) :: InterPolateCylinder
external InterPolateCylinder

   Nas = NProcs - MyRank

   FF = 0.d0

   do i = Nas, N, NProcs

     itype = NBAAType(i)

     dx1 = R(icx1,i)
     dx2 = R(icx2,i)

     dR2 = dx1*dx1 + dx2*dx2

     if(dR2 > Rcylmax2(itype)) cycle

     Rm = Rcylmin2(itype)

     U  =   InterPolateCylinder(dR2,Rm,invdr_size(itype),1,itype)
     Fc = - InterPolateCylinder(dR2,Rm,invdr_size(itype),2,itype)
     Fr =   InterPolateCylinder(dR2,Rm,invdr_size(itype),3,itype)

     F1 = Fc * dx1
     F2 = Fc * dx2

     Frc_OptC(icx1,i) = Frc_OptC(icx1,i) + F1
     Frc_OptC(icx2,i) = Frc_OptC(icx2,i) + F2

     Ene_OptC = Ene_OptC + U

     FF = FF + Fr

   end do

   FCylRs = FF

end subroutine Force_AAcylinder


!######################################################################
!######################################################################


Function InterPolateCylinder(x,drmin,invdr,j,k) Result(rs)
! j : 1 (potential), 2 (force), 3 (force on radius)
! k : nonbond type
use CylParam, only : TabCyl

implicit none

integer :: ii, j, k
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

   y1  = TabCyl( ii  , j, k )
   y2  = TabCyl( ii+1, j, k )
   y3  = TabCyl( ii+2, j, k )
   y4  = TabCyl( ii+3, j, k )

   z3  = y2 - y3 + sb * (y4 - y1)
   z2  = y1 - 2.d0 * y2 + y3
   z1  = -sb2 * y1 - y2 + 2.d0 * y3 - sb * y4

   rs = 0.5d0 * dx * ( z1 + dx * ( z2 + dx * z3 ) ) + y2

   else

   rs = TabCyl( 1, j, k )
   print *, 'Cylinder Bad contact'

   end if

#else

   xx  = (x - drmin) * invdr

   ii  = nint( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = TabCyl( ii  , j, k )
   y2  = TabCyl( ii+1, j, k )
   y3  = TabCyl( ii+2, j, k )

   z2  = y1 - 2.d0 * y2 + y3
   z1  = - y1 + y3
   rs  = 0.5d0 * (z1 + z2 * dx ) * dx + y2

   else

   rs = TabCyl( 1, j, k )
   print *, 'Cylinder Bad contact'

   end if

#endif

end Function InterPolateCylinder


!######################################################################
!######################################################################


subroutine CylPMFsample(istep)

use CylParam
use TimeParam, only : Nstep, lk
use UnitExParam, only : cvol
use CommonBlocks, only : Job_name, QMaster, QFSCyl

implicit none

integer :: istep, nsample
real(8) :: PMFtemp, PMFtempH
character(len=80) :: Filename

   if(istep==lk.and.QMaster) then
     PMFcyl = 0.d0
     if(QFSCyl) PMFcylH = 0.d0
     write(Filename,'(a,a)') trim(adjustl(Job_name)),'_PMF_Cylinder.dat'
     open(76,file=trim(Filename))
   end if

   nsample = istep/lk
   PMFcyl = PMFcyl + FCylRs
   if(QFSCyl) PMFcylH = PMFcylH + FCylH

   if((mod(nsample,100)==0).or.(istep==Nstep)) then
     PMFtemp = PMFcyl / nsample
     call SumPMF_Cyl(PMFtemp)
     if(QFSCyl) then
       PMFtempH = PMFcylH / nsample
       call SumPMF_Cyl(PMFtempH)
       if(QMaster) write(76,'(i10,2e16.8)') nsample, PMFtemp * cvol, PMFtempH * cvol
     else
       if(QMaster) write(76,'(i10,e16.8)') nsample, PMFtemp * cvol
     end if
   end if

   if(istep==Nstep.and.QMaster) close(76)

end subroutine CylPMFsample


!######################################################################
!######################################################################


subroutine Force_FScylinder

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use CylParam, only : icaxis, ngrid_cyl, ngrid_cylZ, icx1, icx2, &
&   Rcylmax2, Zcylmax, TabCylZ, FscylR, FscylZ, FCylRs, FCylH
use CGdata, only : NBAtomType
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, j, j1, k1, Nas, itype
real(8) :: FF1, FF2, F1, F2
real(8) :: dx1, dx2, dx3, dR2, dZ, dR
real(8) :: dRx1, dRx2, InvdR
real(8) :: F3, rslt
real(8) :: x1l, x1u, x2l, x2u
real(8), dimension(4) :: d0,d1,d2,d12
real(8), dimension(5) :: yf
real(8), dimension(ngrid_cyl ) :: pxR
real(8), dimension(ngrid_cylZ) :: pxZ

   Nas = NProcs - MyRank

   FF1 = 0.d0
   FF2 = 0.d0

   do i = Nas, N, NProcs

     itype = NBAtomType(i)

     dx1 = R(icx1,i)
     dx2 = R(icx2,i)
     dx3 = R(icaxis,i)

     dR2 = dx1*dx1 + dx2*dx2
     dZ  = abs(dx3)

     if((dZ > Zcylmax(itype)).or.(dR2 > Rcylmax2(itype))) cycle

     pxR(:) = FscylR(:,itype)
     pxZ(:) = FscylZ(:,itype)

     dR = sqrt(dR2)

     call locate(pxR,ngrid_cyl, dR,j1)
     call locate(pxZ,ngrid_cylZ,dZ,k1)

     x1l = pxR(j1)
     x1u = pxR(j1+1)
     x2l = pxZ(k1)
     x2u = pxZ(k1+1)

     do j = 1, 3

       d0 (1) = TabCylZ(1,j1  ,k1  ,j,itype)
       d1 (1) = TabCylZ(2,j1  ,k1  ,j,itype)
       d2 (1) = TabCylZ(3,j1  ,k1  ,j,itype)
       d12(1) = TabCylZ(4,j1  ,k1  ,j,itype)
       d0 (2) = TabCylZ(1,j1+1,k1  ,j,itype)
       d1 (2) = TabCylZ(2,j1+1,k1  ,j,itype)
       d2 (2) = TabCylZ(3,j1+1,k1  ,j,itype)
       d12(2) = TabCylZ(4,j1+1,k1  ,j,itype)
       d0 (4) = TabCylZ(1,j1  ,k1+1,j,itype)
       d1 (4) = TabCylZ(2,j1  ,k1+1,j,itype)
       d2 (4) = TabCylZ(3,j1  ,k1+1,j,itype)
       d12(4) = TabCylZ(4,j1  ,k1+1,j,itype)
       d0 (3) = TabCylZ(1,j1+1,k1+1,j,itype)
       d1 (3) = TabCylZ(2,j1+1,k1+1,j,itype)
       d2 (3) = TabCylZ(3,j1+1,k1+1,j,itype)
       d12(3) = TabCylZ(4,j1+1,k1+1,j,itype)

       call bicubic(d0,d1,d2,d12,x1l,x1u,x2l,x2u,dR2,dZ,rslt)
       yf(j) = rslt

     end do

     InvdR = 1.d0/dR
     dRx1 = dx1*InvdR
     dRx2 = dx2*InvdR

     F1 = - yf(2) * dRx1
     F2 = - yf(2) * dRx2
     F3 = - yf(3)
     if(dx3 < 0.) F3 = - F3

     Frc_OptC(icx1  ,i) = Frc_OptC(icx1  ,i) + F1
     Frc_OptC(icx2  ,i) = Frc_OptC(icx2  ,i) + F2
     Frc_OptC(icaxis,i) = Frc_OptC(icaxis,i) + F3

     Ene_OptC = Ene_OptC + yf(1)

!     FF1 = FF1 + yf(4)
!     FF2 = FF2 + yf(5)

   end do

!   FCylRs = FF1
!   FCylH  = FF2

end subroutine Force_FScylinder


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! *******************************

subroutine Force_FixonPlane

use OptConstraintParam, only : NumPLC, PLCI, kPLC, rPLC, IaxisPLC
use Configuration, only : R
use CommonMPI
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit NONE

integer :: i, j, Nas
real(8) :: dR, fc

   Nas = NProcs - MyRank

   do j = Nas, NumPLC, NProcs

     i = PLCI(j)

     dR = R(IaxisPLC,i) - rPLC(j)

     fc = kPLC(j) * dR

     Frc_OptC(IaxisPLC,i) = Frc_OptC(IaxisPLC,i) - 2.d0 * fc
     Ene_OptC = Ene_OptC + fc * dR

   end do

end subroutine Force_FixonPlane


!######################################################################
!######################################################################


! *******************************
! ** Optional Constraint Force **
! ** Positional Constraint by  **
! ** Harmonic potential        **
! *******************************

subroutine Force_HamC

use Numbers, only : N
use CommonBlocks, only : QPathInt
use Configuration, only : R
use CommonMPI
use CommonPI
use OptConstraintParam, only : NHam, kHam, HamR, Rrot, &
&   Frc_OptC, Vir_OptC, Ene_OptC

implicit NONE

integer :: i, Nas
real(8) :: R2, R1, dR, Fc
real(8), dimension(3) :: Rij, Fij
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   if( NHam == 0 ) Return

   Nas = NProcsTemp - MyRankTemp

   do i = Nas , N, NProcsTemp

   if(.not.HamR(i)) cycle

     Rij = R(:,i) - Rrot(:,i)
     R2  = dot_product(Rij,Rij)

     R1 = sqrt(R2)

     if(R1==0.) cycle

     dR = R1
     Fc = kHam * dR

     Ene_OptC = Ene_OptC + Fc * dR

     Fij = - 2.d0 * Fc * Rij / R1
     Frc_OptC(:,i) = Frc_OptC(:,i) + Fij

     Vir_OptC(1,1) = Vir_OptC(1,1) + Fij(1) * R(1,i)
     Vir_OptC(1,2) = Vir_OptC(1,2) + Fij(1) * R(2,i)
     Vir_OptC(1,3) = Vir_OptC(1,3) + Fij(1) * R(3,i)
     Vir_OptC(2,2) = Vir_OptC(2,2) + Fij(2) * R(2,i)
     Vir_OptC(2,3) = Vir_OptC(2,3) + Fij(2) * R(3,i)
     Vir_OptC(3,3) = Vir_OptC(3,3) + Fij(3) * R(3,i)

   end do

   Vir_OptC(2,1) = Vir_OptC(1,2)
   Vir_OptC(3,1) = Vir_OptC(1,3)
   Vir_OptC(3,2) = Vir_OptC(2,3)

end subroutine Force_HamC


!#####################################################################
!#####################################################################


subroutine Force_ConstCOM

use OptConstraintParam, only : Nccom, kccom, Idcomf, Idcomt, Rcon, &
&   Xcon, fcom, InvMscom
use AtomParam, only : Mass
use Configuration, only : R
use OptConstraintParam, only : Frc_OptC, Ene_OptC

implicit none

integer :: i, j
real(8), dimension(3) :: dR, Fc

   do i = 1, Nccom
     Rcon(:,i) = 0.d0
     do j = Idcomf(i), Idcomt(i)
       Rcon(:,i) = Rcon(:,i) + Mass(j) * R(:,j)
     end do
     Rcon(:,i) = Rcon(:,i) * InvMscom(i)
     dR(:) = Rcon(:,i) - Xcon(:,i)
     fcom(:,i) = - kccom * dR(:)
     Fc(:) = fcom(:,i) * InvMscom(i)
     do j = Idcomf(i), Idcomt(i)
       Frc_OptC(:,j) = Frc_OptC(:,j) + Fc(:) * Mass(j)
     end do
     Ene_OptC = Ene_OptC + 0.5d0 * kccom * dot_product( dR, dR )
   end do

end subroutine Force_ConstCOM


!#####################################################################
!#####################################################################


subroutine Force_Efield

use Conduct, only : Eforce
use EwaldParam, only : Nel, Nelist, PCh
use OptConstraintParam, only : Frc_OptC

implicit none

integer :: i, j

   do i = 1, Nel
     j = Nelist(i)
     Frc_OptC(:,j) = Frc_OptC(:,j) + Pch(i) * Eforce(:)
   end do

end subroutine Force_Efield


!#####################################################################
!#####################################################################


subroutine Force_Constraint

use Numbers, only : N
use Configuration, only : R
use SHAKEparam

implicit NONE

integer :: i , j , k , ii , jj , kk
real(8), dimension(3) :: dR

! ## Clear
   Frc_Const = 0.d0
   Vir_Const = 0.d0

   do k = 1 , NSHAKEGroup

     do kk = 1 , NCoupleBond(k)

       ii = CouplePair(k,kk,1)
       jj = CouplePair(k,kk,2)

       i  = CoupleAtom(k,ii)
       j  = CoupleAtom(k,jj)

       dR = R(:,i) - R(:,j)

       Frc_Const(:,i) = Frc_Const(:,i) + Lagmultip(k,kk) * dR
       Frc_Const(:,j) = Frc_Const(:,j) - Lagmultip(k,kk) * dR

     end do

   end do

   do i = 1, N

     Vir_Const(1,1) = Vir_Const(1,1) + Frc_Const(1,i) * R(1,i)
     Vir_Const(1,2) = Vir_Const(1,2) + Frc_Const(1,i) * R(2,i)
     Vir_Const(1,3) = Vir_Const(1,3) + Frc_Const(1,i) * R(3,i)
     Vir_Const(2,2) = Vir_Const(2,2) + Frc_Const(2,i) * R(2,i)
     Vir_Const(2,3) = Vir_Const(2,3) + Frc_Const(2,i) * R(3,i)
     Vir_Const(3,3) = Vir_Const(3,3) + Frc_Const(3,i) * R(3,i)

   end do

   Vir_Const(2,1) = Vir_Const(1,2)
   Vir_Const(3,1) = Vir_Const(1,3)
   Vir_Const(3,2) = Vir_Const(2,3)

end subroutine Force_Constraint


!######################################################################
!######################################################################


! ***************************************************
! **  Force and Torque on the COM of the molecule  **
! ***************************************************

subroutine Force_Div_Component(F,FCOM,TORQ)

use Numbers, only : N
use RBparam, only : NumRB, QSingle, RBType, NumRBAtom, Rmolec, Rotation

implicit NONE

integer :: i, j, k, ii, Nc, MyType
real(8), dimension(3)       :: TOR
real(8), dimension(3,NumRB) :: FCOM, TORQ
real(8), dimension(3,N) :: F

   FCOM = 0.d0
   TORQ = 0.d0

   k = 0

   do i = 1 , NumRB

     if(QSingle(i)) then

       k = k + 1
       FCOM(:,i) = FCOM(:,i) + F(:,k)

     else 

       MyType = RBType(i)
       Nc = NumRBAtom(MyType)

       TOR = 0.d0

       do j = 1 , Nc

         ii = j + k

         FCOM(:,i) = FCOM(:,i) + F(:,ii)

         TOR(1) = TOR(1) + Rmolec(2,j,i)*F(3,ii) - Rmolec(3,j,i)*F(2,ii)
         TOR(2) = TOR(2) + Rmolec(3,j,i)*F(1,ii) - Rmolec(1,j,i)*F(3,ii)
         TOR(3) = TOR(3) + Rmolec(1,j,i)*F(2,ii) - Rmolec(2,j,i)*F(1,ii)

       end do

! space fixed ---> body fixed
       TORQ(1,i) = Rotation(1,1,i) * TOR(1) + Rotation(1,2,i) * TOR(2)   &
          &                                 + Rotation(1,3,i) * TOR(3)
       TORQ(2,i) = Rotation(2,1,i) * TOR(1) + Rotation(2,2,i) * TOR(2)   &
          &                                 + Rotation(2,3,i) * TOR(3)
       TORQ(3,i) = Rotation(3,1,i) * TOR(1) + Rotation(3,2,i) * TOR(2)   &
          &                                 + Rotation(3,3,i) * TOR(3)

       k = k + Nc

     end if

   end do

end subroutine Force_Div_Component


!######################################################################
!######################################################################


subroutine SPCF_Intra

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QPathInt
use Configuration, only : R
use CommonPI
use CommonMPI
use UnitExParam, only : reng, pi
use OptConstraintParam, only : Frc_OptC, Ene_OptC
use AtomParam, only : ResidName

implicit NONE

integer :: i, j, Nas, ii, jj, isp
integer :: iOO, iH1, iH2
real(8), parameter :: ScaleParam = 1.d-18 * reng
real(8), parameter :: a_const =  2.361d0*2.361d0*0.708d0*ScaleParam !3.947d0 * ScaleParam
real(8), parameter :: b_const =  1.803d0 * ScaleParam * 0.5d0
real(8), parameter :: c_const = -1.469d0 * ScaleParam
real(8), parameter :: d_const =  0.776d0 * ScaleParam
real(8), parameter :: b_OHeq = 1.d0
real(8) :: b_HHeq
real(8), dimension(3) :: dOH1, dOH2, dHH0
real(8), dimension(3) :: sOH1, sOH2, sHH0
real(8) :: R2OH1, R2OH2, R2HH0, R1OH1, R1OH2, R1HH0
real(8) :: dR_OH1, dR_OH2, dR_HH0
real(8) :: fca1, fca2, fcb0, fcc1, fcc2, fcd1, fcd2
real(8) :: fcoh1, fcoh2, fchh0
integer :: NProcsTemp, MyRankTemp

   b_HHeq = b_OHeq * sin(pi*54.d0/180.d0) * 2.d0

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   jj = 0

   do i = 1, NumSpec

     if(ResidName(jj+1) == 'SPCF') then

       isp = i
       exit

     end if

     jj = jj + NumMol(i) * NumAtm(i)

   end do

   do j = Nas, NumMol(isp), NProcsTemp

     ii = (j-1) * 3 + jj

     iOO = ii + 1
     iH1 = ii + 2
     iH2 = ii + 3

     dOH1(:) = R(:,iOO) - R(:,iH1)
     dOH2(:) = R(:,iOO) - R(:,iH2)
     dHH0(:) = R(:,iH1) - R(:,iH2)

#ifdef PCC
     R2OH1 = dOH1(1) * dOH1(1) + dOH1(2) * dOH1(2) + dOH1(3) * dOH1(3)
     R2OH2 = dOH2(1) * dOH2(1) + dOH2(2) * dOH2(2) + dOH2(3) * dOH2(3)
     R2HH0 = dHH0(1) * dHH0(1) + dHH0(2) * dHH0(2) + dHH0(3) * dHH0(3)
#else
     R2OH1 = dot_product(dOH1, dOH1)
     R2OH2 = dot_product(dOH2, dOH2)
     R2HH0 = dot_product(dHH0, dHH0)
#endif
     R1OH1 = sqrt( R2OH1 )
     R1OH2 = sqrt( R2OH2 )
     R1HH0 = sqrt( R2HH0 )

     dR_OH1 = R1OH1 - b_OHeq
     dR_OH2 = R1OH2 - b_OHeq
     dR_HH0 = R1HH0 - b_HHeq

     sOH1(:) = dOH1(:) / R1OH1
     sOH2(:) = dOH2(:) / R1OH2
     sHH0(:) = dHH0(:) / R1HH0

     Ene_OptC = Ene_OptC + a_const * ( dR_OH1 * dR_OH1 + dR_OH2 * dR_OH2 ) &
     &                   + b_const * dR_HH0 * dR_HH0                       &
     &                   + c_const * ( dR_OH1 + dR_OH2 ) * dR_HH0          &
     &                   + d_const * dR_OH1 * dR_OH2

     fca1 = 2.d0 * a_const * dR_OH1
     fca2 = 2.d0 * a_const * dR_OH2

     fcb0 = 2.d0 * b_const * dR_HH0

     fcc1 = c_const * dR_HH0
     fcc2 = c_const * ( dR_OH1 + dR_OH2 )

     fcd1 = d_const * dR_OH2
     fcd2 = d_const * dR_OH1

     fcoh1 = fca1 + fcc1 + fcd1
     fcoh2 = fca2 + fcc1 + fcd2
     fchh0 = fcb0 + fcc2

     Frc_OptC(:,iOO) = - fcoh1 * sOH1(:) - fcoh2 * sOH2(:)

     Frc_OptC(:,iH1) =   fcoh1 * sOH1(:) - fchh0 * sHH0(:)

     Frc_OptC(:,iH2) =   fcoh2 * sOH2(:) + fchh0 * sHH0(:)

   end do

end subroutine SPCF_Intra


!######################################################################
!######################################################################


subroutine Set_MembrPot

use Numbers, only : N
use CommonBlocks, only: ForceField, QMaster
use AtomParam, only : AtomName, ResidName
use IOparam, only : Extmem_file
use MembrParam
use UnitExParam, only : ExParam

implicit none

integer :: i, ii, j, NumbParam
integer, dimension(:), allocatable :: MembTypeParam
real(8), dimension(:), allocatable :: dEzParam, AmemParam, BmemParam
character(len=4), dimension(:), allocatable :: RNameParam, ANameParam
character(len=80) :: String, String1

! ## Reading the parameter file

   open(99,file=trim(Extmem_file),status='old')

   do
     read(99,'(a)') String1
     String = trim(adjustl(String1))
     if((String(1:1) == '#').or.(String(1:1) == '!')) cycle
     read(String,*) NumbParam
     exit
   end do

   allocate(MembTypeParam(NumbParam))
   allocate(dEzParam(NumbParam))
   allocate(AmemParam(NumbParam))
   allocate(BmemParam(NumbParam))
   allocate(RNameParam(NumbParam))
   allocate(ANameParam(NumbParam))

   i = 0
   do
     read(99,'(a)') String1
     String = trim(adjustl(String1))
     if((String(1:1) == '#').or.(String(1:1) == '!')) cycle
     i = i + 1
     read(String,*) RNameParam(i),ANameParam(i),MembTypeParam(i),dEzParam(i), &
     &              AmemParam(i),BmemParam(i)
     if(i==NumbParam) exit
   end do

   close(99)

! ## assign the parameters

   allocate(Memptype(N))
   allocate(dEz(N))
   allocate(amem(N))
   allocate(bmem(N))

   Memptype(:) = 0

   if(ForceField(1:5)/='CHARM') then
     if(QMaster) write(*,*) 'error : EXP_MEMP is valid only for CHARMM force field'
     call Finalize
   end if

   do i = 1, N
inn: do j = 1, NumbParam
       if(ResidName(i) == RNameParam(j)) then
         if(AtomName(i) == ANameParam(j)) then
           Memptype(i) = MembTypeParam(j)
           dEz(i)  = dEzParam(j) * ExParam
           Amem(i) = 1.d0 / AmemParam(j)
           Bmem(i) = BmemParam(j)
           exit inn
         end if
       end if
     end do inn
   end do

   if(QMaster) then
     do ii = 6, 11, 5
     write(ii,*)
     write(ii,*) 'The membrane potential of mean force will exert on the following atoms'
     write(ii,*)
     write(ii,'(4x,a)') '----------------------------'
     do i = 1, N
     if(Memptype(i)/=0) then
       write(ii,'(4x,a,a4,a,a4)') 'Residue = ', ResidName(i),' : Atom = ',AtomName(i)
     end if
     end do
     write(ii,'(4x,a)') '----------------------------'
     write(ii,*)
     end do
   end if

end subroutine Set_MembrPot


!######################################################################
!######################################################################


subroutine Force_MembrPot

use Numbers, only : N
use CommonBlocks, only: QPathInt
use OptConstraintParam, only : Frc_OptC, Ene_OptC
use Configuration, only : R
use MembrParam
use CommonMPI, only: MyRank, NProcs
use CommonPI, only: MyRankPI, NumProcess
use UnitExParam, only : pi

implicit none

integer :: i, Nas
integer :: NProcsTemp, MyRankTemp
real(8) :: ppi2, x1, xx, pp, fp, aa, ex, pref, zz

   ppi2 = 1.d0 / sqrt( 2.d0 * pi )

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   do i = Nas, N, NProcsTemp

     if(Memptype(i)==0) cycle

     if(Memptype(i)==1) then

       if(R(3,i) >= 0.) then
         zz   =   R(3,i)
         pref =   1.d0
       else
         zz   = - R(3,i)
         pref = - 1.d0
       end if

       x1 = zz * Amem(i)
       xx = 1.d0 / (1.d0 + (x1)**Bmem(i))
       pp = dEz(i) * xx

       Ene_OptC = Ene_OptC + pp

       fp = pref * pp * xx * Amem(i) * Bmem(i) * (x1)**(Bmem(i)-1.d0)

       Frc_OptC(3,i) = Frc_OptC(3,i) + fp

     else if(Memptype(i)==2) then

       if(R(3,i) >= 0.) then
         zz   =   R(3,i)
         pref =   1.d0
       else
         zz   = - R(3,i)
         pref = - 1.d0
       end if

       zz = abs( R(3,i) )
       aa = Amem(i) * Amem(i)
       x1 = zz - Bmem(i)
       xx = - 0.5d0 * x1 * x1 * aa
       ex = exp(xx)
       pp = dEz(i) * Amem(i) * ppi2 * ex

       Ene_OptC = Ene_OptC + pp

       fp = pref * pp * aa * x1

       Frc_OptC(3,i) = Frc_OptC(3,i) + fp

     end if

   end do

end subroutine Force_MembrPot
