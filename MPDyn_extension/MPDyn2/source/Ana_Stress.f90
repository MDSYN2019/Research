

! ###########################################################
! ###########################################################


subroutine Ana_Stress

use Configuration, only : R
use ParamAnalyze, only : NJobs, NtrjStep, QEW, Nbin, dZbin, dRbin, InvdRbin, Zbin, &
&   Ibin, VirProXY, VirProZZ, VirProN, VirTProN, VirTProT, VirProT, Qslab, NComp, &
&   QLC, RcutC2
use CellParam, only : H, InvH, Volume, CellShape
use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster, QCoulomb
use UnitExParam, only : rprs, pi
use AtomParam, only : Mass
use CGdata, only : Rcut_MAX

implicit none

integer :: i, j, k, Nframe, Numa, Numb, ii
real(8) :: Hzz, Hzh, det, AvZbin, AvV, R1
real(8) :: Lz, zz, Pxy, Pzz, RRmax, InvMasG, xmin, x
real(8) :: Rad, vv, PreV, vol, Pn, Pt, R2, InvdZbin
real(8), dimension(3) :: Rg
real(8), dimension(:), allocatable :: VolRad
external det

   if(QMaster) open(1,file='./Analy/StressProfile.dat',status='unknown')

   CellShape = 2
   QLC = .False.

   if(QEW) then
     call ErrorFuncList           !  make a table of error-function
     call RecLatticeList          !  define reciprocal lattice
   end if

   allocate( Zbin(N), Ibin(N) )
   AvZbin = 0.d0
   AvV    = 0.d0

   if(Qslab) then
     allocate( VirProXY(Nbin) )
     allocate( VirProZZ(Nbin) )
     VirProXY(:) = 0.d0
     VirProZZ(:) = 0.d0
   else
     allocate( VirProN(Nbin) )
     allocate( VirProT(Nbin) )
     allocate( VirTProN(Nbin) )
     allocate( VirTProT(Nbin) )
     allocate( VolRad(Nbin) )
     VirProN(:) = 0.d0
     VirProT(:) = 0.d0

     PreV = 0.d0
     do i = 1, Nbin
       Rad = i * dRbin
       vv  = 4.d0/3.d0*pi*Rad*Rad*Rad
       vol = vv - PreV
       VolRad(i) = 1.d0 / vol
       PreV = vv
     end do
   end if

   if(.not.Qslab) then
     RRmax = dRbin*Nbin
     ii = 0
     Numb = 0
     do k = 1, NumSpec
       if(k==NComp) then
         Numa = ii+1
         Numb = ii+NumMol(k)*NumAtm(k)
         exit
       end if
       ii = ii + NumMol(k)*NumAtm(k)
     end do
     InvMasG = 0.d0
     if( Numb /= 0 ) then
       do k = Numa, Numb
         InvMasG = InvMasG + Mass(k)
       end do
       InvMasG = 1.d0 / InvMasG
     end if
   end if

   Nframe = 0
   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif

       Nframe = Nframe + 1

!     ---------------
       call BcastRH
!     ---------------

       if(Nframe == 1) then
         call Set_CGcond ! Cell List
         R2 = Rcut_MAX*Rcut_MAX
         if(QCoulomb.and.(.not.QEW).and.(R2<RcutC2)) QLC = .True.
       else
         call CheckCellList
       end if

       call InversMatrix(H,InvH)

       Volume = det(H)

       if(Qslab) then

         AvV = AvV + Volume

         Hzz = H(3,3)
         dZbin = Hzz / dble(Nbin)
         InvdZbin = 1.d0 / dZbin
         Hzh = Hzz * 0.5d0

         do k = 1, N
           R(3,k) = R(3,k) + Hzh
         end do

         do k = 1, N
           if(R(3,k) > Hzz ) then
             R(3,k) = R(3,k) - Hzz
           else if(R(3,k) < 0.d0) then
             R(3,k) = R(3,k) + Hzz
           end if
         end do

         do k = 1, N
           Ibin(k) = int(R(3,k)*InvdZbin) + 1
         end do

         do k = 1, N
           Zbin(k) = R(3,k) - (Ibin(k)-1) * dZbin
         end do

         AvZbin = AvZbin + dZbin

       else

         xmin = min(H(1,1),H(2,2),H(3,3))
         if(xmin<RRmax) then
           if(QMaster) write(*,*) 'error: too small cell length'
           call Finalize
         end if

         if( Numb /= 0 ) then
           Rg(:) = 0.d0
           do k = Numa, Numb
             Rg(:) = Rg(:) + Mass(k)*R(:,k)
           end do
           Rg(:) = Rg(:) * InvMasG
           do k = 1, N
             R(:,k) = R(:,k) - Rg(:)
           end do
           call PBC
           Rg(:) = 0.d0
           do k = Numa, Numb
             Rg(:) = Rg(:) + Mass(k)*R(:,k)
           end do
           Rg(:) = Rg(:) * InvMasG
           do k = 1, N
             R(:,k) = R(:,k) - Rg(:)
           end do
         end if

         do k = 1, N
           R(1,k) = R(1,k) - nint(InvH(1,1)*R(1,k))*H(1,1)
           R(2,k) = R(2,k) - nint(InvH(2,2)*R(2,k))*H(2,2)
           R(3,k) = R(3,k) - nint(InvH(3,3)*R(3,k))*H(3,3)
         end do

         do k = 1, N
           R1 = sqrt(R(1,k)*R(1,k)+R(2,k)*R(2,k)+R(3,k)*R(3,k))
           Zbin(k) = R1
           Ibin(k) = int(R1*InvdRbin) + 1
         end do

         VirTProN(:) = 0.d0
         VirTProT(:) = 0.d0

       end if

       call TransCellList

       call Cal_StressProfile

       if(.not.Qslab) then
         VirProN(:) = VirProN(:) + VirTProN(:) * VolRad(:)
         VirProT(:) = VirProT(:) + VirTProT(:) * VolRad(:)
       end if

     end do

   end do

   if(Qslab) then

     call Sum_StressProf

     if(QMaster) then
       AvZbin = AvZbin / dble(Nframe)
       AvV    = dble(Nbin) / AvV

       do i = 1, Nbin
         VirProXY(i) = VirProXY(i) * AvV
         VirProZZ(i) = VirProZZ(i) * AvV
       end do
       Lz = AvZbin * Nbin

       x = 1.d-5 / rprs

       do i = 1, Nbin
         zz = - Lz*0.5d0 - AvZbin*0.5d0 + AvZbin*i
         Pxy = (VirProXY(i)) * x * 0.5d0
         Pzz = VirProZZ(i) * x
         write(1,'(f12.6,3e16.8)') zz, Pxy, Pzz, Pxy-Pzz ! [A] [bar]
       end do
     end if

   else

     call Sum_StressProf2

     if(QMaster) then
       do i = 1, Nbin
         VirProN(i) = VirProN(i) / dble(Nframe)
         VirProT(i) = VirProT(i) / dble(Nframe)
       end do

       x = 1.d-5 / rprs

       do i = 1, Nbin
         Rad = (dble(i)-0.5d0)*dRbin
         Pt  = VirProT(i) * x
         Pn  = VirProN(i) * x
         write(1,'(f12.6,3e16.8)') Rad, Pt, Pn, Pt-Pn
       end do
     end if

   end if

   if(QMaster) close(1)

end subroutine Ana_Stress


! ###########################################################
! ###########################################################


subroutine Cal_StressProfile

use CommonBlocks, only : QCoulomb
use ParamAnalyze, only : QEW
use CellParam, only : H

implicit none

real(8) :: clhx, clhy, clhz

   clhx = H(1,1) * 0.5d0
   clhy = H(2,2) * 0.5d0
   clhz = H(3,3) * 0.5d0

   call Stress_Bond_CG(clhx,clhy,clhz)
   call Stress_Angle_CG(clhx,clhy,clhz)
   call Stress_vdW_CG(clhx,clhy,clhz)
   if(QCoulomb.and.QEW) call Stress_Reciprocal

end subroutine Cal_StressProfile


! ###########################################################
! ###########################################################


subroutine Stress_vdW_CG(clhx,clhy,clhz)

use CommonBlocks, only : QCellList
use ParamAnalyze, only : QLC

implicit none

real(8) :: clhx, clhy, clhz

   if(QCellList) then
     call Stress_CELL
     if(QLC) call Stress_LC(clhx,clhy,clhz)
   else
     call Stress_NoCELL(clhx,clhy,clhz)
   end if

end subroutine Stress_vdW_CG


!######################################################################
!######################################################################


subroutine Stress_CELL

use Numbers, only : N
use Configuration, only : R
use CommonMPI, only : MyRank, NProcs
use CellListMethod
use NoLJParam, only : NumNoLJ
use NonbondParam, only : Charge
use CellParam, only : CellShft, CellL, InvCL, H, InvH
use CommonBlocks, only : Qcoulomb, QMaster
use CGdata, only : NBAtomType, Rcut2, NBFuncType, CoefAtype, CoefBtype
use ParamAnalyze, only : QEW, dZbin, Zbin, Ibin, RcutC2, Qslab
use EwaldParam, only : Alpha, msh, EFList, ar2

implicit none

integer :: i, j, k, itype, jtype, Nas, jj, ii
integer :: icell, jcell, jcell0, nabor, ftype
integer :: non, check_bonded, Nz
integer :: IbinI, IbinJ
real(8) :: dZI, dZJ
real(8) :: Rix, Riy, Riz
real(8) :: Rx, Ry, Rz
real(8) :: Risx, Risy, Risz
real(8) :: Fx, Fy, Fz
real(8) :: R2, cf
real(8) :: aij, bij, R1
real(8) :: InvR2, InvR1, InvR3, InvR4, term1, term2, fk, fke
real(8) :: x, dx, yy, y1, y2, y3
real(8) :: z1, z2, xtm, fk1, fk2, ErrorFunc
real(8) :: virxy, virzz
external check_bonded

   Nas = NProcs - MyRank

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

           ftype = NBFuncType(itype,jtype)

           if(ftype==7) then ! LJ12-4

             aij = CoefAtype(itype,jtype)
             bij = CoefBtype(itype,jtype)
             InvR2 = 1.d0 / R2
             InvR4 = InvR2 * InvR2
             term1 = aij * InvR4 * InvR4 * InvR4
             term2 = bij * InvR4
             fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2

           else if(ftype==1) then ! LJ9-6

             aij = CoefAtype(itype,jtype)
             bij = CoefBtype(itype,jtype)
             InvR2 = 1.d0 / R2
             InvR1 = sqrt(InvR2)
             InvR3 = InvR2 * InvR1
             term1 = aij * InvR3 * InvR3 * InvR3
             term2 = bij * InvR3 * InvR3
             fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2

           end if

           if(Qcoulomb) then

             cf = Charge(i) * Charge(j)

             if(cf/=0.) then

               if(QEW) then

                 R1  = sqrt(R2)
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

                 fke  = ( fk1 + fk2 ) * InvR2
                 fk = fk + fke

               else

                 R1  = sqrt(R2)
                 fk1 = cf * R1 * InvR2
                 fke  = fk1 * InvR2
                 fk = fk + fke

               end if

             end if

           end if

           if(Qslab) then

             Fx  = fk * Rx
             Fy  = fk * Ry
             Fz  = fk * Rz

             virxy = Fx * Rx + Fy * Ry
             virzz = Fz * Rz
             IbinI = Ibin(i)
             IbinJ = Ibin(j)

             dZI = Zbin(i)
             dZJ = Zbin(j)
             call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,0)

           else

             call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,j)

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

           if(nabor<=4) then
             Nz = 0
           else if(nabor<=8) then
             Nz = 1
           else
             Nz = -1
           end if

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

               ftype = NBFuncType(itype,jtype)

               if(ftype==7) then ! LJ12-4

                 aij = CoefAtype(itype,jtype)
                 bij = CoefBtype(itype,jtype)
                 InvR2 = 1.d0 / R2
                 InvR4 = InvR2 * InvR2
                 term1 = aij * InvR4 * InvR4 * InvR4
                 term2 = bij * InvR4
                 fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2

               else if(ftype==1) then ! LJ9-6

                 aij = CoefAtype(itype,jtype)
                 bij = CoefBtype(itype,jtype)
                 InvR2 = 1.d0 / R2
                 InvR1 = sqrt(InvR2)
                 InvR3 = InvR2 * InvR1
                 term1 = aij * InvR3 * InvR3 * InvR3
                 term2 = bij * InvR3 * InvR3
                 fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2

               end if

               if(Qcoulomb) then

                 cf = Charge(i) * Charge(j)

                 if(cf/=0.) then

                   if(QEW) then

                     R1  = sqrt(R2)
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

                     fke  = ( fk1 + fk2 ) * InvR2
                     fk = fk + fke

                   else

                     R1  = sqrt(R2)
                     fk1 = cf * R1 * InvR2
                     fke  = fk1 * InvR2
                     fk = fk + fke

                   end if

                 end if

               end if

               if(Qslab) then

                 Fx = fk * Rx
                 Fy = fk * Ry
                 Fz = fk * Rz

                 virxy = Fx * Rx + Fy * Ry
                 virzz = Fz * Rz
                 IbinI = Ibin(i)
                 IbinJ = Ibin(j)

                 dZI = Zbin(i)
                 dZJ = Zbin(j)
                 call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,Nz)

               else

                 call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,j)

               end if

             end if

             j = NextP(j)

           end do

         end if

       end do

       i = NextP(i)

     end do

   end do

end subroutine Stress_CELL


! ###########################################################
! ###########################################################


subroutine Stress_LC(clhx,clhy,clhz)

use Configuration, only : R
use Numbers, only : N
use CellParam, only : CellShft
use ParamAnalyze, only : dZbin, Zbin, Ibin, RcutC2, Qslab
use CommonMPI
use CGdata, only : Rcut_MAX
use EwaldParam, only : Nel, Nelist, PCh

implicit none

integer :: i, j, k, Nas
integer :: Nx, Ny, Nz, IbinI, IbinJ
real(8) :: clhx, clhy, clhz, dZI, dZJ
real(8) :: Rix, Riy, Riz, Rx, Ry, Rz, R1, R2
real(8) :: InvR2
real(8) :: fk, cf, fk1
real(8) :: Fx, Fy, Fz
real(8) :: virxy, virzz
real(8) :: R2_MIN
integer, dimension(:), allocatable :: Kbin
real(8), dimension(:,:), allocatable :: Rel

   Nas = NProcs - MyRank

   allocate( Rel(3,Nel) )
   allocate( Kbin(Nel) )

   R2_MIN = Rcut_MAX*Rcut_MAX

   do i = 1, Nel
     j = Nelist(i)
     Rel(1,i) = R(1,j)
     Rel(2,i) = R(2,j)
     Rel(3,i) = R(3,j)
     Kbin(i) = Ibin(j)
   end do

   do i = Nas, Nel, NProcs

     Rix = Rel(1,i)
     Riy = Rel(2,i)
     Riz = Rel(3,i)

     do j = i-2 , 1, -2

       Rx = Rix - Rel(1,j)
       Ry = Riy - Rel(2,j)
       Rz = Riz - Rel(3,j)
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

       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)

       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if(R2>R2_MIN.and.R2<RcutC2) then

         cf = PCh(i) * PCh(j)
         InvR2 = 1.d0 / R2
         R1  = sqrt(R2)
         fk1 = cf * R1 * InvR2
         fk  = fk1 * InvR2

         if(Qslab) then

           Fx = fk * Rx
           Fy = fk * Ry
           Fz = fk * Rz

           virxy = Fx * Rx + Fy * Ry
           virzz = Fz * Rz
           IbinI = Kbin(i)
           IbinJ = Kbin(j)

           dZI = Zbin(i)
           dZJ = Zbin(j)
           call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,Nz)

         else

           call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,j)

         end if

       end if

     end do

     do j = i+1 , Nel, 2

       Rx = Rix - Rel(1,j)
       Ry = Riy - Rel(2,j)
       Rz = Riz - Rel(3,j)
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

       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)

       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       cf = PCh(i) * PCh(j)

       if(R2>R2_MIN.and.R2<RcutC2) then

         cf = PCh(i) * PCh(j)
         InvR2 = 1.d0 / R2
         R1  = sqrt(R2)
         fk1 = cf * R1 * InvR2
         fk  = fk1 * InvR2

         if(Qslab) then

           Fx = fk * Rx
           Fy = fk * Ry
           Fz = fk * Rz

           virxy = Fx * Rx + Fy * Ry
           virzz = Fz * Rz
           IbinI = Kbin(i)
           IbinJ = Kbin(j)

           dZI = Zbin(i)
           dZJ = Zbin(j)
           call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,Nz)

         else

           call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,j)

         end if

       end if

     end do

   end do

   deallocate(Rel,Kbin)

end subroutine Stress_LC


! ###########################################################
! ###########################################################


subroutine Stress_NoCELL(clhx,clhy,clhz)

use Configuration, only : R
use Numbers, only : N
use CellParam, only : CellShft
use ParamAnalyze, only : QEW, dZbin, Zbin, Ibin, RcutC2, Qslab
use CommonMPI
use CGdata, only : NBAtomType, NBFuncType, CoefAtype, CoefBtype, Rcut2
use EwaldParam, only : Alpha, msh, EFList, ar2
use NonbondParam, only : Charge
use NoLJParam, only : NumNoLJ
use CommonBlocks, only : QCoulomb

implicit none

integer :: i, j, k, ii, Nas, non, check_bonded
integer :: itype, jtype, ftype
integer :: Nx, Ny, Nz, IbinI, IbinJ
real(8) :: clhx, clhy, clhz, dZI, dZJ
real(8) :: Rix, Riy, Riz, Rx, Ry, Rz, R1, R2
real(8) :: aij, bij, InvR1, InvR2, InvR3, InvR4
real(8) :: term1, term2, fk, fke, cf
real(8) :: x, dx, yy, y1, y2, y3
real(8) :: z1, z2, xtm, fk1, fk2, ErrorFunc
real(8) :: Fx, Fy, Fz
real(8) :: virxy, virzz
external check_bonded

   Nas = NProcs - MyRank

   do i = Nas, N, NProcs

     non   = NumNoLJ(i)
     itype = NBAtomType(i)
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     do j = i-2 , 1, -2

       if(check_bonded(i,j,non) == 1) cycle

       jtype = NBAtomType(j)

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

       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)

       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if(R2<Rcut2(itype,jtype)) then

         ftype = NBFuncType(itype,jtype)

         if(ftype==7) then ! LJ12-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR4
           term2 = bij * InvR4
           fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2

         else if(ftype==1) then ! LJ9-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR1 = sqrt(InvR2)
           InvR3 = InvR2 * InvR1
           term1 = aij * InvR3 * InvR3 * InvR3
           term2 = bij * InvR3 * InvR3
           fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2

         end if

       else

         fk = 0.d0

       end if

       if(QCoulomb) then

         cf = Charge(i) * Charge(j)

         if(cf/=0.) then

           InvR2 = 1.d0 / R2

           if(QEW.and.R2<RcutC2) then

             R1  = sqrt(R2)
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

             fke  = ( fk1 + fk2 ) * InvR2
             fk = fk + fke

           else if(R2<RcutC2) then

             R1  = sqrt(R2)
             fk1 = cf * R1 * InvR2
             fke  = fk1 * InvR2
             fk = fk + fke

           end if

         end if

       end if

       if(fk/=0.) then

         if(Qslab) then

           Fx = fk * Rx
           Fy = fk * Ry
           Fz = fk * Rz

           virxy = Fx * Rx + Fy * Ry
           virzz = Fz * Rz
           IbinI = Ibin(i)
           IbinJ = Ibin(j)

           dZI = Zbin(i)
           dZJ = Zbin(j)
           call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,Nz)

         else

           call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,j)

         end if

       end if

     end do

     do j = i+1 , N, 2

       if(check_bonded(i,j,non) == 1) cycle

       jtype = NBAtomType(j)

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

       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)

       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if(R2<Rcut2(itype,jtype)) then

         ftype = NBFuncType(itype,jtype)

         if(ftype==7) then ! LJ12-4

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR4 = InvR2 * InvR2
           term1 = aij * InvR4 * InvR4 * InvR4
           term2 = bij * InvR4
           fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2

         else if(ftype==1) then ! LJ9-6

           aij = CoefAtype(itype,jtype)
           bij = CoefBtype(itype,jtype)
           InvR2 = 1.d0 / R2
           InvR1 = sqrt(InvR2)
           InvR3 = InvR2 * InvR1
           term1 = aij * InvR3 * InvR3 * InvR3
           term2 = bij * InvR3 * InvR3
           fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2

         end if

       else

         fk = 0.d0

       end if

       if(QCoulomb) then

         cf = Charge(i) * Charge(j)

         if(cf/=0.) then

           InvR2 = 1.d0 / R2

           if(QEW.and.R2<RcutC2) then

             R1  = sqrt(R2)
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

             fke  = ( fk1 + fk2 ) * InvR2
             fk = fk + fke

           else if(R2<RcutC2) then

             R1  = sqrt(R2)
             fk1 = cf * R1 * InvR2
             fke  = fk1 * InvR2
             fk = fk + fke

           end if

         end if

       end if

       if(fk/=0.) then

         Fx = fk * Rx
         Fy = fk * Ry
         Fz = fk * Rz

         if(Qslab) then

           virxy = Fx * Rx + Fy * Ry
           virzz = Fz * Rz
           IbinI = Ibin(i)
           IbinJ = Ibin(j)

           dZI = Zbin(i)
           dZJ = Zbin(j)
           call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,Nz)

         else

           call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,j)

         end if

       end if

     end do

   end do

end subroutine Stress_NoCELL

! ###########################################################
! ###########################################################


subroutine Stress_in_bin(virxy,virzz,dZI,dZJ,dZ,Rz,IbinI,IbinJ,Nz)

use ParamAnalyze, only : QIK, VirProXY, VirProZZ, RcutC2, Nbin

implicit none

integer :: i, ii, id
real(8) :: virxy, virzz, virwxy, virwzz, vmidxy, vmidzz
real(8) :: dZI, dZJ
real(8) :: dZ, Rz, W1
integer :: IbinI, IbinJ, Nz

   if(QIK) then

     if(IbinI == IbinJ) then

       VirProXY(IbinI) = VirProXY(IbinI) + virxy
       VirProZZ(IbinI) = VirProZZ(IbinI) + virzz

     else if(Nz == 0) then

       W1 = 1.d0 / abs(Rz)
       virwxy = virxy * W1
       virwzz = virzz * W1
       vmidxy = virwxy * dZ
       vmidzz = virwzz * dZ

       if(IbinI>IbinJ) then

         ii = IbinI - IbinJ - 1
         if(ii>=1) then
           do i = 1, ii
             id = IbinJ + i
             VirProXY(id) = VirProXY(id) + vmidxy
             VirProZZ(id) = VirProZZ(id) + vmidzz
           end do
         end if
         VirProXY(IbinI) = VirProXY(IbinI) + virwxy * dZI
         VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * dZI
         VirProXY(IbinJ) = VirProXY(IbinJ) + virwxy * (dZ - dZJ)
         VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * (dZ - dZJ)

       else

         ii = IbinJ - IbinI - 1
         if(ii>=1) then
           do i = 1, ii
             id = IbinI + i
             VirProXY(id) = VirProXY(id) + vmidxy
             VirProZZ(id) = VirProZZ(id) + vmidzz
           end do
         end if
         VirProXY(IbinI) = VirProXY(IbinI) + virwxy * (dZ - dZI)
         VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * (dZ - dZI)
         VirProXY(IbinJ) = VirProXY(IbinJ) + virwxy * dZJ
         VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * dZJ

       end if

     else if(Nz == 1) then

       W1 = 1.d0 / abs(Rz)
       virwxy = virxy * W1
       virwzz = virzz * W1
       vmidxy = virwxy * dZ
       vmidzz = virwzz * dZ

       if(IbinJ < Nbin) then
         do i = IbinJ+1, Nbin
           VirProXY(i) = VirProXY(i) + vmidxy
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if
       if(IbinI > 1) then
         do i = 1, IbinI-1
           VirProXY(i) = VirProXY(i) + vmidxy
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if
       VirProXY(IbinI) = VirProXY(IbinI) + virwxy * dZI
       VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * dZI
       VirProXY(IbinJ) = VirProXY(IbinJ) + virwxy * (dZ - dZJ)
       VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * (dZ - dZJ)

     else if(Nz == -1) then

       W1 = 1.d0 / abs(Rz)
       virwxy = virxy * W1
       virwzz = virzz * W1
       vmidxy = virwxy * dZ
       vmidzz = virwzz * dZ

       if(IbinI < Nbin) then
         do i = IbinI+1, Nbin
           VirProXY(i) = VirProXY(i) + vmidxy
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if
       if(IbinJ > 1) then
         do i = 1, IbinJ-1
           VirProXY(i) = VirProXY(i) + vmidxy
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if

       VirProXY(IbinI) = VirProXY(IbinI) + virwxy * (dZ - dZI)
       VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * (dZ - dZI)
       VirProXY(IbinJ) = VirProXY(IbinJ) + virwxy * dZJ
       VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * dZJ

     else

       write(*,*) "Error : no bin data is assigned"
       call Finalize

     end if

   else

     if(IbinI == IbinJ) then

       VirProXY(IbinI) = VirProXY(IbinI) + virxy

     else

       VirProXY(IbinI) = VirProXY(IbinI) + virxy*0.5d0
       VirProXY(IbinJ) = VirProXY(IbinJ) + virxy*0.5d0

     end if

     if(IbinI == IbinJ) then

       VirProZZ(IbinI) = VirProZZ(IbinI) + virzz

     else if(Nz == 0) then

       W1 = 1.d0 / abs(Rz)
       virwzz = virzz * W1
       vmidzz = virwzz * dZ

       if(IbinI>IbinJ) then

         ii = IbinI - IbinJ - 1
         do i = 1, ii
           id = IbinJ + i
           VirProZZ(id) = VirProZZ(id) + vmidzz
         end do
         VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * dZI
         VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * (dZ - dZJ)

       else

         ii = IbinJ - IbinI - 1
         do i = 1, ii
           id = IbinI + i
           VirProZZ(id) = VirProZZ(id) + vmidzz
         end do
         VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * (dZ - dZI)
         VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * dZJ

       end if

     else if(Nz == 1) then

       W1 = 1.d0 / abs(Rz)
       virwzz = virzz * W1
       vmidzz = virwzz * dZ

       if(IbinJ < Nbin) then
         do i = IbinJ+1, Nbin
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if
       if(IbinI > 1) then
         do i = 1, IbinI-1
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if
       VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * dZI
       VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * (dZ - dZJ)

     else if(Nz == -1) then

       W1 = 1.d0 / abs(Rz)
       virwzz = virzz * W1
       vmidzz = virwzz * dZ

       if(IbinI < Nbin) then
         do i = IbinI+1, Nbin
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if
       if(IbinJ > 1) then
         do i = 1, IbinJ-1
           VirProZZ(i) = VirProZZ(i) + vmidzz
         end do
       end if

       VirProZZ(IbinI) = VirProZZ(IbinI) + virwzz * (dZ - dZI)
       VirProZZ(IbinJ) = VirProZZ(IbinJ) + virwzz * dZJ

     else

       write(*,*) "Error : no bin data is assigned"
       call Finalize

     end if

   end if


end subroutine Stress_in_bin


!######################################################################
!######################################################################


! ***************************
! ** Bond Stretching Force **
! ***************************

subroutine Stress_Bond_CG(clhx,clhy,clhz)

use Numbers, only : N
use Configuration, only : R
use BondedParam, only : NumBond, BondI, BondJ, kBond, rBond, FTypeBond
use ParamAnalyze, only : QIK, QEW, Nbin, dZbin, Zbin, Ibin, Qslab
use CommonMPI
use CellParam, only : CellShft
use NonbondParam, only : Charge

implicit NONE

integer :: i, j, k, l, Nx, Ny, Nz
integer :: IbinI, IbinJ, Nas
real(8) :: R2, R1, dR, Fc, pref, cf, fke
real(8) :: Rx, Ry, Rz, Fx, Fy, Fz
real(8) :: clhx, clhy, clhz, Ueij
real(8) :: dZI, dZJ
real(8) :: virxy, virzz

   if( NumBond == 0 ) Return

   Nas = NProcs - MyRank

   do k = Nas , NumBond, NProcs

     i = BondI(k)
     j = BondJ(k)

     Rx = R(1,i) - R(1,j)
     Ry = R(2,i) - R(2,j)
     Rz = R(3,i) - R(3,j)

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

     R2  = Rx*Rx + Ry*Ry + Rz*Rz

     if (FTypeBond(k) == 1) then

       R1 = sqrt(R2)
       dR = R1 - rBond(k)
       Fc = kBond(k) * dR

       pref = - 2.d0 * Fc / R1

     else

       write(*,*) 'error: bond type in (Force_Bond_CG)'
       call Finalize

     end if

     if(QEW) then

       cf = Charge(i)*Charge(j)

       if(cf/=0.) then
         call COULOMB_SUB(Ueij,fke,R2,cf)
         pref = pref + fke
       end if

     end if

     Fx = pref * Rx
     Fy = pref * Ry
     Fz = pref * Rz

     if(Qslab) then
       virxy = Fx * Rx + Fy * Ry
       virzz = Fz * Rz
       IbinI = Ibin(i)
       IbinJ = Ibin(j)

       dZI = Zbin(i)
       dZJ = Zbin(j)
       call Stress_in_bin(virxy,virzz,dZI,dZJ,dZbin,Rz,IbinI,IbinJ,Nz)
     else
       call Stress_Sphere(pref,Rx,Ry,Rz,R2,i,j)
     end if

   end do

end subroutine Stress_Bond_CG


!######################################################################
!######################################################################


subroutine Stress_Angle_CG(clhx,clhy,clhz)

use Numbers, only : N
use Configuration, only : R
use UnitExParam, only : pi
use BondedParam, only : NumAngle, AngleI, AngleJ, AngleK, kTheta, &
&   Theta0, FTypeAngle, vdWSubtAng
use ParamAnalyze, only : QIK, QEW, Nbin, dZbin, Zbin, Ibin, Qslab
use CommonMPI
use CellParam, only : CellShft
use NonbondParam, only : Charge
use CGdata

implicit NONE

integer :: i, j, k, l, ll, Nas
real(8) :: clhx, clhy, clhz
real(8) :: Ra2, Rb2, fk, fke, R2, Ueij
real(8) :: rRa2, rRb2, rRab
real(8) :: Cst, Theta, Thpi, Thpi2
real(8) :: cf, dT, Snt, t1
real(8) :: Rijx, Rkjx
real(8) :: Rijy, Rkjy
real(8) :: Rijz, Rkjz
real(8) :: Rx, Ry, Rz
real(8) :: Fix, Fiy, Fiz
real(8) :: Fkx, Fky, Fkz
real(8) :: Fx, Fy, Fz
real(8) :: pref, pref1, pref2
real(8) :: poly_quart
real(8) :: virijxy, virijzz, virkjxy, virkjzz
real(8) :: dZI, dZJ, dZK
real(8) :: aij, bij, InvR2, InvR1, InvR3, InvR4
real(8) :: term1, term2, virxy, virzz
integer :: IbinI, IbinJ, IbinK
integer :: Nx, Ny, Nz, Nijz, Nkjz
integer :: itype, jtype, ftype
external poly_quart

   if( NumAngle == 0 ) Return

   Nas = NProcs - MyRank

   do l = Nas , NumAngle, NProcs

     i = AngleI(l)
     j = AngleJ(l)
     k = AngleK(l)

     Rx = R(1,i) - R(1,j)
     Ry = R(2,i) - R(2,j)
     Rz = R(3,i) - R(3,j)
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
     ll = Nx + Ny + Nz
     Rijx = Rx + CellShft(1,ll)
     Rijy = Ry + CellShft(2,ll)
     Rijz = Rz + CellShft(3,ll)
     Nijz = Nz

     Rx = R(1,k) - R(1,j)
     Ry = R(2,k) - R(2,j)
     Rz = R(3,k) - R(3,j)
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
     ll = Nx + Ny + Nz
     Rkjx = Rx + CellShft(1,ll)
     Rkjy = Ry + CellShft(2,ll)
     Rkjz = Rz + CellShft(3,ll)
     Nkjz = Nz

     Ra2  = Rijx * Rijx + Rijy * Rijy + Rijz * Rijz
     Rb2  = Rkjx * Rkjx + Rkjy * Rkjy + Rkjz * Rkjz

     rRa2 = 1.d0 / Ra2
     rRb2 = 1.d0 / Rb2
     rRab = sqrt( rRa2 * rRb2 )

     Cst = ( Rijx * Rkjx + Rijy * Rkjy + Rijz * Rkjz ) * rRab

     if(FTypeAngle(l) == 1) then ! cosine function

       cf = - kTheta(l)

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

     else if(FTypeAngle(l) == 3) then

       if(Cst <= -1.d0) then

         dT = pi - Theta0(l)
         cf = 0.d0

       else

         Theta = acos(Cst)

         dT = Theta - Theta0(l)
         cf = 2.d0 * kTheta(l) * dT / sin(Theta)

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

     if(Qslab) then
       virijxy = Fix * Rijx + Fiy * Rijy
       virijzz = Fiz * Rijz

       virkjxy = Fkx * Rkjx + Fky * Rkjy
       virkjzz = Fkz * Rkjz

       IbinI = Ibin(i)
       IbinJ = Ibin(j)
       IbinK = Ibin(k)

       dZI = Zbin(i)
       dZJ = Zbin(j)
       dZK = Zbin(k)
       call Stress_in_bin(virijxy,virijzz,dZI,dZJ,dZbin,Rijz,IbinI,IbinJ,Nijz)
       call Stress_in_bin(virkjxy,virkjzz,dZK,dZJ,dZbin,Rkjz,IbinK,IbinJ,Nkjz)
     else
       call Stress_Sphere2(Fix,Fiy,Fiz,rRa2,i,j)
       call Stress_Sphere2(Fkx,Fky,Fkz,rRb2,k,j)
     end if

!#### VDW CORRECTION

     if(vdWSubtAng(l)) cycle

     Rx = R(1,i) - R(1,k)
     Ry = R(2,i) - R(2,k)
     Rz = R(3,i) - R(3,k)
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
     ll = Nx + Ny + Nz
     Rx = Rx + CellShft(1,ll)
     Ry = Ry + CellShft(2,ll)
     Rz = Rz + CellShft(3,ll)

     R2 = Rx*Rx + Ry*Ry + Rz*Rz

     itype = NBAtomType(i)
     jtype = NBAtomType(k)

     fk = 0.d0

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
         fk = ( 9.d0 * term1 - 6.d0 * term2 ) * InvR2

       else if(ftype==7) then ! LJ12-4

         aij = CoefAtype(itype,jtype)
         bij = CoefBtype(itype,jtype)
         InvR2 = 1.d0 / R2
         InvR4 = InvR2 * InvR2
         term1 = aij * InvR4 * InvR4 * InvR4
         term2 = bij * InvR4
         fk = ( 12.d0 * term1 - 4.d0 * term2 ) * InvR2

       end if

     end if

     if(QEW) then

       cf = Charge(i) * Charge(k)
       if(cf/=0.) then
         call COULOMB_SUB(Ueij,fke,R2,cf)
         fk = fk + fke
       end if

     end if

     if(fk/=0.) then

       if(Qslab) then

         Fx = fk * Rx
         Fy = fk * Ry
         Fz = fk * Rz

         virxy = Fx * Rx + Fy * Ry
         virzz = Fz * Rz
         call Stress_in_bin(virxy,virzz,dZI,dZK,dZbin,Rz,IbinI,IbinK,Nz)
       else
         call Stress_Sphere(fk,Rx,Ry,Rz,R2,i,k)
       end if

     end if

   end do

end subroutine Stress_Angle_CG


!#####################################################################
!#####################################################################


subroutine Stress_Reciprocal

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use UnitExParam, only : InvPi, sqpi, pi2
use EwaldParam, only : Nel, Nh, Nelist, ih, alp2, PCh
use CellParam, only : InvH, Volume
use ParamAnalyze, only : Nbin, Ibin, VirProXY, VirProZZ

implicit NONE

integer :: i, j, l
real(8) :: Wkx, Wky, Wkz
real(8) :: CsSum, SnSum, tyo, pref
real(8) :: et, et1, pa2, arkn2
integer :: IbinI
real(8) :: zz
real(8) :: kn2, Epkn, Trp, Invkn2, InvVOL
real(8) :: IHxx, IHxy, IHxz, IHyx, IHyy, IHyz, IHzx, IHzy, IHzz
real(8) :: knx, kny, knz
integer :: ihx, ihy, ihz
integer, dimension(:), allocatable :: Kbin
real(8), dimension(:), allocatable :: Csl, Snl
real(8), dimension(:,:), allocatable :: Rel, Pvir

   pref = -sqpi * alp2
   InvVOL = 1.d0 / Volume

   et  = InvPi * InvVOL
   pa2 = -sqpi * alp2 * 2.d0

   IHxx = InvH(1,1)
   IHxy = InvH(1,2)
   IHxz = InvH(1,3)
   IHyx = InvH(2,1)
   IHyy = InvH(2,2)
   IHyz = InvH(2,3)
   IHzx = InvH(3,1)
   IHzy = InvH(3,2)
   IHzz = InvH(3,3)

   allocate( Rel(3,Nel) )
   allocate( Snl(Nel) )
   allocate( Csl(Nel) )
   allocate( Kbin(Nel) )
   allocate( Pvir(3,Nel) )

   do i = 1, Nel
     j = Nelist(i)
     Rel(1,i) = R(1,j)
     Rel(2,i) = R(2,j)
     Rel(3,i) = R(3,j)
     Kbin(i) = Ibin(j)
   end do

   Pvir(:,:) = 0.d0

!--------------------------------------------------------------------------

   do l = 1 , Nh

     ihx = ih(1,l)
     ihy = ih(2,l)
     ihz = ih(3,l)

     knx = IHxx*ihx + IHyx*ihy + IHzx*ihz
     kny = IHxy*ihx + IHyy*ihy + IHzy*ihz
     knz = IHxz*ihx + IHyz*ihy + IHzz*ihz
     kn2 = knx*knx + kny*kny + knz*knz
     Invkn2 = 1.d0 / kn2
     Epkn  = exp( pref * kn2 ) * Invkn2

     Wkx = knx * pi2
     Wky = kny * pi2
     Wkz = knz * pi2
     CsSum = 0.d0
     SnSum = 0.d0

     do i = 1, Nel

       zz = Pch(i)

       tyo = Wkx*Rel(1,i) + Wky*Rel(2,i) + Wkz*Rel(3,i)

       Csl(i) = zz * cos(tyo)
       Snl(i) = zz * sin(tyo)

       CsSum = CsSum + Csl(i)
       SnSum = SnSum + Snl(i)

     end do

     do i = 1, Nel

       Trp = Csl(i) * CsSum + Snl(i) * SnSum

       et1   =  et * Epkn * Trp
       arkn2 = ( -2.d0 * Invkn2 + pa2 ) * et1

       Pvir(1,i) = Pvir(1,i) + et1 + arkn2 * knx * knx
       Pvir(2,i) = Pvir(2,i) + et1 + arkn2 * kny * kny
       Pvir(3,i) = Pvir(3,i) + et1 + arkn2 * knz * knz

     end do

   end do

!----------------------------------------------------------------------

   do i = 1, Nel

     IbinI = Kbin(i)
     VirProXY(IbinI) = VirProXY(IbinI) + Pvir(1,i) + Pvir(2,i)
     VirProZZ(IbinI) = VirProZZ(IbinI) + Pvir(3,i)

   end do

!----------------------------------------------------------------------
   deallocate( Rel, Snl, Csl, Kbin, Pvir )

end subroutine Stress_Reciprocal


!#####################################################################
!#####################################################################


subroutine Stress_Sphere(fijdivRij,Rx,Ry,Rz,Rij2,i,j)

use Configuration, only : R
use ParamAnalyze, only : Zbin, Ibin, Nbin, dRbin, InvdRbin, VirTProN, VirTProT

implicit none

integer :: L1, L2, L, Nlam1, Nlam2, Nbin2, Nlam, N1, i, j
real(8) :: fijdivRij, Rx, Ry, Rz, Rij2, SR2
real(8) :: R1, R2, Rnx1, Rny1, Rnz1, Rnx2, Rny2, Rnz2
real(8) :: Rx12, Ry12, Rz12, R1R12, c, rr, b
real(8) :: Lamda1, Lamda2, ccsqrt
real(8), dimension(:), allocatable :: Lam1, Lam2
real(8), dimension(:), allocatable :: Lay1, Lay2
real(8) :: b2ac, Invb2ac
real(8), parameter :: epslog=1.d-15
integer :: Jcount, Jlam
real(8) :: Lamda, att, btt
real(8) :: a12, x12, y12, x22, x12_4, b2, bh, xx2
real(8) :: fn, ft

   if((Ibin(i)>Nbin).and.(Ibin(j)>Nbin)) return

   Nbin2 = Nbin*2
   allocate( Lam1(0:Nbin2), Lam2(0:Nbin2), Lay1(0:Nbin2), Lay2(0:Nbin2) )

   SR2 = 1.d0 / Rij2

   if(Zbin(i)<Zbin(j)) then
     R1  = Zbin(j)
     R2  = Zbin(i)
     Rnx1 = R(1,j)
     Rny1 = R(2,j)
     Rnz1 = R(3,j)
     Rnx2 = R(1,i)
     Rny2 = R(2,i)
     Rnz2 = R(3,i)
     Rx12 = Rx
     Ry12 = Ry
     Rz12 = Rz
   else
     R1  = Zbin(i)
     R2  = Zbin(j)
     Rnx1 = R(1,i)
     Rny1 = R(2,i)
     Rnz1 = R(3,i)
     Rnx2 = R(1,j)
     Rny2 = R(2,j)
     Rnz2 = R(3,j)
     Rx12 = - Rx
     Ry12 = - Ry
     Rz12 = - Rz
   end if

   L1 = int(R1*InvdRbin)+1
   L2 = int(R2*InvdRbin)+1
   R1R12 = Rx12 * Rnx1 + Ry12 * Rny1 + Rz12 * Rnz1

   a12 = R1 * R1
   x12 = a12 * SR2
   y12 = a12 * Rij2
   x22 = R2 * R2 * SR2
   x12_4 = 4.d0 * x12
   b  = 1.d0 + x12 - x22
   b2 = b * b
   bh = b * 0.5d0

! Kirkwood
!
! look for the intersects between Line R12 and circle r(=tlayer*dble(l))
! Lam1()=Lamda at intersect (=(intersect-R1)/rij)
! Lay1()=layer number
! Lam2(), Lay2() are used, when there are two intersects.  
! They are included in Lam1(), Lay1(), and later.
! *SR2 is equivalent to /Rij2
   c = 0.d0
   Nlam1 = 0
   Nlam2 = 0
   do l = min(Nbin,L1-1),1,-1
     if (c>=0.d0) then
       rr = dRbin * dble(l)
       xx2 = rr * rr * SR2
       c  = b2 + 4.d0 * xx2 - x12_4
       if (c>=0.d0) then
         ccsqrt = sqrt(c) * 0.5d0
         Lamda1 = bh-ccsqrt
         Lamda2 = bh+ccsqrt
         if (Lamda1>1.d0) c = -1.d0
         if ((Lamda1>=0.d0).and.(Lamda1<=1.d0)) then
            Nlam1 = Nlam1 + 1
            Lam1(Nlam1) = Lamda1
            Lay1(Nlam1) = l
            if (Lamda2<=1.d0) then
              Nlam2 = Nlam2 + 1
              Lam2(Nlam2) = Lamda2
              Lay2(Nlam2) = l
            end if
         end if
       end if
     end if
   end do

!     combine two intersects
   Nlam=Nlam1
   do N1 = Nlam2, 1, -1
     Nlam       = Nlam+1
     Lay1(Nlam) = Lay2(N1)
     Lam1(Nlam) = Lam2(N1)
   end do
   Nlam1       = Nlam + 1
   Lay1(Nlam1) = L2
   Lam1(Nlam1) = 1.d0
   Lay1(0)     = L1
   Lam1(0)     = 0.d0

   b2ac = y12 - R1R12*R1R12

   Lamda  = -R1R12*SR2

   Jcount = 0

   if (b2ac <= epslog) then

     do Nlam = 1, Nlam1

       L = max(Lay1(Nlam),Lay1(Nlam-1))

       if (L <= Nbin) then

         Jcount = Jcount + 1

         if (Lamda < Lam1(Nlam-1)) then
           Jlam =  1
         else if (Lamda < Lam1(Nlam)) then
           Jlam =  0
         else
           Jlam = -1
         end if

         fn =  Rij2 * (Lam1(Nlam)-Lam1(Nlam-1))
         ft =  0.d0

!   fijdivRij=Force/Rij !!!!!!!!!,modified     

         VirTProN(L) = VirTProN(L) + fn * fijdivRij
!         VirTProT(L) = VirTProT(L) + ft * fijdivRij

       end if

     end do

   else

     b2ac=sqrt(b2ac)
     Invb2ac = 1.d0/b2ac

     do Nlam = 1, Nlam1

       L = max(Lay1(Nlam),Lay1(Nlam-1))

       if (L <= Nbin) then

         Jcount = Jcount + 1

         if (Lamda < Lam1(Nlam-1)) then
           Jlam =  1
         else if (Lamda < Lam1(Nlam)) then
           Jlam =  0
         else
           Jlam = -1
         end if

         if(Jcount==1) then
           att=atan((Rij2*Lam1(Nlam-1)+R1R12)*Invb2ac)
         else
           att=btt
         end if
         btt=atan((Rij2*Lam1(Nlam)+R1R12)*Invb2ac)

!    Fn is function that can be found in paper....

         fn =  Rij2 * (Lam1(Nlam)-Lam1(Nlam-1)) - b2ac * (btt-att)
         ft =  b2ac*(btt-att)*0.5d0

!   fijdivRij=Force/Rij !!!!!!!!!,modified     

         VirTProN(L) = VirTProN(L) + fn * fijdivRij
         VirTProT(L) = VirTProT(L) + ft * fijdivRij

       end if

     end do

   end if 

end subroutine Stress_Sphere


!#####################################################################
!#####################################################################


subroutine Stress_Sphere2(Fix,Fiy,Fiz,SR2,i,j)

use Configuration, only : R
use ParamAnalyze, only : Zbin, Ibin, Nbin, dRbin, InvdRbin, VirTProN, VirTProT

implicit none

real(8) :: Fix, Fiy, Fiz
integer :: L1, L2, l, Nlam1, Nlam2, Nbin2, Nlam, N1, i, j
real(8) :: SR2
real(8) :: R1, R2
real(8) :: c, rr, b
real(8) :: Lamda1, Lamda2, ccsqrt
real(8), dimension(:), allocatable :: Lay1, Lay2, Laytmp
integer :: Jcount
real(8) :: a12, x12, x22, x12_4, b2, bh, xx2
real(8) :: fn
real(8) :: aln, bln
integer :: ijflag
real(8) :: ffactor

   if((Ibin(i)>Nbin).and.(Ibin(j)>Nbin)) return

   Nbin2 = Nbin*2
   allocate( Lay1(0:Nbin2), Lay2(0:Nbin2) )

   ijflag = 0
   if(Zbin(i)<Zbin(j)) then
     R1  = Zbin(j)
     R2  = Zbin(i)
   else
     R1  = Zbin(i)
     R2  = Zbin(j)
     ijflag = 1
   end if

   L1 = int(R1*InvdRbin)+1
   L2 = int(R2*InvdRbin)+1

   a12 = R1 * R1
   x12 = a12 * SR2
   x22 = R2 * R2 * SR2
   x12_4 = 4.d0 * x12
   b  = 1.d0 + x12 - x22
   b2 = b * b
   bh = b * 0.5d0

! Kirkwood
!
! look for the intersects between Line R12 and circle r(=tlayer*dble(l))
! Lam1()=Lamda at intersect (=(intersect-R1)/rij)
! Lay1()=layer number
! Lam2(), Lay2() are used, when there are two intersects.  
! They are included in Lam1(), Lay1(), and later.
! *SR2 is equivalent to /Rij2

   c = 0.d0
   Nlam1 = 0
   Nlam2 = 0

   do l = min(Nbin,L1-1),1,-1
     if (c>=0.d0) then
       rr = dRbin * dble(l)
       xx2 = rr * rr * SR2
       c  = b2 + 4.d0 * xx2 - x12_4
      if (c>=0.d0) then
        ccsqrt = sqrt(c) * 0.5d0
        Lamda1 = bh-ccsqrt
        Lamda2 = bh+ccsqrt
        if (Lamda1 > 1.d0) c = -1.d0
        if ((Lamda1 >= 0.d0).and.(Lamda1 <= 1.d0)) then
          Nlam1 = Nlam1 + 1
          Lay1(Nlam1) = l
          if (Lamda2<=1.d0) then
            Nlam2 = Nlam2 + 1
            Lay2(Nlam2) = l
          end if
        end if
      end if
    end if
  end do

  !     combine two intersects
  Nlam = Nlam1
  do N1 = Nlam2, 1, -1
     Nlam       = Nlam + 1
     Lay1(Nlam) = Lay2(N1)
  end do
  Nlam1 = Nlam + 1
  Lay1(Nlam1) = L2
  Lay1(0)     = L1


  ! for the i->j case 
  if ( ijflag == 1) then
     allocate( Laytmp(0:Nlam1))
     do l = 0, Nlam1
        Laytmp(l) = Lay1(Nlam1-l)
     end do
     do l = 0, Nlam1
        Lay1(l) = Laytmp(l)
     end do
     deallocate( Laytmp )
  end if


  ! replace as the case ijflag = 0 

  Jcount = 0

  ! Factor Calculation From Here !
  ffactor = Fix * R(1,j) + Fiy * R(2,j) + Fiz * R(3,j)

  do Nlam = 1, Nlam1

     l=max(Lay1(Nlam),Lay1(Nlam-1))

     if (l<=Nbin) then
        Jcount=Jcount+1

        if(Jcount==1) then
           aln = Lay1(Nlam-1)
        else
           aln = bln
        end if
        bln = Lay1(Nlam)
        fn  = log(bln/aln)

        VirTProN(l) = VirTProN(l) + ffactor * fn
        VirTProT(l) = VirTProT(l) - 0.5d0 * ffactor* fn
     end if

  end do

end subroutine Stress_Sphere2

