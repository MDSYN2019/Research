
subroutine CohesiveEne

use CommonBlocks, only : QMaster, QSwitch, QCorrectCutoff, ForceField
use NonbondParam, only : Ene_Ersp, Ene_LJ
use EwaldParam, only : Ene_Eksp, Ene_Eslf
use CellParam, only : H, InvH, Volume
use TailCorrect
use UnitExParam, only : cvol
use ParamAnalyze, only : NJobs, NtrjStep

implicit NONE

integer :: i, j, istep, TotalStep
real(8), dimension(2) :: Ene, AvE
real(8) :: Ecoh, det
external det

   if(QMaster) then
   open(111,file='cohesive_time.dat')

   write(111,*) '# step, cohesive energy, LJ term, Coulomb term'
   write(111,*) '#       energy unit is kcal/mol'

   if(QCorrectCutoff) then
     write(111,*) '# LJ cutoff correction is taken into account !!'
   end if
   end if

   if(ForceField(1:4)/='OPLS'.and.ForceField(1:5)/='CHARM') then
     if(QMaster) write(11,*) 'error : this analysis is useful for OPLS or CHARMM force field'
     if(QMaster) write(*,*) 'error : this analysis is useful for OPLS or CHARMM force field'
     call Finalize
   end if

   call PrepCohesEne
   if(QCorrectCutoff) call CorrectPotential
   call ErrorFuncList           !  make a table of error-function
   call RecLatticeList          !  define reciprocal lattice
   call Ewald_SelfTerm

! -------------------------------------------------

   AvE(:) = 0.d0

   TotalStep = 0
   istep = 0

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     TotalStep = TotalStep + NTrjStep(i)

     do j = 1 , NTrjStep(i)

       istep = istep + 1

!     ------------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     ------------------------------

!     ---------------
       call BcastRH
!     ---------------

       call InversMatrix(H,InvH)
       Volume = det(H)

!     ---------------------
       call CohEne_Real
       call CohEne_Kspace
!     ---------------------

       Ene(1) = Ene_LJ
       Ene(2) = Ene_Ersp + Ene_Eksp

       call SumEne(Ene,2)

       if(QMaster) then
         Ene(2) = Ene(2) + Ene_Eslf
         if(QCorrectCutoff) then
           Ene_LJ_co = CorrectE / Volume
           Ene(1) = Ene(1) + Ene_LJ_co
         end if
         Ene(:) = Ene(:) * cvol
         Ecoh = Ene(1) + Ene(2)
         write(111,'(i9,3e15.7)') istep,Ecoh,Ene(1),Ene(2)
       end if

       AvE(:) = AvE(:) + Ene(:)

     end do

   end do

   if(QMaster) then

   close(111)

   AvE(:) = AvE(:) / TotalStep

   open(112,file='cohesive_ave.dat')

   if(QCorrectCutoff) then
     write(112,*) '# LJ cutoff correction is taken into account !!'
   end if
   write(112,'(/a,e15.7,a/)') ' Cohesive Energy = ',AvE(1)+AvE(2),'[kcal/mol]'
   write(112,'(a,e15.7,a)')   '   LJ term       = ',AvE(1),'[kcal/mol]'
   write(112,'(a,e15.7,a)')   '   Coulomb term  = ',AvE(2),'[kcal/mol]'

   close(112)

   end if

end subroutine CohesiveEne


subroutine PrepCohesEne

use Numbers, only : N, NumSpec, NumMol, NumAtm
use ParamAnalyze, only : MolID

implicit NONE

integer :: i, j, k, l, ii

   allocate(MolID(N))

   k = 0
   l = 0
   do ii = 1, NumSpec
     do i = 1, NumMol(ii)
       k = k + 1
       do j = 1, NumAtm(ii)
         l = l + 1
         MolID(l) = k
       end do
     end do
   end do

end subroutine PrepCohesEne

subroutine CohEne_Real

use Numbers, only : N
use Configuration, only : R
use CommonBlocks, only : QSwitch, QCorrectCutoff, ForceField
use EwaldParam, only : Alpha, ar2
use CommonMPI
use CellParam, only : H, CellShft
use CutoffParam, only : Ron2, Rcutoff2, swf1
use NonbondParam, only : Charge, SgmLJ, EpsLJ, Ene_Ersp, Ene_LJ, Rminh
use ParamAnalyze, only : MolID

implicit NONE

integer :: i, j, k, Nas
real(8) :: Six, Siy, Siz
real(8) :: Rix, Riy, Riz
real(8) :: Sx, Sy, Sz
real(8) :: Rx, Ry, Rz
integer :: Nx, Ny, Nz
real(8), dimension(3,N) :: ScR
real(8) :: R2, InvR2, R1
real(8) :: Sgm, Sgm2, Eps
real(8) :: SR2, SR6, SR12
real(8) :: ErrorFunc
real(8) :: fk1, ek, cf, x
real(8) :: Error_Function
external Error_Function


   call PBC
   call ScaledCoordinate(ScR)

   Ene_LJ = 0.d0
   Ene_Ersp = 0.d0

   Nas = NProcs - MyRank

if(ForceField(1:4)=='OPLS') then

   do i = Nas, N, NProcs

     Six = ScR(1,i)
     Siy = ScR(2,i)
     Siz = ScR(3,i)
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     do j = i-2 , 1, -2

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
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
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rcutoff2 ) then

         InvR2 = 1.d0 / R2

         if(MolID(i)/=MolID(j)) then

           Eps   = EpsLJ(i) * EpsLJ(j)
           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           ek   = Eps * 4.d0 * ( SR12 - SR6 )
! --------------------------------------------------------
 if(QSwitch) call SwitchF(R2,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------
           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * ErrorFunc * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         else

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * (ErrorFunc - 1.d0) * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         end if

       else if(MolID(i)==MolID(j)) then

         InvR2 = 1.d0 / R2

         if(QCorrectCutoff) then
           Eps   = EpsLJ(i) * EpsLJ(j)
           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           ek   = Eps * 4.d0 * ( SR12 - SR6 )
           Ene_LJ = Ene_LJ - ek
         end if

         cf = Charge(i) * Charge(j)
         R1  = sqrt( R2 )
         fk1 = - cf * R1 * InvR2

         Ene_Ersp = Ene_Ersp + fk1

       end if

     end do

     do j = i+1 , N, 2

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
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
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rcutoff2 ) then

         InvR2 = 1.d0 / R2

         if(MolID(i)/=MolID(j)) then

           Eps   = EpsLJ(i) * EpsLJ(j)
           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           ek   = Eps * 4.d0 * ( SR12 - SR6 )
! --------------------------------------------------------
 if(QSwitch) call SwitchF(R2,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------
           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * ErrorFunc * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         else

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * (ErrorFunc - 1.d0) * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         end if

       else if(MolID(i)==MolID(j)) then

         InvR2 = 1.d0 / R2

         if(QCorrectCutoff) then
           Eps   = EpsLJ(i) * EpsLJ(j)
           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           ek   = Eps * 4.d0 * ( SR12 - SR6 )
           Ene_LJ = Ene_LJ - ek
         end if

         cf = Charge(i) * Charge(j)
         R1  = sqrt( R2 )
         fk1 = - cf * R1 * InvR2

         Ene_Ersp = Ene_Ersp + fk1

       end if

     end do

   end do

else if(ForceField(1:5)=='CHARM') then

   do i = Nas, N, NProcs

     Six = ScR(1,i)
     Siy = ScR(2,i)
     Siz = ScR(3,i)
     Rix = R(1,i)
     Riy = R(2,i)
     Riz = R(3,i)

     do j = i-2 , 1, -2

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
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
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rcutoff2 ) then

         InvR2 = 1.d0 / R2

         if(MolID(i)/=MolID(j)) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)

           SR2  = Sgm2 * InvR2
           SR6  = SR2 * SR2 * SR2
           SR12 = SR6 * SR6
           ek   = Eps * ( SR12 - 2.d0 * SR6 )
! --------------------------------------------------------
 if(QSwitch) call SwitchF(R2,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------
           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * ErrorFunc * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         else

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * (ErrorFunc - 1.d0) * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         end if

       else if(MolID(i)==MolID(j)) then

         InvR2 = 1.d0 / R2

         if(QCorrectCutoff) then
           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)

           SR2  = Sgm2 * InvR2
           SR6  = SR2 * SR2 * SR2
           SR12 = SR6 * SR6
           ek   = Eps * ( SR12 - 2.d0 * SR6 )
           Ene_LJ = Ene_LJ - ek
         end if

         cf = Charge(i) * Charge(j)
         R1  = sqrt( R2 )
         fk1 = - cf * R1 * InvR2

         Ene_Ersp = Ene_Ersp + fk1

       end if

     end do

     do j = i+1 , N, 2

       Sx = Six - ScR(1,j)
       Sy = Siy - ScR(2,j)
       Sz = Siz - ScR(3,j)
       Rx = Rix - R(1,j)
       Ry = Riy - R(2,j)
       Rz = Riz - R(3,j)
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
       k = Nx + Ny + Nz
       Rx = Rx + CellShft(1,k)
       Ry = Ry + CellShft(2,k)
       Rz = Rz + CellShft(3,k)
       R2 = Rx*Rx + Ry*Ry + Rz*Rz

       if ( R2 <= Rcutoff2 ) then

         InvR2 = 1.d0 / R2

         if(MolID(i)/=MolID(j)) then

           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)

           SR2  = Sgm2 * InvR2
           SR6  = SR2 * SR2 * SR2
           SR12 = SR6 * SR6
           ek   = Eps * ( SR12 - 2.d0 * SR6 )
! --------------------------------------------------------
 if(QSwitch) call SwitchF(R2,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------
           Ene_LJ = Ene_LJ + ek

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * ErrorFunc * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         else

           cf = Charge(i) * Charge(j)
           R1  = sqrt( R2 )
           x   = Alpha * R1
           ErrorFunc = Error_Function(x)
           fk1 = cf * (ErrorFunc - 1.d0) * R1 * InvR2

           Ene_Ersp = Ene_Ersp + fk1

         end if

       else if(MolID(i)==MolID(j)) then

         InvR2 = 1.d0 / R2

         if(QCorrectCutoff) then
           Sgm   = Rminh(i) + Rminh(j)
           Sgm2  = Sgm * Sgm
           Eps   = EpsLJ(i) * EpsLJ(j)

           SR2  = Sgm2 * InvR2
           SR6  = SR2 * SR2 * SR2
           SR12 = SR6 * SR6
           ek   = Eps * ( SR12 - 2.d0 * SR6 )
           Ene_LJ = Ene_LJ - ek
         end if

         cf = Charge(i) * Charge(j)
         R1  = sqrt( R2 )
         fk1 = - cf * R1 * InvR2

         Ene_Ersp = Ene_Ersp + fk1

       end if

     end do

   end do

end if

Contains

subroutine SwitchF(R2,ek,Ron2,Rcutoff2,Swpref)

implicit none

real(8) :: R2, ek
real(8) :: Xoff, Xon
real(8) :: Ron2, Rcutoff2, Swpref
real(8) :: Switch, Dswitch

   if(R2 > Ron2) then

     Xoff = R2 - Rcutoff2      ! x-xoff
     Xon  = R2 - Ron2          ! x-xon

     Switch  = ( 3.d0 * Xon - Xoff ) * Xoff * Xoff * Swpref

     ek = ek * Switch

   end if

end subroutine SwitchF

end subroutine CohEne_Real

!######################################################################

subroutine CohEne_Kspace

use Numbers, only : N
use Configuration, only : R
use CommonMPI
use UnitExParam, only : InvPi, sqpi, pi2
use EwaldParam, only : Nel, Nh, Nelist, ih, alp2, Ene_Eksp, PCh
use CellParam, only : InvH, Volume

implicit NONE

integer :: i, j, k, l

real(8) :: Wkx, Wky, Wkz
real(8) :: CsSum, SnSum, tyo, pref
integer :: Nc, MyType, ii
real(8) :: zz
real(8) :: kn2, Epkn, Trp, Invkn2, InvVOL
real(8) :: IHxx, IHxy, IHxz, IHyx, IHyy, IHyz, IHzx, IHzy, IHzz
real(8) :: knx, kny, knz
integer :: ihx, ihy, ihz
real(8), dimension(:), allocatable :: Csl, Snl
real(8), dimension(:,:), allocatable :: Rel, Fel

   pref = -sqpi * alp2
   InvVOL = 1.d0 / Volume

   Ene_Eksp = 0.d0

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
   allocate( Fel(3,Nel) )
   allocate( Snl(Nel) )
   allocate( Csl(Nel) )

   do i = 1, Nel
     j = Nelist(i)
     Rel(1,i) = R(1,j)
     Rel(2,i) = R(2,j)
     Rel(3,i) = R(3,j)
   end do

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

     Trp = CsSum * CsSum + SnSum * SnSum

     Ene_Eksp = Ene_Eksp + Epkn * Trp

   end do

!----------------------------------------------------------------------

   Ene_Eksp = Ene_Eksp * InvVOL * InvPi

!----------------------------------------------------------------------
   deallocate( Rel, Fel, Snl, Csl )

end subroutine CohEne_Kspace
