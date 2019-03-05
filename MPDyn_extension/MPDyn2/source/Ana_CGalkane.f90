! ############################
! ## SUBROUTINE LIST 
! ## -- CGAnalyz (MODULE)
! ## -- CoarseGrainC12 
! ## -- CGC12_data 
! ## -- CalcCOM_Gyr 
! ## -- Force_NonBondIntraCG 
! ## -- ForceIntraCG 
! ## -- EneIntraAssign 
! ## -- ForceInterCG 
! ## -- Print_CGInteraction 
! ############################


module CGAnalyz

   integer :: Kcomp, IRcut
   integer, dimension(:,:), allocatable :: listGR
   integer :: NoAtom, NoMol
! ## Grouping
   integer :: NumGroup
   integer, dimension(:), allocatable :: GroupNumber, NumAtm_Group
   real(8), dimension(:), allocatable :: Mass_Group
   real(8), dimension(:,:,:), allocatable :: Rcom_Group

   integer :: NumTypes
   integer, dimension(:), allocatable :: NtypeGroup
   integer, dimension(:), allocatable :: NinTypeGroup
   integer, dimension(:,:), allocatable :: NonInteractPair

   character(len=50) :: CGdefine_File
   integer :: NumTypePair
   integer, parameter :: ndr = 10
   real(8), parameter :: slabR = 1.d0 / ndr
   integer, dimension(:,:), allocatable :: TypePair

   integer, dimension(:,:), allocatable :: NAtomGroup
   real(8), dimension(:,:), allocatable :: EneNonBond, FrcNonBond

   integer :: NumCGBond, NumCGAngle
   integer, dimension(:), allocatable :: BondCGI, BondCGJ
   integer, dimension(:), allocatable :: AngleCGI, AngleCGJ, AngleCGK

end module CGAnalyz


!######################################################################
!######################################################################


subroutine CoarseGrainC12

use CommonBlocks, only : QMaster
use CellParam, only : H, InvH
use ParamAnalyze, only : NJobs, NtrjStep, Interval
use CGAnalyz

implicit none

integer :: i, j, Numa, totalstep
real(8) :: Vol, det
external det

   call CGC12_data(Numa)

   totalstep = 0
   Vol = 0.d0

   do i = 1 , NJobs

     if(QMaster) call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     -----------------------------
#ifdef MOLFILE
       if(QMaster) call Read_RTraj(i)
#else
       if(QMaster) call Read_RTraj
#endif
!     -----------------------------

       if( mod(j,Interval(i)) == 0 ) then

         totalstep = totalstep + 1

         call InversMatrix(H,InvH)
         call PBC

         call CalcCOMCG(Numa)
         call ForceInterCG(Numa)

         Vol = Vol + det(H)

       end if

     end do

   end do

   if(QMaster) call Print_CGInteraction(Vol,totalstep)

end subroutine CoarseGrainC12


!######################################################################
!######################################################################


subroutine CGC12_data(Numa)

use Numbers, only : NumSpec, NumAtm, NumMol
use CommonBlocks, only : QMaster
use AtomParam , only : Mass
use CutoffParam, only : Rcutoff2
use CGAnalyz

implicit none

integer :: i, ii, Numa, Num
integer :: j, jj, k, MaxNumAtm_Group, nn
integer, dimension(:), allocatable :: Count
character(len=72) :: String1, String

   open(21,file='CGmodel.data',status='unknown')

! ## this number is used to identify the numbering of molecule to be analyzed 

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     select case(nn)

     case(1)
       read(String,*) Kcomp

! ## to get 'Numa'

       ii = 1

       do i = 1 , NumSpec

         if(i == Kcomp) then
           Numa = ii - 1
         end if

         ii = ii + NumMol(i) * NumAtm(i)

       end do

       NoAtom = NumAtm(Kcomp)
       NoMol  = NumMol(Kcomp)

       allocate( GroupNumber(NoAtom) ) ! Group Number ( atom number ) within a molecule 

! ## Number of Group in a molecule, Number of Group Types in a molecule 

     case(2)

       read(String,*) NumGroup, NumTypes

       allocate(NumAtm_Group(NumGroup)) ! ## Number of atoms in a defined group
       allocate(Mass_Group(NumGroup))   ! ## Mass of a defined group 
       allocate(Rcom_Group(3,NumGroup,NoMol)) ! ## Center of mass of a defined group

       allocate( NTypeGroup(NumGroup) )   ! ## Type number of a group in a molecule
       allocate( NinTypeGroup(NumTypes) ) ! ## Number of groups of a type

       exit

     end select

   end do

! ## read group number of an atom in a molecule

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) j, GroupNumber(j)

     if(nn == NoAtom) then
       if(j /= NoAtom) then
         write(*,*) 'ERROR : reading GroupNumber'
         call Finalize
       end if
       exit
     end if

   end do

! ## read type number of a group

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) j, NTypeGroup(j)

     if(nn == NumGroup) then
       if(j /= NumGroup) then
         write(*,*) 'ERROR : reading NTypeGroup'
         call Finalize
       end if
       exit
     end if

   end do

! ## Number of Bond in a molecule

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     read(String,*) NumCGBond

     exit

   end do

   allocate( BondCGI(NumCGBond) )
   allocate( BondCGJ(NumCGBond) )

! ## Bond pair

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) BondCGI(nn), BondCGJ(nn)

     if(nn == NumCGBond) exit

   end do

! ## Number of bending pairs in a molecule

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     read(String,*) NumCGAngle

     exit

   end do

   if(NumCGAngle/=0) then

   allocate( AngleCGI(NumCGAngle) )
   allocate( AngleCGJ(NumCGAngle) )
   allocate( AngleCGK(NumCGAngle) )

! ## Angle pair

   nn = 0

   do

     read(21,'(a72)') String1

     String = trim(adjustl(String1))

     if(String(1:1) == '#' .or. String(1:1) == '!') cycle

     nn = nn + 1

     read(String,*) AngleCGI(nn), AngleCGJ(nn), AngleCGK(nn)

     if(nn == NumCGAngle) exit

   end do

   end if

   close(21)

! ## Number of type pairs

   NumTypePair = ( NumTypes * ( NumTypes - 1 ) ) / 2 + NumTypes

   allocate( TypePair(NumTypes,NumTypes) )

   IRcut = int( sqrt(Rcutoff2) * ndr ) + 1

   allocate( EneNonBond(NumTypePair,IRcut) )
   allocate( FrcNonBond(NumTypePair,IRcut) )
   allocate( listGR(NumTypePair,IRcut) )

   EneNonBond = 0.d0
   FrcNonBond = 0.d0
   listGR = 0

! ## function TypePair yield the type-pair number 

   Num = 0
   do i = 1, NumTypes

     Num = Num + 1
     TypePair(i,i) = Num

   end do

   do i = 1, NumTypes - 1
     do j = i + 1, NumTypes

       Num = Num + 1
       TypePair(i,j) = Num
       TypePair(j,i) = Num

     end do
   end do

   if(Num /= NumTypePair) then
     write(*,*) 'ERROR : NumTypePair'
     call Finalize
   end if

! ## Number of atoms in a group, Mass of a group

   NumAtm_Group = 0
   Mass_Group = 0.d0

   do i = 1, NoAtom

     j = Numa + i
     k = GroupNumber(i)

     NumAtm_Group(k) = NumAtm_Group(k) + 1
     Mass_Group(k)   = Mass_Group(k)   + Mass(j)

   end do

! ## Maximum value of the number of atoms in a group

   MaxNumAtm_Group = 0

   do i = 1, NumGroup

     if(NumAtm_Group(i) > MaxNumAtm_Group) then
       MaxNumAtm_Group = NumAtm_Group(i)
     end if

   end do

! ## atom number of i-th atom in a group 

   allocate( NAtomGroup(MaxNumAtm_Group,NumGroup) )

   NumAtm_Group = 0

   do i = 1, NoAtom

     k = GroupNumber(i)

     NumAtm_Group(k) = NumAtm_Group(k) + 1
     NAtomGroup(NumAtm_Group(k),k) = i

   end do

! ## These CG pairs are not calculated by non-bond interaction 

   allocate( NonInteractPair(7,NumGroup) )

   allocate( Count(NumGroup) ) ! ## for temporary use

   Count = 0

   do i = 1 , NumGroup

     NonInteractPair(:,i) = i

   end do

   do i = 1, NumCGBond

     ii = BondCGI(i)
     jj = BondCGJ(i)

     Count(ii) = Count(ii) + 1
     Count(jj) = Count(jj) + 1

     NonInteractPair(Count(ii),ii) = jj
     NonInteractPair(Count(jj),jj) = ii

   end do

   do i = 1 , NumCGAngle

     ii = AngleCGI(i)
     jj = AngleCGK(i)

     Count(ii) = Count(ii) + 1
     Count(jj) = Count(jj) + 1

     NonInteractPair(Count(ii),ii) = jj
     NonInteractPair(Count(jj),jj) = ii

   end do

   if(QMaster) then

     do i = 1, NumGroup

       print *, i, NonInteractPair(:,i)

     end do

   end if

end subroutine CGC12_data


!######################################################################
!######################################################################


subroutine CalcCOMCG(Numa)

use Configuration, only : R
use AtomParam , only : Mass
use CGAnalyz, only : NoMol, NoAtom, Rcom_Group, Mass_Group, NumGroup, &
& GroupNumber

implicit none

integer :: k, ii, ll, jj, l, Numa

   do k = 1 , NoMol

     ii = Numa + (k-1)*NoAtom

! ##  Center of Mass

     Rcom_Group(:,:,k) = 0.d0

     do l = 1 , NoAtom

       ll = GroupNumber(l)
       jj = ii + l

       Rcom_Group(:,ll,k) = Rcom_Group(:,ll,k) + R(:,jj) * Mass(jj)

     end do

     do l = 1 , NumGroup

       Rcom_Group(:,l,k) = Rcom_Group(:,l,k) / Mass_Group(l)

     end do

   end do

end subroutine CalcCOMCG


!######################################################################
!######################################################################


! ********************************
! ** Intramolecular interaction **
! ********************************

subroutine ForceInterCG(Numa)

use Configuration, only : R
use CutoffParam, only : Ron2, Rcutoff2, swf1
use CommonBlocks, only : QSwitch
use CellParam, only : H, InvH
use NonbondParam, only : Rminh, EpsLJ, Charge
use CGAnalyz

implicit none

integer :: i, j, Numa, Natm
integer :: imol, igrp, ii, jj, jmol, jgrp
integer :: itype, jtype, kk, iatm, jatm
real(8), dimension(3,NumGroup,NoMol) :: ScGR
real(8), dimension(3) :: SGij, RGij, Rtra
integer, dimension(3) :: NGij, Rij
real(8) :: RG2, RG1, R2, Fsm, Frc, Ene
real(8) :: Eps, Sgm, Sgm2, InvR2, SR2, SR6, SR12, fkLJ, ek
real(8) :: cf, R1, fk1, fk, Ene_Ersp, Ene_LJ
integer :: IRG

   Natm = NoAtom

   do imol = 1, NoMol

     do igrp = 1, NumGroup

       ScGR(:,igrp,imol) = matmul( InvH, Rcom_Group(:,igrp,imol) )

     end do

   end do

! ## INTRAMOLECULE
!
!   do imol = 1, NoMol
!
!     ii = (imol-1) * Natm + Numa
!
!     do igrp = 1, NumGroup
!
!       itype = NTypeGroup(igrp)
!
!       do jgrp = 1, NumGroup
!
!         if((jgrp==NonInteractPair(1,igrp)).or. &
!         &  (jgrp==NonInteractPair(2,igrp)).or. &
!         &  (jgrp==NonInteractPair(3,igrp)).or. &
!         &  (jgrp==NonInteractPair(4,igrp)).or. &
!         &  (jgrp==NonInteractPair(5,igrp)).or. &
!         &  (jgrp==NonInteractPair(6,igrp)).or. &
!         &  (jgrp==NonInteractPair(7,igrp))) cycle
!
!         jtype = NTypeGroup(jgrp)
!
!         SGij = ScGR(:,igrp,imol) - ScGR(:,jgrp,jmol)
!         NGij = - nint( SGij )
!         SGij = SGij + NGij
!         Rtra = matmul( H, NGij )
!         RGij = matmul( H, SGij )
!
!         RG2 = dot_product(RGij, RGij)
!         RG1 = sqrt(RG2)
!         IRG = int( RG1 * ndr ) + 1
!
!         RGij = RGij / RG1
!
!         kk = TypePair(itype,jtype)
!
!         if(RG2 < Rcutoff2) then
!
!           listGR(kk,IRG) = listGR(kk,IRG) + 1
!
!           Ene = 0.d0
!           Fsm = 0.d0
!
!           do iatm = 1, NumAtm_Group(igrp)
!
!             i = ii + NAtomGroup(iatm,igrp)
!
!             do jatm = 1, NumAtm_Group(jgrp)
!
!               j = ii + NAtomGroup(jatm,jgrp)
!
!               Rij = R(:,i) - R(:,j) + Rtra
!               R2  = dot_product( Rij, Rij )
!
!               Sgm   = Rminh(i) + Rminh(j)
!               Sgm2  = Sgm * Sgm
!               Eps   = EpsLJ(i) * EpsLJ(j)
!               InvR2 = 1.d0 / R2
!
!               SR2  = Sgm2 * InvR2                    !(sigma/r)^2
!               SR6  = SR2 * SR2 * SR2                 !         ^6
!               SR12 = SR6 * SR6                       !         ^12
!               fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
!               ek   = Eps * ( SR12 - 2.d0 * SR6 )
!
!               if(QSwitch) then
!                 if ( R2 > Rcutoff2) then
!                   fkLJ = 0.d0
!                   ek   = 0.d0
!                 else
!! --------------------------------------------------------
!  call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
!! --------------------------------------------------------
!                 end if
!               end if
!
!               Ene_LJ = ek
!
!               cf = Charge(i) * Charge(j)
!               R1  = sqrt( R2 )
!               fk1 = cf / R1
!               fk  = fk1 * InvR2
!
!               Ene_Ersp = fk1
!               Frc = (fkLJ + fk) * dot_product( Rij, RGij )
!
!               Ene = Ene + Ene_LJ + Ene_Ersp
!               Fsm = Fsm + Frc
!
!             end do
!
!           end do
!
!           EneNonBond(kk,IRG) = EneNonBond(kk,IRG) + Ene
!           FrcNonBond(kk,IRG) = FrcNonBond(kk,IRG) + Fsm
!
!         end if
!
!       end do
!
!     end do
!
!   end do
!

! ## INTERMOLECULE

   do imol = 1, NoMol-1

     ii = (imol-1) * Natm + Numa

     do jmol = imol+1, NoMol

       jj = (jmol-1) * Natm + Numa

       do igrp = 1, NumGroup

         itype = NTypeGroup(igrp)

         do jgrp = 1, NumGroup

           jtype = NTypeGroup(jgrp)

           SGij = ScGR(:,igrp,imol) - ScGR(:,jgrp,jmol)
           NGij = - nint( SGij )
           SGij = SGij + NGij
           Rtra = matmul( H, NGij )
           RGij = matmul( H, SGij )

           RG2 = dot_product(RGij, RGij)
           RG1 = sqrt(RG2)
           IRG = int( RG1 * ndr ) + 1

           RGij = RGij / RG1

           kk = TypePair(itype,jtype)

           if(RG2 < Rcutoff2) then

             listGR(kk,IRG) = listGR(kk,IRG) + 1

             Ene = 0.d0
             Fsm = 0.d0

             do iatm = 1, NumAtm_Group(igrp)

               i = ii + NAtomGroup(iatm,igrp)

               do jatm = 1, NumAtm_Group(jgrp)

                 j = jj + NAtomGroup(jatm,jgrp)

                 Rij = R(:,i) - R(:,j) + Rtra
                 R2  = dot_product( Rij, Rij )

                 Sgm   = Rminh(i) + Rminh(j)
                 Sgm2  = Sgm * Sgm
                 Eps   = EpsLJ(i) * EpsLJ(j)
                 InvR2 = 1.d0 / R2

                 SR2  = Sgm2 * InvR2                    !(sigma/r)^2
                 SR6  = SR2 * SR2 * SR2                 !         ^6
                 SR12 = SR6 * SR6                       !         ^12
                 fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
                 ek   = Eps * ( SR12 - 2.d0 * SR6 )

                 if(QSwitch) then
                   if ( R2 > Rcutoff2) then
                     fkLJ = 0.d0
                     ek   = 0.d0
                   else
! --------------------------------------------------------
  call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! --------------------------------------------------------
                   end if
                 end if

                 Ene_LJ = ek

                 cf = Charge(i) * Charge(j)

                 R1  = sqrt( R2 )

                 fk1 = cf / R1
                 fk  = fk1 * InvR2

                 Ene_Ersp = fk1

                 Frc = (fkLJ + fk) * dot_product( Rij, RGij )

                 Ene = Ene + Ene_LJ + Ene_Ersp
                 Fsm = Fsm + Frc

               end do

             end do

             EneNonBond(kk,IRG) = EneNonBond(kk,IRG) + Ene
             FrcNonBond(kk,IRG) = FrcNonBond(kk,IRG) + Fsm

           end if

         end do

       end do

     end do

   end do


end subroutine ForceInterCG


!######################################################################
!######################################################################


subroutine Print_CGInteraction(Vol,totalstep)

use UnitExParam, only : pi, cvol
use CGAnalyz

implicit none

character(len=72) :: FileName
integer :: i, j, k, ii, jj, kk, NumSample, totalstep
real(8) :: xx, en, fc, Vol, drr, qk

   Vol = Vol / dble(totalstep)

   NinTypeGroup(:) = 0
   do i = 1, NumGroup
     j = NTypeGroup(i)
     NinTypeGroup(j) = NinTypeGroup(j) + 1
   end do

   do i = 1, NumTypes

     ii = NinTypeGroup(i)

     do j = i, NumTypes

       jj = NinTypeGroup(j)

       write(FileName,'(a,i1,a,i1,a)') './Analy/NonBond_Type',i,'--Type',j,'.dat'

       open(1,file = trim(FileName), status='unknown')

       kk = TypePair(i,j)

       NumSample = 0
       do k = 1 , IRcut
         NumSample = NumSample + listGR(kk,k)
       end do

       do k = 1 , IRcut

         drr = (k-0.5d0) * slabR
         qk = Vol / (NoMol * ii * 4.d0 * pi * drr * drr * slabR)

         if(listGR(kk,k) == 0) cycle
         xx = dble( listGR(kk,k) ) * qk / dble(totalstep*NoMol*jj)
         if(i==j) xx = xx * 2.d0
         en = EneNonBond(kk,k) / dble( listGR(kk,k) ) * cvol
         fc = FrcNonBond(kk,k) / dble( listGR(kk,k) ) * cvol

         write(1,'(f8.3,f15.10,2e15.7)') drr, xx, en, fc

       end do

       close(1)

     end do

   end do


end subroutine Print_CGInteraction
