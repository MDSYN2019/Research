! ############################
! ## SUBROUTINE LIST 
! ## -- Read_Charmm_Parameter 
! ## -- Read_Charmm_Topology 
! ## -- Read_PSF 
! ############################


!######################################################################
!######################################################################


! *****************************
! ** Charmm Parameter File   **
! *****************************

subroutine Read_Charmm_Parameter

use CommonBlocks, only : QMaster, Qstdout
use FFParameters
use IOparam, only : Parameter_file
use UnitExParam, only : pi, ExParam

implicit NONE

integer :: i, j, NDummy
character(len=90) :: String, String1
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4
character(len=90) :: Dummy
character(len=1), dimension(90) :: cline
integer :: ich, lentxt

! Parameters

open(1,file=trim(Parameter_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

   do

     read(1,'(90a)') String
     if(String(1:5)=='BONDS') exit

   end do

! ----------------------------------------------------------------------
! ## Bond parameters
!
!    V(bond) = Kb(b - b0)**2
!
!    Kb: kcal/mole/A**2
!    b0: A
!
!    atom type    Kb          b0

   NumBondParam=0

   do

     read(1,'(90a)') String1

     String = adjustl(String1)

     if(String(1:5)=='ANGLE') exit !until "ANGLE" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumBondParam = NumBondParam + 1

     read(String,*) BondPairAtoms(1,NumBondParam),  &
     &              BondPairAtoms(2,NumBondParam),  &
     &              kBondParam(NumBondParam),       &
     &              rBondParam(NumBondParam)

   end do

   if(QMaster.and.Qstdout) write(*,*) 'NumBondParam=',NumBondParam

! ----
! Unit
! ----
   do i = 1, NumBondParam

     kBondParam(i) = kBondParam(i) * ExParam

   end do

! ----------------------------------------------------------------------
! ## Angle paramters
!
!    V(angle) = Ktheta(Theta - Theta0)**2
!
!    V(Urey-Bradley) = Kub(S - S0)**2
!
!    Ktheta: kcal/mole/rad**2
!    Theta0: degrees
!    Kub: kcal/mole/A**2 (Urey-Bradley)
!    S0: A
!
!  atom types      Ktheta    Theta0    Kub      S0

   NumAngleParam=0

   do

     read(1,'(90a)') String1

     String = adjustl(String1)

     if(String(1:9)=='DIHEDRALS') exit !until "DIHEDRALS" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     do i = 1, 90
       cline(i) = String(i:i)
     end do

re1: do i = 1, 90
       if(cline(i) == Wk2) then
         lentxt = i-1
         exit re1
       end if
     end do re1

     ich = 1
     do i = 2, lentxt
       if((cline(i)/=Wk3).and.(cline(i-1)==Wk3)) then
         ich = ich + 1
       end if
     end do

     do i = 1, lentxt
       String1(i:i) = cline(i)
     end do
     do i = lentxt+1,90
       String1(i:i) = Wk3
     end do

     NumAngleParam = NumAngleParam + 1

     KubParam(NumAngleParam) = 0.d0
     S0Param(NumAngleParam)  = 0.d0

     if(ich == 5) then

       read(String1,*) (AnglePairAtoms(j,NumAngleParam),j=1,3), &
       &              kThetaParam(NumAngleParam),              &
       &              Theta0Param(NumAngleParam)

     else

       read(String1,*) (AnglePairAtoms(j,NumAngleParam),j=1,3),  &
       &               kThetaParam(NumAngleParam),              &
       &               Theta0Param(NumAngleParam),              &
       &               KubParam(NumAngleParam),                 &
       &               S0Param(NumAngleParam)

     end if

   end do

   if(QMaster.and.Qstdout) write(*,*) 'NumAngleParam=',NumAngleParam

! ----
! Unit
! ----

   do i = 1, NumAngleParam

     kThetaParam(i) = kThetaParam(i) * ExParam
     KubParam(i)    = KubParam(i)    * ExParam
     Theta0Param(i) = Theta0Param(i) / 180.d0 * pi

   end do

! ----------------------------------------------------------------------
! ## Dihedral parameters
!
!    V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!    Kchi: kcal/mole
!    n: multiplicity
!    delta: degrees
!
!    atom types             Kchi    n   delta

   NumDihedralParam=0

   do

     read(1,'(90a)') String1

     String = adjustl(String1)

     if(String(1:8)=='IMPROPER') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumDihedralParam = NumDihedralParam + 1

     read(String,*) (DihedralPairAtoms(j,NumDihedralParam),j=1,4), &
     &               kChiParam(NumDihedralParam),                  &
     &               NDihParam(NumDihedralParam),                  &
     &               DeltaDihParam(NumDihedralParam)

   end do

   if(QMaster.and.Qstdout) write(*,*) 'NumDihedralParam=',NumDihedralParam

! ----
! Unit
! ----

   do i = 1, NumDihedralParam

     kChiParam(i)     = kChiParam(i)     * ExParam
     DeltaDihParam(i) = DeltaDihParam(i) / 180.d0 * pi

   end do

! ----------------------------------------------------------------------
! ## Improper parameters
!
!    V(improper) = Kpsi(psi - psi0)**2
!
!    Kpsi: kcal/mole/rad**2
!    psi0: degrees
!    note that the second column of numbers (0) is ignored
!
!    atom types           Kpsi                   psi0

   NumImproperParam=0

   do

     read(1,'(90a)') String1

     String = adjustl(String1)

     if(String(1:7)=='NONBOND') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumImproperParam = NumImproperParam + 1

     read(String,*) (ImproperPairAtoms(j,NumImproperParam),j=1,4), &
     &               kPsiParam(NumImproperParam),NDummy,           &
     &               PsiImpParam(NumImproperParam)

   end do

   if(QMaster.and.Qstdout) write(*,*) 'NumImproperParam=',NumImproperParam

! ----
! Unit
! ----

   do i = 1, NumImproperParam

     kPsiParam(i)   = kPsiParam(i)   * ExParam
     PsiImpParam(i) = PsiImpParam(i) / 180.d0 * pi

   end do

! ----------------------------------------------------------------------
! ## LJ parameters
!
!    V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
!
!    epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
!    Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
!
!    atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
!

   NumLJParam = 0

   do

     read(1,'(90a)') String1

     String = adjustl(String1)

     if(String(1:5)=='cutnb') cycle
     if(String(1:5)=='CUTNB') cycle
     if(String(1:5)=='HBOND') exit
     if(String(1:3)=='END') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     do i = 1, 90
       cline(i) = String(i:i)
     end do

re2: do i = 1, 90
       if(cline(i) == Wk2) then
         lentxt = i-1
         exit re2
       end if
     end do re2

     ich = 1
     do i = 2, lentxt
       if((cline(i)/=Wk3).and.(cline(i-1)==Wk3)) then
         ich = ich + 1
       end if
     end do

     do i = 1, lentxt
       String1(i:i) = cline(i)
     end do
     do i = lentxt+1,90
       String1(i:i) = Wk3
     end do

     NumLJParam = NumLJParam + 1

     if(ich == 4) then

       read(String1,*) LJAtoms(NumLJParam),        &
       &              ignoredParam(NumLJParam),   &
       &              EpsLJParam(NumLJParam),     &
       &              RminhParam(NumLJParam)

       ignored14Param(NumLJParam)= ignoredParam(NumLJParam)
       Eps14Param(NumLJParam)    = EpsLJParam(NumLJParam)
       Rminh14Param(NumLJParam)  = RminhParam(NumLJParam)

     else

       read(String,*) LJAtoms(NumLJParam),       &
       &              ignoredParam(NumLJParam),  &
       &              EpsLJParam(NumLJParam),    &
       &              RminhParam(NumLJParam),    &
       &              ignored14Param(NumLJParam),&
       &              Eps14Param(NumLJParam),    &
       &              Rminh14Param(NumLJParam)

     end if

   end do

   if(QMaster.and.Qstdout) write(*,*) 'NumLJParam=',NumLJParam


! ----
! Unit
! ----

   do i = 1, NumLJParam

     EpsLJParam(i) = EpsLJParam(i) * ExParam
     Eps14Param(i) = Eps14Param(i) * ExParam

   end do

!do i = 1 , NumImproperParam
!write(51,'(4(a4,x),f10.4,i3,f10.4)') (ImproperPairAtoms(j,i),j=1,4),  &
!     &          kPsiParam(i),PsiImpParam(i)
!end do

! ----------------------------------------------------------------------

close(1)

end subroutine Read_Charmm_Parameter


!######################################################################
!######################################################################


! ****************************
! ** Charmm Topology File   **
! ****************************

subroutine Read_Charmm_Topology

use CommonBlocks, only : QMaster, Qstdout
use FFParameters
use IOparam, only : Topology_file

implicit NONE

character(len=4) :: Dummy1
character(len=6) :: Dummy2

real(8), dimension(200) :: ResiChargeParam

integer :: i,count
character(len=90) :: Ch, Ch1
character(len=1) :: Wk2
character(len=1) :: Wk3

! Topology

open(1,file=trim(Topology_file),status='old')

   Wk2 = '!'
   Wk3 = ' '

! ---------------------------------------------------------------------

   NumAtomTypeParam=0
   NumResidueParam=0

   do

     read(1,'(90a)') Ch1

     Ch = adjustl(Ch1)

     if(Ch(1:3)=='END') exit

     if(( Ch(1:5) == 'GROUP'    ) .or. ( Ch(1:5) == 'DONOR' ) .or.  &
     &  ( Ch(1:8) == 'ACCEPTOR' ) .or. ( Ch(1:1) == Wk3     ) .or.  &
     &  ( Ch(1:2) == 'IC'       ) .or. ( Ch(1:5) == 'Group' ) ) cycle

     if( Ch(1:4) == 'MASS' ) then

       NumAtomTypeParam = NumAtomTypeParam + 1
       read(Ch,*) Dummy1,i,AtomTypeParam(NumAtomTypeParam),&
       &          MassParam(NumAtomTypeParam)

     end if

!     if(Ch(1:4)=='DECL') then
!       NumAtomTypeParam = NumAtomTypeParam + 1
!       read(Ch,*) Dummy1,AtomTypeParam(NumAtomTypeParam)
!     end if

!     if(Ch(1:4)=='AUTO') exit

! ---------------------------------------------------------------------

     if( ( Ch(1:4) == 'RESI' ) .or. ( Ch(1:4) == 'PRES' ) .or. &
     &   ( Ch(1:4) == 'Resi' ) ) then

       NumResidueParam = NumResidueParam + 1

       read(Ch,*) Dummy1,ResiNameParam(NumResidueParam),&
       &                 ResiChargeParam(NumResidueParam)

       NumAtom_inResi(NumResidueParam) = 0
       NumBond_inResi(NumResidueParam) = 0
       NumDoub_inResi(NumResidueParam) = 0
       NumIMPR_inResi(NumResidueParam) = 0
       NumUNIT_inResi(NumResidueParam) = 0

       cycle

     end if

     if( Ch(1:4) == 'UNIT' ) then
       NumUNIT_inResi(NumResidueParam) = NumUNIT_inResi(NumResidueParam) + 1

       read(Ch,*) Dummy1, &
       &   UNITNameParam(NumUNIT_inResi(NumResidueParam),NumResidueParam), &
       &   NumAtom_inUNITParam(NumUNIT_inResi(NumResidueParam),NumResidueParam)

       cycle

     end if

     if( ( Ch(1:4) == 'ATOM' ).or.( Ch(1:4) == 'Atom' ).or.( Ch(1:4) == 'atom' ) ) then

       NumAtom_inResi(NumResidueParam) = NumAtom_inResi(NumResidueParam) + 1

       read(Ch,*) Dummy1,                                                        &
       &          AtomNameParam(NumAtom_inResi(NumResidueParam),NumResidueParam),&
       &          AtomType(NumAtom_inResi(NumResidueParam),NumResidueParam),     &
       &          ChargeParam(NumAtom_inResi(NumResidueParam),NumResidueParam)

       cycle

     end if

     if( ( Ch(1:4) == 'BOND' ).or.( Ch(1:4) == 'Bond' ) ) then

       count=0

       do i = 6 , 90

         if(  Ch(i:i) == Wk2 ) exit
         if(( Ch(i:i) /= Wk3 ) .and. ( Ch(i-1:i-1) == Wk3 ) ) count = count + 1

       end do

       if( mod(count,2) /= 0 ) then

         if(QMaster) write(*,*) 'ERROR: Bond',count,&
         & ResiNameParam(NumResidueParam),NumAtom_inResi(NumResidueParam)

         call Finalize

       end if

       NumBond_inResi(NumResidueParam) = NumBond_inResi(NumResidueParam) + count/2

       read(Ch,*) Dummy1, (BondPair_inResi(1,i,NumResidueParam),     &
       &                   BondPair_inResi(2,i,NumResidueParam),     &
       &          i = NumBond_inResi(NumResidueParam) - count/2 + 1, &
       &              NumBond_inResi(NumResidueParam) )

       cycle

     end if

     if( ( Ch(1:6) == 'DOUBLE' ).or.( Ch(1:6) == 'Double' ).or. &
     &   ( Ch(1:6) == 'TRIPLE' ).or.( Ch(1:6) == 'Triple' ) ) then

       count = 0

       do i = 8 , 90

         if(  Ch(i:i) == Wk2 ) exit

         if(( Ch(i:i) /= Wk3 ) .and. ( Ch(i-1:i-1) == Wk3 ) ) &
         &         count = count + 1

       end do

       if( mod(count,2) /= 0 ) then

         if(QMaster) write(*,*) 'ERROR: Double'

         call Finalize

       end if

       NumDoub_inResi(NumResidueParam) = NumDoub_inResi(NumResidueParam) + count / 2

       read(Ch,*) Dummy2, ( DoubPair_inResi(1,i,NumResidueParam),      &
       &                    DoubPair_inResi(2,i,NumResidueParam),      &
       &          i = NumDoub_inResi(NumResidueParam) - count / 2 + 1, &
       &              NumDoub_inResi(NumResidueParam) )

       cycle

     end if

     if( Ch(1:4) == 'IMPR' ) then

       count = 0

       do i = 6 , 90

         if(  Ch(i:i) == Wk2 ) exit

         if(( Ch(i:i) /= Wk3 ) .and. ( Ch(i-1:i-1) == Wk3 ) ) &
         &    count = count + 1

       end do

       if( mod(count,4) /= 0 ) then

         if(QMaster) write(*,*) 'ERROR: IMPR'

         call Finalize

       end if

       NumIMPR_inResi(NumResidueParam) = NumIMPR_inResi(NumResidueParam) + count / 4

       read(Ch,*) Dummy1,( ImprPair_inResi(1,i,NumResidueParam),        &
       &                   ImprPair_inResi(2,i,NumResidueParam),        &
       &                   ImprPair_inResi(3,i,NumResidueParam),        &
       &                   ImprPair_inResi(4,i,NumResidueParam),        &
       &          i = NumIMPR_inResi(NumResidueParam) - count / 4 + 1,  &
       &              NumIMPR_inResi(NumResidueParam) )

       cycle

     end if

   end do

! ---------------------------------------------------------------------

close(1)

   if(QMaster.and.Qstdout) write(*,*) 'NumAtomTypeParam = ',NumAtomTypeParam
   if(QMaster.and.Qstdout) write(*,*) 'NumResidueParam  = ',NumResidueParam

end subroutine Read_Charmm_Topology


!######################################################################
!######################################################################


! *********************
! **  Read PSF File  **
! *********************

subroutine Read_PSF

use Numbers, only : N, NumSpec, NumMol, NumAtm
use IOparam, only : PSF_file
use CommonBlocks, only : QMaster, QPBC, QSHAKE, QRigidBody, Qstdout, QCoulomb, &
&   QCyl
use FFParameters
use RBparam, only : NumRB, AtomUnitNum, AtomUnitName
use SHAKEparam, only : NSHAKEGroup
use UnitExParam, only : Avogadro, ec, kb
use NonbondParam, only : Charge, Rminh, EpsLJ, Eps14, Rminh14
use EwaldParam, only : Nel, Nelist, PCh
use BondedParam, only : NumBond, NumAngle, NumUB, NumDihedral, NumImproper, &
&   BondI, BondJ, kBond, rBond, AngleI, AngleJ, AngleK, kTheta, Theta0,     &
&   UB_I, UB_J, Kub, S0, DihedI, DihedJ, DihedK, DihedL, vdWSubtDih,        &
&   kChi, DeltaDih, NDih, DupFlag, DupkChi, DupDeltaDih, DupNDih, Ts,       &
&   ImproI, ImproJ, ImproK, ImproL, kPsi, PsiImp
use AtomParam, only : AtomName, ResidName, ResidNum, Mass, InvMass, TypeName, &
&   NAAType, AATypeList, NBAAType
use CylParam, only : CylFile, invdr_size, TabCyl, Rcylmin2, Rcylmax2, FCylRs, &
&   ngrid_cyl

implicit none

integer :: i, j, k, l, ii, jj, eofile
integer :: iatom, iunit, NResP, num
integer :: NTitle, NAtom
integer :: i1, i2, j1, j2
integer :: count
character(len=3) :: FileConfirm
character(len=4) :: TmpMolName
character(len=4) :: TmpResidName, TmpAtomName
character(len=6), dimension(4) :: Name, PName
character(len=4), dimension(4) :: TmpRName
character(len=80) :: String1, String
integer, dimension(:), allocatable :: CalUB, ParmUB
real :: cdip, cdip14
real(8) :: TotalCharge
integer :: iTotalC

integer, dimension(:), allocatable :: TmpBondI, TmpBondJ
integer, dimension(:), allocatable :: TmpAngleI, TmpAngleJ, TmpAngleK
integer, dimension(:), allocatable :: TmpDihedI, TmpDihedJ, TmpDihedK, TmpDihedL
integer, dimension(:), allocatable :: TmpImproI, TmpImproJ, TmpImproK, TmpImproL
logical, dimension(:), allocatable :: Qnewtype
character(len=6), dimension(:), allocatable :: Tmpname
character(len=80), dimension(:), allocatable :: TmpFile
real(8), dimension(:), allocatable :: TmpRmin, TmpRmax
real(8) :: x, xx, vp, fp, fr, gridsize, Rmx

open(1,file=trim(PSF_file),status='old')

   read(1,'(a3/)') FileConfirm

   if( FileConfirm /= 'PSF' ) then

#ifndef BMONI
     if(QMaster) write(11,*) 'PSF file error!'
#endif
     if(QMaster) write( 6,*) 'PSF file error!'

     call Finalize

   end if

   read(1,'(i8)') NTitle

   do i = 1, NTitle

     read(1,*)

   end do

   read(1,'(/i8)') NAtom

   if( NAtom /= N ) then  ! File Check

#ifndef BMONI
     if(QMaster) write(11,*) 'FILE ERROR! Number of Atom'
#endif
     if(QMaster) write( 6,*) 'FILE ERROR! Number of Atom'

     call Finalize

   end if

! --------------------------------------------------------------

   allocate( Charge  (N) )
   allocate( Mass    (N) )
   allocate( InvMass (N) )
   allocate( TypeName(N) )

   do i = 1 , N

     read(1,*)                    &
     &    ii, TmpMolName, jj, TmpResidName, TmpAtomName, TypeName(i), &
     &    Charge(i), Mass(i)

       if( TmpResidName /= ResidName(i) ) then  ! File Check

#ifndef BMONI
         if(QMaster) write(11,*) 'FILE ERROR! ResidName : Atom = ',i
#endif
         if(QMaster) write( 6,*) 'FILE ERROR! ResidName : Atom = ',i

         call Finalize

       end if

       if( TmpAtomName /= AtomName(i) ) then  ! File Check

#ifndef BMONI
         if(QMaster) write(11,*) 'FILE ERROR! AtomName : Atom = ',i
#endif
         if(QMaster) write( 6,*) 'FILE ERROR! AtomName : Atom = ',i

         call Finalize

       end if

   end do

   allocate( Qnewtype(N) )

   Qnewtype(:) = .False.
   Qnewtype(1) = .True.
   do i = 2, N
u1:  do j = 1, i-1
       if(TypeName(j)==TypeName(i)) exit u1
       if(j==(i-1)) Qnewtype(i) = .True.
     end do u1
   end do
   ii = 0
   do i = 1, N
     if(Qnewtype(i)) then
       ii = ii + 1
     end if
   end do

   NAAType = ii

   allocate( AATypeList(NAAType) )
   ii = 0
   do i = 1, N
     if(Qnewtype(i)) then
       ii = ii + 1
       AATypeList(ii) = TypeName(i)
     end if
   end do

   if(QMaster.and.Qstdout) then
     write(*,*) "Number of atomtypes = ",NAAType
     do i = 1, NAAType
       print *, 'Type ',i,' = ',AATypeList(i)
     end do
   end if

   allocate( NBAAType(N) )
   do i = 1, N
     do j = 1, NAAType
       if(AATypeList(j)==TypeName(i)) then
         NBAAType(i) = j
         exit
       end if
       if(j==NAAType) then
         write(*,*) "ERROR: no atomtype for ",TypeName(i)
         call Finalize
       end if
     end do
   end do

! ----
! Unit
! ----
   Mass = Mass * 1.d-3 / Avogadro

   InvMass = 1.d0 / Mass

   QCoulomb = .False.

   do i = 1, N
     if(Charge(i)/=0.) then
       QCoulomb = .True.
       exit
     end if
   end do

   if(QCoulomb) then
     TotalCharge = Sum( Charge )

     ItotalC = nint( TotalCharge * 1.d+5 )
     TotalCharge = dble(ItotalC) * 1.d-5

     Charge = Charge * sqrt(ec)

     Nel = 0
     do i = 1, N
       if(Charge(i)/=0.) Nel = Nel + 1
     end do

     if(Nel/=0) then
       allocate( Nelist(Nel) )
       allocate( PCh(Nel) )
       ii = 0
       do i = 1, N
         if(Charge(i)/=0.) then
           ii = ii + 1
           Nelist(ii) = i
           PCh(ii) = Charge(i)
         end if
       end do
     end if
   end if

   if(QMaster.and.Qstdout) then
     if(QCoulomb) then
       write(*,*) 'Checking charge neutrality: total charge =',TotalCharge
     else
       write(*,*) 'No Coulomb interaction'
     end if
   end if

! ----------------------------------------------------------------------
! ## Rigid Body
! NOTE! ## In this version, molecule that composed by just one rb-unit is
! NOTE! ## able to treated as a partial rigid-body molecule.

   if(QRigidBody) then

     allocate( AtomUnitNum(N) )
     allocate( AtomUnitName(N) )

     iatom = 0
     iunit = 0

     do i = 1 , NumSpec
       if(i == 1) then
         ii = 1
       else
         ii = ii + NumMol(i-1)*NumAtm(i-1)
       end if

       TmpRName(1) = ResidName(ii)

       do k = 1 , NumResidueParam

         if(ResiNameParam(k) == TmpRName(1)) then

           NResP = k

           exit

         end if

         if(k == NumResidueParam) then

           write(*,*) 'ERROR: residue param'
           write(*,*) 'residue = ',TmpRName(1)
           call Finalize

         end if

       end do

       do j = 1 , NumMol(i)

         do k = 1 , NumUnit_inResi(NResP)

           iunit = iunit + 1

           do l = 1 , NumAtom_inUNITParam(k,NResP)

             iatom = iatom + 1
             AtomUnitName(iatom) = UNITNameParam(k,NResP)
             AtomUnitNum (iatom) = iunit

           end do

         end do

       end do

     end do

     NumRB = AtomUnitNum(N)

   end if

! ----------------------------------------------------------------------
! ## Lennard-Jones parameter

   allocate( EpsLJ(N) )
   allocate( Rminh(N) )
   allocate( Eps14(N) )
   allocate( Rminh14(N) )

   do i = 1 , N

!   -----------------------
     Name(1) = TypeName(i)
!   -----------------------

     do j = 1 , NumLJParam

       if( Name(1) == LJAtoms(j) ) then

         cdip       = ignoredParam(j)
         if( EpsLJParam(j) /= 0. ) then
           EpsLJ(i) = sqrt( abs( EpsLJParam(j) ) )
         else
           EpsLJ(i) = 0.d0
         end if
         Rminh(i)   = RminhParam  (j)

         cdip14     = ignored14Param(j)
         if( Eps14Param(j) /= 0. ) then
           Eps14(i) = sqrt( abs( Eps14Param(j) ) )
         else
           Eps14(i) = 0.d0
         end if
         Rminh14(i) = Rminh14Param  (j)

         if( cdip  /= 0. ) then

           if(QMaster) write(*,*) 'induced dipole must be considered!'
           if(QMaster) write(*,*) i , Name(1)

           call Finalize

         end if

         if( cdip14 /= 0. ) then

           if(QMaster) write(*,*) 'induced dipole must be considered! 1-4'
           if(QMaster) write(*,*) i , Name(1)

           call Finalize

         end if

         exit

       end if

       if( j == NumLJParam ) then

         if(QMaster) write(*,*) 'No LJ parameter was found for',Name(1)
         if(QMaster) write(*,*) 'Atom Number =' , i

         call Finalize

       end if

     end do

   end do

! ----------------------------------------------------------------------
! ## Bond Parameters

   read(1,'(/i8)') NumBond

   if(NumBond /= 0) then

     allocate( BondI(NumBond) )
     allocate( BondJ(NumBond) )

     ii = NumBond / 4

     do i = 1 , ii

       k = (i-1) * 4
       read(1,'(8i8)') ( BondI(k+j) , BondJ(k+j) , j = 1 , 4 )

     end do

     ii = ii * 4
     jj = NumBond - ii

     if( jj /= 0 ) read(1,'(8i8)') ( BondI(ii+j) , BondJ(ii+j) , j = 1 , jj )

     if(QRigidBody) then

       allocate( TmpBondI(NumBond) )
       allocate( TmpBondJ(NumBond) )

       num = 0

       do ii = 1 , NumBond

         i = BondI(ii)
         j = BondJ(ii)

         if(AtomUnitNum(i)==AtomUnitNum(j)) cycle

         num = num + 1

         TmpBondI(num) = i
         TmpBondJ(num) = j

       end do

       deallocate( BondI, BondJ )

       NumBond = num

       allocate( BondI(NumBond) )
       allocate( BondJ(NumBond) )

       do ii = 1, NumBond
         BondI(ii) = TmpBondI(ii)
         BondJ(ii) = TmpBondJ(ii)
       end do

       deallocate( TmpBondI, TmpBondJ )

     end if

     allocate( kBond(NumBond) )
     allocate( rBond(NumBond) )

     do i = 1, NumBond

!   --------------------------------
       Name(1) = TypeName(BondI(i))
       Name(2) = TypeName(BondJ(i))
!   --------------------------------

       do j = 1, NumBondParam

         if(( (BondPairAtoms(1,j) == Name(1) )  &
      & .and. (BondPairAtoms(2,j) == Name(2)) ) &
      & .or.( (BondPairAtoms(1,j) == Name(2) )  &
      & .and. (BondPairAtoms(2,j) == Name(1)) )) then

           kBond(i) = kBondParam(j)
           rBond(i) = rBondParam(j)
           exit

         end if

         if( j == NumBondParam ) then

           if(QMaster) then

#ifndef BMONI
             write(11,*) 'lost parameter : Bond'
#endif
             write( 6,*) 'lost parameter : Bond'
             write( 6,*) 'Bond Number=',i
             write( 6,*) 'TypeNames=',Name(1),Name(2)

           end if

           call Finalize

         end if

       end do

     end do

     if(.not.QPBC) call NoLJList(1)

! ----------------------------------------
     if(QSHAKE) then

       call MakeSHAKEList

     else

       NSHAKEGroup = 0

     end if
! ----------------------------------------

   else if(.not.QPBC) then

     call NoLJList(0)

   end if

! ----------------------------------------------------------------------
! ## Angle & UB Parameters

   read(1,'(/i8)') NumAngle

   if(NumAngle /= 0) then

     allocate( AngleI(NumAngle) )
     allocate( AngleJ(NumAngle) )
     allocate( AngleK(NumAngle) )

     ii = NumAngle / 3

     do i = 1 , ii

       k = (i-1) * 3
       read(1,'(9i8)') ( AngleI(k+j) , AngleJ(k+j) , AngleK(k+j) , j = 1 , 3 )

     end do

     ii = ii * 3
     jj = NumAngle - ii

     if( jj /= 0 ) &
     &  read(1,'(8i8)') ( AngleI(ii+j) , AngleJ(ii+j) , AngleK(ii+j) , j = 1 , jj )

     if(QRigidBody) then

       allocate( TmpAngleI(NumAngle) )
       allocate( TmpAngleJ(NumAngle) )
       allocate( TmpAngleK(NumAngle) )

       num = 0

       do ii = 1 , NumAngle

         i = AngleI(ii)
         j = AngleJ(ii)
         k = AngleK(ii)

         if((AtomUnitNum(i)==AtomUnitNum(j)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(k))) cycle

         num = num + 1

         TmpAngleI(num) = i
         TmpAngleJ(num) = j
         TmpAngleK(num) = k

       end do

       deallocate( AngleI, AngleJ, AngleK )

       NumAngle = num

       allocate( AngleI(NumAngle) )
       allocate( AngleJ(NumAngle) )
       allocate( AngleK(NumAngle) )

       do ii = 1, NumAngle
         AngleI(ii) = TmpAngleI(ii)
         AngleJ(ii) = TmpAngleJ(ii)
         AngleK(ii) = TmpAngleK(ii)
       end do

       deallocate( TmpAngleI, TmpAngleJ, TmpAngleK )

     end if

     if(.not.QPBC) call NoLJList(2)

     if(QSHAKE) then

       allocate( TmpAngleI(NumAngle) )
       allocate( TmpAngleJ(NumAngle) )
       allocate( TmpAngleK(NumAngle) )

       num = 0

       do ii = 1 , NumAngle

         i = AngleI(ii)
         j = AngleJ(ii)
         k = AngleK(ii)

         if((ResidName(i) == 'TIP3').or. &
         &  (ResidName(i) == 'SPC' ).or. &
         &  (ResidName(i) == 'SPCE')) cycle

         num = num + 1

         TmpAngleI(num) = i
         TmpAngleJ(num) = j
         TmpAngleK(num) = k

       end do

       deallocate( AngleI, AngleJ, AngleK )

       NumAngle = num

       allocate( AngleI(NumAngle) )
       allocate( AngleJ(NumAngle) )
       allocate( AngleK(NumAngle) )

       do ii = 1, NumAngle
         AngleI(ii) = TmpAngleI(ii)
         AngleJ(ii) = TmpAngleJ(ii)
         AngleK(ii) = TmpAngleK(ii)
       end do

       deallocate( TmpAngleI, TmpAngleJ, TmpAngleK )

     end if

     allocate( kTheta(NumAngle) )
     allocate( Theta0(NumAngle) )

     allocate( CalUB(NumAngle) )
     allocate( ParmUB(NumAngle) )

     NumUB = 0

     do i = 1, NumAngle

!   ----------------------------------
       Name(1) = TypeName(AngleI(i))
       Name(2) = TypeName(AngleJ(i))
       Name(3) = TypeName(AngleK(i))
!   ----------------------------------

       do j = 1 , NumAngleParam

         if( AnglePairAtoms(2,j) == Name(2) ) then

           if(((AnglePairAtoms(1,j) == Name(1))  &
       &  .and.(AnglePairAtoms(3,j) == Name(3))) &
       &  .or.((AnglePairAtoms(1,j) == Name(3))  &
       &  .and.(AnglePairAtoms(3,j) == Name(1)))) then

             kTheta(i) = kThetaParam(j)
             Theta0(i) = Theta0Param(j)

             if( KubParam(j) /= 0.d0 ) then

               NumUB         = NumUB + 1
               CalUB(NumUB)  = i
               ParmUB(NumUB) = j

             end if

             exit

           end if

         end if

         if( j == NumAngleParam ) then

           if(QMaster) then

#ifndef BMONI
             write(11,*) 'lost parameter : Angle',i
#endif
             write( 6,*) 'lost parameter : Angle',i
             write( 6,*) AngleI(i),AngleJ(i),AngleK(i)
             write( 6,*) ( Name    (k) , k = 1 , 3 )

           end if

           call Finalize

         end if

       end do

     end do

     allocate( UB_I(NumUB) )
     allocate( UB_J(NumUB) )
     allocate( Kub(NumUB) )
     allocate( S0(NumUB) )

     do i = 1, NumUB

       j  = ParmUB(i)
       ii = CalUB(i)

       UB_I(i) = AngleI(ii)
       UB_J(i) = AngleK(ii)

       Kub(i) = KubParam(j)
       S0(i)  = S0Param(j)

     end do

   end if

! ----------------------------------------------------------------------
! ## Dihedral Parameters

   read(1,'(/i8)') NumDihedral

   if(NumDihedral /= 0) then

     allocate( DihedI(NumDihedral) )
     allocate( DihedJ(NumDihedral) )
     allocate( DihedK(NumDihedral) )
     allocate( DihedL(NumDihedral) )

     ii = NumDihedral / 2

     do i = 1 , ii

       k = (i-1) * 2
       read(1,'(8i8)') ( DihedI(k+j) , DihedJ(k+j) , &
       &                 DihedK(k+j) , DihedL(k+j) , j = 1 , 2 )

     end do

     ii = ii * 2
     jj = NumDihedral - ii

     if( jj /= 0 )                                     &
     & read(1,'(8i8)') ( DihedI(ii+j) , DihedJ(ii+j) , &
     &                   DihedK(ii+j) , DihedL(ii+j) , j = 1 , jj )

!-------------------------------------------
     if(QRigidBody) then

       allocate( TmpDihedI(NumDihedral) )
       allocate( TmpDihedJ(NumDihedral) )
       allocate( TmpDihedK(NumDihedral) )
       allocate( TmpDihedL(NumDihedral) )

       num = 0

       do ii = 1 , NumDihedral

         i = DihedI(ii)
         j = DihedJ(ii)
         k = DihedK(ii)
         l = DihedL(ii)

         if((AtomUnitNum(i)==AtomUnitNum(j)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(k)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(l))) cycle

         num = num + 1

         TmpDihedI(num) = i
         TmpDihedJ(num) = j
         TmpDihedK(num) = k
         TmpDihedL(num) = l

       end do

       deallocate( DihedI, DihedJ, DihedK, DihedL )

       NumImproper = num

       allocate( DihedI(NumDihedral) )
       allocate( DihedJ(NumDihedral) )
       allocate( DihedK(NumDihedral) )
       allocate( DihedL(NumDihedral) )

       do ii = 1 , NumImproper

         DihedI(ii) = TmpDihedI(ii)
         DihedJ(ii) = TmpDihedJ(ii)
         DihedK(ii) = TmpDihedK(ii)
         DihedL(ii) = TmpDihedL(ii)

       end do

       deallocate( TmpDihedI, TmpDihedJ, TmpDihedK, TmpDihedL )

     end if

     allocate( vdWSubtDih(NumDihedral) )

     allocate( kChi(NumDihedral) )
     allocate( DeltaDih(NumDihedral) )
     allocate( NDih(NumDihedral) )

     allocate( DupFlag(NumDihedral) )
     allocate( DupkChi(6,NumDihedral) )
     allocate( DupDeltaDih(6,NumDihedral) )
     allocate( DupNDih(6,NumDihedral) )

     allocate( Ts(NumDihedral) )

     vdWSubtDih = .False.

     DupFlag = 0

     do i = 1 , NumDihedral

!   ----------------------------------
       Name(1) = TypeName(DihedI(i))
       Name(2) = TypeName(DihedJ(i))
       Name(3) = TypeName(DihedK(i))
       Name(4) = TypeName(DihedL(i))
!   ----------------------------------

       count = 0

       do j = 1, NumDihedralParam

         PName(1) = DihedralPairAtoms(1,j)
         PName(2) = DihedralPairAtoms(2,j)
         PName(3) = DihedralPairAtoms(3,j)
         PName(4) = DihedralPairAtoms(4,j)

         if((PName(2) == Name(2)) .and. (PName(3) == Name(3)) .and. &
         &  (PName(1) == Name(1)) .and. (PName(4) == Name(4))) then

             count = count + 1

             if( count == 1 ) then

               kChi(i)     = kChiParam(j)
               DeltaDih(i) = DeltaDihParam(j)
               NDih(i)     = NdihParam(j)

             else

               DupFlag(i)                = count - 1
               DupkChi(DupFlag(i),i)     = kChiParam(j)
               DupDeltaDih(DupFlag(i),i) = DeltaDihParam(j)
               DupNDih(DupFlag(i),i)     = NdihParam(j)

             end if

         else if((PName(2) == Name(3)) .and. (PName(3)==Name(2)) .and.&
         &       (PName(1) == Name(4)) .and. (PName(4)==Name(1))) then

             count = count + 1

             if(count==1) then

               kChi(i)     = kChiParam(j)
               DeltaDih(i) = DeltaDihParam(j)
               NDih(i)     = NdihParam(j)

             else

               DupFlag(i)                = count - 1
               DupkChi(DupFlag(i),i)     = kChiParam(j)
               DupDeltaDih(DupFlag(i),i) = DeltaDihParam(j)
               DupNDih(DupFlag(i),i)     = NdihParam(j)

             end if

         end if

         if( (count == 0) .and. (j == NumDihedralParam) ) then

           do jj = 1 , NumDihedralParam

             PName(1) = DihedralPairAtoms(1,jj)
             PName(2) = DihedralPairAtoms(2,jj)
             PName(3) = DihedralPairAtoms(3,jj)
             PName(4) = DihedralPairAtoms(4,jj)

             if((PName(2)==Name(2)).and.(PName(3)==Name(3)).and.&
             &  (PName(1)=='X')    .and.(PName(4)=='X')   ) then

                 count = count + 1

                 if( count == 1 ) then

                   kChi(i)     = kChiParam(jj)
                   DeltaDih(i) = DeltaDihParam(jj)
                   NDih(i)     = NdihParam(jj)

                 else

                   DupFlag(i)                = count - 1
                   DupkChi(DupFlag(i),i)     = kChiParam(jj)
                   DupDeltaDih(DupFlag(i),i) = DeltaDihParam(jj)
                   DupNDih(DupFlag(i),i)     = NdihParam(jj)

                 end if

             else if((PName(2) == Name(3)) .and. (PName(3) == Name(2)) .and.&
             &       (PName(1) == 'X')     .and. (PName(4) == 'X')    ) then

                 count = count + 1

                 if( count == 1 ) then

                   kChi(i)     = kChiParam(jj)
                   DeltaDih(i) = DeltaDihParam(jj)
                   NDih(i)     = NdihParam(jj)

                 else

                   DupFlag(i)                = count - 1
                   DupkChi(DupFlag(i),i)     = kChiParam(jj)
                   DupDeltaDih(DupFlag(i),i) = DeltaDihParam(jj)
                   DupNDih(DupFlag(i),i)     = NdihParam(jj)

                 end if

             end if

             if( count > 7 ) then

               if(QMaster) then
                 write(*,*) 'too many dihedral parameters for the pair'
                 write(*,*)  ( Name    (l) , l = 1 , 4 )
               end if

               call Finalize

             end if

             if( ( count == 0 ) .and. ( jj == NumDihedralParam ) ) then

               if(QMaster) then

#ifndef BMONI
                 write(11,*) 'lost parameter : Dihedral'
#endif
                 write( 6,*) 'lost parameter : Dihedral'
                 write( 6,*) 'Dihedral Number=' , i
                 write( 6,*) 'Name    =' , ( Name    (l) , l = 1 , 4 )

               end if

               call Finalize

             end if

           end do

         end if

       end do

     end do

ex:  do i = 1 , NumDihedral

       i1 = DihedI(i)
       i2 = DihedL(i)

inn:   do j = 1 , NumAngle

         j1 = AngleI(j)
         j2 = AngleK(j)

         if( ( ( i1 == j1 ) .and. (i2 == j2 ) ) .or. &
         &   ( ( i1 == j2 ) .and. (i2 == j1 ) ) ) then

           vdWSubtDih(i) = .True.
!          if(QMaster) write(*,*) 'DUPLICATE 1-3, 1-4',i

           exit inn

         end if

       end do inn

       if(i == NumDihedral) exit ex

jnn:   do j = i+1 , NumDihedral

         if(j==i) cycle

         j1 = DihedI(j)
         j2 = DihedL(j)

         if( ( ( i1 == j1 ) .and. (i2 == j2 ) ) .or. &
         &   ( ( i1 == j2 ) .and. (i2 == j1 ) ) ) then

           vdWSubtDih(i) = .True.
!          if(QMaster) write(*,*) 'DUPLICATE 1-4',i
           exit jnn

         end if

       end do jnn

     end do ex

   end if

! ----------------------------------------------------------------------
! ## Improper torsion parameters

   read(1,'(/i8)') NumImproper

   if(NumImproper /= 0) then

     allocate( ImproI(NumImproper) )
     allocate( ImproJ(NumImproper) )
     allocate( ImproK(NumImproper) )
     allocate( ImproL(NumImproper) )

     ii = NumImproper / 2

     do i = 1 , ii

       k = (i-1) * 2
       read(1,'(8i8)') ( ImproI(k+j) , ImproJ(k+j) , &
       &                 ImproK(k+j) , ImproL(k+j) , j = 1 , 2 )

     end do

     ii = ii * 2
     jj = NumImproper - ii

     if( jj /= 0 )                                     &
     & read(1,'(8i8)') ( ImproI(ii+j) , ImproJ(ii+j) , &
     &                   ImproK(ii+j) , ImproL(ii+j) , j = 1 , jj )

     if(QRigidBody) then

       allocate( TmpImproI(NumImproper) )
       allocate( TmpImproJ(NumImproper) )
       allocate( TmpImproK(NumImproper) )
       allocate( TmpImproL(NumImproper) )

       num = 0

       do ii = 1 , NumImproper

         i = ImproI(ii)
         j = ImproJ(ii)
         k = ImproK(ii)
         l = ImproL(ii)

         if((AtomUnitNum(i)==AtomUnitNum(j)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(k)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(l))) cycle

         num = num + 1

         TmpImproI(num) = i
         TmpImproJ(num) = j
         TmpImproK(num) = k
         TmpImproL(num) = l

       end do

       deallocate( ImproI, ImproJ, ImproK, ImproL )

       NumImproper = num

       allocate( ImproI(NumImproper) )
       allocate( ImproJ(NumImproper) )
       allocate( ImproK(NumImproper) )
       allocate( ImproL(NumImproper) )

       do ii = 1 , NumImproper

         ImproI(ii) = TmpImproI(ii)
         ImproJ(ii) = TmpImproJ(ii)
         ImproK(ii) = TmpImproK(ii)
         ImproL(ii) = TmpImproL(ii)

       end do

       deallocate( TmpImproI, TmpImproJ, TmpImproK, TmpImproL )

     end if

     allocate( kPsi(NumImproper) )
     allocate( PsiImp(NumImproper) )

     do i = 1, NumImproper

!   ----------------------------------
       Name(1) = TypeName(ImproI(i))
       Name(2) = TypeName(ImproJ(i))
       Name(3) = TypeName(ImproK(i))
       Name(4) = TypeName(ImproL(i))
!   ----------------------------------

       do j = 1, NumImproperParam

         PName(1) = ImproperPairAtoms(1,j)
         PName(2) = ImproperPairAtoms(2,j)
         PName(3) = ImproperPairAtoms(3,j)
         PName(4) = ImproperPairAtoms(4,j)

         if((PName(1) == Name(1)) .and. (PName(4) == Name(4)) .and.&
         &  (PName(2) == Name(2)) .and. (PName(3) == Name(3))) then

             kPsi(i)   = kPsiParam(j)
             PsiImp(i) = PsiImpParam(j)
             exit

         else if((PName(1) == Name(4)) .and. (PName(4) == Name(1)) .and.&
         &       (PName(2) == Name(3)) .and. (PName(3) == Name(2))) then

             kPsi(i)   = kPsiParam(j)
             PsiImp(i) = PsiImpParam(j)
             exit

         end if

         if( j == NumImproperParam ) then

           do jj = 1 , NumImproperParam

             PName(1) = ImproperPairAtoms(1,jj)
             PName(2) = ImproperPairAtoms(2,jj)
             PName(3) = ImproperPairAtoms(3,jj)
             PName(4) = ImproperPairAtoms(4,jj)

             if(( PName(1) == Name(1) ) .and. ( PName(4) == Name(4) ) .and.&
             &  ( PName(2) == 'X'     ) .and. ( PName(3) == 'X'     ) ) then

                 kPsi(i)   = kPsiParam(jj)
                 PsiImp(i) = PsiImpParam(jj)
                 exit

             else if(( PName(1) == Name(4) ) .and. ( PName(4) == Name(1) ) .and.&
             &       ( PName(2) == 'X'     ) .and. ( PName(3) == 'X'     ) ) then

                 kPsi(i)   = kPsiParam(jj)
                 PsiImp(i) = PsiImpParam(jj)
                 exit

             end if

             if( jj == NumImproperParam ) then

               if(QMaster) then

#ifndef BMONI
                 write(11,*) 'lost parameter : Improper'
#endif
                 write( 6,*) 'lost parameter : Improper'
                 write( 6,*) 'Improper Number=',i
                 write( 6,*) 'Name=',(Name(l),l=1,4)

               end if

               call Finalize

             end if

           end do

         end if

       end do

     end do

   end if


   if(QMaster.and.Qstdout) then

     write(*,*) 'NumBOND=',NumBond
     write(*,*) 'NumANGLE=',NumAngle
     write(*,*) 'NumDIHED=',NumDihedral
     write(*,*) 'NumIMPRO=',NumImproper

   end if

close(1)

! ----------------------------------------------
! ## Cylinder potential
   if(QCyl) then

     open(19,file=trim(CylFile),status='old')

     allocate( Tmpname(100) )
     allocate( Tmpfile(100) )
     allocate( TmpRmin(100) )
     allocate( TmpRmax(100) )

     num = 0
     do
       read(19,'(a80)',iostat=eofile) String1
       if(eofile == -1) exit
       String = trim(adjustl(String1))
       if(String(1:1)=='#'.or.String(1:1)=='!') cycle
       num = num + 1
       read(String,*) Tmpname(num), TmpRmin(num), TmpRmax(num), Tmpfile(num)
     end do

     close(19)

     if(num<NAAType) then
       write(*,*) 'error : the number of wall files is fewer than needed'
       call Finalize
     end if

     allocate( invdr_size(NAAType) )
     allocate( TabCyl(ngrid_cyl,3,NAAType) )
     allocate( Rcylmin2(NAAType) )
     allocate( Rcylmax2(NAAType) )

     do i = 1, NAAType
irn:   do j = 1, num
         if(Tmpname(j)==AATypeList(i)) then
           open(99,file=trim(Tmpfile(j)),status='old')
           Rcylmin2(i) = TmpRmin(j)*TmpRmin(j)
           Rcylmax2(i) = TmpRmax(j)*TmpRmax(j)
           l = 0
           do k = 1, ngrid_cyl
             read(99,*) x, vp, fp, fr
             xx = x * x
             if(xx>(Rcylmin2(i)-1.d-10)) then
               l = l + 1
               TabCyl(l,1,i) = vp * kb !ExParam ! U
               TabCyl(l,2,i) = fp / x * kb !ExParam ! dU/dR / R
               TabCyl(l,3,i) = fr * kb !ExParam ! dU/dRadius
               if(l==1) then
                 Rcylmin2(i) = xx
               end if
               Rmx = xx
             end if
           end do
           if(Rmx < Rcylmax2(i)) then
             Rcylmax2(i) = Rmx
             write(*,*) 'WARNING: Rcylmax for Atomtype',i,' has been changed to ',&
             &          sqrt(Rcylmax2(i))
           end if
           gridsize = (Rmx - Rcylmin2(i)) / (dble(l) - 1.d0)
           invdr_size(i) = 1.d0 / gridsize
           exit irn
         end if
         if(j==num) then
           write(*,*) 'error : cannot find a cylinder file'
           call Finalize
         end if
       end do irn
     end do

     deallocate( Tmpname, Tmpfile, TmpRmin, TmpRmax )

     FCylRs = 0.d0

   end if

end subroutine Read_PSF
