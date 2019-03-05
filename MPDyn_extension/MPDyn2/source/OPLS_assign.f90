! ############################
! ## SUBROUTINE LIST 
! ## -- Read_OPLS_Parameter 
! ## -- Read_OPLS_Topology 
! ## -- MakeTopologicalData 
! ## -- AllocateOPLS 
! ############################


!######################################################################
!######################################################################


! ***************************
! ** OPLS Parameter File   **
! ***************************

subroutine Read_OPLS_Parameter

use CommonBlocks, only : QMaster, Qdebug
use FFParameters
use IOparam, only : Parameter_file
use UnitExParam, only : pi, ExParam

implicit NONE

integer :: i, j
character(len=72) :: String, String1
character(len=1) :: Wk1
character(len=1) :: Wk2
character(len=1) :: Wk3
character(len=1) :: Wk4

! Parameters

open(1,file=trim(Parameter_file),status='old')

   Wk1 = '*'
   Wk2 = '!'
   Wk3 = ' '
   Wk4 = '	'

   do

     read(1,'(72a)') String

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

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:5)=='ANGLE') exit !until "ANGLE" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumBondParam = NumBondParam + 1

     read(String,*) BondPairAtoms(1,NumBondParam), &
     &              BondPairAtoms(2,NumBondParam), &
     &              kBondParam(NumBondParam),      &
     &              rBondParam(NumBondParam)

   end do

   if(QMaster.and.Qdebug) write(*,*) 'NumBondParam=',NumBondParam

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
!    Ktheta: kcal/mole/rad**2
!    Theta0: degrees
!
!  atom types      Ktheta    Theta0

   NumAngleParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:9)=='DIHEDRALS') exit !until "DIHEDRALS" is found

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumAngleParam = NumAngleParam + 1

     read(String,*) (AnglePairAtoms(j,NumAngleParam),j=1,3),  &
     &               kThetaParam(NumAngleParam),              &
     &               Theta0Param(NumAngleParam)

   end do

   if(QMaster.and.Qdebug) write(*,*) 'NumAngleParam=',NumAngleParam

! ----
! Unit
! ----

   do i = 1, NumAngleParam

     kThetaParam(i) = kThetaParam(i) * ExParam
     Theta0Param(i) = Theta0Param(i) / 180.d0 * pi

   end do

! ----------------------------------------------------------------------
! ## Dihedral parameters
!
!    V(dihedral) = Kchi(1 + cos(n(chi) + delta))
!
!    Kchi: kcal/mole
!    n: multiplicity
!    delta: degrees
!
!    atom types             n       Kchi    delta

   NumDihedralParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:8)=='IMPROPER') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumDihedralParam = NumDihedralParam + 1

     read(String,*) (DihedralPairAtoms(j,NumDihedralParam),j=1,4), &
     &          NDihParam(NumDihedralParam),                       &
     &          kChiParam(NumDihedralParam),                       &
     &          DeltaDihParam(NumDihedralParam)

   end do

   if(QMaster.and.Qdebug) write(*,*) 'NumDihedralParam=',NumDihedralParam

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
!    V(improper) = Kpsi(1 + cos(n(psi) - delta))
!
!    Kchi: kcal/mole
!    n: multiplicity
!    delta: degrees
!
!    atom types             n       Kpsi    delta

   NumImproperParam=0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:7)=='NONBOND') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumImproperParam = NumImproperParam + 1

     read(String,*) (ImproperPairAtoms(j,NumImproperParam),j=1,4),  &
     &               NImpParam(NumImproperParam),                    &
     &               kImpParam(NumImproperParam),                    &
     &               DeltaImpParam(NumImproperParam)

   end do

   if(QMaster.and.Qdebug) write(*,*) 'NumImproperParam=',NumImproperParam

! ----
! Unit
! ----

   do i = 1, NumImproperParam

     kImpParam(i)     = kImpParam(i)   * ExParam
     DeltaImpParam(i) = DeltaImpParam(i) / 180.d0 * pi

   end do

! ----------------------------------------------------------------------
! ## LJ parameters
!
!    V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
!
!    epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
!    sigma: A, Sigma,i,j = sqrt(sigma,i * sigma,j)
!
!    atom   sigma   epsilon
!

   NumLJParam = 0

   do

     read(1,'(72a)') String1

     String = adjustl(String1)

     if(String(1:3)=='END') exit

     if( (String(1:1)==Wk1).or.(String(1:1)==Wk2).or. &
     &   (String(1:1)==Wk3).or.(String(1:1)==Wk4) ) cycle

     NumLJParam = NumLJParam + 1

     read(String,*) LJAtoms(NumLJParam),      &
     &              SigmaLJParam(NumLJParam), &
     &              EpsLJParam(NumLJParam)

   end do

   if(QMaster.and.Qdebug) write(*,*) 'NumLJParam=',NumLJParam

! ----
! Unit
! ----

   do i = 1, NumLJParam

     EpsLJParam(i) = EpsLJParam(i) * ExParam

   end do

! ----------------------------------------------------------------------

close(1)

end subroutine Read_OPLS_Parameter


!######################################################################
!######################################################################


! **************************
! ** OPLS Topology File   **
! **************************

subroutine Read_OPLS_Topology

use CommonBlocks, only : QMaster, Qdebug
use FFParameters
use IOparam, only : Topology_file

implicit NONE

character(len=4) :: Dummy1

real(8), dimension(200) :: ResiChargeParam

integer :: i,count
character(len=80) :: String
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

     read(1,'(80a)') String

     if(String(1:3)=='END') exit

     if(( String(1:5) == 'GROUP' ) .or. ( String(1:1) == Wk3 ) ) cycle

     if( String(1:4) == 'MASS' ) then

       NumAtomTypeParam = NumAtomTypeParam + 1
       read(String,*) Dummy1,i,AtomTypeParam(NumAtomTypeParam),&
       &          MassParam(NumAtomTypeParam)

     end if

! ---------------------------------------------------------------------

     if( ( String(1:4) == 'RESI' ) .or. ( String(1:4) == 'PRES' ) ) then

       NumResidueParam = NumResidueParam + 1

       read(String,*) Dummy1,ResiNameParam(NumResidueParam),&
       &                 ResiChargeParam(NumResidueParam)

       NumAtom_inResi(NumResidueParam) = 0
       NumBond_inResi(NumResidueParam) = 0
       NumIMPR_inResi(NumResidueParam) = 0
       NumUNIT_inResi(NumResidueParam) = 0

       cycle

     end if

     if( String(1:4) == 'UNIT' ) then
       NumUNIT_inResi(NumResidueParam) = NumUNIT_inResi(NumResidueParam) + 1

       read(String,*) Dummy1, &
       &   UNITNameParam(NumUNIT_inResi(NumResidueParam),NumResidueParam), &
       &   NumAtom_inUNITParam(NumUNIT_inResi(NumResidueParam),NumResidueParam)

       cycle

     end if

     if( String(1:4) == 'ATOM' ) then

       NumAtom_inResi(NumResidueParam) = NumAtom_inResi(NumResidueParam) + 1

       read(String,*) Dummy1,                                                    &
       &          AtomNameParam(NumAtom_inResi(NumResidueParam),NumResidueParam),&
       &          AtomType(NumAtom_inResi(NumResidueParam),NumResidueParam),     &
       &          ChargeParam(NumAtom_inResi(NumResidueParam),NumResidueParam)

       cycle

     end if

     if( String(1:4) == 'BOND' ) then

       count=0

       do i = 6 , 80

         if(  String(i:i) == Wk2 ) exit
         if(( String(i:i) /= Wk3 ) .and. ( String(i-1:i-1) == Wk3 ) ) count = count + 1

       end do

       if( mod(count,2) /= 0 ) then

         if(QMaster) write(*,*) 'ERROR: Bond',count,&
         & ResiNameParam(NumResidueParam),NumAtom_inResi(NumResidueParam)

         call Finalize

       end if

       NumBond_inResi(NumResidueParam) = NumBond_inResi(NumResidueParam) + count/2

       read(String,*) Dummy1, (BondPair_inResi(1,i,NumResidueParam),     &
       &                       BondPair_inResi(2,i,NumResidueParam),     &
       &              i = NumBond_inResi(NumResidueParam) - count/2 + 1, &
       &                  NumBond_inResi(NumResidueParam) )

       cycle

     end if

     if( String(1:4) == 'IMPR' ) then

       count = 0

       do i = 6 , 80

         if(  String(i:i) == Wk2 ) exit

         if(( String(i:i) /= Wk3 ) .and. ( String(i-1:i-1) == Wk3 ) ) &
         &    count = count + 1

       end do

       if( mod(count,4) /= 0 ) then

         if(QMaster) write(*,*) 'ERROR: IMPR'

         call Finalize

       end if

       NumIMPR_inResi(NumResidueParam) = NumIMPR_inResi(NumResidueParam) + count / 4

       read(String,*) Dummy1,( ImprPair_inResi(1,i,NumResidueParam),        &
       &                       ImprPair_inResi(2,i,NumResidueParam),        &
       &                       ImprPair_inResi(3,i,NumResidueParam),        &
       &                       ImprPair_inResi(4,i,NumResidueParam),        &
       &              i = NumIMPR_inResi(NumResidueParam) - count / 4 + 1,  &
       &                  NumIMPR_inResi(NumResidueParam) )

       cycle

     end if

   end do

! ---------------------------------------------------------------------

close(1)

   if(QMaster.and.Qdebug) then
     write(*,*) 'NumAtomTypeParam = ',NumAtomTypeParam
     write(*,*) 'NumResidueParam  = ',NumResidueParam
     do i = 1, NumResidueParam
       print *, ResiNameParam(i)
     end do
   end if

end subroutine Read_OPLS_Topology


!######################################################################
!######################################################################


subroutine MakeTopologicalData

use Numbers, only : N, NumSpec, NumMol, NumAtm, NumMer
use CommonBlocks, only : QMaster, Qstdout, QPBC, PolyFlag, QRigidBody, QMacro, &
&   QCoulomb
use FFParameters
use RBparam, only : NumRB, AtomUnitNum, AtomUnitName
use UnitExParam, only : ec
use NonbondParam, only : Charge
use EwaldParam, only : Nel, Nelist, PCh
use BondedParam, only : NumBond, NumAngle, NumDihedral, NumImproper, &
&   BondI, BondJ, AngleI, AngleJ, AngleK, DihedI, DihedJ, DihedK, DihedL, &
&   ImproI, ImproJ, ImproK, ImproL
use AtomParam, only : MolName, AtomName

implicit NONE

integer, parameter :: tempdimension = 500000
integer :: i, j, k, l
integer :: NumAtomMol, NumBondMol, NumIMPRMol, NumUnitMol
integer :: Nbond, Nimpr, NResP, NumUnitAcc, Nadd
real(8), dimension(tempdimension) :: ChargeMol
integer, dimension(2,tempdimension) :: BondPairMol
integer, dimension(4,tempdimension) :: ImprPairMol
integer, dimension(tempdimension) :: TmpBondI
integer, dimension(tempdimension) :: TmpBondJ
integer, dimension(tempdimension) :: TmpAngleI
integer, dimension(tempdimension) :: TmpAngleJ
integer, dimension(tempdimension) :: TmpAngleK
integer, dimension(tempdimension) :: TmpDihedI
integer, dimension(tempdimension) :: TmpDihedJ
integer, dimension(tempdimension) :: TmpDihedK
integer, dimension(tempdimension) :: TmpDihedL
integer, dimension(tempdimension) :: TmpImproI
integer, dimension(tempdimension) :: TmpImproJ
integer, dimension(tempdimension) :: TmpImproK
integer, dimension(tempdimension) :: TmpImproL

character(len=8) :: String, String1
character(len=4) :: RName, RNameT
character(len=4) :: AtomI, AtomJ, AtomK, AtomL, NameA

integer :: ISearch, FSearch
integer :: IJSearch, FJSearch
integer :: ISSearch, FSSearch
integer :: IPSearch, FPSearch
integer :: IMSearch, FMSearch
integer :: JAtom, IAtom, NAtom
logical :: ConfAllocI, ConfAllocJ, ConfAllocK, ConfAllocL

integer, dimension(N) :: AtomHands
integer, dimension(N,6) :: AtomHandBond
integer :: Ipair, Jpair
integer :: i1, i2, j1, j2, NbI, NbJ, ii, jj
character(len=6), dimension(N) :: UName
integer, dimension(N) :: UNum

real(8) :: TotalCharge
integer :: ItotalC
integer :: num

   NumBond = 0
   NumImproper = 0
   NumUnitAcc = 0

   allocate( Charge(N) )
   allocate( AtomUnitNum(N) )
   allocate( AtomUnitName(N) )

   do i = 1 , NumSpec

     if(i == 1) then
       ISearch = 0
     else
       ISearch = ISearch + NumMol(i-1)*NumAtm(i-1)
       NumUnitAcc = NumUnitAcc + NumMol(i-1)*NumUnitMol
     end if

     FSearch = ISearch + NumAtm(i)

     if(PolyFlag(i)) then
       String1 = MolName(i)(5:8)
     else
       String1 = MolName(i)
     end if

     String = trim(adjustl(String1))

     RName = String(1:4)

     NumAtomMol = 0
     NumBondMol = 0
     NumIMPRMol = 0
     NumUnitMol = 0

     if(QMacro) then
       if(MolName(i)=='CGSP') then

         NumAtomMol = 1
         NumUnitMol = 1

         NAtom = ISearch
         do j = 1 , NumMol(i)
           Nadd = NumAtm(i) * (j-1)
           do k = 1 , NumAtomMol
             Charge(NAtom+k)       = 0.d0
             AtomUnitName(NAtom+k) = 'SNGL'
             AtomUnitNum(NAtom+k)  = 1 + NumUnitAcc + (j-1)*NumUnitMol
           end do
           NAtom = NAtom + NumAtm(i)
         end do

         cycle

       end if
     end if

! ----------------------------------------------------------------------
! ## polymers
! ----------------------------------------------------------------------

     if(PolyFlag(i)) then

       JAtom = 0
       IJSearch = ISearch

       do j = 1 , NumMer(i)

         IJSearch = IJSearch + JAtom

         if( j == 1 ) then

           write(RName,'(a,a)') trim(String),'I'

         else if( j == NumMer(i) ) then

           write(RName,'(a,a)') trim(String),'F'

         else

           write(RName,'(a)') trim(String)

         end if

! ##################

         do k = 1 , NumResidueParam

           if(ResiNameParam(k) == RName) then

             Nbond = NumBond_inResi(k)
             Nimpr = NumIMPR_inResi(k)
             IAtom = NumAtom_inResi(k)
             NResP = k

             exit

           end if

           if(k == NumResidueParam) then

             write(*,*) 'ERROR: polymer resid.'
             write(*,*) 'residue = ',RName
             call Finalize

           end if

         end do

         FJSearch = IJSearch + IAtom

! ## Charge ##

         do l = IJSearch + 1, FJSearch

           NumAtomMol = NumAtomMol + 1

           AtomI = AtomName(l)

           ConfAllocI = .False.

           do k = 1 , IAtom

             if(AtomI == AtomNameParam(k,NResP)) then

               ChargeMol(NumAtomMol) = ChargeParam(k,NResP)
               ConfAllocI = .True.

             end if

             if((k == IAtom).and.(.not.ConfAllocI)) then

               write(*,*) 'ERROR : charge allocation'
               call Finalize

             end if

           end do

         end do

! ## Unit Name ##

         NumAtomMol = NumAtomMol - IAtom

         do k = 1 , NumUnit_inResi(NResP)

           NumUnitMol = NumUnitMol + 1

           do l = 1 , NumAtom_inUNITParam(k,NResP)

             NumAtomMol = NumAtomMol + 1
             UName(NumAtomMol) = UNITNameParam(k,NResP)
             UNum(NumAtomMol)  = NumUnitMol

           end do

         end do

! ## Atoms in the neighbour residues
! ## IPSearch, FPSearch
! ## IMSearch, FMSearch

         IPSearch = IJSearch + IAtom

         if(j==(NumMer(i)-1)) then

           write(RNameT,'(a,a)') trim(RName),'F'

           do l = 1 , NumResidueParam

             if(ResiNameParam(l) == RNameT) then

               NAtom = NumAtom_inResi(l)

               exit

             end if

             if(l == NumResidueParam) then

               write(*,*) 'ERROR: polymer resid.'
               write(*,*) 'residue = ',RNameT
               call Finalize

             end if

           end do

         else if(j==1) then

           l = len( trim(RName) )
           write(RNameT,'(a)') RName(1:l-1)

           do l = 1 , NumResidueParam

             if(ResiNameParam(l) == RNameT) then

               NAtom = NumAtom_inResi(l)

               exit

             end if

             if(l == NumResidueParam) then

               write(*,*) 'ERROR: polymer resid.'
               write(*,*) 'residue = ',RNameT
               call Finalize

             end if

           end do

         else

           NAtom = IAtom

         end if

         FPSearch = FJSearch + NAtom

         IMSearch = IJSearch - JAtom
         FMSearch = FJSearch - IAtom


! ## Bond ##

         do k = 1 , Nbond

           NumBondMol = NumBondMol + 1

           AtomI = BondPair_inResi(1,k,NResP)
           AtomJ = BondPair_inResi(2,k,NResP)

           ConfAllocI = .False.
           ConfAllocJ = .False.

           if((AtomI(1:1)=='+').or.(AtomJ(1:1)=='+')) then

             if(j==NumMer(i)) then

               write(*,*) 'ERROR : bond atoms in the topology file.'
               call Finalize

             end if

           else if( (AtomI(1:1) == '-').or.(AtomJ(1:1) == '-') ) then

             if( j == 1 ) then

               write(*,*) 'ERROR : bond atoms in the topology file.'
               call Finalize

             end if

           end if

           if( (AtomI(1:1) == '+').or.(AtomI(1:1) == '-') ) then

             NameA = AtomI(2:)

             if( AtomI(1:1) == '+' ) then
               ISSearch = IPSearch
               FSSearch = FPSearch
             else
               ISSearch = IMSearch
               FSSearch = FMSearch
             end if

             do l = ISSearch+1, FSSearch

               if(AtomName(l) == NameA) then

                 BondPairMol(1,NumBondMol) = l
                 ConfAllocI = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomJ) then

                 BondPairMol(2,NumBondMol) = l
                 ConfAllocJ = .True.

                 exit

               end if

             end do

           else if((AtomJ(1:1)=='+').or.(AtomJ(1:1)=='-')) then

             NameA = AtomJ(2:)

             if( AtomJ(1:1) == '+' ) then
               ISSearch = IPSearch
               FSSearch = FPSearch
             else
               ISSearch = IMSearch
               FSSearch = FMSearch
             end if

             do l = ISSearch+1, FSSearch

               if(AtomName(l) == NameA) then

                 BondPairMol(2,NumBondMol) = l
                 ConfAllocJ = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomI) then

                 BondPairMol(1,NumBondMol) = l
                 ConfAllocI = .True.

                 exit

               end if

             end do

           else

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomI) then

                 BondPairMol(1,NumBondMol) = l
                 ConfAllocI = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomJ) then

                 BondPairMol(2,NumBondMol) = l
                 ConfAllocJ = .True.

                 exit

               end if

             end do

           end if

           if((.not.ConfAllocI).or.(.not.ConfAllocJ)) then

             write(*,*) 'ERROR : missing bond'
             write(*,*) AtomI,AtomJ,' in residue ',RName
             call Finalize

           end if

         end do

! ## Improper ##

         do k = 1 , Nimpr

           NumImprMol = NumImprMol + 1

           AtomI = ImprPair_inResi(1,k,NResP)
           AtomJ = ImprPair_inResi(2,k,NResP)
           AtomK = ImprPair_inResi(3,k,NResP)
           AtomL = ImprPair_inResi(4,k,NResP)

           ConfAllocI = .False.
           ConfAllocJ = .False.
           ConfAllocK = .False.
           ConfAllocL = .False.

           if((AtomI(1:1)=='+').or.(AtomJ(1:1)=='+').or. &
           &  (AtomK(1:1)=='+').or.(AtomL(1:1)=='+')) then

             if(j==NumMer(i)) then

               write(*,*) 'ERROR : bond atoms in the topology file.'
               call Finalize

             end if

           else if((AtomI(1:1) == '-').or.(AtomJ(1:1) == '-').or. &
           &       (AtomK(1:1) == '-').or.(AtomL(1:1) == '-')) then

             if( j == 1 ) then

               write(*,*) 'ERROR : bond atoms in the topology file.'
               call Finalize

             end if

           end if

           if((AtomI(1:1) == '+').or.(AtomI(1:1) == '-')) then

             NameA = AtomI(2:)

             if( AtomI(1:1) == '+' ) then
               ISSearch = IPSearch
               FSSearch = FPSearch
             else
               ISSearch = IMSearch
               FSSearch = FMSearch
             end if

             do l = ISSearch+1, FSSearch

               if(AtomName(l) == NameA) then

                 ImprPairMol(1,NumImprMol) = l
                 ConfAllocI = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomJ) then

                 ImprPairMol(2,NumImprMol) = l
                 ConfAllocJ = .True.

               else if(AtomName(l) == AtomK) then

                 ImprPairMol(3,NumImprMol) = l
                 ConfAllocK = .True.

               else if(AtomName(l) == AtomL) then

                 ImprPairMol(4,NumImprMol) = l
                 ConfAllocL = .True.

               end if

             end do

           else if((AtomJ(1:1)=='+').or.(AtomJ(1:1)=='-')) then

             NameA = AtomJ(2:)

             if( AtomJ(1:1) == '+' ) then
               ISSearch = IPSearch
               FSSearch = FPSearch
             else
               ISSearch = IMSearch
               FSSearch = FMSearch
             end if

             do l = ISSearch+1, FSSearch

               if(AtomName(l) == NameA) then

                 ImprPairMol(2,NumImprMol) = l
                 ConfAllocJ = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomI) then

                 ImprPairMol(1,NumImprMol) = l
                 ConfAllocI = .True.

               else if(AtomName(l) == AtomK) then

                 ImprPairMol(3,NumImprMol) = l
                 ConfAllocK = .True.

               else if(AtomName(l) == AtomL) then

                 ImprPairMol(4,NumImprMol) = l
                 ConfAllocL = .True.

               end if

             end do

           else if((AtomK(1:1)=='+').or.(AtomK(1:1)=='-')) then

             if( AtomK(1:1) == '+' ) then
               ISSearch = IPSearch
               FSSearch = FPSearch
             else
               ISSearch = IMSearch
               FSSearch = FMSearch
             end if

             NameA = AtomK(2:)

             do l = ISSearch+1, FSSearch

               if(AtomName(l) == NameA) then

                 ImprPairMol(3,NumImprMol) = l
                 ConfAllocK = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomI) then

                 ImprPairMol(1,NumImprMol) = l
                 ConfAllocI = .True.

               else if(AtomName(l) == AtomJ) then

                 ImprPairMol(2,NumImprMol) = l
                 ConfAllocJ = .True.

               else if(AtomName(l) == AtomL) then

                 ImprPairMol(4,NumImprMol) = l
                 ConfAllocL = .True.

               end if

             end do

           else if((AtomL(1:1)=='+').or.(AtomL(1:1)=='-')) then

             NameA = AtomL(2:)

             if( AtomL(1:1) == '+' ) then
               ISSearch = IPSearch
               FSSearch = FPSearch
             else
               ISSearch = IMSearch
               FSSearch = FMSearch
             end if

             do l = ISSearch+1, FSSearch

               if(AtomName(l) == NameA) then

                 ImprPairMol(4,NumImprMol) = l
                 ConfAllocL = .True.

                 exit

               end if

             end do

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomI) then

                 ImprPairMol(1,NumImprMol) = l
                 ConfAllocI = .True.

               else if(AtomName(l) == AtomJ) then

                 ImprPairMol(2,NumImprMol) = l
                 ConfAllocJ = .True.

               else if(AtomName(l) == AtomK) then

                 ImprPairMol(3,NumImprMol) = l
                 ConfAllocK = .True.

               end if

             end do

           else

             do l = IJSearch+1, FJSearch

               if(AtomName(l) == AtomI) then

                 ImprPairMol(1,NumImprMol) = l
                 ConfAllocI = .True.

               else if(AtomName(l) == AtomJ) then

                 ImprPairMol(2,NumImprMol) = l
                 ConfAllocJ = .True.

               else if(AtomName(l) == AtomK) then

                 ImprPairMol(3,NumImprMol) = l
                 ConfAllocK = .True.

               else if(AtomName(l) == AtomL) then

                 ImprPairMol(4,NumImprMol) = l
                 ConfAllocL = .True.

               end if

             end do

           end if

           if((.not.ConfAllocI).or.(.not.ConfAllocJ).or. &
           &  (.not.ConfAllocK).or.(.not.ConfAllocL)) then

             write(*,*) 'ERROR : missing impr'
             write(*,*) AtomI,AtomJ,AtomK,AtomL,' in residue ',RName
             call Finalize

           end if

         end do

         JAtom = IAtom

       end do

     else

! ----------------------------------------------------------------------
! ## non-polymer
! ----------------------------------------------------------------------

       do k = 1 , NumResidueParam

         if(ResiNameParam(k) == RName) then

           Nbond = NumBond_inResi(k)
           Nimpr = NumIMPR_inResi(k)
           IAtom = NumAtom_inResi(k)
           NResP = k

           exit

         end if

         if(k == NumResidueParam) then

           write(*,*) 'ERROR: polymer resid.'
           write(*,*) 'residue = ',RName
           call Finalize

         end if

       end do

! ## Charge ##

       do l = ISearch + 1, FSearch

         NumAtomMol = NumAtomMol + 1

         AtomI = AtomName(l)

         ConfAllocI = .False.

         do k = 1 , IAtom

           if(AtomI == AtomNameParam(k,NResP)) then

             ChargeMol(NumAtomMol) = ChargeParam(k,NResP)
             ConfAllocI = .True.

           end if

           if((k == IAtom).and.(.not.ConfAllocI)) then

             write(*,*) 'ERROR : charge allocation'
             call Finalize

           end if

         end do

       end do

! ## Unit Name ##

       NumAtomMol = NumAtomMol - IAtom

       do k = 1 , NumUnit_inResi(NResP)

         NumUnitMol = NumUnitMol + 1

         do l = 1 , NumAtom_inUNITParam(k,NResP)

           NumAtomMol = NumAtomMol + 1
           UName(NumAtomMol) = UNITNameParam(k,NResP)
           UNum(NumAtomMol)  = NumUnitMol

         end do

       end do

! ## Bond ##

       do k = 1 , Nbond

         NumBondMol = NumBondMol + 1

         AtomI = BondPair_inResi(1,k,NResP)
         AtomJ = BondPair_inResi(2,k,NResP)

         ConfAllocI = .False.
         ConfAllocJ = .False.

         do l = ISearch+1, FSearch

           if(AtomName(l) == AtomI) then

             BondPairMol(1,NumBondMol) = l
             ConfAllocI = .True.

             exit

           end if

         end do

         do l = ISearch+1, FSearch

           if(AtomName(l) == AtomJ) then

             BondPairMol(2,NumBondMol) = l
             ConfAllocJ = .True.

             exit

           end if

         end do

         if((.not.ConfAllocI).or.(.not.ConfAllocJ)) then

           write(*,*) 'ERROR : missing bond'
           write(*,*) AtomI,AtomJ,' in residue ',RName
           call Finalize

         end if

       end do

! ## Improper ##

       do k = 1 , Nimpr

         NumImprMol = NumImprMol + 1

         AtomI = ImprPair_inResi(1,k,NResP)
         AtomJ = ImprPair_inResi(2,k,NResP)
         AtomK = ImprPair_inResi(3,k,NResP)
         AtomL = ImprPair_inResi(4,k,NResP)

         ConfAllocI = .False.
         ConfAllocJ = .False.
         ConfAllocK = .False.
         ConfAllocL = .False.

         do l = ISearch+1, FSearch

           if(AtomName(l) == AtomI) then

             ImprPairMol(1,NumImprMol) = l
             ConfAllocI = .True.

           else if(AtomName(l) == AtomJ) then

             ImprPairMol(2,NumImprMol) = l
             ConfAllocJ = .True.

           else if(AtomName(l) == AtomK) then

             ImprPairMol(3,NumImprMol) = l
             ConfAllocK = .True.

           else if(AtomName(l) == AtomL) then

             ImprPairMol(4,NumImprMol) = l
             ConfAllocL = .True.

           end if

         end do

         if((.not.ConfAllocI).or.(.not.ConfAllocJ).or. &
         &  (.not.ConfAllocK).or.(.not.ConfAllocL)) then

           write(*,*) 'ERROR : missing impr'
           write(*,*) AtomI,AtomJ,AtomK,AtomL,' in residue ',RName
           call Finalize

         end if

       end do

     end if

! ## check
     if(NumAtomMol /= NumAtm(i)) then
       write(*,*) 'ERROR : NumAtomMol in (MakeTopologicalData) '
       call Finalize
     end if

     NAtom = ISearch

     do j = 1 , NumMol(i)

       Nadd = NumAtm(i) * (j-1)

       do k = 1 , NumAtomMol

         Charge(NAtom+k)       = ChargeMol(k)
         AtomUnitName(NAtom+k) = UName(k)
         AtomUnitNum(NAtom+k)  = UNum(k) + NumUnitAcc + (j-1)*NumUnitMol

       end do

       do k = 1 , NumBondMol

         NumBond = NumBond + 1

         TmpBondI(NumBond) = BondPairMol(1,k) + Nadd
         TmpBondJ(NumBond) = BondPairMol(2,k) + Nadd

       end do

       do k = 1 , NumImprMol

         NumImproper = NumImproper + 1

         TmpImproI(NumImproper) = ImprPairMol(1,k) + Nadd
         TmpImproJ(NumImproper) = ImprPairMol(2,k) + Nadd
         TmpImproK(NumImproper) = ImprPairMol(3,k) + Nadd
         TmpImproL(NumImproper) = ImprPairMol(4,k) + Nadd

       end do

       NAtom = NAtom + NumAtm(i)

     end do

   end do

   NumRB = AtomUnitNum(N)
   if(QMaster) print *,'NumRB=',NumRB

!  ################ Angle and Dihedrals ###################

   AtomHands = 0
   AtomHandBond = 0

   do k = 1, NumBond

     i = TmpBondI(k)
     j = TmpBondJ(k)

     AtomHands(i) = AtomHands(i) + 1
     AtomHands(j) = AtomHands(j) + 1

     AtomHandBond(i,AtomHands(i)) = k
     AtomHandBond(j,AtomHands(j)) = k

   end do

   NumAngle = 0

   do i = 1 , N

     if(AtomHands(i)<2) cycle

     do k = 1, AtomHands(i) - 1

       do l = k + 1, AtomHands(i)

         Ipair = AtomHandBond(i,k)
         Jpair = AtomHandBond(i,l)

         i1 = TmpBondI(Ipair)
         j1 = TmpBondJ(Ipair)
         i2 = TmpBondI(Jpair)
         j2 = TmpBondJ(Jpair)

         NumAngle = NumAngle + 1

         TmpAngleJ(NumAngle) = i

         if(i1==i) then

           TmpAngleI(NumAngle) = j1

         else if(j1==i) then

           TmpAngleI(NumAngle) = i1

         else

           write(*,*) 'ERROR : angle'
           call Finalize

         end if

         if(i2==i) then

           TmpAngleK(NumAngle) = j2

         else if(j2==i) then

           TmpAngleK(NumAngle) = i2

         else

           write(*,*) 'ERROR : angle'
           call Finalize

         end if

       end do

     end do

   end do


   NumDihedral = 0

   do k = 1 , NumBond

     i = TmpBondI(k)
     j = TmpBondJ(k)

     if((AtomHands(i)<2).or.(AtomHands(j)<2)) cycle

     NbI = AtomHands(i)
     NbJ = AtomHands(j)

     do ii = 1 , NbI

       if(AtomHandBond(i,ii) == k) cycle

       do jj = 1 , NbJ

         if(AtomHandBond(j,jj) == k) cycle

         Ipair = AtomHandBond(i,ii)
         Jpair = AtomHandBond(j,jj)

         i1 = TmpBondI(Ipair)
         j1 = TmpBondJ(Ipair)

         i2 = TmpBondI(Jpair)
         j2 = TmpBondJ(Jpair)

         NumDihedral = NumDihedral + 1

         TmpDihedJ(NumDihedral) = i
         TmpDihedK(NumDihedral) = j

         if(i1==i) then

           TmpDihedI(NumDihedral) = j1

         else if(j1==i) then

           TmpDihedI(NumDihedral) = i1

         else

           write(*,*) 'ERROR : dihedral'
           write(*,*) 'pair I, j = ',i,j
           call Finalize

         end if

         if(i2==j) then

           TmpDihedL(NumDihedral) = j2

         else if(j2==j) then

           TmpDihedL(NumDihedral) = i2

         else

           write(*,*) 'ERROR : dihedral'
           write(*,*) 'pair i, J = ',i,j
           call Finalize

         end if

       end do

     end do

   end do

! ## checking charge neutrality

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

! ## rigid body ##

   if(QRigidBody) then

     if(NumBond /= 0) then

       allocate( BondI(NumBond) )
       allocate( BondJ(NumBond) )

       num = 0

       do ii = 1, NumBond

         i = TmpBondI(ii)
         j = TmpBondJ(ii)

         if(AtomUnitNum(i)==AtomUnitNum(j)) cycle

         num = num + 1

         BondI(num) = i
         BondJ(num) = j

       end do

       NumBond = num

       do i = 1 , NumBond
         TmpBondI(i) = BondI(i)
         TmpBondJ(i) = BondJ(i)
       end do

       deallocate( BondI, BondJ )

     end if

     if(NumAngle /= 0) then

       allocate( AngleI(NumAngle) )
       allocate( AngleJ(NumAngle) )
       allocate( AngleK(NumAngle) )

       num = 0

       do ii = 1 , NumAngle

         i = TmpAngleI(ii)
         j = TmpAngleJ(ii)
         k = TmpAngleK(ii)

         if((AtomUnitNum(i)==AtomUnitNum(j)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(k))) cycle

         num = num + 1

         AngleI(num) = i
         AngleJ(num) = j
         AngleK(num) = k

       end do

       NumAngle = num

       do i = 1 , NumAngle
         TmpAngleI(i) = AngleI(i)
         TmpAngleJ(i) = AngleJ(i)
         TmpAngleK(i) = AngleK(i)
       end do

       deallocate( AngleI, AngleJ, AngleK )

     end if

     if(NumDihedral /= 0) then

       allocate( DihedI(NumDihedral) )
       allocate( DihedJ(NumDihedral) )
       allocate( DihedK(NumDihedral) )
       allocate( DihedL(NumDihedral) )

       num = 0

       do ii = 1 , NumDihedral

         i = TmpDihedI(ii)
         j = TmpDihedJ(ii)
         k = TmpDihedK(ii)
         l = TmpDihedL(ii)

         if((AtomUnitNum(i)==AtomUnitNum(j)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(k)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(l))) cycle

         num = num + 1

         DihedI(num) = i
         DihedJ(num) = j
         DihedK(num) = k
         DihedL(num) = l

       end do

       NumDihedral = num

       do i = 1 , NumDihedral
         TmpDihedI(i) = DihedI(i)
         TmpDihedJ(i) = DihedJ(i)
         TmpDihedK(i) = DihedK(i)
         TmpDihedL(i) = DihedL(i)
       end do

       deallocate( DihedI, DihedJ, DihedK, DihedL )

     end if

     if(NumImproper /= 0) then

       allocate( ImproI(NumImproper) )
       allocate( ImproJ(NumImproper) )
       allocate( ImproK(NumImproper) )
       allocate( ImproL(NumImproper) )

       num = 0

       do ii = 1 , NumImproper

         i = TmpImproI(ii)
         j = TmpImproJ(ii)
         k = TmpImproK(ii)
         l = TmpImproL(ii)

         if((AtomUnitNum(i)==AtomUnitNum(j)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(k)).and. &
         &  (AtomUnitNum(i)==AtomUnitNum(l))) cycle

         num = num + 1

         ImproI(num) = i
         ImproJ(num) = j
         ImproK(num) = k
         ImproL(num) = l

       end do

       NumImproper = num

       do i = 1 , NumImproper
         TmpImproI(i) = ImproI(i)
         TmpImproJ(i) = ImproJ(i)
         TmpImproK(i) = ImproK(i)
         TmpImproL(i) = ImproL(i)
       end do

       deallocate( ImproI, ImproJ, ImproK, ImproL )

     end if

   end if

   if(NumBond /= 0) then

     allocate( BondI(NumBond) )
     allocate( BondJ(NumBond) )

     do i = 1 , NumBond
       BondI(i) = TmpBondI(i)
       BondJ(i) = TmpBondJ(i)
     end do

     if(.not.QPBC) call NoLJList(1)

   else if(.not.QPBC) then

     call NoLJList(0)

   end if

   if(NumAngle /= 0) then

     allocate( AngleI(NumAngle) )
     allocate( AngleJ(NumAngle) )
     allocate( AngleK(NumAngle) )

     do i = 1 , NumAngle
       AngleI(i) = TmpAngleI(i)
       AngleJ(i) = TmpAngleJ(i)
       AngleK(i) = TmpAngleK(i)
     end do

     if(.not.QPBC) call NoLJList(2)

   end if

   if(NumDihedral /= 0) then

     allocate( DihedI(NumDihedral) )
     allocate( DihedJ(NumDihedral) )
     allocate( DihedK(NumDihedral) )
     allocate( DihedL(NumDihedral) )

     do i = 1 , NumDihedral
       DihedI(i) = TmpDihedI(i)
       DihedJ(i) = TmpDihedJ(i)
       DihedK(i) = TmpDihedK(i)
       DihedL(i) = TmpDihedL(i)
     end do

   end if

   if(NumImproper /= 0) then

     allocate( ImproI(NumImproper) )
     allocate( ImproJ(NumImproper) )
     allocate( ImproK(NumImproper) )
     allocate( ImproL(NumImproper) )

     do i = 1 , NumImproper
       ImproI(i) = TmpImproI(i)
       ImproJ(i) = TmpImproJ(i)
       ImproK(i) = TmpImproK(i)
       ImproL(i) = TmpImproL(i)
     end do

   end if

!   print *, 'AtomUnitNum',N
!   do i = 1 , N
!     print *, i,AtomUnitNum(i),ResidNum(i),ResidName(i)
!   end do

end subroutine MakeTopologicalData


!######################################################################
!######################################################################


subroutine AllocateOPLS

use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster, QSHAKE, Qstdout, PolyFlag, QMacro
use FFParameters
use SHAKEparam, only : NSHAKEGroup
use UnitExParam, only : Avogadro
use NonbondParam, only : EpsLJ, SgmLJ
use BondedParam, only : NumBond, NumAngle, NumDihedral, NumImproper, &
&   BondI, BondJ, kBond, rBond, AngleI, AngleJ, AngleK, kTheta, Theta0,     &
&   DihedI, DihedJ, DihedK, DihedL, vdWSubtDih, CsDelDih, DupCsDelDih,      &
&   kChi, DeltaDih, NDih, DupFlag, DupkChi, DupDeltaDih, DupNDih, Ts,       &
&   ImproI, ImproJ, ImproK, ImproL, kImp, DeltaImp, NImp, CsDelImp
use AtomParam, only : AtomName, ResidName, ResidNum, Mass, InvMass

implicit none

integer :: ii, i, j, k, l, jj
integer :: i1, i2, j1, j2
integer, dimension(:), allocatable :: MinResid, MaxResid
character(len=4), dimension(4) :: Name, PName
character(len=4), dimension(4) :: TmpAName, TmpRName
integer, dimension(4) :: TmpResNum
integer :: count, num
integer, dimension(:), allocatable :: TmpAngleI, TmpAngleJ, TmpAngleK

   ii = 0

   do i = 1 , NumSpec

     do j = 1 , NumMol(i)

       ii = ii + 1

     end do

   end do

   allocate( MinResid(ii) )
   allocate( MaxResid(ii) )

   MinResid = 100000
   MaxResid = 0

   ii = 0
   l  = 0

   do i = 1 , NumSpec

     do j = 1 , NumMol(i)

       ii = ii + 1

       do k = 1 , NumAtm(i)

         l = l + 1
         MinResid(ii) = Min( ResidNum(l) , MinResid(ii) )
         MaxResid(ii) = Max( ResidNum(l) , MaxResid(ii) )

       end do

     end do

   end do

   ii = 0
   l  = 0

   do i = 1 , NumSpec

     do j = 1 , NumMol(i)

       ii = ii + 1

       do k = 1 , NumAtm(i)

         l = l + 1

         if(PolyFlag(i)) then

           if( ResidNum(l) == MinResid(ii) ) then
             write(TmpRName(1),'(a,a)') Trim(adjustl(ResidName(l))),'I'
             ResidName(l) = TmpRName(1)
           else if( ResidNum(l) == MaxResid(ii) ) then
             write(TmpRName(1),'(a,a)') Trim(adjustl(ResidName(l))),'F'
             ResidName(l) = TmpRName(1)
           end if

         end if

       end do

     end do

   end do

! --------------------------------------------------------------

   allocate( Mass   (N) )
   allocate( InvMass(N) )
   allocate( EpsLJ(N) )
   allocate( SgmLJ(N) )

   do i = 1 , N

     TmpAName (1) = AtomName (i)
     TmpRName (1) = ResidName(i)
     TmpResNum(1) = ResidNum (i)

     if(QMacro) then
       if(ResidName(i)=='CGSP') cycle
     end if

!   ------------------
     call NameExch(1)
!   ------------------

     do j = 1 , NumLJParam

       if( Name(1) == LJAtoms(j) ) then

         if( EpsLJParam(j) /= 0. ) then
           EpsLJ(i) = sqrt( EpsLJParam(j) )
         else
           EpsLJ(i) = 0.d0
         end if
         SgmLJ(i)   = SigmaLJParam(j)

         exit

       end if

       if( j == NumLJParam ) then

         if(QMaster) write(*,*) 'No LJ parameter was found for',Name(1)
         if(QMaster) write(*,*) 'Atom Number =' , i

         call Finalize

       end if

     end do

     do j = 1 , NumAtomTypeParam

       if( Name(1) == AtomTypeParam(j) ) then

         Mass(i) = MassParam(j)

         exit

       end if

       if( j == NumAtomTypeParam ) then

         if(QMaster) write(*,*) 'No Mass parameter was found for',Name(1)
         if(QMaster) write(*,*) 'Atom Number =' , i

         call Finalize

       end if

     end do

   end do

!## Macroparticle
   if(QMacro) then
     call Set_Sphere(1)
   end if

! ----
! Unit
! ----
   Mass = Mass * 1.d-3 / Avogadro

   InvMass = 1.d0 / Mass


! ----------------------------------------------------------------------
! ## Bond Parameters

   if(NumBond /= 0) then

     allocate( kBond(NumBond) )
     allocate( rBond(NumBond) )

     do i = 1, NumBond

       TmpAName(1) = AtomName(BondI(i))
       TmpAName(2) = AtomName(BondJ(i))

       TmpRName(1) = ResidName(BondI(i))
       TmpRName(2) = ResidName(BondJ(i))

       TmpResNum(1) = ResidNum(BondI(i))
       TmpResNum(2) = ResidNum(BondJ(i))

!     ------------------
       call NameExch(2)
!     ------------------

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
             write( 6,*) 'Bond Number=',i,' Pair=',BondI(i),BondJ(i)
             write( 6,*) 'TmpRName=',TmpRName(1),TmpRName(2)
             write( 6,*) 'TmpAName=',TmpAName(1),TmpAName(2)
             write( 6,*) 'Name=',Name(1),Name(2)

           end if

           call Finalize

         end if

       end do

     end do

! ----------------------------------------
     if(QSHAKE) then

       call MakeSHAKEList

     else

       NSHAKEGroup = 0

     end if
! ----------------------------------------

   end if

! ----------------------------------------------------------------------
! ## Angle Parameters

   if(NumAngle /= 0) then

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

     do i = 1, NumAngle

       TmpAName(1) = AtomName(AngleI(i))
       TmpAName(2) = AtomName(AngleJ(i))
       TmpAName(3) = AtomName(AngleK(i))

       TmpRName(1) = ResidName(AngleI(i))
       TmpRName(2) = ResidName(AngleJ(i))
       TmpRName(3) = ResidName(AngleK(i))

       TmpResNum(1) = ResidNum(AngleI(i))
       TmpResNum(2) = ResidNum(AngleJ(i))
       TmpResNum(3) = ResidNum(AngleK(i))

!     ------------------
       call NameExch(3)
!     ------------------

       do j = 1 , NumAngleParam

         if( AnglePairAtoms(2,j) == Name(2) ) then

           if(((AnglePairAtoms(1,j) == Name(1))  &
       &  .and.(AnglePairAtoms(3,j) == Name(3))) &
       &  .or.((AnglePairAtoms(1,j) == Name(3))  &
       &  .and.(AnglePairAtoms(3,j) == Name(1)))) then

             kTheta(i) = kThetaParam(j)
             Theta0(i) = Theta0Param(j)

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
             write( 6,*) ( TmpRName(k) , k = 1 , 3 )
             write( 6,*) ( TmpAName(k) , k = 1 , 3 )
             write( 6,*) ( Name    (k) , k = 1 , 3 )

           end if

           call Finalize

         end if

       end do

     end do

   end if

! ----------------------------------------------------------------------
! ## Dihedral Parameters

   if(NumDihedral /= 0) then

     allocate( vdWSubtDih(NumDihedral) )

     allocate( kChi(NumDihedral) )
     allocate( DeltaDih(NumDihedral) )
     allocate( NDih(NumDihedral) )

     allocate( DupFlag(NumDihedral) )
     allocate( DupkChi(3,NumDihedral) )
     allocate( DupDeltaDih(3,NumDihedral) )
     allocate( DupNDih(3,NumDihedral) )

     allocate( CsDelDih(NumDihedral) )
     allocate( DupCsDelDih(3,NumDihedral) )

     allocate( Ts(NumDihedral) )

     vdWSubtDih = .False.

     DupFlag = 0

     do i = 1 , NumDihedral

       TmpAName(1) = AtomName(DihedI(i))
       TmpAName(2) = AtomName(DihedJ(i))
       TmpAName(3) = AtomName(DihedK(i))
       TmpAName(4) = AtomName(DihedL(i))

       TmpRName(1) = ResidName(DihedI(i))
       TmpRName(2) = ResidName(DihedJ(i))
       TmpRName(3) = ResidName(DihedK(i))
       TmpRName(4) = ResidName(DihedL(i))

       TmpResNum(1) = ResidNum(DihedI(i))
       TmpResNum(2) = ResidNum(DihedJ(i))
       TmpResNum(3) = ResidNum(DihedK(i))
       TmpResNum(4) = ResidNum(DihedL(i))

!     ------------------
       call NameExch(4)
!     ------------------

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

             if( ( count == 0 ) .and. ( jj == NumDihedralParam ) ) then

               if(QMaster) then

#ifndef BMONI
                 write(11,*) 'lost parameter : Dihedral'
#endif
                 write( 6,*) 'lost parameter : Dihedral'
                 write( 6,*) 'Dihedral Number=' , i
                 write( 6,*) 'TmpRName=' , ( TmpRName(l) , l = 1 , 4 )
                 write( 6,*) 'TmpAName=' , ( TmpAName(l) , l = 1 , 4 )
                 write( 6,*) 'Name    =' , ( Name    (l) , l = 1 , 4 )

               end if

               call Finalize

             end if

           end do

         end if

       end do

     end do

     do i = 1 , NumDihedral

       CsDelDih(i) = cos( DeltaDih(i) )

       do j = 1, DupFlag(i)

         DupCsDelDih(j,i) = cos( DupDeltaDih(j,i) )

       end do

     end do

     do i = 1 , NumDihedral

       i1 = DihedI(i)
       i2 = DihedL(i)

       do j = 1 , NumAngle

         j1 = AngleI(j)
         j2 = AngleK(j)

         if( ( ( i1 == j1 ) .and. (i2 == j2 ) ) .or. &
         &   ( ( i1 == j2 ) .and. (i2 == j1 ) ) ) then

           vdWSubtDih(i) = .True.
!          if(QMaster) write(*,*) 'DUPLICATE 1-3, 1-4',i

           exit

         end if

       end do

       do j = i+1 , NumDihedral

         j1 = DihedI(j)
         j2 = DihedL(j)

         if( ( ( i1 == j1 ) .and. (i2 == j2 ) ) .or. &
         &   ( ( i1 == j2 ) .and. (i2 == j1 ) ) ) then

           vdWSubtDih(i) = .True.
!          if(QMaster) write(*,*) 'DUPLICATE 1-4',i
           exit

         end if

       end do

     end do

   end if

! ----------------------------------------------------------------------
! ## Improper torsion parameters

   if(NumImproper /= 0) then

     allocate( kImp(NumImproper) )
     allocate( DeltaImp(NumImproper) )
     allocate( NImp(NumImproper) )
     allocate( CsDelImp(NumImproper) )

     do i = 1, NumImproper

       TmpAName(1) = AtomName(ImproI(i))
       TmpAName(2) = AtomName(ImproJ(i))
       TmpAName(3) = AtomName(ImproK(i))
       TmpAName(4) = AtomName(ImproL(i))

       TmpRName(1) = ResidName(ImproI(i))
       TmpRName(2) = ResidName(ImproJ(i))
       TmpRName(3) = ResidName(ImproK(i))
       TmpRName(4) = ResidName(ImproL(i))

       TmpResNum(1) = ResidNum(ImproI(i))
       TmpResNum(2) = ResidNum(ImproJ(i))
       TmpResNum(3) = ResidNum(ImproK(i))
       TmpResNum(4) = ResidNum(ImproL(i))

!     ------------------
       call NameExch(4)
!     ------------------

       do j = 1, NumImproperParam

         PName(1) = ImproperPairAtoms(1,j)
         PName(2) = ImproperPairAtoms(2,j)
         PName(3) = ImproperPairAtoms(3,j)
         PName(4) = ImproperPairAtoms(4,j)

         if((PName(3) == Name(3)) .and. &
         &(((PName(1) == Name(1)) .and. (PName(2) == Name(2)) .and. &
         &  (PName(4) == Name(4))) .or. &
         & ((PName(1) == Name(2)) .and. (PName(2) == Name(1)) .and. &
         &  (PName(4) == Name(4))) .or. &
         & ((PName(1) == Name(4)) .and. (PName(2) == Name(2)) .and. &
         &  (PName(4) == Name(1))) .or. &
         & ((PName(1) == Name(1)) .and. (PName(2) == Name(4)) .and. &
         &  (PName(4) == Name(2))) .or. &
         & ((PName(1) == Name(2)) .and. (PName(2) == Name(4)) .and. &
         &  (PName(4) == Name(1))) .or. &
         & ((PName(1) == Name(4)) .and. (PName(2) == Name(1)) .and. &
         &  (PName(4) == Name(2))))) then

             kImp(i)     = kImpParam(j)
             NImp(i)     = NImpParam(j)
             DeltaImp(i) = DeltaImpParam(j)
             exit

         end if

         if( j == NumImproperParam ) then

           do jj = 1 , NumImproperParam

             PName(1) = ImproperPairAtoms(1,jj)
             PName(2) = ImproperPairAtoms(2,jj)
             PName(3) = ImproperPairAtoms(3,jj)
             PName(4) = ImproperPairAtoms(4,jj)

             if(( PName(3) == Name(3) ) .and. &
             &((( PName(2) == Name(2) ) .and. ( PName(1) == 'X' ) .and. &
             &  ( PName(4) == 'X'     )) .or. &
             & (( PName(1) == 'X'     ) .and. ( PName(1) == 'X' ) .and. &
             &  ( PName(4) == Name(4) )))) then

                 kImp(i)     = kImpParam(jj)
                 NImp(i)     = NImpParam(jj)
                 DeltaImp(i) = DeltaImpParam(jj)
                 exit

             end if

             if( jj == NumImproperParam ) then

               if(QMaster) then

#ifndef BMONI
                 write(11,*) 'lost parameter : Improper'
#endif
                 write( 6,*) 'lost parameter : Improper'
                 write( 6,*) 'Improper Number=',i
                 write( 6,*) 'TmpRName=',(TmpRName(l),l=1,4)
                 write( 6,*) 'TmpAName=',(TmpAName(l),l=1,4)
                 write( 6,*) 'Name=',(Name(l),l=1,4)

               end if

               call Finalize

             end if

           end do

         end if

       end do

     end do

     do i = 1, NumImproper

       CsDelImp(i) = cos( DeltaImp(i) )

     end do

   end if

! ----------------------------------------------------------------------

   if(QMaster.and.Qstdout) then

     write(*,*) 'NumBOND=',NumBond
     write(*,*) 'NumANGLE=',NumAngle
     write(*,*) 'NumDIHED=',NumDihedral
     write(*,*) 'NumIMPRO=',NumImproper

   end if

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

Contains

   subroutine NameExch(Na)

   integer :: j, k, l, Na, count1

     do j = 1, Na

       count1 = 0

       do k = 1, NumResidueParam

         if( TmpRName(j) == ResiNameParam(k) ) then

           count1 = count1 + 1

           do l = 1, NumAtom_inResi(k)

             if(TmpAName(j) == AtomNameParam(l,k)) then

               Name(j) = AtomType(l,k)

             end if

           end do

         end if

       end do

       if( count1 /= 1 ) then

         if(QMaster) then

           if( Na == 1 ) then

             write(*,'(a/a,a/a,a/a,i6)')             &
             &          'ERROR: no residue (ATOM)',  &
             &          'Residue = ',TmpRName(1),    &
             &          'Atom    = ',TmpAName(1),    &
             &          'count   = ',count1

           else if( Na == 2 ) then

             write(*,'(a/a,a/a,a/a,i6)')             &
             &          'ERROR: no residue (BOND)',  &
             &          'Residue = ',TmpRName(j),    &
             &          'Atom    = ',TmpAName(j),    &
             &          'count   = ',count1

           else if( Na == 3 ) then

             write(*,'(a/a,a/a,a/a,i6)')              &
             &          'ERROR: no residue (ANGLE)',  &
             &          'Residue = ',TmpRName(j),     &
             &          'Atom    = ',TmpAName(j),     &
             &          'count   = ',count1

           else if( Na == 4 ) then

             write(*,'(a/a,a/a,a/a,i6)')                           &
             &          'ERROR: no residue (DIHEDRAL or IMPROPER)',&
             &          'Residue = ',TmpRName(j),                  &
             &          'Atom    = ',TmpAName(j),                  &
             &          'count   = ',count1

           end if

         end if

         call Finalize

       end if

     end do

   end subroutine NameExch

! ----------------------------------------------------------------------

end subroutine AllocateOPLS
