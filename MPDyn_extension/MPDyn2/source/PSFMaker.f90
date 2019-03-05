! ############################
! ## SUBROUTINE LIST 
! ## -- PSFMaker 
! ## -- SetCondition 
! ## -- MakeCharmPSF 
! ## -- Write_PSF 
! ############################


! ######################################################################
! ######################################################################


Program PSFMaker

implicit none

   call SetCondition ! local file

   call Read_Charmm_Parameter
   call Read_Charmm_Topology

   call MakeCharmPSF

   call Write_PSF

end program PSFMaker


! ######################################################################
! ######################################################################


subroutine SetCondition

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster, QRigidBody
use FFParameters
use IOparam, only : PSF_file, Topology_file, Parameter_file
use AtomParam, only : MolName

implicit none

integer :: i, j, ii, jj, Count
character(len=80) :: String, String1
logical, dimension(80) :: empty
integer,dimension(10) :: first, last
logical :: Qdone

   QRigidBody = .False.
   QMaster = .True.
   Qdone = .False.

   Count = 0

   do

     if(Qdone) exit

     select case(Count)

     case(0)
       print *, 'Name of PARAMETER file = ?'
       print *, 'Put "D" for a default choice, "param/par_charmm27.prm"'
     case(1)
       print *, 'Name of TOPOLOGY file = ?'
       print *, 'Put "D" for a default choice, "param/top_charmm27.prm"'
     case(2)
       print *, 'Name of PSF file = ?'
     case(3)
       print *, 'Number of molecular species = ?'
     case(4)
       print *, 'Molecular names ? Give them in one line.'
     case(5)
       print *, 'Number of Molecules ? Give them in one line.'

     end select

     read( 5,'(a80)') String1

     String = trim(adjustl(String1))

     if((String(1:1) == '!').or.(String(1:1) == '#')) cycle

     if((String(1:5) == '<end>').or.(String(1:5) == '<END>')) exit

     Count = Count + 1

     select case(Count)

     case(1)

       if(String=='D') then
         write(Parameter_file,*) 'param/par_charmm27.prm'
       else
         read(String,*) Parameter_file
       end if

     case(2)

       if(String=='D') then
         write(Topology_file,*) 'param/top_charmm27.prm'
       else
         read(String,*) Topology_file
       end if

     case(3)

       read(String,*) PSF_file


     case(4)

       read(String,*) NumSpec

       allocate( MolName(NumSpec) )
       allocate( NumMol(NumSpec) )
       allocate( NumAtm(NumSpec) )

     case(5)

! ## Name of molecules

       do i = 1, 80
         if(String(i:i) == ' ') then
           empty(i) = .True.
         else
           empty(i) = .False.
         end if
       end do

       ii = 0
       j  = 1
       first(j) = 1
       do i = 2, 80
         if(empty(i).and.(.not.empty(i-1))) then
           ii = ii + 1
           last(ii) = i-1
         end if
         if(empty(i-1).and.(.not.empty(i))) then
           j = j + 1
           first(j) = i
         end if
         if(ii==NumSpec) then
           jj = i
           exit
         end if
       end do

       do i = 1, NumSpec
         MolName(i) = String(first(i):last(i))
       end do

     case(6)

! ## Number of Molecules

       do i = 1, 80
         if(String(i:i) == ' ') then
           empty(i) = .True.
         else
           empty(i) = .False.
         end if
       end do

       ii = 0
       j  = 1
       first(j) = 1
       do i = 2, 80
         if(empty(i).and.(.not.empty(i-1))) then
           ii = ii + 1
           last(ii) = i-1
         end if
         if(empty(i-1).and.(.not.empty(i))) then
           j = j + 1
           first(j) = i
         end if
         if(ii==NumSpec) then
           jj = i
           exit
         end if
       end do

       do i = 1, NumSpec
         read(String(first(i):last(i)),*) NumMol(i)
       end do

       Qdone = .True.

     end select

   end do

   write(*,*)
   write(*,*) '***********************************************************'
   write(*,*) ' Parameter file = ',trim(Parameter_file)
   write(*,*) ' Topology  file = ',trim(Topology_file)
   write(*,*) ' PSF       file = ',trim(PSF_file)
   write(*,*) ' The number of molecule types = ',NumSpec
   do i = 1, NumSpec
     write(*,'(a,i2,a,a,a,i6)') ' Molecule ',i,' = ',MolName(i),'NumMol =',NumMol(i)
   end do 
   write(*,*) '***********************************************************'


end subroutine SetCondition


! ######################################################################
! ######################################################################


subroutine MakeCharmPSF

use Numbers, only : N, NumSpec, NumMol, NumAtm
use FFParameters
use NonbondParam, only : Charge
use BondedParam, only : NumBond, NumAngle, NumDihedral, NumImproper, &
&   BondI, BondJ, AngleI, AngleJ, AngleK, DihedI, DihedJ, DihedK, DihedL, &
&   ImproI, ImproJ, ImproK, ImproL
use AtomParam, only : MolName, AtomName, ResidName, ResidNum, Mass, TypeName
use CommonBlocks, only : QMaster

implicit none

integer :: i, j, k, l, ii, jj, ia, ib, im
integer :: i1, i2, j1, j2, NbI, NbJ

integer, parameter :: MolDIM = 10000

real(8), dimension(MolDIM,NumSpec) :: MolCharge
real(8), dimension(MolDIM,NumSpec) :: MolMass

integer, dimension(MolDIM,NumSpec) :: MolBondI
integer, dimension(MolDIM,NumSpec) :: MolBondJ
integer, dimension(MolDIM,NumSpec) :: MolImproI
integer, dimension(MolDIM,NumSpec) :: MolImproJ
integer, dimension(MolDIM,NumSpec) :: MolImproK
integer, dimension(MolDIM,NumSpec) :: MolImproL

integer, parameter :: tempdimension = 100000
integer, dimension(tempdimension) :: TmpAngleI
integer, dimension(tempdimension) :: TmpAngleJ
integer, dimension(tempdimension) :: TmpAngleK
integer, dimension(tempdimension) :: TmpDihedI
integer, dimension(tempdimension) :: TmpDihedJ
integer, dimension(tempdimension) :: TmpDihedK
integer, dimension(tempdimension) :: TmpDihedL

character(len=4) :: AtomI, AtomJ, AtomK, AtomL, RName
integer, dimension(NumSpec) :: NResP
integer, dimension(NumSpec) :: NumBondMol, NumImproMol

integer, dimension(:),   allocatable :: AtomHands
integer, dimension(:,:), allocatable :: AtomHandBond
integer :: Ipair, Jpair

character(len=4), dimension(1) :: TmpAName, TmpRName
integer, dimension(1) :: TmpResNum
character(len=4), dimension(1) :: Name
character(len=4) :: Temporary

   do i = 1 , NumSpec

     RName = trim(adjustl(MolName(i)))

     do j = 1, NumResidueParam

       if( ResiNameParam(j) == RName ) then
         NumAtm(i) = NumAtom_inResi(j)
         NResP(i) = j
         exit
       end if
       if( j == NumResidueParam ) then
         write(*,*) 'ERROR : missing parameter'
         write(*,*) RName
         call Finalize
       end if

     end do

   end do

   N = 0
   do i = 1 , NumSpec
     N = N + NumMol(i) * NumAtm(i)
   end do

   do i = 1 , NumSpec

     ii = NResP(i) ! ii-th parameter should be read 

     do j = 1 , NumAtm(i)

       MolCharge(j,i) = ChargeParam(j,ii)

       AtomI = AtomType(j,ii)

       do k = 1 , NumAtomTypeParam
         if(AtomTypeParam(k) == AtomI) then
           MolMass(j,i) = MassParam(k)
           exit
         end if
         if(k == NumAtomTypeParam) then
           write(*,*) 'ERROR : missing atomtype'
           stop
         end if
       end do

     end do

     do j = 1 , NumBond_inResi(ii)

       AtomI = BondPair_inResi(1,j,ii)
       AtomJ = BondPair_inResi(2,j,ii)

       l = 0
       do k = 1 , NumAtm(i)

         if(AtomI == AtomNameParam(k,ii)) then
           MolBondI(j,i) = k
           l = l + 1
         else if(AtomJ == AtomNameParam(k,ii)) then
           MolBondJ(j,i) = k
           l = l + 1
         end if

       end do

       if(l/=2) then
         write(*,*) 'ERROR : bond'
         write(*,*) 'l=',l,AtomI,AtomJ
         stop
       end if

     end do

     do j = 1 , NumDoub_inResi(ii)

       AtomI = DoubPair_inResi(1,j,ii)
       AtomJ = DoubPair_inResi(2,j,ii)

       l = 0
       do k = 1 , NumAtm(i)

         if(AtomI == AtomNameParam(k,ii)) then
           MolBondI(j+NumBond_inResi(ii),i) = k
           l = l + 1
         else if(AtomJ == AtomNameParam(k,ii)) then
           MolBondJ(j+NumBond_inResi(ii),i) = k
           l = l + 1
         end if

       end do

       if(l/=2) then
         write(*,*) 'ERROR : double'
         stop
       end if

     end do

     NumBondMol(i) = NumBond_inResi(ii) + NumDoub_inResi(ii)

     do j = 1 , NumIMPR_inResi(ii)

       AtomI = ImprPair_inResi(1,j,ii)
       AtomJ = ImprPair_inResi(2,j,ii)
       AtomK = ImprPair_inResi(3,j,ii)
       AtomL = ImprPair_inResi(4,j,ii)

       l = 0
       do k = 1 , NumAtm(i)

         if(AtomI == AtomNameParam(k,ii)) then
           MolImproI(j,i) = k
           l = l + 1
         else if(AtomJ == AtomNameParam(k,ii)) then
           MolImproJ(j,i) = k
           l = l + 1
         else if(AtomK == AtomNameParam(k,ii)) then
           MolImproK(j,i) = k
           l = l + 1
         else if(AtomL == AtomNameParam(k,ii)) then
           MolImproL(j,i) = k
           l = l + 1
         end if

       end do

       if(l/=4) then
         write(*,*) 'ERROR : impro'
         stop
       end if

     end do

     NumImproMol(i) = NumIMPR_inResi(ii)

   end do

   NumBond     = 0
   NumImproper = 0

   do i = 1 , NumSpec
     NumBond = NumBond + NumBondMol(i) * NumMol(i)
     NumImproper = NumImproper + NumImproMol(i) * NumMol(i)
   end do


! #######################

   allocate( Charge   (N) )
   allocate( Mass     (N) )
   allocate( AtomName (N) )
   allocate( ResidName(N) )
   allocate( ResidNum (N) )
   allocate( TypeName (N) )

   if(NumBond /= 0) then
     allocate( BondI(NumBond) )
     allocate( BondJ(NumBond) )
   end if

   if(NumImproper /= 0) then
     allocate( ImproI(NumImproper) )
     allocate( ImproJ(NumImproper) )
     allocate( ImproK(NumImproper) )
     allocate( ImproL(NumImproper) )
   end if

   l  = 0
   ia = 0
   ib = 0
   im = 0

   do i = 1, NumSpec

     ii = NResP(i)

     do j = 1 , NumMol(i)

       im = im + 1

       do k = 1 , NumBondMol(i)

         ib = ib + 1

         BondI(ib) = MolBondI(k,i) + ia
         BondJ(ib) = MolBondJ(k,i) + ia

       end do

       do k = 1 , NumImproMol(i)

         l = l + 1

         ImproI(l) = MolImproI(k,i) + ia
         ImproJ(l) = MolImproJ(k,i) + ia
         ImproK(l) = MolImproK(k,i) + ia
         ImproL(l) = MolImproL(k,i) + ia

       end do

       do k = 1 , NumAtm(i)

         ia = ia + 1

         Charge   (ia) = MolCharge(k,i)
         Mass     (ia) = MolMass  (k,i)
         AtomName (ia) = AtomNameParam(k,ii)
         ResidName(ia) = MolName(i)(1:4)
         ResidNum (ia) = im

       end do

     end do

   end do

   do i = 1 , N

     TmpAName (1) = AtomName (i)
     TmpRName (1) = ResidName(i)
     TmpResNum(1) = ResidNum (i)

!   ------------------
     call NameExch(1)
!   ------------------

     TypeName(i) = Name(1)

   end do

   allocate( AtomHands(N) )
   allocate( AtomHandBond(N,6) )


!  ################ Angle and Dihedrals ###################

print *, 'NumBond=',NumBond
   AtomHands = 0
   AtomHandBond = 0

   do k = 1, NumBond

     i = BondI(k)
     j = BondJ(k)
     print *, k, i, j

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

         i1 = BondI(Ipair)
         j1 = BondJ(Ipair)
         i2 = BondI(Jpair)
         j2 = BondJ(Jpair)

         NumAngle = NumAngle + 1

         TmpAngleJ(NumAngle) = i

         if(i1==i) then

           TmpAngleI(NumAngle) = j1

         else if(j1==i) then

           TmpAngleI(NumAngle) = i1

         else

           write(*,*) 'ERROR : angle'
           stop

         end if

         if(i2==i) then

           TmpAngleK(NumAngle) = j2

         else if(j2==i) then

           TmpAngleK(NumAngle) = i2

         else

           write(*,*) 'ERROR : angle'

         end if

       end do

     end do

   end do

   NumDihedral = 0

   do k = 1 , NumBond

     i = BondI(k)
     j = BondJ(k)

     if((AtomHands(i)<2).or.(AtomHands(j)<2)) cycle

     NbI = AtomHands(i)
     NbJ = AtomHands(j)

     do ii = 1 , NbI

       if(AtomHandBond(i,ii) == k) cycle

       do jj = 1 , NbJ

         if(AtomHandBond(j,jj) == k) cycle

         Ipair = AtomHandBond(i,ii)
         Jpair = AtomHandBond(j,jj)

         i1 = BondI(Ipair)
         j1 = BondJ(Ipair)

         i2 = BondI(Jpair)
         j2 = BondJ(Jpair)

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
           stop

         end if

         if(i2==j) then

           TmpDihedL(NumDihedral) = j2

         else if(j2==j) then

           TmpDihedL(NumDihedral) = i2

         else

           write(*,*) 'ERROR : dihedral'
           write(*,*) 'pair i, J = ',i,j
           stop

         end if

       end do

     end do

   end do

   allocate( AngleI(NumAngle) )
   allocate( AngleJ(NumAngle) )
   allocate( AngleK(NumAngle) )

   allocate( DihedI(NumDihedral) )
   allocate( DihedJ(NumDihedral) )
   allocate( DihedK(NumDihedral) )
   allocate( DihedL(NumDihedral) )

   do i = 1 , NumAngle
     AngleI(i) = TmpAngleI(i)
     AngleJ(i) = TmpAngleJ(i)
     AngleK(i) = TmpAngleK(i)
   end do

   do i = 1 , NumDihedral
     DihedI(i) = TmpDihedI(i)
     DihedJ(i) = TmpDihedJ(i)
     DihedK(i) = TmpDihedK(i)
     DihedL(i) = TmpDihedL(i)
   end do


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

Contains

   subroutine NameExch(Na)

   integer :: j, k, l, Na, count1
   integer, dimension(4) :: Occ

     Occ = 0

     do j = 1, Na

       count1 = 0

       do k = 1, NumResidueParam

         if( TmpRName(j) == ResiNameParam(k) ) then

           count1 = count1 + 1

           do l = 1, NumAtom_inResi(k)

             if(TmpAName(j) == AtomNameParam(l,k)) then

               Name(j) = AtomType(l,k)
               Occ (j) = 1

             end if

           end do

         end if

       end do

       if(count1 == 0) then

         Temporary = TmpRName(j)

!         if( Temporary == 'HIS' ) then
!
!           ii = 0
!
!           do k = 1, NumSpec
!
!           do l = 1, NumMol(k)
!
!             ii = ii + 1
!
!             if( ( TmpResNum(j) == MinResid(ii) ) .or. &
!             &   ( TmpResNum(j) == MaxResid(ii) ) ) then
!
!               Temporary = 'HSE'
!
!             end if
!
!           end do
!
!           end do
!
!         end if

         do k = 1, NumResidueParam

           if( Temporary(1:3) == ResiNameParam(k) ) then

             count1 = count1 + 1

             do l = 1 , NumAtom_inResi(k)

               if( TmpAName(j) == AtomNameParam(l,k) ) then

                 Name(j) = AtomType(l,k)
                 Occ (j) = 1

               end if

             end do

           end if

         end do

       end if

       if( count1 /= 1 ) then

         if(QMaster) then

           if( Na == 1 ) then

             write(*,*) 'ERROR: no residue (ATOM)',&
             &          'Residue=',TmpRName(1),    &
             &          'Atom=',   TmpAName(1),    &
             &          'count=',  count1

           else if( Na == 2 ) then

             write(*,*) 'ERROR: no residue (BOND)',&
             &          'Residue=',TmpRName(j),    &
             &          'Atom=',   TmpAName(j),    &
             &          'count=',  count1

           else if( Na == 3 ) then

             write(*,*) 'ERROR: no residue (ANGLE)',&
             &          'Residue=',TmpRName(j),     &
             &          'Atom=',   TmpAName(j),     &
             &          'count=',  count1

           else if( Na == 4 ) then

             write(*,*) 'ERROR: no residue (DIHEDRAL or IMPROPER)',&
             &          'Residue=',TmpRName(j),                    &
             &          'Atom=',   TmpAName(j),                    &
             &          'count=',  count1

           end if

         end if

         call Finalize

       end if

       if( Occ(j) == 0 ) then

         write( *,'(a/2a/2a)') 'ERROR: no parameter ',     &
         &                     'Atom    = ' , TmpAName(j), &
         &                     'Residue = ' , TmpRName(j)
#ifndef BMONI
         write(11,'(a/2a/2a)') 'ERROR: no parameter ',     &
         &                     'Atom    = ' , TmpAName(j), &
         &                     'Residue = ' , TmpRName(j)
#endif

         call Finalize

       end if

     end do

   end subroutine NameExch

end subroutine MakeCharmPSF


! ######################################################################
! ######################################################################


subroutine Write_PSF

use Numbers, only : N
use IOparam, only : PSF_file
use NonbondParam, only : Charge
use BondedParam, only : NumBond, NumAngle, NumDihedral, NumImproper, &
&   BondI, BondJ, AngleI, AngleJ, AngleK, DihedI, DihedJ, DihedK, DihedL, &
&   ImproI, ImproJ, ImproK, ImproL
use AtomParam, only : AtomName, ResidName, ResidNum, Mass, TypeName

implicit none

integer :: i, j, ii, jj, k
character(len=20) :: username
character(len=17) :: ptime

   open(1,file=trim(PSF_file),status='unknown')

   write(1,'(a/)') 'PSF '
   write(1,'(a)')  '       2 !NTITLE'
   write(1,'(a)')  &
   & '* SCRIPT FILE PRODUCED BY PSFMaker version 2                                    '

   call system('whoami > MPDyn_temporary_file')
   open(18,file='MPDyn_temporary_file')
   read(18,*) username
   close(18)
   call system('rm -f MPDyn_temporary_file')

   call AcTime(ptime)

   write(1,'(4a/)') &
   & '*  DATE:    ',ptime,'      CREATED BY USER: ',username

   write(1,'(i8,x,a)') N,'!NATOM'

   j = 0
   do i = 1 , N
     write(1,'(i8,x,a4,x,i4,3(x,a4),e15.6,f10.4,i12)') &
     &  i,ResidName(i),ResidNum(i),ResidName(i),AtomName(i),TypeName(i),Charge(i),Mass(i),j
   end do

   write(1,'(/i8,x,a)') NumBond,'!NBOND: bonds'

   ii = NumBond / 4

   do i = 1 , ii

     k = (i-1) * 4
     write(1,'(8i8)') ( BondI(k+j) , BondJ(k+j) , j = 1 , 4 )

   end do

   ii = ii * 4
   jj = NumBond - ii

   if( jj /= 0 ) write(1,'(8i8)') ( BondI(ii+j) , BondJ(ii+j) , j = 1 , jj )

   write(1,'(/i8,x,a)') NumAngle,'!NTHETA: angles'

   ii = NumAngle / 3

   do i = 1 , ii

     k = (i-1) * 3
     write(1,'(9i8)') ( AngleI(k+j) , AngleJ(k+j) , AngleK(k+j) , j = 1 , 3 )

   end do

   ii = ii * 3
   jj = NumAngle - ii

   if( jj /= 0 ) &
   &  write(1,'(8i8)') ( AngleI(ii+j) , AngleJ(ii+j) , AngleK(ii+j) , j = 1 , jj )

   write(1,'(/i8,x,a)') NumDihedral,'!NPHI: dihedrals'

   ii = NumDihedral / 2

   do i = 1 , ii

     k = (i-1) * 2
     write(1,'(8i8)') ( DihedI(k+j) , DihedJ(k+j) , &
     &                  DihedK(k+j) , DihedL(k+j) , j = 1 , 2 )

   end do

   ii = ii * 2
   jj = NumDihedral - ii

   if( jj /= 0 )                                      &
   & write(1,'(8i8)') ( DihedI(ii+j) , DihedJ(ii+j) , &
   &                    DihedK(ii+j) , DihedL(ii+j) , j = 1 , jj )

   write(1,'(/i8,x,a)') NumImproper,'!NIMPHI: impropers'

   ii = NumImproper / 2

   do i = 1 , ii

     k = (i-1) * 2
     write(1,'(8i8)') ( ImproI(k+j) , ImproJ(k+j) , &
     &                  ImproK(k+j) , ImproL(k+j) , j = 1 , 2 )

   end do

   ii = ii * 2
   jj = NumImproper - ii

   if( jj /= 0 )                                     &
   & write(1,'(8i8)') ( ImproI(ii+j) , ImproJ(ii+j) , &
   &                    ImproK(ii+j) , ImproL(ii+j) , j = 1 , jj )

   j = 0

   write(1,'(/i8,x,a)') j,'!NDON: donors'

   write(1,'(/i8,x,a)') j,'!NACC: acceptors'

   close(1)

end subroutine Write_PSF
