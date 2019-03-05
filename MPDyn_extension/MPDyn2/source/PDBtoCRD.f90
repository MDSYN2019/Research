! ############################
! ## SUBROUTINE LIST 
! ## -- PDBtoCRD 
! ## -- SetCondition 
! ## -- Read_PDB 
! ## -- RenameAtom 
! ## -- Read_CRD 
! ## -- Write_CRD 
! ## -- Write_RB 
! ## -- Renumbering_to_topology 
! ## -- RB_assign 
! ## -- Def_Quaternion 
! ## -- Rotation_Quaternion 
! ## -- MinimEuler 
! ############################


! ######################################################################
! ######################################################################


program PDBtoCRD

use Numbers, only : N, Number
use CommonBlocks, only : QRigidBody, QInitial, QThermostat, QBarostat, ForceField
use IOparam, only : PSF_file
use BathParam, only : NHchain
use AtomParam, only : AtomName, ResidName, DefaultAtomName, DefaultResidName

implicit none

logical :: QCell, QPDB
character(len= 7) :: Flag
character(len= 3) :: cThermostat
character(len= 3) :: cBarostat
real(8), dimension(:,:), allocatable :: tesR
real(8) :: snt, cst, Theta
integer :: i
character(len=80) :: String, String1

   call SetCondition(QPDB)

   if( ForceField(1:4) == 'OPLS' ) then

     call Read_OPLS_Parameter
     write(*,*) 'OK:Parameter'
     call Read_OPLS_Topology
     write(*,*) 'OK:Topology'

   else if( ForceField(1:5) == 'CHARM' ) then

     call Read_Charmm_Parameter
     write(*,*) 'OK:Parameter'
     call Read_Charmm_Topology
     write(*,*) 'OK:Topology'

   else if( ForceField(1:3) == 'BKS' ) then

     call Read_BKS_Parameter
     write(*,*) 'OK:Parameter'
     call Read_OPLS_Topology
     write(*,*) 'OK:Topology'

   else if( ForceField(1:2) == 'CG' ) then

     call Read_CG_Parameter
     write(*,*) 'OK:Parameter'
     call Read_CG_Topology
     write(*,*) 'OK:Topology'

   end if

   call MemoryAllocation
   write(*,*) 'OK:Memory Allocation'

   N = Number

   if(ForceField(1:4) == 'OPLS' ) then

     call Read_PDB(QCell)
     write(*,*) 'OK:Read Conf.'

     call RenameAtom
     write(*,*) 'OK:Rename'

     call Renumbering_to_topology(1)
     write(*,*) 'OK:renumbering'

     call Write_CRD(QCell)

   else if(ForceField(1:3) == 'BKS') then

     call Read_PDB(QCell)
     write(*,*) 'OK:Read Conf.'

     call RenameAtom
     write(*,*) 'OK:Rename'

     call Renumbering_to_topology(1)
     write(*,*) 'OK:renumbering'

     call Write_CRD(QCell)

   else if(ForceField(1:5) == 'CHARM') then

     print *, 'CHARMM force field was selected! Is this correct?'
     read(5,*)
     print *, 'STATUS of your configuration file to be generated? [Initial or Restart]'
     read(5,*) Flag
     if(Flag=='Initial') then
       QInitial=.True.
     else if(Flag=='Restart') then
       QInitial=.False.
     else
       write(*,*) 'ERROR : Wrong flag on initial velocity'
     end if

     print *, 'Use thermostat? [ON or OFF] In case of ON, put the NH chain length in the same line'
     read(5,'(a80)') String1
     String = adjustl(trim(String1))
     read(String,*) cThermostat
     if(cThermostat=='ON') then
       QThermostat=.True.
       read(String,*) cThermostat, NHchain
     else if(cThermostat=='OFF') then
       QThermostat=.False.
     else
       write(*,*) 'ERROR : Thermostat condition'
     end if

     print *, 'Use barostat? [ON or OFF]'
     read(5,'(a)') cBarostat
     if(cBarostat=='ON') then
       QBarostat=.True.
     else if(cBarostat=='OFF') then
       QBarostat=.False.
     else
       write(*,*) 'ERROR : Barostat condition'
     end if

     print *, 'Periodic boundary condition? [ .T. or .F. ]'
     read(5,*) QCell
     write(*,*) QCell

     if(QPDB) then

       call Read_PDB(QCell)
       write(*,*) 'OK:Read Conf.'

       call RenameAtom
       write(*,*) 'OK:Rename'

       call Write_CRD(QCell)

     else

       call Read_CRD(QCell)
       write(*,*) 'OK:read crd'

     end if

     QThermostat=.False.

   else if(ForceField(1:2) == 'CG') then

     call Read_PDB(QCell)
     write(*,*) 'OK:Read Conf.'

     call RenameAtom
     write(*,*) 'OK:Rename'

     call Write_CRD(QCell)

   end if

   if(QRigidBody) then

     if(ForceField(1:4) == 'OPLS' ) then

       call MakeTopologicalData
       write(*,*) 'OK:TopologyData'

       call AllocateOPLS
       write(*,*) 'OK:OPLS'

       call Renumbering_to_topology(2)
       write(*,*) 'OK:renumbering2'

     else if(ForceField(1:5) == 'CHARM') then

       print *, 'put PSF filename'
       read(*,'(a)') PSF_File
       call Read_PSF
       write(*,*) 'OK:CHARMM'

     end if

     call RB_data
     write(*,*) 'OK:RB data'

     call RB_assign
     write(*,*) 'OK:RB assign'

     call Write_RB(QCell)
     write(*,*) 'OK:Write RB'

     call IntraMolVec
     write(*,*) 'OK:Intramolecular coordinates'

     call Write_PDB
     write(*,*) 'OK:write pdb file'

     AtomName  = DefaultAtomName
     ResidName = DefaultResidName

     call Write_CRD(QCell)

   end if

end program PDBtoCRD


!######################################################################
!######################################################################


subroutine SetCondition(QPDB)

use Numbers, only : NumSpec, NumMol, NumAtm, NumMer
use CommonBlocks, only : QMaster, QPBC, QThermostat, QSHAKE, QRigidBody, &
&   ForceField, PolyFlag
use IOparam, only : Parameter_file, Topology_file
use AtomParam, only : MolName

implicit none

integer :: i, ii, j, jj
character(len=80) :: String1, String
logical :: QPDB
logical, dimension(80) :: empty
integer,dimension(10) :: first, last

   QPBC = .False.
   QThermostat = .False.
   QSHAKE = .False.
   QMaster = .True.

   print *, 'Which force field do you use ?'
   read(5,*) ForceField

   print *, 'Name of PARAMETER file ? if you use a default file, write "D"'
   read(5,'(a80)') String1
   String = adjustl(String1)

   if(String(1:1)=='D') then
     if(ForceField(1:5)=='CHARM') then
     write(Parameter_file,*) 'param/par_charmm27.prm'
     else if(ForceField(1:4)=='OPLS') then
     write(Parameter_file,*) 'param/par_oplsaa.prm'
     else if(ForceField(1:3)=='BKS') then
     write(Parameter_file,*) 'param/par_BKS.prm'
     else if(ForceField(1:2)=='CG') then
     write(Parameter_file,*) '~/CGparam/par_CG.prm'
     end if
   else
     read(String,*) Parameter_file
   end if

   print *, 'Name of TOPOLOGY file ? if you use a default file, write "D"'
   read(5,'(a80)') String1
   String = adjustl(String1)

   if(String(1:1)=='D') then
     if(ForceField(1:5)=='CHARM') then
     write(Topology_file,*) 'param/top_charmm27.prm'
     else if(ForceField(1:4)=='OPLS') then
     write(Topology_file,*) 'param/top_oplsaa.prm'
     else if(ForceField(1:3)=='BKS') then
     write(Topology_file,*) 'param/top_BKS.prm'
     else if(ForceField(1:2)=='CG') then
     write(Topology_file,*) '~/CGparam/top_CG.prm'
     end if
   else
     read(String,*) Topology_file
   end if

   print *, 'rigid-body model ? [ .T. or .F. ]'
   read(5,*) QRigidBody

   QPDB = .False.
   if((ForceField(1:5) == 'CHARM').and.(.not.QRigidBody)) then

     print *, 'Convert : PDB --> CRD file'
     QPDB = .True.

   end if

   print *, 'number of molecular species ?'
   read(5,*) NumSpec

   allocate( MolName(NumSpec) )
   allocate( NumMol(NumSpec) )
   allocate( NumAtm(NumSpec) )

   allocate( NumMer(NumSpec) )
   allocate( PolyFlag(NumSpec) )

   print *, 'name of molecules ? write all in one line with space inbetween'
   read(5,'(a80)') String1
   String = trim(adjustl(String1))

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

   do i = 1 , NumSpec

     if( ForceField(1:4) == 'OPLS' ) then

       if(MolName(i)(1:4) == 'Poly') then

         PolyFlag(i) = .True.

       else

         PolyFlag(i) = .False.

       end if

     else

       PolyFlag(i) = .False.

     end if

   end do

   print *, 'number of molecules ? write all in one line.'
   read(5,'(a80)') String1
   String = trim(adjustl(String1))

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


   print *, 'number of atoms per molecule ? write all in one line.'
   read(5,'(a80)') String1
   String = trim(adjustl(String1))

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
     if(( ForceField(1:4) == 'OPLS' ).and.(PolyFlag(i))) then
       read(String(first(i):last(i)),*) NumMer(i)
     else
       read(String(first(i):last(i)),*) NumAtm(i)
     end if
   end do

end subroutine SetCondition


!######################################################################
!######################################################################


! **************************
! **  Read Configration   **
! **************************

subroutine Read_PDB(QCell)

use Numbers, only : N
use Configuration, only : R
use UnitExParam, only : pi
use CellParam, only : H, InvH
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer :: i, ii
character(len=80) :: String
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca
logical :: QCell
integer :: eofile

   QCell = .False.

   open( 7, file='Initial.pdb',status='old')

   ii = 0

   allocate( AtomName (N) )
   allocate( ResidName(N) )
   allocate( ResidNum (N) )

   print *, 'N=', N

   do

     read(7,'(a80)',iostat=eofile) String

     if( eofile == -1 ) exit
     if( String(1:3) == 'END' ) exit

     if((String(1:6) == 'HETATM') .or. (String(1:4) == 'ATOM')) then

       ii = ii + 1

       read(String,'(12x,a4,x,a4,x,i4,4x,3f8.3)') &
       &       AtomName(ii), ResidName(ii), ResidNum(ii), R(:,ii)

     end if

     if(String(1:6)=='CRYST1') then
       read(String,'(6x,3f9.3,3f7.2)') LLa,LLb,LLc,Abc,Aca,Aab
       QCell = .True.
     end if

   end do

   close(7)

   if(QCell) then

     H = 0.

     if((Aab/=90.).or.(Abc/=90.).or.(Aca/=90.)) then

       Aab = pi * Aab / 180.d0
       Abc = pi * Abc / 180.d0
       Aca = pi * Aca / 180.d0

       H(1,1) = LLa
       H(1,2) = LLb * cos(Aab)
       H(2,2) = LLb * sin(Aab)
       H(1,3) = LLc * cos(Aca)
       H(2,3) = LLc * cos(Abc)
       H(3,3) = LLc * sin(Aca) * sin(Abc)

     else

       H(1,1) = LLa
       H(2,2) = LLb
       H(3,3) = LLc

     end if

     call InversMatrix(H,InvH)

   end if

end subroutine Read_PDB


!######################################################################
!######################################################################


subroutine RenameAtom

use Numbers, only : N
use AtomParam, only : AtomName, ResidName, DefaultAtomName, DefaultResidName

implicit none

integer :: i
character(len=4) :: String, String1
character(len=1) :: Wk, Ch

   Wk = ' '

   do i = 1, N

     String1 = adjustl(AtomName(i))

     if((String1(2:2) == 'H').and.((String1(1:1) == '0').or.(String1(1:1) == '1') &
     &.or.(String1(1:1) == '2').or.(String1(1:1) == '3').or.(String1(1:1) == '4') &
     &.or.(String1(1:1) == '5').or.(String1(1:1) == '6').or.(String1(1:1) == '7') &
     &.or.(String1(1:1) == '8').or.(String1(1:1) == '9'))) then

       if(String1(4:4) == Wk) then
         String(1:2) = String1(2:3)
         String(3:3) = String1(1:1)
       else
         String(1:3) = String1(2:4)
         String(4:4) = String1(1:1)
       end if

     else

       String = String1

     end if

     AtomName(i) = String

   end do

! ## for check

   do i = 1 , N

     String = ResidName(i)

     Ch = String(1:1)

     if((Ch=='0').or.(Ch=='1').or.(Ch=='2').or.(Ch=='3').or.(Ch=='4').or. &
     &  (Ch=='5').or.(Ch=='6').or.(Ch=='7').or.(Ch=='8').or.(Ch=='9').or. &
     &  (Ch==Wk)) then

       write(*,*) 'ERROR : mistakes in the defined residue names.'
       write(*,*) 'particle = ',i
       stop

     end if

   end do

   allocate( DefaultAtomName (N) )
   allocate( DefaultResidName(N) )

   DefaultAtomName  = AtomName
   DefaultResidName = ResidName

end subroutine RenameAtom


!######################################################################
!######################################################################


! ***************************
! **  Write Configration   **
! ***************************

subroutine Read_CRD(QCell)

use Numbers, only : N
use CommonBlocks, only : QInitial, QThermostat
use Configuration, only : R
use BathParam, only : NHchain
use CellParam, only : H
use AtomParam, only : AtomName, ResidName, ResidNum, DefaultAtomName, DefaultResidName

implicit none

integer :: i, ii
character(len=4) :: ModelName
real(8), dimension(N) :: VV
real(8) :: DumR, DumV
logical :: QCell

     allocate( AtomName (N) )
     allocate( ResidName(N) )
     allocate( ResidNum (N) )

     if(QInitial) then

     open(7,file='initial.crd')

     read(7,*)
     read(7,*)
     read(7,*) N

     do i = 1, N

        read(7,'(2i5,2(x,a4),3f10.5,x,a4)')  &
       & ii  , ResidNum(i) , ResidName(i) ,  &
       & AtomName(i) , R(:,i) , ModelName

     end do

     if(QCell) then

       read(7,*) (H(1,i), i = 1, 3)
       read(7,*) (H(2,i), i = 1, 3)
       read(7,*) (H(3,i), i = 1, 3)

     end if

     else

     open(7,file='restart.dat',form='unformatted')

     read(7) N

     read(7) ResidNum
     read(7) ResidName
     read(7) AtomName
     read(7) R
     read(7) VV

     if(QThermostat) then
       do i = 1 , NHChain
         read(7) DumR, DumV
       end do
     end if

     if(QCell) then

       read(7) H

     end if

     end if

     close(7)

   allocate( DefaultAtomName (N) )
   allocate( DefaultResidName(N) )

   DefaultAtomName  = AtomName
   DefaultResidName = ResidName

end subroutine Read_CRD


!######################################################################
!######################################################################


! ***************************
! **  Write Configration   **
! ***************************

subroutine Write_CRD(QCell)

use Numbers, only : N
use Configuration, only : R
use CellParam, only : H
use AtomParam, only : AtomName, ResidName, ResidNum
use IOparam, only : DirectoryName
use CommonBlocks, only : QPBC

implicit none

integer :: i
character(len=4) :: ModelName
logical :: QCell

     ModelName = ResidName(1)

     open(7,file='initial.crd')

     write(7,'(a)') '! CRD file'
     write(7,'(a)') '!'
     write(7,'(i10)') N

     if(N>100000) then
       do i = 1, N
         write(7,'(2i9,2(x,a4),3f12.5)')  &
         & i  , ResidNum(i) , ResidName(i) ,  &
         & AtomName(i) , R(:,i)
       end do
     else
       do i = 1, N
         write(7,'(2i5,2(x,a4),3f10.5)')  &
         & i  , ResidNum(i) , ResidName(i) ,  &
         & AtomName(i) , R(:,i)
       end do
     end if

     if(QCell) then

       write(7,'(3f12.4)') (H(1,i), i = 1, 3)
       write(7,'(3f12.4)') (H(2,i), i = 1, 3)
       write(7,'(3f12.4)') (H(3,i), i = 1, 3)

     end if

     close(7)

     write(DirectoryName,'(a)') './'
     QPBC = QCell

     call Write_PDB


end subroutine Write_CRD


!######################################################################
!######################################################################


! ***************************
! **  Write Configration   **
! ***************************

subroutine Write_RB(QCell)

use RBparam, only : NumRB, R_RB, Quaternion
use CellParam, only : H

implicit none

integer :: i, j
logical :: QCell

open( 7, file='initialRB.dat', status='unknown')

   do i = 1 , NumRB

     write(7,'(3d24.16)') ( R_RB(j,i) , j = 1 , 3 )

   end do

   do i = 1 , NumRB

     write(7,'(4d24.16)') ( Quaternion(j,i) , j = 1 , 4 )

   end do

   if(QCell) then

     write(7,'(3f12.4)') (H(1,i), i = 1, 3)
     write(7,'(3f12.4)') (H(2,i), i = 1, 3)
     write(7,'(3f12.4)') (H(3,i), i = 1, 3)

   end if

close(7)

end subroutine Write_RB


!######################################################################
!######################################################################


subroutine Renumbering_to_topology(iflag)

use Numbers, only : N, NumSpec, NumMol, NumAtm, NumMer
use Configuration, only : R
use CommonBlocks, only : PolyFlag
use FFParameters
use AtomParam, only : MolName, AtomName, ResidNum, Mass

implicit NONE

integer :: i, j, k, l, kk, iflag
integer :: NResP

character(len=8) :: String, String1
character(len=4) :: AtomI, RName

integer :: ISearch, FSearch
integer :: IJSearch, FJSearch
integer :: JAtom, IAtom
logical :: ConfAllocI

character(len=4), dimension(N) :: TmpAtomName
real(8), dimension(3,N) :: TmpR
real(8), dimension(N) :: TmpMass

integer, dimension(:), allocatable :: OriginalOrder
save OriginalOrder

if(iflag == 1) then

   allocate( OriginalOrder(N) )

   do i = 1 , NumSpec

     if(i == 1) then
       ISearch = 0
     else
       ISearch = ISearch + NumMol(i-1)*NumAtm(i-1)
     end if

     FSearch = ISearch + NumAtm(i)

     if(PolyFlag(i)) then
       String1 = MolName(i)(5:8)
     else
       String1 = MolName(i)
     end if

     String = trim(adjustl(String1))

     RName = String(1:4)

! ## polymers

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

             IAtom = NumAtom_inResi(k)
             NResP = k

             exit

           end if

           if(k == NumResidueParam) then

             write(*,*) 'ERROR: polymer resid.'
             write(*,*) 'residue : ',RName
             stop

           end if

         end do

         FJSearch = IJSearch + IAtom

! ## Consistency check! ##

         do l = IJSearch + 2, FJSearch

           if(ResidNum(IJSearch+1) /= ResidNum(l)) then

             write(*,*) 'ERROR : the defined structure in the pdb file is not consistent'
             write(*,*) '        with the condition.'
             write(*,*) 'number : ',l,ResidNum(l)
             stop

           end if

         end do

! ## Name ##

         do k = 1 , IAtom

           kk = IJSearch + k

           AtomI = AtomNameParam(k,NResP)

           ConfAllocI = .False.

           do l = IJSearch + 1, FJSearch

             if(AtomI == AtomName(l)) then

               OriginalOrder(kk) = l

               ConfAllocI = .True.

               exit

             end if

             if((l == FJSearch).and.(.not.ConfAllocI)) then

               write(*,*) 'ERROR : allocation'
               write(*,*) 'Atom : ',AtomI,' Residue : ',RName
               stop

             end if

           end do

         end do

         JAtom = IAtom

       end do

     else

! ## non-polymer

       do k = 1 , NumResidueParam

         if(ResiNameParam(k) == RName) then

           IAtom = NumAtom_inResi(k)
           NResP = k

           exit

         end if

         if(k == NumResidueParam) then

           write(*,'(a,a)') 'ERROR: no residue / ',RName
           stop

         end if

       end do

! ## Consistency check! ##

       do l = ISearch + 2, FSearch

         if(ResidNum(ISearch+1) /= ResidNum(l)) then

           write(*,*) 'ERROR : the defined structure in the pdb file is not consistent'
           write(*,*) '        with the condition.'
           write(*,*) 'number : ',l,ResidNum(l)
           stop

         end if

       end do

! ## Name ##

       do k = 1 , IAtom

         kk = ISearch + k

         AtomI = AtomNameParam(k,NResP)

         ConfAllocI = .False.

         do l = ISearch + 1, FSearch

           if(AtomI == AtomName(l)) then

             OriginalOrder(kk) = l

             ConfAllocI = .True.

             exit

           end if

           if((l == FSearch).and.(.not.ConfAllocI)) then

             write(*,*) 'ERROR : allocation'
             write(*,*) 'Atom : ',AtomI,' Residue : ',RName
             stop

           end if

         end do

       end do

     end if


     do j = 2 , NumMol(i)

       k = (j-1)*NumAtm(i) + ISearch

       do l = 1 , NumAtm(i)

         kk = k + l
         OriginalOrder(kk) = OriginalOrder(ISearch+l) + (j-1)*NumAtm(i)

       end do

     end do


   end do


   do i = 1 , N

     j = OriginalOrder(i)

     TmpAtomName(i) = AtomName(j)
     TmpR(:,i) = R(:,j)

   end do

   R = TmpR
   AtomName = TmpAtomName

else if(iflag==2) then

   do i = 1 , N

     j = OriginalOrder(i)

     TmpMass(i) = Mass(j)

   end do

   Mass = TmpMass

end if

end subroutine Renumbering_to_topology


!######################################################################
!######################################################################


subroutine RB_assign

use Configuration, only : R
use RBparam, only : NumRB, R_RB, Quaternion, RBType, NumRBAtom, InvMassRB, &
&   R_onMol, QLinear, MaxNumAtomRB
use AtomParam, only : Mass

implicit NONE

integer :: i, j, k, Nc, MyType
real(8), dimension(3) :: Rg
real(8), dimension(3,MaxNumAtomRB) :: RonU, RonM
real(8), dimension(4) :: qq

   allocate( R_RB(3,NumRB) )
   allocate( Quaternion(4,NumRB) )

   open(10,file='MSDcheck.dat',status='unknown')

   k = 0

   do i = 1 , NumRB

     MyType = RBType(i)
     Nc = NumRBAtom(MyType)

     if(Nc == 1) then

       R_RB(:,i) = R(:,k + 1)
       Quaternion(:,i) = 0.d0

     else

       Rg = 0.d0

       do j = 1 , Nc

         Rg = Rg + Mass(k + j) * R(:,k + j)

       end do

       Rg = Rg * InvMassRB(MyType)
       R_RB(:,i) = Rg

       do j = 1 , Nc

         RonU(:,j) = R(:,k + j) - Rg
         RonM(:,j) = R_onMol(:,j,MyType)

       end do

!##       call MinimEuler(RonU,RonM,phi,theta,psi,Nc,QLinear(i))
!##
!##       Quaternion(1,i) = cos( 0.5*theta ) * cos( 0.5 * (phi + psi) )
!##       Quaternion(2,i) = sin( 0.5*theta ) * cos( 0.5 * (phi - psi) )
!##       Quaternion(3,i) = sin( 0.5*theta ) * sin( 0.5 * (phi - psi) )
!##       Quaternion(4,i) = cos( 0.5*theta ) * sin( 0.5 * (phi + psi) )
!##
!##      call Def_Quaternion(RonU,RonM,qq,Nc,QLinear(i))

       call Rotation_Quaternion(RonU,RonM,qq,Nc,QLinear(i),i,k+1)

       Quaternion(:,i) = qq

     end if

     k = k + Nc

   end do

!   do i = 1 , NumRB
!     if(QSingle(i)) cycle
!     q1 = sqrt( dot_product(Quaternion(:,i),Quaternion(:,i)) )
!     Quaternion(:,i) = Quaternion(:,i) / q1
!     write(*,*) i,'QQ=', q1
!   end do

   close(10)

end subroutine RB_assign


!######################################################################
!######################################################################


subroutine Def_Quaternion(RonU,RonM,qq,Nc,QLimit)

use RBparam

implicit none

integer :: i, j, k
logical :: QLimit
real(8), dimension(3,MaxNumAtomRB) :: RonU ! space-fixed
real(8), dimension(3,MaxNumAtomRB) :: RonM ! body-fixed
real(8), dimension(4) :: qq
real(8), dimension(MaxNumAtomRB,MaxNumAtomRB) :: dist
integer, dimension(3) :: MaximumDist
integer :: Nc
real(8), dimension(3,MaxNumAtomRB) :: Rtmp
integer :: MinN, MaxN
real(8) :: MinX, MaxX, sng
real(8) :: phi,theta,psi
real(8) :: csph,snph
real(8) :: csth,snth
real(8), dimension(3) :: Rvec, dR
real(8), dimension(3,3) :: Es, InvEs, Rot
real(8), dimension(3) :: Ebx, Eby, Ebz
real(8) :: S, TrRot, MaxR, dis, Mdis
integer :: MaxI

   if(QLimit) then

     MaxN = 1
     MinN = 1

     MaxX = RonM(1,1)
     MinX = RonM(1,1)

     do i = 2, Nc

       if(RonM(1,i)>MaxX) then

         MaxX = RonM(1,i)
         MaxN = i

       else if(RonM(1,i)<MinX) then

         MinX = RonM(1,i)
         MinN = i

       end if

     end do

     Rvec = RonU(:,MaxN) - RonU(:,MinN)
     Rvec = Rvec / sqrt(dot_product(Rvec, Rvec) )

     csth = Rvec(3)
     snth = sqrt( 1.d0 - csth*csth )

     if(snth/=0.) then
       csph = Rvec(1) / snth
       snph = Rvec(2) / snth
     else
       csph = 1.d0
       snph = 0.d0
     end if

! ## phi

     phi = acos( csph )
     sng = sngl(sin( phi ))

     if(sng /= sngl(snph)) then
       if(sng == -sngl(snph)) then
         phi  = - phi
         snph = - snph
       else
         write(*,*) 'WARNING : rotation error'
       end if
     end if

     Rot = 0.d0

     Rot(1,1) =   csph
     Rot(1,2) =   snph
     Rot(2,1) = - snph
     Rot(2,2) =   csph
     Rot(3,3) =   1.d0

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

! ## theta

     theta = acos( csth )
     sng   = sngl( snth )

     if( sng /= sngl(snth) ) then
       if(sng == -sngl(snth)) then
         theta = - theta
         snth  = - snth
       else
         write(*,*) 'WARNING : rotation error'
       end if
     end if

     Rot = 0.d0

     Rot(1,1) =   1.d0
     Rot(2,2) =   csth
     Rot(2,3) =   snth
     Rot(3,2) = - snth
     Rot(3,3) =   csth

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

     do i = 1 , Nc
       print *, 'body',RonM(:,i)
       print *, 'RotS',RonU(:,i)
     end do

     psi = 0.d0

     qq(1) = cos( 0.5*theta ) * cos( 0.5 * (phi + psi) )
     qq(2) = sin( 0.5*theta ) * cos( 0.5 * (phi - psi) )
     qq(3) = sin( 0.5*theta ) * sin( 0.5 * (phi - psi) )
     qq(4) = cos( 0.5*theta ) * sin( 0.5 * (phi + psi) )

   else

     if(Nc > 3) then

       do i = 1 , Nc - 1
         do j = i+1 , Nc
           dR = RonM(:,i) - RonM(:,j)
           dist(i,j) = dot_product(dR,dR)
         end do
       end do

       Mdis = 0.
       do i = 1 , Nc - 2
         do j = i + 1 , Nc - 1
           do k = j + 1 , Nc
             dis = dist(i,j) + dist(j,k) + dist(i,k)
             if(dis > Mdis) then
               Mdis = dis
               MaximumDist(1) = i
               MaximumDist(2) = j
               MaximumDist(3) = k
             end if
           end do
         end do
       end do

     else

       MaximumDist(1) = 1
       MaximumDist(2) = 2
       MaximumDist(3) = 3

     end if

     Ebx(:) = RonM(1,MaximumDist(:))
     Eby(:) = RonM(2,MaximumDist(:))
     Ebz(:) = RonM(3,MaximumDist(:))

     Es(:,1) = RonU(:,MaximumDist(1))
     Es(:,2) = RonU(:,MaximumDist(2))
     Es(:,3) = RonU(:,MaximumDist(3))

     call InversMatrix(Es,InvEs)

!     print *, Es
!     print *, InvEs

     Rot(1,:) = matmul(InvEs,Ebx)
     Rot(2,:) = matmul(InvEs,Eby)
     Rot(3,:) = matmul(InvEs,Ebz)

     TrRot = Rot(1,1) + Rot(2,2) + Rot(3,3) + 1.d0

     if(TrRot > 0.d0) then

       S = 0.5d0 / sqrt(TrRot)

       qq(1) = 0.25d0 / S
       qq(2) = (Rot(2,3) - Rot(3,2)) * S
       qq(3) = (Rot(3,1) - Rot(1,3)) * S
       qq(4) = (Rot(1,2) - Rot(2,1)) * S

     else

       MaxR = Rot(1,1)
       MaxI = 1

       do i = 2, 3

         if(Rot(i,i) > MaxR) then

           MaxR = Rot(i,i)
           MaxI = i

         end if

       end do

       select case(MaxI)

       case(1)

         S = sqrt( 1.d0 + Rot(1,1) - Rot(2,2) - Rot(3,3) ) * 2.d0

         qq(1) = (Rot(2,3) - Rot(3,2)) / S
         qq(2) = 0.25 * S
         qq(3) = (Rot(1,2) + Rot(2,1)) / S
         qq(4) = (Rot(1,3) + Rot(3,1)) / S

       case(2)

         S = sqrt( 1.d0 - Rot(1,1) + Rot(2,2) - Rot(3,3) ) * 2.d0

         qq(1) = (Rot(3,1) - Rot(1,3)) / S
         qq(2) = (Rot(1,2) + Rot(2,1)) / S
         qq(3) = 0.25 * S
         qq(4) = (Rot(2,3) + Rot(3,2)) / S

       case(3)

         S = sqrt( 1.d0 - Rot(1,1) - Rot(2,2) + Rot(3,3) ) * 2.d0

         qq(1) = (Rot(1,2) - Rot(2,1)) / S
         qq(2) = (Rot(1,3) + Rot(3,1)) / S
         qq(3) = (Rot(2,3) + Rot(3,2)) / S
         qq(4) = 0.25 * S

       end select

     end if

   end if


end subroutine Def_Quaternion


!######################################################################
!######################################################################


subroutine Rotation_Quaternion(RonU,RonM,qq,Nc,QLimit,ii,jj)

use Configuration, only : R
use RBparam
use UnitExParam, only : pi

implicit none

integer :: i, j, k, l, Count, ii, jj
logical :: QLimit
real(8), dimension(3,MaxNumAtomRB) :: RonU ! space-fixed
real(8), dimension(3,MaxNumAtomRB) :: RonM ! body-fixed
real(8), dimension(3,MaxNumAtomRB) :: Rtmp
real(8), dimension(4) :: qq
integer, dimension(3) :: MinDis
integer :: Nc
real(8) :: MSD, MinMSD
real(8) :: phi,theta,psi
real(8), dimension(3) :: dR, Rsave
integer :: IIa, JJa, KKa, IIb, JJb, KKb, IIc, JJc, KKc

   if(QLimit) then

     Count = 0

     IIa = -180
     JJa =  180
     KKa =   30

     call SearchAngle2(IIa,JJa,KKa,IIa,JJa,KKa)

     IIa = MinDis(1)-30
     JJa = MinDis(1)+30
     KKa = 10
     IIb = MinDis(2)-30
     JJb = MinDis(2)+30
     KKb = 10

     call SearchAngle2(IIa,JJa,KKa,IIb,JJb,KKb)

     IIa = MinDis(1)-10
     JJa = MinDis(1)+10
     KKa = 2
     IIb = MinDis(2)-10
     JJb = MinDis(2)+10
     KKb = 2

     call SearchAngle2(IIa,JJa,KKa,IIb,JJb,KKb)

     IIa = MinDis(1)-2
     JJa = MinDis(1)+2
     KKa = 1
     IIb = MinDis(2)-2
     JJb = MinDis(2)+2
     KKb = 1

     call SearchAngle2(IIa,JJa,KKa,IIb,JJb,KKb)

! ## check

     call RotCalc(MinDis(1),MinDis(2),0)

     call MSDc(MSD)
     MSD = MSD / dble(Nc)
     write(10,*) ii, sqrt(MSD)

     phi   = dble(MinDis(1)) * pi / 180.
     theta = dble(MinDis(2)) * pi / 180.
     psi = 0.d0

     qq(1) = cos( 0.5*theta ) * cos( 0.5 * (phi + psi) )
     qq(2) = sin( 0.5*theta ) * cos( 0.5 * (phi - psi) )
     qq(3) = sin( 0.5*theta ) * sin( 0.5 * (phi - psi) )
     qq(4) = cos( 0.5*theta ) * sin( 0.5 * (phi + psi) )

   else

     Count = 0

     IIa = -180
     JJa =  180
     KKa =   30

     call SearchAngle(IIa,JJa,KKa,IIa,JJa,KKa,IIa,JJa,KKa)

     IIa = MinDis(1)-30
     JJa = MinDis(1)+30
     KKa = 10
     IIb = MinDis(2)-30
     JJb = MinDis(2)+30
     KKb = 10
     IIc = MinDis(3)-30
     JJc = MinDis(3)+30
     KKc = 10

     call SearchAngle(IIa,JJa,KKa,IIb,JJb,KKb,IIc,JJc,KKc)

     IIa = MinDis(1)-10
     JJa = MinDis(1)+10
     KKa = 2
     IIb = MinDis(2)-10
     JJb = MinDis(2)+10
     KKb = 2
     IIc = MinDis(3)-10
     JJc = MinDis(3)+10
     KKc = 2

     call SearchAngle(IIa,JJa,KKa,IIb,JJb,KKb,IIc,JJc,KKc)

     IIa = MinDis(1)-2
     JJa = MinDis(1)+2
     KKa = 1
     IIb = MinDis(2)-2
     JJb = MinDis(2)+2
     KKb = 1
     IIc = MinDis(3)-2
     JJc = MinDis(3)+2
     KKc = 1

     call SearchAngle(IIa,JJa,KKa,IIb,JJb,KKb,IIc,JJc,KKc)

! ## check

     MinMSD = MinMSD / dble(Nc)

     if( MinMSD > 0.05 ) then

       write(*,*) 're-calculating!',ii, AtomUnitName(jj)

       IIa = -180
       JJa =  180
       KKa =    3

       Count = 0

       call SearchAngle(IIa,JJa,KKa,IIa,JJa,KKa,IIa,JJa,KKa)

       MinMSD = MinMSD / dble(Nc)

     end if

     if( ( MinMSD > 0.05 ).and.( AtomUnitName(jj) == 'METL' ) ) then

       Rsave     = R(:,jj+2)
       R(:,jj+2) = R(:,jj+3)
       R(:,jj+3) = Rsave

       Rsave     = RonU(:,3)
       RonU(:,3) = RonU(:,4)
       RonU(:,4) = Rsave

       write(*,*) 'exchange the tags for H1 and H2 in the unit METL'

       IIa = -180
       JJa =  180
       KKa =    3

       Count = 0

       call SearchAngle(IIa,JJa,KKa,IIa,JJa,KKa,IIa,JJa,KKa)

       MinMSD = MinMSD / dble(Nc)

     end if

     if( MinMSD > 0.05 ) then

       write(*,*) ' Euler angles have not converged yet ! '
       write(*,*) ' Rigid-Body = ', ii, MinMSD
       write(*,*) ' Check the configuration ! '

       write(*,*) ' MSD optimization           : ',MinDis(1),MinDis(2),MinDis(3)

       call RotCalc(MinDis(1),MinDis(2),MinDis(3))

       do l = 1 , Nc
         write(*,'(i3,3f8.3)') l,Rtmp(:,l)
       end do

       write(*,*) ' Original Cartesian coordinates '

       do l = 1 , Nc
         write(*,'(i3,3f8.3)') l,RonU(:,l)
       end do

       write(*,*) ' Original body-fixed coordinates '

       do l = 1 , Nc
         write(*,'(i3,3f8.3)') l,RonM(:,l)
       end do

       stop

     end if

! ## check

     call RotCalc(MinDis(1),MinDis(2),MinDis(3))

     call MSDc(MSD)

     MSD = MSD / dble(Nc)
     write(10,*) ii, sqrt(MSD)

     phi   = dble(MinDis(1)) * pi / 180.
     theta = dble(MinDis(2)) * pi / 180.
     psi   = dble(MinDis(3)) * pi / 180.

     qq(1) = cos( 0.5*theta ) * cos( 0.5 * (phi + psi) )
     qq(2) = sin( 0.5*theta ) * cos( 0.5 * (phi - psi) )
     qq(3) = sin( 0.5*theta ) * sin( 0.5 * (phi - psi) )
     qq(4) = cos( 0.5*theta ) * sin( 0.5 * (phi + psi) )

   end if

Contains

   subroutine SearchAngle2(MinA,MaxA,dA,MinB,MaxB,dB)

   integer :: MinA, MaxA, dA
   integer :: MinB, MaxB, dB

     do i = MinA, MaxA, dA

       do j = MinB, MaxB, dB

         call RotCalc(i,j,0)

         Count = Count + 1

         call MSDc(MSD)

         if(Count == 1) then
           MinMSD = MSD
           MinDis(1) = i
           MinDis(2) = j
         else if(MSD < MinMSD) then
           MinMSD = MSD
           MinDis(1) = i
           MinDis(2) = j
         end if

       end do

     end do

   end subroutine SearchAngle2


   subroutine SearchAngle(MinA,MaxA,dA,MinB,MaxB,dB,MinC,MaxC,dC)

   integer :: MinA, MaxA, dA
   integer :: MinB, MaxB, dB
   integer :: MinC, MaxC, dC

     do i = MinA, MaxA, dA

       do j = MinB , MaxB, dB

         do k = MinC , MaxC, dC

           call RotCalc(i,j,k)

           Count = Count + 1

           call MSDc(MSD)

           if(Count == 1) then
             MinMSD = MSD
             MinDis(1) = i
             MinDis(2) = j
             MinDis(3) = k
           else if(MSD < MinMSD) then
             MinMSD = MSD
             MinDis(1) = i
             MinDis(2) = j
             MinDis(3) = k
           end if

         end do

       end do

     end do

   end subroutine SearchAngle

   subroutine RotCalc(i1,j1,k1)

   integer :: i1, j1, k1
   real(8), dimension(10) :: qt2
   real(8), dimension(3,3) :: Rot

     phi   = dble(i1) * pi / 180.
     theta = dble(j1) * pi / 180.
     psi   = dble(k1) * pi / 180.

     qq(1) = cos( 0.5*theta ) * cos( 0.5 * (phi + psi) )
     qq(2) = sin( 0.5*theta ) * cos( 0.5 * (phi - psi) )
     qq(3) = sin( 0.5*theta ) * sin( 0.5 * (phi - psi) )
     qq(4) = cos( 0.5*theta ) * sin( 0.5 * (phi + psi) )

     qt2( 1) = qq(1) * qq(1)
     qt2( 2) = qq(2) * qq(2)
     qt2( 3) = qq(3) * qq(3)
     qt2( 4) = qq(4) * qq(4)
     qt2( 5) = qq(1) * qq(2)
     qt2( 6) = qq(1) * qq(3)
     qt2( 7) = qq(1) * qq(4)
     qt2( 8) = qq(2) * qq(3)
     qt2( 9) = qq(2) * qq(4)
     qt2(10) = qq(3) * qq(4)

     Rot(1,1) = qt2(1) + qt2(2) - qt2(3) - qt2(4)
     Rot(2,2) = qt2(1) - qt2(2) + qt2(3) - qt2(4)
     Rot(3,3) = qt2(1) - qt2(2) - qt2(3) + qt2(4)
     Rot(1,2) = 2.d0 * ( qt2( 8) + qt2( 7) )
     Rot(2,1) = 2.d0 * ( qt2( 8) - qt2( 7) )
     Rot(1,3) = 2.d0 * ( qt2( 9) - qt2( 6) )
     Rot(3,1) = 2.d0 * ( qt2( 9) + qt2( 6) )
     Rot(2,3) = 2.d0 * ( qt2(10) + qt2( 5) )
     Rot(3,2) = 2.d0 * ( qt2(10) - qt2( 5) )

     do l = 1 , Nc

       Rtmp(1,l) = dot_product( Rot(:,1) , RonM(:,l) )
       Rtmp(2,l) = dot_product( Rot(:,2) , RonM(:,l) )
       Rtmp(3,l) = dot_product( Rot(:,3) , RonM(:,l) )

     end do

   end subroutine RotCalc


   subroutine MSDc(MSD)

   real(8) :: MSD

     MSD = 0.d0
     do l = 1 , Nc
       dR = RonU(:,l) - Rtmp(:,l)
       MSD = MSD + dot_product(dR,dR)
     end do

   end subroutine

end subroutine Rotation_Quaternion


!######################################################################
!######################################################################


subroutine MinimEuler(RonU,RonM,phi,theta,psi,Nc,QLimit)

use RBparam
use UnitExParam, only : pi

implicit none

integer :: i, Ir, Nc
real(8) :: fb, fc
real(8) :: phi,theta,psi
real(8) :: csph,snph,tnph
real(8) :: csth,snth,tnth
real(8) :: csps,snps,tnps
real(8), dimension(3,MaxNumAtomRB) :: RonU ! space-fixed
real(8), dimension(3,MaxNumAtomRB) :: RonM ! body-fixed
real(8), dimension(3,3) :: Rot
real(8), dimension(3,MaxNumAtomRB) :: Rtmp
real(8) :: Z1, Z2
logical :: QLimit
integer :: MinN, MaxN
real(8) :: MinX, MaxX, sng
real(8), dimension(3) :: Rvec

   if(QLimit) then

     MaxN = 1
     MinN = 1

     MaxX = RonM(1,1)
     MinX = RonM(1,1)

     do i = 2, Nc

       if(RonM(1,i)>MaxX) then

         MaxX = RonM(1,i)
         MaxN = i

       else if(RonM(1,i)<MinX) then

         MinX = RonM(1,i)
         MinN = i

       end if

     end do

     Rvec = RonU(:,MaxN) - RonU(:,MinN)
     Rvec = Rvec / sqrt(dot_product(Rvec, Rvec) )

     csth = Rvec(3)
     snth = sqrt( 1.d0 - csth*csth )

     if(snth/=0.) then
       csph = Rvec(1) / snth
       snph = Rvec(2) / snth
     else
       csph = 1.d0
       snph = 0.d0
     end if

! ## phi

     phi = acos( csph )
     sng = sngl(sin( phi ))

     if(sng /= sngl(snph)) then
       if(sng == -sngl(snph)) then
         phi  = - phi
         snph = - snph
       else
         write(*,*) 'WARNING : rotation error'
       end if
     end if

     Rot = 0.d0

     Rot(1,1) =   csph
     Rot(1,2) =   snph
     Rot(2,1) = - snph
     Rot(2,2) =   csph
     Rot(3,3) =   1.d0

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

! ## theta

     theta = acos( csth )
     sng   = sngl( snth )

     if( sng /= sngl(snth) ) then
       if(sng == -sngl(snth)) then
         theta = - theta
         snth  = - snth
       else
         write(*,*) 'WARNING : rotation error'
       end if
     end if

     Rot = 0.d0

     Rot(1,1) =   1.d0
     Rot(2,2) =   csth
     Rot(2,3) =   snth
     Rot(3,2) = - snth
     Rot(3,3) =   csth

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

     do i = 1 , Nc
       print *, 'body',RonM(:,i)
       print *, 'RotS',RonU(:,i)
     end do

     psi = 0.d0

   else

! ## phi

     fb = 0.d0
     fc = 0.d0
     do i = 1, Nc
       Ir = int(dot_product(RonU(:,i),RonU(:,i))*10.)
       if(Ir==0) cycle
       fb = fb + (   RonM(2,i)*RonU(1,i) - RonM(1,i)*RonU(2,i) )
       fc = fc + ( - RonM(1,i)*RonU(1,i) - RonM(2,i)*RonU(2,i) )
     end do

     tnph = fb / fc
     phi  = atan( tnph )

     snph = sin( phi )
     csph = cos( phi )

     Z1 =   fb * snph + fc * csph
     Z2 = - fb * snph - fc * csph

     if(Z2 < Z1) then
       snph = - snph
       csph = - csph
       phi  = phi + pi
     end if

     Rot = 0.d0

     Rot(1,1) =   csph
     Rot(1,2) =   snph
     Rot(2,1) = - snph
     Rot(2,2) =   csph
     Rot(3,3) =   1.d0

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

! ## theta

     fb = 0.d0
     fc = 0.d0
     do i = 1, Nc
       fb = fb + (   RonM(3,i)*RonU(2,i) - RonM(2,i)*RonU(3,i) )
       fc = fc + ( - RonM(2,i)*RonU(2,i) - RonM(3,i)*RonU(3,i) )
     end do

     tnth  = fb / fc
     theta = atan( tnth )

     snth = sin( theta )
     csth = cos( theta )

     Z1 =   fb * snth + fc * csth
     Z2 = - fb * snth - fc * csth

     if(Z2 < Z1) then
       snth = - snth
       csth = - csth
       theta = theta + pi
     end if

     Rot = 0.d0

     Rot(1,1) =   1.d0
     Rot(2,2) =   csth
     Rot(2,3) =   snth
     Rot(3,2) = - snth
     Rot(3,3) =   csth

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

! ## psi

     fb = 0.d0
     fc = 0.d0
     do i = 1, Nc
       fb = fb + (   RonM(2,i)*RonU(1,i) - RonM(1,i)*RonU(2,i) )
       fc = fc + ( - RonM(1,i)*RonU(1,i) - RonM(2,i)*RonU(2,i) )
     end do

     tnps = fb / fc
     psi  = atan( tnps )

     snps = sin( psi )
     csps = cos( psi )

     Z1 =   fb * snps + fc * csps
     Z2 = - fb * snps - fc * csps

     if(Z2 < Z1) then
       snps = - snps
       csps = - csps
       psi  = psi + pi
     end if

     Rot = 0.d0

     Rot(1,1) =   csps
     Rot(1,2) =   snps
     Rot(2,1) = - snps
     Rot(2,2) =   csps
     Rot(3,3) =   1.d0

     do i = 1 , Nc
       Rtmp(:,i) = matmul( Rot, RonU(:,i) )
     end do

     do i = 1 , Nc
       RonU(:,i) = Rtmp(:,i)
     end do

   end if

end subroutine MinimEuler
