! ############################
! ## SUBROUTINE LIST 
! ## --  Gene_Config_CG 
! ## --  Read_Inputs 
! ## --  Read_Parameters 
! ## --  NameAssign 
! ## --  BuildMolecule 
! ## --  GeneSingleMol_Vacuume 
! ## --  Read_Library 
! ## --  ReadLib 
! ## --  MakeBox 
! ## --  Write_CRD 
! ############################


!######################################################################
!######################################################################


module MolLib
   integer :: NMolLib
   logical :: QLibrary, QPartSt, Qnopart
   integer :: NumSpecies, NumStSpec, idir
   real(8) :: Xmin, Xmax
   real(8), dimension(:,:,:,:), allocatable :: Rlib
   integer, dimension(:), allocatable :: molparnum
   character(len=80) :: StructureFile
   integer :: fileflag
end module MolLib


!######################################################################
!######################################################################


program Gene_Config_CG

use MolLib

implicit none

real(8) :: Rho

   call PreRandomNumber
   call Read_Inputs(Rho)
   call Read_Parameters
   if(QLibrary) then
     call BuildMolecule
   else
     call Read_Library
     call MakeBox(Rho)
     call PBC
     call Write_CRD
     call Write_PDB
   end if

end program Gene_Config_CG


!######################################################################
!######################################################################


subroutine Read_Inputs(Rho)

use Numbers, only : NumSpec, NumMol, NumAtm
use MolLib
use CommonBlocks, only : QMaster, QGeneconf, QRigidBody, QSHAKE, &
     &   QPathInt, ForceField, QPBC, PolyFlag
use IOparam, only : Topology_file, Parameter_file
use CellParam, only : CellL
use AtomParam, only : MolName

implicit none

integer :: i, Count, filelen, j, k
character(len=80) :: String, String1
logical :: Qdone
character(len=1) :: Ch, Cdir
character(len=1), dimension(80) :: Stbr
real(8) :: Rho

   QGeneconf = .True.

   QRigidBody = .False.
   QMaster = .True.
   QSHAKE = .False.
   Qdone = .False.
   QPathInt = .False.
   ForceField = 'CG'
   QPBC = .False.
   QLibrary= .False.
   QPartSt = .False.
   NumStSpec = 0

   Count = 0

   do

     if(Qdone) exit

     select case(Count)

     case(0)
       print '(/a/)', 'Do you want to make a library for a molecule? [Y or N]'
     case(1)
       print '(/a)', 'How many molecular structures do you (want to) have in the library?'
       print '(a/)', "When you DON'T need the molecular library, put any POSITIVE number"
     case(2)
       print '(/a/)', 'Do you have have a reference structure for some components of the system? [Y or N]'
     case(3)
       print '(/a)', 'Name of PARAMETER file = ?'
       print '(a/)', 'Put "D" for a default choice, "~/CGparam/par_CG.prm"'
     case(4)
       print '(/a)', 'Name of TOPOLOGY file = ?'
       print '(a/)', 'Put "D" for a default choice, "~/CGparam/top_CG.prm"'
     case(5)
       print '(/a/)', 'Number of molecular species = ?'
     case(6)
       print '(/a/)', 'Molecular names ? Give them in one line.'
     case(7)
       print '(/a/)', 'Number of Molecules ? Give them in one line.'
     case(8)
       print '(/a/)', 'Number of Atoms in a Molecule ? Give them in one line.'
     case(9)
       print '(/a/)', 'Periodic Boundary Condition is used or not? [Y or N]'
     case(10)
       print '(/a)', 'Density { rho [g/cm^3] } or the size of the simulation box{ Lx, Ly, Lz [A] }?'
       print '(a/)', '**Giving a single value, it will be understood as density**'
     case(11)
       print '(/a)', '########## Solvation option ###########'
       print '(a)' , '    A reference structure is needed    '
       print '(a/)', '#######################################'
       print '(/a/)', 'Number of the molecular components you have structure ?'
       print '(/a)', 'Note: if you put "2", the first TWO molecular species are'
       print '(a/)', '     considered to have strctural data in a prepared file'
     case(12)
       print '(/a)', 'The filename containing the reference structural data?'
       print '(a/)', 'It should be in PDB or CRD formated'
     case(13)
       print '(/a)', 'Do you have a region to have no particle? [Y or N]'
       print '(a/)', 'e.g., no water in the lipid membrane core.'
     case(14)
       print '(/a)', 'Specify the region where no particle should be.'
       print '(a)' , 'e.g., "Z -20 20" means the region is Z=-20to20A'
       print '(a/)', 'Use "X" or "Y" or "Z" to specify the axis.'

     end select

     read( 5,'(a80)') String1

     String = trim(adjustl(String1))

     if((String(1:1) == '!').or.(String(1:1) == '#')) cycle

     if((String(1:5) == '<end>').or.(String(1:5) == '<END>')) exit

     Count = Count + 1

     select case(Count)

     case(1)

       Ch = String(1:1)
       if(Ch=='y'.or.Ch=='Y') then
         QLibrary = .True.
       end if

     case(2)

       read(String,*) NMolLib

     case(3)

       Ch = String(1:1)
       if(Ch=='y'.or.Ch=='Y') then
          QPartSt = .True.
       end if

     case(4)

       if(String=='D') then
         write(Parameter_file,*) '~/CGparam/par_CG.prm'
       else
         read(String,*) Parameter_file
       end if

     case(5)

       if(String=='D') then
         write(Topology_file,*) '~/CGparam/top_CG.prm'
       else
         read(String,*) Topology_file
       end if

     case(6)

       read(String,*) NumSpec

       allocate( MolName(NumSpec) )
       allocate( NumMol(NumSpec) )
       allocate( NumAtm(NumSpec) )
       allocate( PolyFlag(NumSpec) )

     case(7)

! ## Name of molecules

       read(String,*) (MolName(i),i=1,NumSpec)
       if(QLibrary) Count = Count + 1

     case(8)

! ## Number of Molecules

       read(String,*) (NumMol(i),i=1,NumSpec)

     case(9)

! ## Number of Molecules

       read(String,*) (NumAtm(i),i=1,NumSpec)
       if(QLibrary) Qdone = .True.

     case(10)

       Ch = String(1:1)
       if(Ch=='y'.or.Ch=='Y') then
          QPBC = .True.
       else if(QPartSt) then
          Count = Count + 1
       else
          Qdone = .True.
       end if

     case(11)

       if(QPBC) then
         
         do j = 1, 80
           Stbr(j) =  String(j:j)
         end do
         k = 0
         do j = 2, 80
           if(Stbr(j)==' '.and.Stbr(j-1)/=' ') then
             k = k + 1
           end if
         end do
         if(k>1) then
           Rho = 0.d0
           read(String,*) CellL(:)
         else
           read(String,*) Rho
           CellL(:) = 0.
         end if
       end if
       if(.not.QPartSt) exit

     case(12)

       read(String,*) NumStSpec
       if(NumStSpec>=NumSpec) then
         print *, "The number doesn't make sense!"
         Count = Count - 1
       end if

     case(13)

       read(String,*) StructureFile
       filelen = len_trim(StructureFile)
       if(StructureFile(filelen-2:filelen)=='PDB'.or. &
       &  StructureFile(filelen-2:filelen)=='pdb') then
         fileflag = 1
       else if(StructureFile(filelen-2:filelen)=='CRD'.or. &
       &       StructureFile(filelen-2:filelen)=='crd') then
         fileflag = 2
       else
         print *, 'only PDB or CRD file is acceptable'
         Count = Count - 1
       end if

     case(14)

       Ch = String(1:1)
       if(Ch=='y'.or.Ch=='Y') then
          Qnopart = .True.
       else
         Qdone = .True.
       end if

     case(15)

       read(String,*) Ch, Xmin, Xmax
       Cdir = Ch
       if(Ch=='X'.or.Ch=='x') then
         idir = 1
         Qdone = .True.
       else if(Ch=='Y'.or.Ch=='y') then
         idir = 2
         Qdone = .True.
       else if(Ch=='Z'.or.Ch=='z') then
         idir = 3
         Qdone = .True.
       else
         write(*,*) 'X or Y or Z should be used to give the prohibited region'
         Count = Count - 1
       end if

     end select

   end do

   write(*,*) 'Parameter file = ',trim(Parameter_file)
   write(*,*) 'Topology  file = ',trim(Topology_file)
   write(*,*) 'The number of molecule types = ',NumSpec
   if(QLibrary) then
     do i = 1, NumSpec
       write(*,'(a,i2,a,a)') ' Molecule ',i,' = ',MolName(i)
     end do
     NumMol(:) = 1
   else
     do i = 1, NumSpec
       write(*,'(a,i2,a,a,2i6)') ' Molecule ',i,' = ',MolName(i),NumMol(i),NumAtm(i)
     end do
   end if

   if(QPBC.and.Rho==0.) write(*,'(a,3f10.2)') ' Cell size = ', CellL(:)

   if(QPartSt) then
     write(*,'(2a,i1,a)') trim(StructureFile),' shuold contain the first ',NumStSpec,&
     &                    ' molecular-components configurations'
   end if

   if(Qnopart) then
     write(*,'(3a,f5.1,a,f5.1)') &
     & 'No particle will be placed in ',Cdir,', from ',Xmin,' to ',Xmax
   end if

   PolyFlag = .False.

end subroutine Read_Inputs


!######################################################################
!######################################################################


subroutine Read_Parameters

   call Read_CG_Parameter
   print *, 'OK:Parameter'
   call Read_CG_Topology 
   print *, 'OK:Topology'
   call MemoryAllocation
   print *, 'OK:Memory'
   call NameAssign
   print *, 'OK:Name'
   call AllocateCGdata
   print *, 'OK:ParameterAssign'

end subroutine Read_Parameters


!######################################################################
!######################################################################


subroutine NameAssign

use Numbers, only : N, NumSpec, NumMol, NumAtm
use MolLib
use CommonBlocks, only : QMaster
use CGParameters
use AtomParam, only : MolName, AtomName, ResidName, ResidNum

implicit none

integer :: i, j, k, ii, jj, NAtom
integer :: NResd, l

   allocate( molparnum(NumSpec) )

   N = 0
   do i = 1, NumSpec
     N = N + NumMol(i) * NumAtm(i)
   end do

   allocate( AtomName (N) )
   allocate( ResidName(N) )
   allocate( ResidNum (N) )

   NAtom = 0
   NResd = 0

   do i = 1, NumSpec

inn: do j = 1, NumMolParam
       if(MolName(i)==MolNameParam(j)) then
         molparnum(i) = j
         exit inn
       end if
       if(j==NumMolParam) then
         write(*,'(a,a)') 'error : no parameter for ',MolName(i)
         write(*,'(a)') 'topology file provides; '
         do k = 1, NumMolParam
           write(*,'(a)') MolNameParam(k)
         end do
         call Finalize
       end if
     end do inn

     if(NumAtm(i)/=NumAtom_inMol(molparnum(i))) then
       if(QMaster) then
         write(*,*) 'error : consistency of NumAtm (configuration vs. topology data)'
         write(*,*) 'configuration =',NumAtm(i)
         write(*,*) 'topology file =',NumAtom_inMol(molparnum(i))
         call Finalize
       end if
     end if

     do j = 1, NumMol(i)

       do k = 1, NumAtm(i)

         NAtom = NAtom + 1

         AtomName(NAtom) = AtomNameParam(k,molparnum(i))

         ii = 0
inn1:    do l = 1, NumResid_inMol(molparnum(i))
           ii = ii + NumAtom_inResidParam(l,molparnum(i))
           if(ii>=k) then
             jj = l
             exit inn1
           end if
         end do inn1

         ResidName(NAtom) = ResidNameParam(jj,molparnum(i))

         ResidNum(NAtom) = NResd + jj

       end do

       NResd = NResd + jj

     end do

   end do

end subroutine NameAssign


!######################################################################
!######################################################################


subroutine BuildMolecule

use Numbers, only : NumSpec, NumMol, NumAtm
use MolLib

implicit none

integer :: i, Num

   Num = 0
   do i = 1, NumSpec
     if(NumAtm(i)/=1) then
       call GeneSingleMol_Vacuume(i,Num)
     end if
     Num = Num + NumMol(i)*NumAtm(i)
   end do

end subroutine BuildMolecule


!######################################################################
!######################################################################


subroutine GeneSingleMol_Vacuume(ispec,Num)

use Numbers, only : NumMol, NumAtm
use CommonBlocks, only : QMaster
use MolLib
use CGParameters
use CGdata
use UnitExParam, only : pi
use BondedParam, only : NumBond, NumAngle, BondI, BondJ, rBond, &
&   AngleI, AngleJ, AngleK, Theta0
use AtomParam, only : AtomName, ResidName, ResidNum, MolName, Mass

implicit none

integer :: i, j, k, l, Num, Ibond, Jbond
integer :: icit, Ipair, Jpair, i1, j1, i2, j2
integer :: ispec, ii, jj, kk, ll, icc, ini, imol
integer :: curr, next, prev, side
integer :: Natm, Nbon, Nang, iang
real(8) :: dist, theta, theta1, theta2
real(8), dimension(3) :: delR
character(len=30) :: Filename
real(8), dimension(:,:), allocatable :: Ratm
real(8), dimension(:), allocatable :: Mss
logical, dimension(:), allocatable :: DoneAtom,DoneBond,DoneAngle
integer, dimension(:), allocatable :: BI,BJ,AI,AJ,AK,Hands,TypeID,BondID,AngleID
integer, dimension(:,:), allocatable :: HandB
real(8), dimension(3) :: Rg
real(8) :: Tmass
logical :: Qfail

! ## getting mol info

   icit = molparnum(ispec)

   ini = 0
   if(ispec /= 1) then
     do i = 1, ispec-1
       ini = ini + NumAtm(i)*NumMol(i)
     end do
   end if

   Natm = NumAtom_inMol(icit)
   Nbon = NumBond_inMol(icit)

   allocate( Ratm(3,Natm) )
   allocate( Mss(Natm) )
   allocate( DoneAtom(Natm) )
   allocate( BI(Nbon) )
   allocate( BJ(Nbon) )
   allocate( BondID(Nbon) )
   allocate( DoneBond(Nbon) )
   allocate( Hands(Natm) )
   allocate( HandB(Natm,4) )
   allocate( TypeID(Natm) )

! ## mass
   Tmass = 0.d0
   do i = 1, Natm
     Mss(i) = Mass(i+ini)
     Tmass = Tmass + Mss(i)
   end do

   do i = 1, Natm
     TypeID(i) = NBAtomType(Num+i)
   end do

   do i = 1, Nbon
     BI(i) = BondPair_inMol(1,i,icit)
     BJ(i) = BondPair_inMol(2,i,icit)
     Ibond = BondPair_inMol(1,i,icit) + Num
     Jbond = BondPair_inMol(2,i,icit) + Num
L01: do j = 1, NumBond
       if((Ibond==BondI(j).and.Jbond==BondJ(j)).or. &
       &  (Ibond==BondJ(j).and.Jbond==BondI(j))) then
         BondID(i) = j
         exit L01
       end if
       if(j==NumBond) then
         if(QMaster) write(*,*) 'error : no bond in (GeneSingleMol_Vacuume)'
         call Finalize
       end if
     end do L01
   end do

   Hands = 0
   HandB = 0

   do k = 1, Nbon
     i = BI(k)
     j = BJ(k)
     Hands(i) = Hands(i) + 1
     Hands(j) = Hands(j) + 1
     HandB(i,Hands(i)) = j
     HandB(j,Hands(j)) = i
   end do

   Nang = 0

   do i = 1, Natm
     k = 0
     do j = 1, Hands(i) - 1
       k = k + 1
     end do
     Nang = Nang + k
   end do

   if(Nang/=0) then

     allocate(AI(Nang))
     allocate(AJ(Nang))
     allocate(AK(Nang))
     allocate(AngleID(Nang))
     allocate(DoneAngle(Nang))

     iang = 0

     do i = 1, Natm
       if(Hands(i)<2) cycle

       do k = 1, Hands(i) - 1
         do l = k+1, Hands(i)

           Ipair = HandB(i,k)
           Jpair = HandB(i,l)
           i1 = BI(Ipair)
           j1 = BJ(Ipair)
           i2 = BI(Jpair)
           j2 = BJ(Jpair)

           iang = iang + 1

           AI(iang) = HandB(i,k)
           AJ(iang) = i
           AK(iang) = HandB(i,l)

         end do
       end do

     end do

   end if

   do l = 1, Nang
     i = AI(l) + Num
     j = AJ(l) + Num
     k = AK(l) + Num
L02: do ll = 1, NumAngle
       ii = AngleI(ll)
       jj = AngleJ(ll)
       kk = AngleK(ll)
       if(j == jj) then
         if(((i==ii).and.(k==kk)).or.((i==kk).and.(k==ii))) then
           AngleID(l) = ll
           exit L02
         end if
       end if
     end do L02
   end do

   icc = 0
   do i = 1, Natm
     if(Hands(i)>=3) then
       icc = icc + 1
       ini = i
     end if
   end do
   if(icc > 1) then
     write(*,*) 'error : the molecule is too complex'
     call Finalize
   end if

   if(icc==0) then
     do i = 1, Natm
       if(Hands(i)==1) then
         ini = i
         exit
       end if
     end do
   end if

!   print *, 'ini=',ini
   print *, 'Start generating molecular structures'

   do imol = 1, NMolLib

   print *, 'Lib No =', imol

   if(imol<10) then
     write(Filename,'(3a,i1)') './MolLib/',trim(MolName(ispec)),'.000',imol
   else if(imol<100) then
     write(Filename,'(3a,i2)') './MolLib/',trim(MolName(ispec)),'.00',imol
   else if(imol<1000) then
     write(Filename,'(3a,i3)') './MolLib/',trim(MolName(ispec)),'.0',imol
   else if(imol<10000) then
     write(Filename,'(3a,i4)') './MolLib/',trim(MolName(ispec)),'.',imol
   else
     write(*,*) 'error : the number of mol. library cannot exceed 10000'
     call Finalize
   end if

!   print *, trim(FileName)

   open(7,file=trim(FileName))

11 Qfail = .False.
   DoneAtom(:) = .False.
   DoneBond(:) = .False.
   if(Nang /= 0) DoneAngle(:) = .False.

   Ratm(:,ini) = 0.d0
   DoneAtom(ini) = .True.

   curr = ini
   next = HandB(ini,1)
   call GetBondLength(dist)
   call RandomSph(delR,dist)
   Ratm(:,next) = Ratm(:,curr)+delR(:)
   DoneAtom(next) = .True.

   do

     prev = curr
     curr = next

     if(Hands(curr)==1) exit

     if(DoneAtom(HandB(curr,1))) next = HandB(curr,2)
     if(DoneAtom(HandB(curr,2))) next = HandB(curr,1)

     call GetAngleT(theta)
     call GetBondLength(dist)
     call GetNewPos(theta,dist)
     if(Qfail) go to 11

   end do

   if(icc/=0) then

     prev = HandB(ini,1)
     curr = ini
     next = HandB(ini,2)

     call GetAngleT(theta)
     call GetBondLength(dist)
     call GetNewPos(theta,dist)
     if(Qfail) go to 11

     do

       prev = curr
       curr = next

       if(Hands(curr)==1) exit

       if(DoneAtom(HandB(curr,1))) next = HandB(curr,2)
       if(DoneAtom(HandB(curr,2))) next = HandB(curr,1)

       call GetAngleT(theta)
       call GetBondLength(dist)
       call GetNewPos(theta,dist)
       if(Qfail) go to 11

     end do

     prev = HandB(ini,1)
     curr = ini
     next = HandB(ini,3)
     call GetAngleT(theta1)
     side = prev
     prev = HandB(ini,2)
     call GetAngleT(theta2)
     call GetBondLength(dist)
     call GetNewPos2(theta1,theta2,dist)
     if(Qfail) go to 11

     do

       prev = curr
       curr = next

       if(Hands(curr)==1) exit

       if(DoneAtom(HandB(curr,1))) next = HandB(curr,2)
       if(DoneAtom(HandB(curr,2))) next = HandB(curr,1)

       call GetAngleT(theta)
       call GetBondLength(dist)
       call GetNewPos(theta,dist)
       if(Qfail) go to 11

     end do

   end if

   do i = 1, Natm
     if(.not.DoneAtom(i)) then
       write(*,*) 'error : atom',i
       call Finalize
     end if
   end do

   do i = 1, Nbon
     if(.not.DoneBond(i)) then
       write(*,*) 'error : bond', i
       call Finalize
     end if
   end do

   do i = 1, Nang
     if(.not.DoneAngle(i)) then
       write(*,*) 'error : angle', i
       call Finalize
     end if
   end do

   Rg = 0.d0
   do i = 1, Natm
     Rg(:) = Rg(:) + Mss(i)*Ratm(:,i)
   end do
   Rg(:) = Rg(:) / Tmass
   do i = 1, Natm
     Ratm(:,i) = Ratm(:,i) - Rg(:)
   end do

   call Write_Lib

   close(7)

   end do

Contains

   subroutine GetBondLength(dd)

   real(8) :: dd

   do k = 1, Nbon
     if(DoneBond(k)) cycle
     if(((curr==BI(k)).and.(next==BJ(k))).or.&
     &  ((curr==BJ(k)).and.(next==BI(k)))) then
       DoneBond(k) = .True.
       dd = rBond(BondID(k))
     end if
   end do

   end subroutine GetBondLength

   subroutine GetAngleT(thet)

   real(8) :: thet

   do k = 1, Nang
     if(DoneAngle(k)) cycle
     if(curr==AJ(k)) then
       if(((prev==AI(k)).and.(next==AK(k))).or.&
       &  ((prev==AK(k)).and.(next==AI(k)))) then
         DoneAngle(k) = .True.
         thet = Theta0(AngleID(k))
       end if
     end if
   end do

   end subroutine GetAngleT

   subroutine GetNewPos(thet,dis)

   real(8) :: thet, dis
   real(8), dimension(3) :: Rnew, Rpre
   real(8), dimension(3,3) :: Rot
   real(8) :: R1, R2, cst, snt, psi, csp, snp
   real(8) :: ranf
   integer :: ll
   integer, parameter :: maxloop = 1000
   real(8), parameter :: rcrit2 = 9.d0
   logical :: Qdone
   external ranf

   Qdone = .False.

   ll = 0

   do while(.not.Qdone)

     ll = ll + 1
     if(ll>maxloop) then
!       write(*,*) 'error in generating a new site; too many trials'
!       write(*,'(a,i3,3f10.3)') 'previous = ',prev,Ratm(:,prev)
!       write(*,'(a,i3,3f10.3)') 'current  = ',curr,Ratm(:,curr)
!       do i = 1, Natm
!         if(DoneAtom(i)) then
!           write(*,'(a,3f12.5)') AtomName(i+Num),Ratm(:,i)
!         end if
!       end do
       Qfail = .True.
       exit
     end if

     cst = cos(thet)
     snt = sin(thet)
     psi = pi*(ranf()*2.d0-1.d0)
     csp = cos(psi)
     snp = sin(psi)
     Rnew(1) = - cst * dis
     Rnew(2) = csp * snt * dis
     Rnew(3) = snp * snt * dis

     Rpre(:) = Ratm(:,curr) - Ratm(:,prev)
     R1      = sqrt(dot_product(Rpre,Rpre))
     Rpre(:) = Rpre(:) / R1

     cst = Rpre(3)
     snt = sqrt(1.d0 - cst*cst)
     csp = Rpre(1) / snt
     snp = Rpre(2) / snt
     Rot(1,1) =  csp*snt
     Rot(1,2) = -snp
     Rot(1,3) = -csp*cst
     Rot(2,1) =  snp*snt
     Rot(2,2) =  csp
     Rot(2,3) = -snp*cst
     Rot(3,1) =  cst
     Rot(3,2) =  0.d0
     Rot(3,3) =  snt

     Rnew(1) = Rnew(1) + R1

     delR = matmul( Rot, Rnew )

     Ratm(:,next) = Ratm(:,prev) + delR(:)
     DoneAtom(next) = .True.

     Qdone = .True.

     do i = 1, Natm
       if(DoneAtom(i).and.i/=prev.and.i/=curr.and.i/=next) then
         delR = Ratm(:,i) - Ratm(:,next)
         R2 = dot_product( delR, delR )
         if(R2 < rcrit2) Qdone = .False.
       end if
     end do

   end do

   end subroutine GetNewPos


   subroutine GetNewPos2(thet1,thet2,dis)

   real(8) :: thet1, thet2, dis
   real(8), dimension(3) :: Rnew, Rpre, Rsid
   real(8), dimension(3,3) :: Rot, Rot_trans
   real(8) :: R1, R2, cst1, snt1, cst2, snt2, cst, snt, psi, csp, snp
   integer :: ll
   integer, parameter :: maxloop = 1000
   real(8), parameter :: rcrit2 = 9.d0
   real(8) :: fmin, dmin
   real(8) :: ranf, a1, b1, c1, fp, ff, xxx
   real(8), dimension(2) :: kai
   logical :: Qdone
   external ranf

   Qdone = .False.
   ll = 0

   do while(.not.Qdone)

     ll = ll + 1
     if(ll>maxloop) then
!       write(*,*) 'error in generating a new site; too many trials'
!       write(*,'(a,i3,3f10.3)') 'previous = ',prev,Ratm(:,prev)
!       write(*,'(a,i3,3f10.3)') 'current  = ',curr,Ratm(:,curr)
!       do i = 1, Natm
!         if(DoneAtom(i)) then
!           write(*,'(a,3f12.5)') AtomName(i+Num),Ratm(:,i)
!         end if
!       end do
       Qfail = .True.
       exit
     end if

     Rpre(:) = Ratm(:,curr) - Ratm(:,prev)
     R1      = sqrt(dot_product(Rpre,Rpre))
     Rpre(:) = Rpre(:) / R1

     cst = Rpre(3)
     snt = sqrt(1.d0 - cst*cst)
     csp = Rpre(1) / snt
     snp = Rpre(2) / snt
     Rot(1,1) =  csp*snt
     Rot(1,2) = -snp
     Rot(1,3) = -csp*cst
     Rot(2,1) =  snp*snt
     Rot(2,2) =  csp
     Rot(2,3) = -snp*cst
     Rot(3,1) =  cst
     Rot(3,2) =  0.d0
     Rot(3,3) =  snt

     Rot_trans = transpose( Rot )

     Rsid(:) = Ratm(:,side) - Ratm(:,curr)
     R1      = sqrt(dot_product(Rsid,Rsid))
     Rsid(:) = Rsid(:) / R1

     cst1 = cos(thet1)
     snt1 = sin(thet1)
     cst2 = cos(thet2)
     snt2 = sin(thet2)

     a1  = cst2 - cst1*Rsid(1)
     b1  = snt1*Rsid(2)
     c1  = snt1*Rsid(3)

     fp = a1 - b1*cos(-pi) -c1*sin(-pi)

     ii = 0

     fmin = fp
     dmin = -pi

     do i = -175, 180, 5

      psi = i/180.d0*pi
      ff  = a1 - b1*cos(psi) -c1*sin(psi)
      if(ff<fmin) then
        fmin = ff
        dmin = psi
      end if

      xxx = ff * fp

      if(xxx < 0.) then
        ii = ii + 1
        kai(ii) = dble(i)-2.5d0
      end if

      fp  = ff

     end do

     if(ii==1) then
       psi = kai(1)
     else if(ii==2) then
       if(ranf()<0.5) then
         psi = kai(1)
       else
         psi = kai(2)
       end if
     else if(ii==0) then
       psi = dmin
     else
       write(*,*) 'error:psi'
       call Finalize
     end if

      psi = psi/180.d0*pi

     csp = cos(psi)
     snp = sin(psi)

     Rnew(1) = - cst1 * dis
     Rnew(2) =   csp  * snt1 * dis
     Rnew(3) =   snp  * snt1 * dis

     Rnew(1) = Rnew(1) + R1

     delR = matmul( Rot, Rnew )

     Ratm(:,next) = Ratm(:,prev) + delR(:)
     DoneAtom(next) = .True.

     Qdone = .True.

     do i = 1, Natm
       if(DoneAtom(i).and.i/=prev.and.i/=curr.and.i/=next.and.i/=side) then
         delR = Ratm(:,i) - Ratm(:,next)
         R2 = dot_product( delR, delR )
         if(R2 < rcrit2) Qdone = .False.
       end if
     end do

   end do

   end subroutine GetNewPos2

   subroutine RandomSph(xx,dis)

   real(8) :: pa, pb, pc, tm, dis
   real(8) :: ranf
   real(8), dimension(3) :: xx
   external ranf

   pc = 1.1d0

   do while(pc > 1.0)
     pa = 2.d0 * ranf() - 1.d0
     pb = 2.d0 * ranf() - 1.d0
     pc = pa*pa + pb*pb
   end do

   tm = 2.d0*sqrt(1.d0-pc)

   xx(1) = (pa*tm)*dis
   xx(2) = (pb*tm)*dis
   xx(3) = (1.d0-2.d0*pc)*dis

   end subroutine RandomSph

   subroutine Write_Lib

     write(7,'(a)') '* CRD file'
     write(7,'(a)') '*'
     write(7,*) Natm

     do i = 1, Natm

       write(7,'(2i5,2(x,a4),3f10.5,x,a4)')      &
       & i , ResidNum(i+Num), ResidName(i+Num),  &
       & AtomName(i+Num), Ratm(:,i), MolName(ispec)(1:4)

     end do

   end subroutine Write_Lib

end subroutine GeneSingleMol_Vacuume


!#####################################################################
!#####################################################################


subroutine Read_Library

use Numbers, only : NumSpec, NumAtm
use MolLib
use AtomParam, only : MolName

implicit none

integer :: MaxAtm, i, imol
character(len=20) :: Filename

   MaxAtm = 1
   do i = NumStSpec+1, NumSpec
     if(NumAtm(i)>MaxAtm) MaxAtm = NumAtm(i)
   end do

   allocate( Rlib(3,MaxAtm,NMolLib,NumSpec) )

   do i = NumStSpec+1, NumSpec

     if(NumAtm(i)/=1) then

     do imol = 1, NMolLib

       if(imol<10) then
         write(Filename,'(3a,i1)') './MolLib/',trim(MolName(i)),'.000',imol
       else if(imol<100) then
         write(Filename,'(3a,i2)') './MolLib/',trim(MolName(i)),'.00',imol
       else if(imol<1000) then
         write(Filename,'(3a,i3)') './MolLib/',trim(MolName(i)),'.0',imol
       else if(imol<10000) then
         write(Filename,'(3a,i4)') './MolLib/',trim(MolName(i)),'.',imol
       else
         write(*,*) 'error : the number of mol. library cannot exceed 10000'
         call Finalize
       end if

       open(7,file=trim(FileName))

       call ReadLib(i,imol)

       close(7)

     end do

     else

       Rlib(:,:,:,i) = 0.d0

     end if

   end do

end subroutine Read_Library


!#####################################################################
!#####################################################################


subroutine ReadLib(ispec,llib)

use MolLib
use Numbers, only : NumAtm

implicit none

integer :: Natm, i, ii, kk, ispec, llib
character(len=4) :: Dum, Dum1, AName

   read(7,'(/)')
   read(7,*) Natm

   if(Natm/=NumAtm(ispec)) then
     write(*,*) ' errors in mol library : Number of atoms '
     call Finalize
   end if

   do i = 1, Natm

     read(7,'(2i5,2(x,a4),3f10.5,x,a4)')      &
     & ii , kk, Dum, AName, Rlib(:,i,llib,ispec), Dum1

   end do

end subroutine ReadLib


!#####################################################################
!#####################################################################


subroutine MakeBox(Rho)

use Numbers, only : N, NumSpec, NumMol, NumAtm
use MolLib
use CellParam, only : CellL
use Configuration, only : R
use AtomParam, only : Mass

implicit none

integer :: Num, ispec, i, NumPres
logical, dimension(:), allocatable :: Qassign
integer :: Nix, Niy, Niz
real(8) :: CellLx, CellLy, CellLz
real(8) :: InvCLx, InvCLy, InvCLz
real(8) :: Rho, V, totalm
#ifdef CHA
integer, dimension(:), allocatable :: Occupied
integer :: Ibx, Iby, Ibz, Ibxy
integer :: Nbox, Ubox
integer :: ILx, ILy, ILz
integer :: icl
real(8) :: CLhx, CLhy, CLhz
real(8) :: DivCLx, DivCLy, DivCLz
real(8) :: Rix, Riy, Riz
#endif

   if(Rho == 0.) then
     CellLx = CellL(1)
     CellLy = CellL(2)
     CellLz = CellL(3)
   else
     Rho = Rho * 1.d-27 ![kg/A^3]
     totalm = 0.d0
     do i = 1, N
       totalm = totalm + Mass(i)
     end do
     V = totalm/Rho
     CellLx = V**(1.d0/3.d0)
     CellLy = CellLx
     CellLz = CellLx
     CellL(1) = CellLx
     CellL(2) = CellLy
     CellL(3) = CellLz
     write(*,'(a,3f10.2)') ' Cell size = ', CellLx,CellLy,CellLz
   end if
   InvCLx = 1.d0/CellLx
   InvCLy = 1.d0/CellLy
   InvCLz = 1.d0/CellLz

#ifdef CHA
   Ibx = int(0.3333333333333333d0 * CellLx)
   Iby = int(0.3333333333333333d0 * CellLy)
   Ibz = int(0.3333333333333333d0 * CellLz)
   DivCLx = CellLx / Ibx
   DivCLy = CellLy / Iby
   DivCLz = CellLz / Ibz

   CLhx = 0.5d0 * CellLx
   CLhy = 0.5d0 * CellLy
   CLhz = 0.5d0 * CellLz

   Nbox = Ibx*Iby*Ibz
   Ibxy = Ibx*Iby

   allocate( Occupied(Nbox) )
   Occupied = 0

   Ubox = Nbox
#endif

   allocate( Qassign(N) )
   Qassign(:) = .False.

   Num = 0
   NumPres = 0

   if(QPartSt) then

     do ispec = 1, NumStSpec
       Num = Num + NumAtm(ispec)*NumMol(ispec)
     end do

     NumPres = Num

     if(fileflag==1) then
       call Read_PDB(Num)
     else if(fileflag==2) then
       call Read_CRD(Num)
     end if

     do i = 1, Num
       Qassign(i) = .True.
#ifdef CHA
       Rix = R(1,i) + CLhx
       Riy = R(2,i) + CLhy
       Riz = R(3,i) + CLhz
       if(Rix<0.d0) then
         Rix = Rix + CellLz
       else if(Rix>CellLz) then
         Rix = Rix - CellLz
       end if
       if(Riy<0.d0) then
         Riy = Riy + CellLy
       else if(Riy>CellLy) then
         Riy = Riy - CellLy
       end if
       if(Riz<0.d0) then
         Riz = Riz + CellLz
       else if(Riz>CellLz) then
         Riz = Riz - CellLz
       end if
       Rix = Rix * InvCLx * Ibx
       Riy = Riy * InvCLy * Iby
       Riz = Riz * InvCLz * Ibz
       ILx = int(Rix) + 1
       ILy = int(Riy) + 1
       ILz = int(Riz) + 1
       icl = ILx + ILy*Ibx + ILz*Ibxy
       if(Occupied(icl)/=1) Ubox = Ubox - 1
       Occupied(icl) = 1
#endif
     end do

   end if

   do ispec = NumStSpec+1, NumSpec

     if(NumAtm(ispec)/=1) then

       call GetInsert

     end if

     Num = Num + NumAtm(ispec)*NumMol(ispec)

   end do

   Num = NumPres

   do ispec = NumStSpec+1, NumSpec

     if(NumAtm(ispec)==1) then

       call GetInsert

     end if

     Num = Num + NumAtm(ispec)*NumMol(ispec)

   end do

   deallocate( Qassign )
#ifdef CHA
   deallocate( Occupied )
#endif

Contains

#ifndef CHA
   subroutine GetInsert

   integer, parameter :: maxtry = 100000
   integer :: Nmol, Natm
   real(8) :: Rcx, Rcy, Rcz, Rijx, Rijy, Rijz
   real(8) :: Rtx, Rty, Rtz
   integer :: j, k, imol, itry, ii
   logical :: Qmolass, Qprof
   real(8), dimension(3,N) :: Rtry
   real(8), parameter :: rcrit2 = 9.d0
   real(8) :: R2, ranf
   external ranf

     Nmol = NumMol(ispec)
     Natm = NumAtm(ispec)

     do i = 1, Nmol

       Qmolass = .False.
       itry = 0

       do while(.not.Qmolass)
         itry = itry + 1

         if(itry>maxtry) then
           write(*,*) 'error : insertion trial exceeds threshold'
           write(*,*) 'reconsider the acceptance distance or density'
           call Finalize
         end if

         imol = int(ranf() * NMolLib) + 1
         if(Qnopart) then
           Qprof=.False.
           if(idir==1) then
             do while(.not.Qprof)
               Rcx = CellLx * (ranf() - 0.5d0)
               if((Rcx<Xmin).or.(Rcx>Xmax)) Qprof=.True.
             end do
             Rcy = CellLy * (ranf() - 0.5d0)
             Rcz = CellLz * (ranf() - 0.5d0)
           else if(idir==2) then
             do while(.not.Qprof)
               Rcy = CellLy * (ranf() - 0.5d0)
               if((Rcy<Xmin).or.(Rcy>Xmax)) Qprof=.True.
             end do
             Rcx = CellLx * (ranf() - 0.5d0)
             Rcz = CellLz * (ranf() - 0.5d0)
           else if(idir==3) then
             do while(.not.Qprof)
               Rcz = CellLz * (ranf() - 0.5d0)
               if((Rcz<Xmin).or.(Rcz>Xmax)) Qprof=.True.
             end do
             Rcx = CellLx * (ranf() - 0.5d0)
             Rcy = CellLy * (ranf() - 0.5d0)
           end if
         else
           Rcx = CellLx * (ranf() - 0.5d0)
           Rcy = CellLy * (ranf() - 0.5d0)
           Rcz = CellLz * (ranf() - 0.5d0)
         end if

         Qmolass = .True.

inn0:    do j = 1, Natm

           Rtx = Rcx + Rlib(1,j,imol,ispec)
           Rty = Rcy + Rlib(2,j,imol,ispec)
           Rtz = Rcz + Rlib(3,j,imol,ispec)

           do k = 1, N

             if(Qassign(k)) then
               Rijx = R(1,k) - Rtx
               Rijy = R(2,k) - Rty
               Rijz = R(3,k) - Rtz
               Nix  = nint(Rijx * InvCLx)
               Niy  = nint(Rijy * InvCLy)
               Niz  = nint(Rijz * InvCLz)
               Rijx = Rijx - Nix * CellLx
               Rijy = Rijy - Niy * CellLy
               Rijz = Rijz - Niz * CellLz
               R2   = Rijx*Rijx + Rijy*Rijy + Rijz*Rijz

               if(R2 < rcrit2) then
                 Qmolass = .False.
                 exit inn0
               end if

             end if

           end do

           Rtry(1,j) = Rtx
           Rtry(2,j) = Rty
           Rtry(3,j) = Rtz

         end do inn0

       end do

       do j = 1, Natm
         ii = Num + (i-1)*Natm + j
         R(:,ii) = Rtry(:,j)
         Qassign(ii) = .True.
       end do

     end do

   end subroutine GetInsert
#else
   subroutine GetInsert

   integer, parameter :: maxtry = 100000
   real(8) :: Rcx, Rcy, Rcz
   real(8) :: Rtx, Rty, Rtz
   real(8) :: Rijx, Rijy, Rijz
   real(8) :: Orix, Oriy, Oriz
   integer :: j, k, imol, itry, ii, count, Ixy
   integer :: Natm, Nmol, TryBox
   logical :: Qmolass
   real(8), dimension(3,N) :: Rtry
   real(8), parameter :: rcrit2 = 9.d0
   real(8) :: R2, ranf
   external ranf

     Nmol = NumMol(ispec)
     Natm = NumAtm(ispec)

     do i = 1, Nmol

       Qmolass = .False.
       itry = 0

       do while(.not.Qmolass)
         itry = itry + 1

         if(itry>maxtry) then
           write(*,*) 'error : insertion trial exceeds threshold'
           write(*,*) 'reconsider the acceptance distance or density'
           call Finalize
         end if

         TryBox = int(ranf()*Ubox) + 1

         count = 0
         do j = 1 , Nbox
           if(Occupied(j)==0) then
             count = count + 1
             if(count==TryBox) then
               TryBox = j
               exit
             end if
             if(count==Ubox) then
               write(*,*) 'error : bug in the selection of box'
               call Finalize
             end if
           end if
         end do

         ILz = (TryBox-1) / Ibxy + 1
         Ixy   = mod(TryBox-1,Ibxy)
         ILy = Ixy / Ibx + 1
         ILx = mod(Ixy,Ibx) + 1

         Orix = (ILx - 1.d0) * DivCLx - CLhx
         Oriy = (ILy - 1.d0) * DivCLy - CLhy
         Oriz = (ILz - 1.d0) * DivCLz - CLhz

         imol = int(ranf() * NMolLib) + 1
         Rcx = DivCLx * ranf() + Orix
         Rcy = DivCLy * ranf() + Oriy
         Rcz = DivCLz * ranf() + Oriz
         Qmolass = .True.

inn1:    do j = 1, Natm

           Rtx = Rcx + Rlib(1,j,imol,ispec)
           Rty = Rcy + Rlib(2,j,imol,ispec)
           Rtz = Rcz + Rlib(3,j,imol,ispec)

           do k = 1, N

             if(Qassign(k)) then

               Rijx = R(1,k) - Rtx
               Rijy = R(2,k) - Rty
               Rijz = R(3,k) - Rtz
               Nix  = nint(Rijx * InvCLx)
               Niy  = nint(Rijy * InvCLy)
               Niz  = nint(Rijz * InvCLz)
               Rijx = Rijx - Nix * CellLx
               Rijy = Rijy - Niy * CellLy
               Rijz = Rijz - Niz * CellLz
               R2   = Rijx*Rijx + Rijy*Rijy + Rijz*Rijz

               if(R2 < rcrit2) then
                 Qmolass = .False.
                 exit inn1
               end if

             end if

           end do

           Rtry(1,j) = Rtx
           Rtry(2,j) = Rty
           Rtry(3,j) = Rtz

         end do inn1

       end do

       do j = 1, Natm

         ii = Num + (i-1) * Natm + j
         R(:,ii) = Rtry(:,j)
         Qassign(ii) = .True.

         Rtx = Rtry(1,j) + CLhx
         Rty = Rtry(2,j) + CLhy
         Rtz = Rtry(3,j) + CLhz
         if(Rtx<0.d0) then
           Rtx = Rtx + CellLz
         else if(Rtx>CellLz) then
           Rtx = Rtx - CellLz
         end if
         if(Rty<0.d0) then
           Rty = Rty + CellLy
         else if(Rty>CellLy) then
           Rty = Rty - CellLy
         end if
         if(Rtz<0.d0) then
           Rtz = Rtz + CellLz
         else if(Rtz>CellLz) then
           Rtz = Rtz - CellLz
         end if
         Rtx = Rtx * InvCLx * Ibx
         Rty = Rty * InvCLy * Iby
         Rtz = Rtz * InvCLz * Ibz
         ILx = int(Rtx) + 1
         ILy = int(Rty) + 1
         ILz = int(Rtz) + 1
         icl = ILx + ILy*Ibx + ILz*Ibxy
         if(Occupied(icl)/=1) Ubox = Ubox - 1
         Occupied(icl) = 1

       end do

     end do

   end subroutine GetInsert
#endif

end subroutine MakeBox


!######################################################################
!######################################################################


! ***************************
! **  Write Configration   **
! ***************************

subroutine Write_CRD

use Numbers, only : N
use CommonBlocks, only : QPBC
use Configuration, only : R
use CellParam, only : H, CellL
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer :: i
character(len=4) :: ModelName

     ModelName = ResidName(1)

     open(7,file='initial.crd')

     write(7,'(a)') '* CRD file'
     write(7,'(a)') '*'
     write(7,'(i10)') N

     do i = 1, N

       write(7,'(2i5,2(x,a4),3f10.5,x,a4)')  &
       & i  , ResidNum(i) , ResidName(i) ,  &
       & AtomName(i) , R(:,i) , ModelName

     end do

     if(QPBC) then

       H = 0.d0
       H(1,1) = CellL(1)
       H(2,2) = CellL(2)
       H(3,3) = CellL(3)

       write(7,*) (H(1,i), i = 1, 3)
       write(7,*) (H(2,i), i = 1, 3)
       write(7,*) (H(3,i), i = 1, 3)

     end if

     close(7)

end subroutine Write_CRD


!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine Write_PDB

use Numbers, only : N, NumSpec, NumMol
use CommonBlocks, only : QPBC
use Configuration, only : R
use UnitExParam, only : pi
use CellParam, only : H
use AtomParam, only : MolName, AtomName, ResidName, ResidNum
use TimeParam, only : Timeps

implicit NONE

integer :: i
! character(len=4) :: ModelName
character(len=4) :: NameA, RName
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0

real(8), dimension(3) :: va, vb, vc
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca

open(2,file='final.pdb')

   write(2,'(a)') &
   & 'TITLE     The Generated Structure by GenCGconf'

   write(2,'(a)') &
   & 'REMARK   1 The System is composed with'

   do i = 1 , NumSpec

     write(2,'(a21,i1,a3,a10,i5)') &
     &     'REMARK   1 Component ',i,' : ',MolName(i),NumMol(i)

   end do

   write(2,'(a)')         'REMARK   2 '
   write(2,'(a,f11.3,a)') 'REMARK   2  Time = ',Timeps,' ps'
   write(2,'(a)')         'REMARK   2 '

   if(QPBC) then

     va = H(:,1)
     vb = H(:,2)
     vc = H(:,3)

     LLa = sqrt( dot_product(va,va) )
     LLb = sqrt( dot_product(vb,vb) )
     LLc = sqrt( dot_product(vc,vc) )
     Aab = acos( dot_product(va,vb) / LLa / LLb ) * 180. / pi
     Abc = acos( dot_product(vb,vc) / LLb / LLc ) * 180. / pi
     Aca = acos( dot_product(vc,va) / LLc / LLa ) * 180. / pi

     write(2,'(a,3f9.3,3f7.2)') 'CRYST1',LLa,LLb,LLc,Abc,Aca,Aab

   end if

   do i = 1, N

     NameA = AtomName(i)
     RName = ResidName(i)

     if(ResidName(i)=='HSD') RName='HIS'
     if(NameA(4:4)/=' ') then

       write(2,'(a4,i7,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,AtomName(i),RName,ResidNum(i),  &
       &  R(:,i),zero,zero,NameA(1:1)

     else

       write(2,'(a4,i7,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,AtomName(i),RName,ResidNum(i), &
       & R(:,i),zero,zero,NameA(1:1)

     end if

   end do

   write(2,'(a)') 'END'

close(2)

end subroutine Write_PDB


!#####################################################################
!#####################################################################


subroutine Read_CRD(Num)

use MolLib
use Configuration, only : R
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer :: i, ii, Num, NN
character(len=4) :: TmpAtomName, TmpResidName
integer :: TmpResidNum

     open(7,file=trim(StructureFile),status='old')

     read(7,*)
     read(7,*)
     read(7,*) NN

     if(Num/=NN) then
       write(*,*) 'The number of particles in the file,',trim(StructureFile),&
       &          ', is not consistent with the provided condition'
       stop
     end if

     do i = 1, NN

        read(7,'(2i5,2(x,a4),3f10.5)')  &
       & ii  , TmpResidNum , TmpResidName ,  &
       & TmpAtomName, R(:,i)

       if(TmpResidNum/=ResidNum(i)) then
         write(*,*) 'WARNING: ResidNum is not consistent for atom ',i
         write(*,*) TmpResidNum,' .vs. ',ResidNum(i)
       end if
       if(adjustl(TmpResidName)/=ResidName(i)) then
         write(*,*) 'WARNING: ResidName is not consistent for atom ',i
         write(*,*) TmpResidName,' .vs. ',ResidName(i)
       end if
       if(adjustl(TmpAtomName)/=AtomName(i)) then
         write(*,*) 'WARNING: AtomName is not consistent for atom ',i
         write(*,*) TmpAtomName,' .vs. ',AtomName(i)
       end if

     end do

     close(7)

end subroutine Read_CRD


!######################################################################
!######################################################################


! **************************
! **  Read Configration   **
! **************************

subroutine Read_PDB(Num)

use MolLib
use Configuration, only : R
use AtomParam, only : AtomName, ResidName, ResidNum

implicit none

integer :: i, ii, Num
character(len=80) :: String
character(len=4) :: TmpAtomName, TmpResidName
integer :: TmpResidNum

   open( 7, file=trim(StructureFile),status='old')

   i = 0

   do

     read(7,'(a80)') String

     if((i==Num).or.(String(1:3)=='END')) exit

     if((String(1:6) == 'HETATM') .or. (String(1:4) == 'ATOM')) then

       i = i + 1

       read(String,'(6x,i5,x,a4,x,a4,x,i5,4x,3f8.3)') &
       &    ii, TmpAtomName, TmpResidName, TmpResidNum, R(:,i)

       if(TmpResidNum/=ResidNum(i)) then
         write(*,*) 'WARNING: ResidNum is not consistent for atom ',i
         write(*,*) TmpResidNum,' .vs. ',ResidNum(i)
       end if
       if(adjustl(TmpResidName)/=ResidName(i)) then
         write(*,*) 'WARNING: ResidName is not consistent for atom ',i
         write(*,*) TmpResidName,' .vs. ',ResidName(i)
       end if
       if(adjustl(TmpAtomName)/=AtomName(i)) then
         write(*,*) 'WARNING: AtomName is not consistent for atom ',i
         write(*,*) TmpAtomName,' .vs. ',AtomName(i)
       end if

     end if

   end do

   close(7)

end subroutine Read_PDB
