! ############################
! ## SUBROUTINE LIST 
! ## -- PSFMaker 
! ## -- SetCondition 
! ## -- MakeXPLORPSF 
! ## -- Write_PSF 
! ############################


! ######################################################################
! ######################################################################


Program PSFMakerCG

implicit none

   call SetCondition1 ! local file

   call Read_CG_Parameter
   call Read_CG_Topology

   call MakeXPLORPSF

   call Write_XPLORPSF

end program PSFMakerCG


! ######################################################################
! ######################################################################


subroutine SetCondition1

use Numbers, only : NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster, QRigidBody
use FFParameters
use IOparam, only : PSF_file, Topology_file, Parameter_file
use AtomParam, only : MolName

implicit none

integer :: i, j, ii, jj, Count
character(len=130) :: String, String1
logical, dimension(130) :: empty
integer,dimension(50) :: first, last
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
       print *, 'Put "D" for a default choice, "~/CGparam/par_CG.prm"'
     case(1)
       print *, 'Name of TOPOLOGY file = ?'
       print *, 'Put "D" for a default choice, "~/CGparam/top_CG.prm"'
     case(2)
       print *, 'Name of PSF file = ?'
     case(3)
       print *, 'Number of molecular species = ?'
     case(4)
       print *, 'Molecular names ? Give them in one line.'
     case(5)
       print *, 'Number of Molecules ? Give them in one line.'

     end select

     read( 5,'(a130)') String1

     String = trim(adjustl(String1))

     if((String(1:1) == '!').or.(String(1:1) == '#')) cycle

     if((String(1:5) == '<end>').or.(String(1:5) == '<END>')) exit

     Count = Count + 1

     select case(Count)

     case(1)

       if(String=='D') then
         write(Parameter_file,*) '~/CGparam/par_CG.prm'
       else
         read(String,*) Parameter_file
       end if

     case(2)

       if(String=='D') then
         write(Topology_file,*) '~/CGparam/top_CG.prm'
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

       do i = 1, 130
         if(String(i:i) == ' ') then
           empty(i) = .True.
         else
           empty(i) = .False.
         end if
       end do

       ii = 0
       j  = 1
       first(j) = 1
       do i = 2, 130
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

       do i = 1, 130
         if(String(i:i) == ' ') then
           empty(i) = .True.
         else
           empty(i) = .False.
         end if
       end do

       ii = 0
       j  = 1
       first(j) = 1
       do i = 2, 130
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


end subroutine SetCondition1


! ######################################################################
! ######################################################################


subroutine MakeXPLORPSF

use Numbers, only : N, NumSpec, NumMol, NumAtm
use NonbondParam, only : Charge
use BondedParam, only : NumBond, NumAngle, NumDihedral, NumImproper, &
&   BondI, BondJ
use AtomParam, only : MolName, AtomName, ResidName, ResidNum, Mass
use CGParameters
use CommonBlocks, only : QMacro

implicit none

integer :: i, j, k, l, ii, jj, ibond, Imol, Natom
integer, dimension(NumSpec) :: imolparam

   QMacro = .False.
   do i = 1, NumSpec
     if(MolName(i)=='CGSP') then
       imolparam(i) = 0
       NumAtm(i) = 1
       QMacro = .True.
       cycle
     end if
inn: do j = 1, NumMolParam
       if(MolName(i)==MolNameParam(j)) then
         imolparam(i) = j
         NumAtm(i) = NumAtom_inMol(j)
         exit inn
       end if
       if(j==NumMolParam) then
         write(*,*) 'error : no data for your molname'
         write(*,*) 'mol. name = ',MolName(i)
         call Finalize
       end if
     end do inn
   end do
   ii = 0
   do i = 1, NumSpec
     ii = ii + NumMol(i)*NumAtm(i)
   end do
   N = ii

   allocate( Charge  (N) )
   allocate( Mass    (N) )
   allocate( ResidNum(N) )
   allocate( AtomName(N) )
   allocate( ResidName(N) )

   NAtom = 0
   Imol  = 0
   do i = 1, NumSpec

     if(imolparam(i) == 0) then
       do j = 1, NumMol(i)
         Imol = Imol + 1
         do k = 1, NumAtm(i)
           NAtom = NAtom + 1
           AtomName(NAtom) = 'CGSP'
           ResidName(NAtom) = 'CGSP'
           Charge(NAtom) = 0.d0
           Mass  (NAtom) = 0.d0
           ResidNum(Natom) = Imol
         end do
       end do
     else
       do j = 1, NumMol(i)
         Imol = Imol + 1
         do k = 1, NumAtm(i)
           NAtom = NAtom + 1
           Charge(NAtom) = ChargeParam(k,imolparam(i))
           Mass  (NAtom) = AtomMass   (k,imolparam(i))
           ResidName(NAtom) = MolNameParam(imolparam(i))(1:4)
          AtomName(NAtom) = AtomNameParam(k,imolparam(i))
           ResidNum(Natom) = Imol
         end do
       end do
     end if

   end do

   if(QMacro) call Set_Sphere(1)
! ----------------------------

   NumBond = 0
   do i = 1 , NumSpec
     NumBond = NumBond + NumMol(i) * NumBond_inMol(imolparam(i))
   end do
   allocate( BondI(NumBond) )
   allocate( BondJ(NumBond) )

   ii = 0
   ibond = 0
   do i = 1 , NumSpec
     if(i>=2) ii = ii + NumAtm(i-1) * NumMol(i-1)
     do j = 1 , NumMol(i)
       jj = (j-1)*NumAtm(i) + ii
       do k = 1 , NumBond_inMol(imolparam(i))
         ibond = ibond + 1
         BondI(ibond) = jj+BondPair_inMol(1,k,imolparam(i))
         BondJ(ibond) = jj+BondPair_inMol(2,k,imolparam(i))
       end do
     end do
   end do

   NumAngle = 0
   NumImproper = 0
   NumDihedral = 0

end subroutine MakeXPLORPSF


! ######################################################################
! ######################################################################


subroutine Write_XPLORPSF

use Numbers, only : N
use IOparam, only : PSF_file
use NonbondParam, only : Charge
use BondedParam, only : NumBond, NumAngle, NumDihedral, NumImproper, &
&   BondI, BondJ, AngleI, AngleJ, AngleK, DihedI, DihedJ, DihedK, DihedL, &
&   ImproI, ImproJ, ImproK, ImproL
use AtomParam, only : AtomName, ResidName, ResidNum, Mass

implicit none

integer :: i, j, ii, jj, k
character(len=20) :: username
character(len=17) :: ptime
character(len=4) :: Dummy

   open(1,file=trim(PSF_file),status='unknown')

   write(1,'(a/)') 'PSF '
   write(1,'(a)')  '       2 !NTITLE'
   write(1,'(a)')  &
   & '* SCRIPT FILE PRODUCED BY PSFMaker version 1                                    '

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
     write(1,'(i8,x,a4,x,i4,2(x,a4),x,a1,3x,f11.6,f14.4)') &
     &  i,ResidName(i),ResidNum(i),ResidName(i),AtomName(i),&
     &  AtomName(i)(1:1),Charge(i),Mass(i)
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

   j = 0
   write(1,'(/i8,x,a)') j,'!NTHETA: angles'

   write(1,'(/i8,x,a)') j,'!NPHI: dihedrals'

   write(1,'(/i8,x,a)') j,'!NIMPHI: impropers'

   write(1,'(/i8,x,a)') j,'!NDON: donors'

   write(1,'(/i8,x,a)') j,'!NACC: acceptors'

   close(1)

end subroutine Write_XPLORPSF
