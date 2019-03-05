! ############################
! ## SUBROUTINE LIST 
! ## -- Read_CG_Parameter 
! ## -- Read_CG_Topology 
! ## -- Read_Top_Lines
! ## -- AllocateCGdata
! ## -- Read_TableFunctions
! ############################


!######################################################################
!######################################################################


subroutine Read_CG_Parameter

use CommonBlocks, only : QMaster, QGeneconf, Qdebug, Qstdout
use CGParameters
use IOparam, only : Parameter_file
use UnitExParam, only : pi, reng, ExParam

implicit NONE

integer :: ibond, iangl, itors, iimpr, ivdws, idkey, i, ilen
character(len=80) :: String, String1
character(len=10) :: Keyword
integer :: eofile
real(8) :: Rmin, Rmax
real(8) :: Eps, Sig, coeb, coea, coet, coei, coev
character(len=12) :: cUnitbond, cUnitangle, cUnittors, cUnitimpr, cUnitvdw
character(len=10) :: PotType
logical :: Qunitbond, Qunitangle, Qunittors, Qunitimpr, Qunitvdw


open(1,file=trim(Parameter_file),status='old')

! ## getting the number of parameters

   ibond = 0
   iangl = 0
   itors = 0
   iimpr = 0
   ivdws = 0

   idkey = 0

   do

     read(1,'(a)',iostat=eofile) String1
     if(eofile == -1) exit
     String = trim(adjustl(String1))

     if(String(1:1)=='!'.or.String(1:1)=='#'.or.String(1:1)==' ') cycle
     if(String(1:5) == 'UNIT=') cycle

     if(String(1:2) == '>>') then
       read(String(3:),*) Keyword
       if(trim(adjustl(Keyword))=='BOND')     idkey = 1
       if(trim(adjustl(Keyword))=='ANGLE')    idkey = 2
       if(trim(adjustl(Keyword))=='TORSION')  idkey = 3
       if(trim(adjustl(Keyword))=='IMPROPER') idkey = 4
       if(trim(adjustl(Keyword))=='NONBOND')  idkey = 5
     else if(String(1:2) == '<<') then
       idkey = 0
     else
       select case(idkey)
       case(0)
         cycle
       case(1)
         ibond = ibond + 1
       case(2)
         iangl = iangl + 1
       case(3)
         itors = itors + 1
       case(4)
         iimpr = iimpr + 1
       case(5)
         ivdws = ivdws + 1
       end select
     end if

   end do

   NumBondParam    = ibond
   NumAngleParam   = iangl
   NumDihedParam   = itors
   NumImproParam   = iimpr
   NumNonBondParam = ivdws

   if(Qstdout) then
     print *, 'NumBondParam =',NumBondParam
     print *, 'NumAngleParam =',NumAngleParam
     print *, 'NumDihedParam =',NumDihedParam
     print *, 'NumImproParam =',NumImproParam
     print *, 'NumNonBondParam =',NumNonBondParam
   end if

! ## memory allocations

   allocate( BondPairAtoms(2,ibond) )
   allocate( BondType(ibond) )
   allocate( kBondParam(ibond) )
   allocate( rBondParam(ibond) )

   allocate( AnglePairAtoms(3,iangl) )
   allocate( AngleType(iangl) )
   allocate( kThetaParam(iangl) )
   allocate( Theta0Param(iangl) )

   if(itors/=0) then
     allocate( DihedPairAtoms(4,itors) )
     allocate( DihedType(itors) )
     allocate( kChiParam(itors) )
     allocate( DeltaDihParam(itors) )
     allocate( NdihParam(itors) )
   end if

   if(iimpr/=0) then
     allocate( ImproPairAtoms(4,iimpr) )
     allocate( ImproType(iimpr) )
     allocate( kPsiParam(iimpr) )
     allocate( NimpParam(iimpr) )
     allocate( Psi0Param(iimpr) )
   end if

   allocate( NonBondPairAtoms(2,ivdws) )
   allocate( NonBondType(ivdws) )
   allocate( NonBondFileHeader(ivdws) )
   allocate( AijParam(ivdws) )
   allocate( BijParam(ivdws) )
   allocate( CijParam(ivdws) )
   allocate( R13cutParam(ivdws) )
   allocate( EpsParam(ivdws) )
   allocate( RminParam(ivdws) )
   allocate( RmaxParam(ivdws) )

   REWIND(1)

! ## reading parameters

   ibond = 0
   iangl = 0
   itors = 0
   iimpr = 0
   ivdws = 0

   idkey = 0

   do

     read(1,'(a)',iostat=eofile) String
     if(eofile == -1) exit
     String1 = trim(adjustl(String))

     ilen = len(String1)
cuto:do i = 1, ilen
       if(String1(i:i)=='!') then
         ilen = i-1
         exit cuto
       end if
     end do cuto
     write(String,'(a)') String1(1:ilen)

     if(String(1:1)=='!'.or.String(1:1)=='#'.or.String(1:1)==' ') cycle

     if(String(1:2) == '>>') then
       read(String(3:),*) Keyword
       if(trim(adjustl(Keyword))=='BOND')     idkey = 1
       if(trim(adjustl(Keyword))=='ANGLE')    idkey = 2
       if(trim(adjustl(Keyword))=='TORSION')  idkey = 3
       if(trim(adjustl(Keyword))=='IMPROPER') idkey = 4
       if(trim(adjustl(Keyword))=='NONBOND')  idkey = 5
     else if(String(1:2) == '<<') then
       idkey = 0
     else
       select case(idkey)
       case(0)
         cycle
       case(1)
         if(String(1:5) == 'UNIT=') then
           read(String(6:),*) cUnitbond
           Qunitbond = .True.
         else
           ibond = ibond + 1
           read(String,*) BondPairAtoms(:,ibond), PotType, &
           &              kBondParam(ibond), rBondParam(ibond)
           if(PotType=='harm') then
             BondType(ibond) = 1
           else
             print '(/a/)', 'ERROR : function type for bond',ibond
             call Finalize
           end if
         end if
       case(2)
         if(String(1:5) == 'UNIT=') then
           read(String(6:),*) cUnitangle
           Qunitangle = .True.
         else
           iangl = iangl + 1
           read(String,*) AnglePairAtoms(:,iangl), PotType, &
           &              kThetaParam(iangl), Theta0Param(iangl)
           if(PotType=='cosine') then
             AngleType(iangl) = 1
           else if(PotType=='quartic') then
             AngleType(iangl) = 2
           else if(PotType=='harm') then
             AngleType(iangl) = 3
           else
             print '(/a/)', 'ERROR : function type for angle',iangl
             call Finalize
           end if
         end if
       case(3)
         if(String(1:5) == 'UNIT=') then
           read(String(6:),*) cUnittors
           Qunittors = .True.
         else
           itors = itors + 1
           read(String,*) DihedPairAtoms(:,itors), PotType
           if(PotType=='cos') then
             DihedType(itors) = 1
             read(String,*) DihedPairAtoms(:,itors), PotType, &
             &              kChiParam(itors), NdihParam(itors), DeltaDihParam(itors)
           else if(PotType=='harm') then
             DihedType(itors) = 2
             read(String,*) DihedPairAtoms(:,itors), PotType, &
             &              kChiParam(itors), DeltaDihParam(itors)
           else
             print '(/a/)', 'ERROR : function type for torsion',itors
             call Finalize
           end if
         end if
       case(4)
         if(String(1:5) == 'UNIT=') then
           read(String(6:),*) cUnitimpr
           Qunitimpr = .True.
         else
           iimpr = iimpr + 1
           read(String,*) ImproPairAtoms(:,iimpr), PotType
           if(PotType=='cos') then
             ImproType(iimpr) = 1
             read(String,*) ImproPairAtoms(:,iimpr), PotType, &
             &              kPsiParam(iimpr), NimpParam(iimpr), Psi0Param(iimpr)
           else if(PotType=='harm') then
             ImproType(iimpr) = 2
             read(String,*) ImproPairAtoms(:,iimpr), PotType, &
             &              kPsiParam(iimpr), Psi0Param(iimpr)
           else
             print '(/a/)', 'ERROR : function type for improper tosion',iimpr
             call Finalize
           end if
         end if
       case(5)
         if(String(1:5) == 'UNIT=') then
           read(String(6:),*) cUnitvdw
           Qunitvdw = .True.
         else
           ivdws = ivdws + 1
           read(String,*) NonBondPairAtoms(:,ivdws), PotType
           if(PotType=='LJ9-6') then
             NonBondType(ivdws) = 1
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              Eps, Sig, Rmin, Rmax
             AijParam(ivdws) = 27.d0 / 4.d0 * Eps * (Sig ** 9)
             BijParam(ivdws) = 27.d0 / 4.d0 * Eps * (Sig ** 6)
             R13cutParam(ivdws) = (1.5d0 ** (1.d0/3.d0) ) * Sig
             EpsParam(ivdws) = Eps
           else if(PotType=='LJ6-4') then
             NonBondType(ivdws) = 2
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              Eps, Sig, Rmin, Rmax
             AijParam(ivdws) = 27.d0 / 4.d0 * Eps * (Sig ** 6)
             BijParam(ivdws) = 27.d0 / 4.d0 * Eps * (Sig ** 4)
             R13cutParam(ivdws) = sqrt(1.5d0) * Sig
             EpsParam(ivdws) = Eps
           else if(PotType=='LJ8-4') then
             NonBondType(ivdws) = 5
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              Eps, Sig, Rmin, Rmax
             AijParam(ivdws) = 4.d0 * Eps * (Sig ** 8)
             BijParam(ivdws) = 4.d0 * Eps * (Sig ** 4)
             R13cutParam(ivdws) = (2.d0**(0.25d0)) * Sig
             EpsParam(ivdws) = Eps
           else if(PotType=='LJ10-4') then
             NonBondType(ivdws) = 6
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              Eps, Sig, Rmin, Rmax
             AijParam(ivdws) = 5.d0**(5.d0/3.d0)/(2.d0**(2.d0/3.d0)*5.d0-2.d0**(5.d0/3.d0))*&
             &                 Eps * (Sig **(10) )
             BijParam(ivdws) = 5.d0**(5.d0/3.d0)/(2.d0**(2.d0/3.d0)*5.d0-2.d0**(5.d0/3.d0))*&
             &                 Eps * (Sig ** 4)
             R13cutParam(ivdws) = (2.5d0**(1.d0/6.d0)) * Sig
             EpsParam(ivdws) = Eps
           else if(PotType=='LJ12-4') then
             NonBondType(ivdws) = 7
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              Eps, Sig, Rmin, Rmax
             AijParam(ivdws) = 3.d0*sqrt(3.d0)*0.5d0 * Eps * (Sig **(12) )
             BijParam(ivdws) = 3.d0*sqrt(3.d0)*0.5d0 * Eps * (Sig ** 4)
             R13cutParam(ivdws) = (3.d0**(1.d0/8.d0)) * Sig
             EpsParam(ivdws) = Eps
           else if(PotType=='LJ12-6') then
             NonBondType(ivdws) = 8
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              Eps, Sig, Rmin, Rmax
             AijParam(ivdws) = 4.d0 * Eps * (Sig **(12) )
             BijParam(ivdws) = 4.d0 * Eps * (Sig ** 6)
             R13cutParam(ivdws) = (2.d0**(1.d0/6.d0)) * Sig
             EpsParam(ivdws) = Eps
           else if(PotType=='morse') then ! morse : A=epsilon, B=A(beta), C=r_min
             NonBondType(ivdws) = 4
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              AijParam(ivdws), BijParam(ivdws), CijParam(ivdws), Rmin, Rmax
             R13cutParam(ivdws) = CijParam(ivdws)
             EpsParam(ivdws) = Eps
           else if(PotType=='table') then
             NonBondType(ivdws) = 3
             read(String,*) NonBondPairAtoms(:,ivdws), PotType, &
             &              NonBondFileHeader(ivdws), Rmin, Rmax
             R13cutParam(ivdws) = Rmin + 2.d0
             EpsParam(ivdws) = 0.
           end if
           RminParam(ivdws) = Rmin
           RmaxParam(ivdws) = Rmax
         end if
       end select
     end if

   end do


   do i = 1, NumAngleParam
     if((AngleType(i)==1).and.(Theta0Param(i)/=180.)) then ! cosine
       if(QMaster) then
         write(*,*) 'error : when angle func type 1 is chosen, '
         write(*,*) '        theta_0 should be 180.00 in the current implementation'
       end if
       call Finalize
     end if
   end do

! ## CHECK

   if(Qdebug) then
   do i = 1, NumNonBondParam
     if(NonBondType(i)==3) then
     write(*,'(2a,i3,a,2e15.5)') NonBondPairAtoms(:,i), NonBondType(i), &
     &                         NonBondFileHeader(i), RminParam(i), RmaxParam(i)
     else
     write(*,'(2a,i3,4e15.5)') NonBondPairAtoms(:,i), NonBondType(i), &
     &                         AijParam(i), BijParam(i), RminParam(i), RmaxParam(i)
     end if
   end do
   end if

! ----
! Unit
! ----
   if(Qunitbond) then
     if(trim(cUnitbond)=='kcal_per_mol') then
       coeb = ExParam
     else if(trim(cUnitbond)=='kJ_per_mol') then
       coeb = ExParam/4.184d0
     else if(trim(cUnitbond)=='K') then
       coeb = 1.98719137d-03*ExParam
     else
       print '(/a/)', 'ERROR : unit for bond'
       call Finalize
     end if
   else
     coeb = ExParam
   end if

   if(Qunitangle) then
     if(trim(cUnitangle)=='kcal_per_mol') then
       coea = ExParam
     else if(trim(cUnitangle)=='kJ_per_mol') then
       coea = ExParam/4.184d0
     else if(trim(cUnitangle)=='K') then
       coea = 1.98719137d-03*ExParam
     else
       print '(/a/)', 'ERROR : unit for angle'
       call Finalize
     end if
   else
     coea = ExParam
   end if

   if(Qunittors.and.NumDihedParam/=0) then
     if(trim(cUnittors)=='kcal_per_mol') then
       coet = ExParam
     else if(trim(cUnittors)=='kJ_per_mol') then
       coet = ExParam/4.184d0
     else if(trim(cUnittors)=='K') then
       coet = 1.98719137d-03*ExParam
     else
       print '(/a/)', 'ERROR : unit for torsion'
       call Finalize
     end if
   else
     coet = ExParam
   end if

   if(Qunitimpr.and.NumImproParam/=0) then
     if(trim(cUnitimpr)=='kcal_per_mol') then
       coei = ExParam
     else if(trim(cUnitimpr)=='kJ_per_mol') then
       coei = ExParam/4.184d0
     else if(trim(cUnitimpr)=='K') then
       coei = 1.98719137d-03*ExParam
     else
       print '(/a/)', 'ERROR : unit for improper torsion'
       call Finalize
     end if
   else
     coei = ExParam
   end if

   if(Qunitvdw) then
     if(trim(cUnitvdw)=='kcal_per_mol') then
       coev = ExParam
     else if(trim(cUnitvdw)=='kJ_per_mol') then
       coev = ExParam/4.184d0
     else if(trim(cUnitvdw)=='K') then
       coev = 1.98719137d-03*ExParam
     else
       print '(/a/)', 'ERROR : unit for vdW'
       call Finalize
     end if
   else
     coev = ExParam
   end if


   kBondParam (:) = kBondParam (:) * coeb
   kThetaParam(:) = kThetaParam(:) * coea
   Theta0Param(:) = Theta0Param(:) / 180.d0 * pi
   AijParam(:) = AijParam(:) * coev
   EpsParam(:) = EpsParam(:) * coev

   if(NumDihedParam/=0) then
     kChiParam(:) = kChiParam(:) * coet
     DeltaDihParam(:) = DeltaDihParam(:) / 180.d0 * pi
   end if

   if(NumImproParam/=0) then
     kPsiParam(:) = kPsiParam(:) * coei
     Psi0Param(:) = Psi0Param(:) / 180.d0 * pi
   end if

   do i = 1, NumNonBondParam
     if(NonBondType(i)/=4) then ! not for morse
       BijParam(i) = BijParam(i) * coev
     end if
   end do

   if(.not.QGeneconf) then
     do i = 1, NumAngleParam
       if(AngleType(i)==2) then ! Quartic
         Theta0Param(i) = ( Theta0Param(i) - pi ) ** 2
         kThetaParam(i) = kThetaParam(i) / (2.d0 * Theta0Param(i))
       end if
     end do
   end if

close(1)

end subroutine Read_CG_Parameter


!######################################################################
!######################################################################


! ****************************
! ** Charmm Topology File   **
! ****************************

subroutine Read_CG_Topology

use CommonBlocks, only : QMaster, Qstdout, Qdebug
use CGParameters
use IOparam, only : Topology_file

implicit NONE

integer :: eofile
character(len=80) :: String, String1

integer :: MaxAtom, MaxBond, MaxImpr, MaxUnit, MaxResid
integer :: iatom, ibond, iimpr, iunit, iresid, Nmol
integer :: i, j
integer, dimension(1000) :: icheck
logical :: Qscript
character(len=1), parameter :: Blanc = ' '
integer :: ichar
integer :: id, k

! Topology

open(1,file=trim(Topology_file),status='old')

! ## Counting for memory allocations

   NumMolParam = 0
   MaxAtom = 0
   MaxBond = 0
   MaxImpr = 0
   MaxUnit = 0
   MaxResid = 0

   do

     read(1,'(a)',iostat=eofile) String1
     if(eofile == -1) exit
     String = trim(adjustl(String1))

     if(String(1:1)=='!'.or.String(1:1)=='#') cycle

     if(String(1:2) == '>>') NumMolParam = NumMolParam + 1
     if(String(1:8) == 'NUMATOM=') then
       read(String(9:),*) iatom
       if(iatom>MaxAtom) MaxAtom = iatom
     else if(String(1:8) == 'NUMBOND=') then
       read(String(9:),*) ibond
       if(ibond>MaxBond) MaxBond = ibond
     else if(String(1:8) == 'NUMIMPR=') then
       read(String(9:),*) iimpr
       if(iimpr>MaxImpr) MaxImpr = iimpr
     else if(String(1:9) == 'NUMRIGID=') then
       read(String(10:),*) iunit
       if(iunit>MaxUnit) MaxUnit = iunit
     else if(String(1:9) == 'NUMRESID=') then
       read(String(10:),*) iresid
       if(iresid>MaxResid) MaxResid = iresid
     end if

   end do

   allocate( NumAtom_inMol(NumMolParam) )
   allocate( NumBond_inMol(NumMolParam) )
   allocate( NumIMPR_inMol(NumMolParam) )
   allocate( NumUnit_inMol(NumMolParam) )
   allocate( NumResid_inMol(NumMolParam) )

   allocate( MolNameParam(NumMolParam) )
   allocate( AtomNameParam(MaxAtom,NumMolParam) )
   allocate( AtomType(MaxAtom,NumMolParam) )
   allocate( AtomMass(MaxAtom,NumMolParam) )
   allocate( ChargeParam(MaxAtom,NumMolParam) )
   allocate( ScrnParam(MaxAtom,NumMolParam) )

   if(MaxBond/=0) then
     allocate( BondPair_inMol(2,MaxBond,NumMolParam) )
   end if
   if(MaxImpr/=0) then
     allocate( ImprPair_inMol(4,MaxImpr,NumMolParam) )
   end if

   if(MaxUnit == 0) MaxUnit = 1
   allocate( UNITNameParam(MaxUnit,NumMolParam) )
   allocate( NumAtom_inUNITParam(MaxUnit,NumMolParam) )

   if(MaxResid == 0) MaxResid = 1
   allocate( ResidNameParam(MaxResid,NumMolParam) )
   allocate( NumAtom_inResidParam(MaxResid,NumMolParam) )

! ## 

   if(Qstdout) then

     if(Qdebug) then
      print *, 'MaxAtom=',MaxAtom
      print *, 'MaxBond=',MaxBond
      print *, 'MaxImpr=',MaxImpr
      print *, 'MaxUnit=',MaxUnit
      print *, 'MaxResid=',MaxResid
     end if

     print *, 'NumMolParam = ',NumMolParam
   end if

   if(NumMolParam==0) then
     if(QMaster) write(*,*) 'Error : No topological data is given (in Read_CG_Topology)'
     call Finalize
   end if

! ## Reading details 

   REWIND(1)

   Nmol = 0
   Qscript = .False.

   do

     read(1,'(a)',iostat=eofile) String1
     if(eofile == -1) exit
     String = trim(adjustl(String1))

     if(String(1:1)=='#'.or.String(1:1)=='!'.or.String(1:1)==' ') cycle

     if(String(1:2) == '<<') then
       Qscript = .False.

! ## Checking the data

       if(NumAtom_inMol(Nmol)==0) then
         if(QMaster) write(*,*) 'Reading error : NumAtom_inMol(in Read_CG_Topology)'
         call Finalize
       end if
       if(iatom /= NumAtom_inMol(Nmol)) then
         if(QMaster) write(*,*) 'Error : NumAtom_inMol (in Read_CG_Topology)'
         if(QMaster) write(*,*) 'NumAtom_inMol = ',NumAtom_inMol(Nmol)
         if(QMaster) write(*,*) 'Counted atoms = ',iatom
         call Finalize
       end if
       if(ibond /= NumBond_inMol(Nmol)) then
         if(QMaster) write(*,*) 'Error : NumBond_inMol (in Read_CG_Topology)'
         call Finalize
       end if
       if(iimpr /= NumIMPR_inMol(Nmol)) then
         if(QMaster) write(*,*) 'Error : NumIMPR_inMol (in Read_CG_Topology)'
         call Finalize
       end if
       if(iunit /= NumUnit_inMol(Nmol)) then
         if(QMaster) write(*,*) 'Error : NumUnit_inMol (in Read_CG_Topology)'
         call Finalize
       end if
       if(iresid /= NumResid_inMol(Nmol)) then
         if(QMaster) write(*,*) 'Error : NumResid_inMol (in Read_CG_Topology)'
         call Finalize
       end if

       j = 0
       do i = 1, iatom
         j = j + icheck(i)
       end do
       if(j /= NumAtom_inMol(Nmol)) then
         if(QMaster) write(*,*) 'Error : ATOM numbers (in Read_CG_Topology)'
         call Finalize
       end if

       if(NumUnit_inMol(Nmol)==0) then
         NumUnit_inMol(Nmol) = 1
         UNITNameParam(1,Nmol) = MolNameParam(Nmol)(1:6)
         NumAtom_inUNITParam(1,Nmol) = NumAtom_inMol(Nmol)
       else
         j = 0
         do i = 1, NumUnit_inMol(Nmol)
           j = j + NumAtom_inUNITParam(i,Nmol)
         end do
         if(j/=NumAtom_inMol(Nmol)) then
           if(QMaster) write(*,*) 'Error : UNIT numbers (in Read_CG_Topology)'
           call Finalize
         end if
       end if

       if(NumResid_inMol(Nmol)==0) then
         NumResid_inMol(Nmol) = 1
         ResidNameParam(1,Nmol) = MolNameParam(Nmol)(1:4)
         NumAtom_inResidParam(1,Nmol) = NumAtom_inMol(Nmol)
       else
         j = 0
         do i = 1, NumResid_inMol(Nmol)
           j = j + NumAtom_inResidParam(i,Nmol)
         end do
         if(j/=NumAtom_inMol(Nmol)) then
           if(QMaster) write(*,*) 'Error : RESID numbers (in Read_CG_Topology)'
           call Finalize
         end if
       end if

       if(Qdebug) then
       print *, 'NumAtom_inMol(Nmol)', NumAtom_inMol(Nmol)
       print *, 'NumBond_inMol(Nmol)', NumBond_inMol(Nmol)
       print *, 'NumIMPR_inMol(Nmol)', NumIMPR_inMol(Nmol)
       print *, 'NumUnit_inMol(Nmol)', NumUnit_inMol(Nmol)
       print *, 'NumResid_inMol(Nmol)',NumResid_inMol(Nmol)

       do i = 1, NumAtom_inMol(Nmol)
         print *, 'AtomName and Type  :', AtomNameParam(i,Nmol), AtomType(i,Nmol)
       end do
       do i = 1, NumAtom_inMol(Nmol)
         print *, 'Mass and Charges  :', AtomMass(i,Nmol),ChargeParam(i,Nmol)
       end do
       do i = 1, NumAtom_inMol(Nmol)
         print *, 'Screening Coulomb parameter :', ScrnParam(i,Nmol)
       end do
       end if


       if(Nmol==NumMolParam) exit

     end if

     if(Qscript) then

       if(String(1:8) == 'NUMATOM=') then

         read(String(9:),*) NumAtom_inMol(Nmol)

       else if(String(1:8) == 'NUMBOND=') then

         read(String(9:),*) NumBond_inMol(Nmol)

       else if(String(1:8) == 'NUMIMPR=') then

         read(String(9:),*) NumIMPR_inMol(Nmol)

       else if(String(1:9) == 'NUMRIGID=') then

         read(String(10:),*) NumUnit_inMol(Nmol)

       else if(String(1:9) == 'NUMRESID=') then

         read(String(10:),*) NumResid_inMol(Nmol)

       else if(String(1:4) == 'ATOM') then

         iatom = iatom + 1
         read(String(5:),*) id,AtomNameParam(id,Nmol), &
         &  AtomType(id,Nmol),AtomMass(id,Nmol),ChargeParam(id,Nmol),ScrnParam(id,Nmol)
         icheck(id) = 1

       else if(String(1:4) == 'BOND') then

         ichar = 0
         do j = 6, 80
           if( (String(j  :j  ) == Blanc) .and. &
           &   (String(j-1:j-1) /= Blanc) ) then
             ichar = ichar + 1
           end if
         end do
         if(String(80:80) /= Blanc) ichar = ichar + 1

         if(mod(ichar,2)/=0) then
           if(QMaster) write(*,*) 'error : BOND in topology file'
           call Finalize
         end if

         ichar = ichar / 2
         read(String(5:),*) ((BondPair_inMol(j,k,Nmol),j=1,2),k=ibond+1,ibond+ichar)
         ibond = ibond + ichar

       else if(String(1:4) == 'IMPR') then

         ichar = 0
         do j = 6, 80
           if( (String(j  :j  ) == Blanc) .and. &
           &   (String(j-1:j-1) /= Blanc) ) then
             ichar = ichar + 1
           end if
         end do
         if(String(80:80) /= Blanc) ichar = ichar + 1

         if(mod(ichar,4)/=0) then
           if(QMaster) write(*,*) 'error : IMPR in topology file'
           call Finalize
         end if

         ichar = ichar / 4
         read(String(5:),*) ((ImprPair_inMol(j,k,Nmol),j=1,4),k=iimpr+1,iimpr+ichar)
         iimpr = iimpr + ichar

       else if(String(1:4) == 'UNIT') then

         iunit = iunit + 1
         read(String(5:),*) id, UNITNameParam(id,Nmol),NumAtom_inUNITParam(id,Nmol)

       else if(String(1:5) == 'RESID') then

         iresid = iresid + 1
         read(String(6:),*) id, ResidNameParam(id,Nmol),NumAtom_inResidParam(id,Nmol)

       else

         if(QMaster) then
           write(*,'(a)') 'Error : unknown script in topology file'
           write(*,'(a)') trim(String)
         end if
         call Finalize

       end if

     end if

     if(String(1:2) == '>>') then

       Nmol = Nmol + 1

       read(String(3:),*) MolNameParam(Nmol)
       if(QMaster.and.Qstdout) then
         print *, 'Molecule',Nmol,' = ', MolNameParam(Nmol)
       end if
       Qscript = .True.

! ## for check
       iatom = 0
       ibond = 0
       iimpr = 0
       iunit = 0
       iresid = 0
       icheck = 0

! ##
       NumAtom_inMol(Nmol) = 0
       NumBond_inMol(Nmol) = 0
       NumIMPR_inMol(Nmol) = 0
       NumUnit_inMol(Nmol) = 0
       NumResid_inMol(Nmol) = 0

     end if

   end do

close(1)

end subroutine Read_CG_Topology


!######################################################################
!######################################################################


subroutine AllocateCGdata

use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QMaster, QPBC, QSHAKE, QRigidBody, QCoulomb, &
&   Qstdout, QCGWALL, QMacro, QCyl, QFSCyl
use CGParameters
use CGdata
use EwaldParam, only : Nel, Nelist, PCh
use TableFuncs
use RBparam, only : NumRB, AtomUnitNum, AtomUnitName
use SHAKEparam, only : NSHAKEGroup
use UnitExParam, only : ec, Avogadro, ExParam, kb
use NonbondParam, only : Charge, ScrnCR
use BondedParam, only : NumBond, NumAngle, NumUB, NumDihedral, NumImproper, &
&   BondI, BondJ, FTypeBond, kBond, rBond, AngleI, AngleJ, AngleK, FTypeAngle, &
&   kTheta, Theta0, vdwSubtAng, DihedI, DihedJ, DihedK, DihedL, FTypeDihed, &
&   kChi, DeltaDih, NDih, DupFlag, DupkChi, DupDeltaDih, DupNDih, ImproI, &
&   ImproJ, ImproK, ImproL, FTypeImpro, kImp, NImp, DeltaImp
use CutoffParam, only : Rcutoff2, Rbook2, Rskin
use AtomParam, only : MolName, AtomName, ResidName, Mass, InvMass
use WallParam, only : ngrid_wall, invdz_size, WallFile, Rwmin, Rwmax, TabWall
use CylParam, only : CylFile, invdr_size, TabCyl, Rcylmin2, Rcylmax2, &
&   TabCyl, ngrid_cyl, FCylRs, TabCylZ, FscylR, FscylZ, Zcylmax, ngrid_cylZ, FCylH

implicit NONE

integer :: eofile
integer :: i, j, k, l, i1, i2, j1, j2, ii, jj, kk, ll
integer :: it, jt, kt, lt, m
integer :: ichtype, NAtom, Ipair, Jpair, num
integer, dimension(NumSpec) :: imolparam
character(len=6), dimension(:), allocatable :: Tmpname
logical, dimension(:), allocatable :: chflg
character(len=80), dimension(:), allocatable :: TmpFile
real(8), dimension(:), allocatable :: TmpRmin, TmpRmax, TmpZmax
character(len=80) :: String, String1
integer, dimension(:), allocatable :: AtomHands
integer, dimension(:,:), allocatable :: AtomHandBond
integer, dimension(:), allocatable :: TmpBondI, TmpBondJ
integer, dimension(:), allocatable :: TmpAngleI, TmpAngleJ, TmpAngleK
integer, dimension(:), allocatable :: TmpDihedI, TmpDihedJ, TmpDihedK, TmpDihedL
integer, dimension(:), allocatable :: TmpImproI, TmpImproJ, TmpImproK, TmpImproL
integer :: ItotalC, id, iatom, iunit, ibond, iangle, idihed, iimpr
real(8) :: TotalCharge, Rcutmin2, xx
real(8) :: Rcutoff, Rbook, Rmx, x, vp, fp, fpx, fpz, gridsize, fr, fh, z
real(8) :: Upot, dUdr, dUdz, d2Udr2, d2Udz2, d2Udrdz, d3Ud2rdz, d3Udrd2z
integer :: ixp, ixm, iyp, iym, ntmp
real(8) :: invx, invy, dx12, dx23, dy12, dy23
real(8) :: f11, f12, f21, f22
integer :: NbI, NbJ
integer, dimension(:), allocatable :: iFlagDihed

   allocate( Charge (N) )
   allocate( ScrnCR (N) )
   allocate( Mass   (N) )
   allocate( InvMass(N) )

   allocate( AtomHands(N) )
   allocate( AtomHandBond(N,6) )

   NAtom = 0

   do i = 1, NumSpec

     if(QMacro) then
       if(MolName(i)=='CGSP') then

         do j = 1, NumMol(i)
           do k = 1, NumAtm(i)
             NAtom = NAtom + 1
             Charge(NAtom) = 0.d0
             ScrnCR(NAtom) = 0.d0
             Mass  (NAtom) = 0.d0
           end do
         end do

         cycle

       end if
     end if

inn: do j = 1, NumMolParam
       if(MolName(i)==MolNameParam(j)) then
         imolparam(i) = j
         exit inn
       end if
       if(j==NumMolParam) then
         if(QMaster) then
           write(*,*) 'error : no data for your molname'
           write(*,*) 'mol. name = ',MolName(i)
         end if
         call Finalize
       end if
     end do inn

     if(NumAtm(i)/=NumAtom_inMol(imolparam(i))) then
       if(QMaster) write(*,*) 'error : consistency of NumAtm (configuration vs. topology data)'
       call Finalize
     end if

     do j = 1, NumMol(i)

       do k = 1, NumAtm(i)

         NAtom = NAtom + 1
         Charge(NAtom) = ChargeParam(k,imolparam(i))
         ScrnCR(NAtom) = ScrnParam  (k,imolparam(i))
         Mass  (NAtom) = AtomMass   (k,imolparam(i))

         if(AtomName(NAtom)/=AtomNameParam(k,imolparam(i))) then
           if(QMaster) write(*,*) 'error : consistency of atomname (cofig. vs topology)'
           call Finalize
         end if

         ii = 0
inn1:    do l = 1, NumResid_inMol(imolparam(i))
           ii = ii + NumAtom_inResidParam(l,imolparam(i))
           if(ii>=k) then
             jj = l
             exit inn1
           end if
         end do inn1

         if(ResidName(Natom)/=ResidNameParam(jj,imolparam(i))) then
           if(QMaster) write(*,*) 'error : consistency of residname (cofig. vs topology)'
           call Finalize
         end if

       end do
     end do

   end do

   if(NAtom/=N) then
     if(QMaster) write(*,*) 'error : check number of atoms'
     call Finalize
   end if

   if(QMacro) then
     call Set_Sphere(1)
   end if

! ----
! Unit
! ----
   Mass(:) = Mass(:) * 1.d-3 / Avogadro

   InvMass(:) = 1.d0 / Mass(:)

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

     xx = sqrt(ec/Eps_relative)
     Charge(:) = Charge(:) * xx
     ScrnCR(:) = ScrnCR(:) * 0.5d0

     Nel = 0
     do i = 1, N
       if(Charge(i)/=0.) Nel = Nel + 1
     end do

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

   if(QMaster.and.Qstdout) then
     if(QCoulomb) then
       write(*,*) 'Checking charge neutrality: total charge =',TotalCharge
     else
       write(*,*) 'No Coulomb interaction'
     end if
   end if

! ----------------------------------------------
! ## check number of atom types in the system

   ichtype = 0
   do i = 1, NumSpec
     if(QMacro.and.MolName(i)=='CGSP') cycle
     ichtype = ichtype + NumAtm(i)
   end do

   allocate( Tmpname(ichtype) )
   allocate( chflg(ichtype) )

   ii = 0
   do i = 1, NumSpec
     if(QMacro.and.MolName(i)=='CGSP') cycle
     do j = 1, NumAtom_inMol(imolparam(i))
       ii = ii + 1
       Tmpname(ii) = AtomType(j,imolparam(i))
     end do
   end do

   chflg(:) = .True.
   do i = 2, ichtype
     do j = 1, i-1
       if(Tmpname(i)==Tmpname(j)) then
         chflg(i) = .False.
       end if
     end do
   end do

   NumAtype = 0
   do i = 1, ichtype
     if(chflg(i)) NumAtype = NumAtype + 1
   end do

   if(QMaster.and.Qstdout) print *, 'NumAtype =', NumAtype

! ## 

   allocate( AtomtypeList(NumAtype) )

   ii = 0
   do i = 1, ichtype
     if(chflg(i)) then
       ii = ii + 1
       AtomtypeList(ii) = Tmpname(i)
     end if
   end do

   deallocate( Tmpname, chflg )

   if(QMaster.and.Qstdout) then
     do i = 1, ii
       print *, 'Type',i,' = ',AtomTypeList(i)
     end do
   end if

   allocate( NBAtomType(N) )

   ii = 0
   do i = 1, NumSpec
     if(QMacro.and.MolName(i)=='CGSP') then
       do j = 1, NumMol(i)
         do k = 1, NumAtm(i)
           ii = ii + 1  
           NBAtomType(ii) = 0
         end do
       end do
     else
       do j = 1, NumMol(i)
         do k = 1, NumAtm(i)
           ii = ii + 1
inn2:      do l = 1, NumAtype
             if(AtomType(k,imolparam(i))==AtomtypeList(l)) then
               id = l
               exit inn2
             end if
             if(l==NumAtype) then
               write(*,*) 'No atom type is available for component',i,' atom',k
               call Finalize
             end if
           end do inn2
           NBAtomType(ii) = id
         end do
       end do
     end if
   end do

! ## Macrosphere
   if(QMacro) then
     call Set_Sphere(2)
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

       do j = 1 , NumMol(i)

         do k = 1 , NumUnit_inMol(imolparam(i))

           iunit = iunit + 1

           do l = 1 , NumAtom_inUNITParam(k,imolparam(i))

             iatom = iatom + 1
             AtomUnitName(iatom) = UNITNameParam(k,imolparam(i))
             AtomUnitNum (iatom) = iunit

           end do

         end do

       end do

     end do

     NumRB = AtomUnitNum(N)

   end if

!
! ## Complete topology data
! ----------------------------------------------------------------------
! ## Bond 

   NumBond = 0

   do i = 1, NumSpec
     if(QMacro.and.MolName(i)=='CGSP') cycle
     NumBond = NumBond + NumMol(i)*NumBond_inMol(imolparam(i))
   end do

   if(NumBond/=0) then

     allocate( BondI(NumBond) )
     allocate( BondJ(NumBond) )

     ibond = 0

     do i = 1, NumSpec

       if(i==1) then
         ii = 0
       else
         ii = ii + NumMol(i-1)*NumAtm(i-1)
       end if

       if(QMacro.and.MolName(i)=='CGSP') cycle

       do j = 1, NumMol(i)

         jj = (j-1) * NumAtom_inMol(imolparam(i))

         do k = 1, NumBond_inMol(imolparam(i))

           ibond = ibond + 1
           BondI(ibond) = BondPair_inMol(1,k,imolparam(i)) + jj + ii
           BondJ(ibond) = BondPair_inMol(2,k,imolparam(i)) + jj + ii

         end do

       end do

     end do

   end if

! ## Improper Tosion 

   NumImproper = 0

   do i = 1, NumSpec
     if(QMacro.and.MolName(i)=='CGSP') cycle
     NumImproper = NumImproper + NumMol(i)*NumIMPR_inMol(imolparam(i))
   end do

   if(NumImproper/=0) then

     allocate( ImproI(NumImproper) )
     allocate( ImproJ(NumImproper) )
     allocate( ImproK(NumImproper) )
     allocate( ImproL(NumImproper) )

     iimpr = 0

     do i = 1, NumSpec

       if(i==1) then
         ii = 0
       else
         ii = ii + NumMol(i-1)*NumAtm(i-1)
       end if

       if(QMacro.and.MolName(i)=='CGSP') cycle

       do j = 1, NumMol(i)

         jj = (j-1) * NumAtom_inMol(imolparam(i))

         do k = 1, NumIMPR_inMol(imolparam(i))

           iimpr = iimpr + 1
           ImproI(iimpr) = ImprPair_inMol(1,k,imolparam(i)) + jj + ii
           ImproJ(iimpr) = ImprPair_inMol(2,k,imolparam(i)) + jj + ii
           ImproK(iimpr) = ImprPair_inMol(3,k,imolparam(i)) + jj + ii
           ImproL(iimpr) = ImprPair_inMol(4,k,imolparam(i)) + jj + ii

         end do

       end do

     end do

   end if

! ################ Angle and Dihedrals ###################
! ## Number of Bond per atom 

   AtomHands = 0
   AtomHandBond = 0

   do k = 1, NumBond

     i = BondI(k)
     j = BondJ(k)

     AtomHands(i) = AtomHands(i) + 1
     AtomHands(j) = AtomHands(j) + 1

     AtomHandBond(i,AtomHands(i)) = k
     AtomHandBond(j,AtomHands(j)) = k

   end do

! ## number of angles 

   NumAngle = 0

   do i = 1 , N

     k = 0
     do j = 1, AtomHands(i)-1
       k = k + j
     end do

     NumAngle = NumAngle + k

   end do

   if(NumAngle /= 0) then

     allocate( AngleI(NumAngle) )
     allocate( AngleJ(NumAngle) )
     allocate( AngleK(NumAngle) )

! ## atom numbers for angles 

     iangle = 0

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

           iangle = iangle + 1

           AngleJ(iangle) = i

           if(i1==i) then
             AngleI(iangle) = j1
           else if(j1==i) then
             AngleI(iangle) = i1
           end if

           if(i2==i) then
             AngleK(iangle) = j2
           else if(j2==i) then
             AngleK(iangle) = i2
           end if

         end do

       end do

     end do

   end if

! ## number of dihedrals

   NumDihedral = 0

   if(NumDihedParam/=0) then

     do k = 1 , NumBond

       i = BondI(k)
       j = BondJ(k)

       NbI = AtomHands(i)-1
       NbJ = AtomHands(j)-1
       ii  = NbI * NbJ

       NumDihedral = NumDihedral + ii

     end do

     if(NumDihedral/=0) then
       allocate( DihedI(NumDihedral) )
       allocate( DihedJ(NumDihedral) )
       allocate( DihedK(NumDihedral) )
       allocate( DihedL(NumDihedral) )

       idihed = 0

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

             idihed = idihed + 1

             DihedJ(idihed) = i
             DihedK(idihed) = j

             if(i1==i) then
               DihedI(idihed) = j1
             else if(j1==i) then
               DihedI(idihed) = i1
             end if

             if(i2==j) then
               DihedL(idihed) = j2
             else if(j2==j) then
               DihedL(idihed) = i2
             end if

           end do

         end do

       end do

     end if

   end if

   NumUB       = 0

   deallocate( AtomHands, AtomHandBond )

! ----------------------------------------------------------------------
! ## nonbonded parameters

   allocate( NBFuncType(NumAtype,NumAtype) )
   allocate( CoefAtype(NumAtype,NumAtype) )
   allocate( CoefBtype(NumAtype,NumAtype) )
   allocate( CoefCtype(NumAtype,NumAtype) )
   allocate( IDtable(NumAtype,NumAtype) )
   allocate( Rmin2(NumAtype,NumAtype) )
   allocate( Rcut2(NumAtype,NumAtype) )
   allocate( Rbk2(NumAtype,NumAtype) )
   allocate( Rsw2(NumAtype,NumAtype) )
   allocate( Swch(NumAtype,NumAtype) )
   allocate( Rc13(NumAtype,NumAtype) )
   allocate( Eps13(NumAtype,NumAtype) )

   kk = 0
   do i = 1, NumAtype
     do j = i, NumAtype

inn3:  do k = 1, NumNonBondParam

         if(((AtomTypeList(i)==NonBondPairAtoms(1,k)).and.   &
         &   (AtomTypeList(j)==NonBondPairAtoms(2,k))) .or.  &
         &  ((AtomTypeList(i)==NonBondPairAtoms(2,k)).and.   &
         &   (AtomTypeList(j)==NonBondPairAtoms(1,k)))) then

           ii = NonBondType(k)
           NBFuncType(i,j) = ii

           Rcut2(i,j)     = RmaxParam(k)**2
           Rmin2(i,j)     = RminParam(k)**2
           Rbk2 (i,j)     = (RmaxParam(k)+Rskin)**2
           Rc13 (i,j)     = R13cutParam(k)**2
           Eps13(i,j)     = EpsParam(k)

           if((ii==1).or.(ii==2).or.(ii==5).or.(ii==6).or.&
           &  (ii==7).or.(ii==8)) then ! LJ 
             CoefAtype(i,j) = AijParam(k)
             CoefBtype(i,j) = BijParam(k)
             CoefCtype(i,j) = 0.d0
             IDtable(i,j)   = 0
             Rsw2 (i,j)     = (RmaxParam(k)-2.d0)**2
             Swch (i,j)     = 1.d0 / (Rcut2(i,j) - Rsw2(i,j)) ** 3
           else if(ii==3) then ! tabulated 
             kk = kk + 1
             CoefAtype(i,j) = 0.d0
             CoefBtype(i,j) = 0.d0
             CoefCtype(i,j) = 0.d0
             IDtable(i,j)   = kk
           else if(ii==4) then ! morse
             CoefAtype(i,j) = AijParam(k)
             CoefBtype(i,j) = BijParam(k)
             CoefCtype(i,j) = CijParam(k)
             IDtable(i,j)   = 0
             Rsw2 (i,j)     = CijParam(k)**2
             Swch (i,j)     = 1.d0 / (Rcut2(i,j) - Rsw2(i,j)) ** 3
           else
             write(*,*) 'error : NonBondType, missing pair ', &
             &          AtomTypeList(i),AtomTypeList(j)
             call Finalize
           end if

           exit inn3

         end if

         if(k==NumNonBondParam) then
           write(*,*) 'error : NonBondType, missing pair ', &
           &          AtomTypeList(i),AtomTypeList(j)
           write(*,*) 'Check your parameter file: NONBOND'
           call Finalize
         end if

       end do inn3

     end do
   end do

   allocate( TabFunc(2000,2,kk) )
   allocate( Invgs(kk) )
   allocate( TabFileHead(kk) )

   do i = 1, NumAtype
     do j = i, NumAtype
       ii = IDtable(i,j)

       if(ii/=0) then

inn6:    do k = 1, NumNonBondParam
           if(((AtomTypeList(i)==NonBondPairAtoms(1,k)).and.   &
           &   (AtomTypeList(j)==NonBondPairAtoms(2,k))) .or.  &
           &  ((AtomTypeList(i)==NonBondPairAtoms(2,k)).and.   &
           &   (AtomTypeList(j)==NonBondPairAtoms(1,k)))) then
             TabFileHead(ii)= NonBondFileHeader(k)
             exit inn6
           end if
           if(k==NumNonBondParam) then
             if(QMaster) then 
               write(*,*) 'Error : FileName for table func. is not found'
               write(*,*) '        Pair: ',AtomTypeList(i),' - ',AtomTypeList(j)
             end if
             call Finalize
           end if
         end do inn6

         call Read_TableFunctions(TabFileHead(ii),Rmin2(i,j),ii)

       end if

     end do
   end do

   do j = 1, NumAtype-1
     do i = j+1, NumAtype
       NBFuncType(i,j) = NBFuncType(j,i)
       CoefAtype(i,j) = CoefAtype(j,i)
       CoefBtype(i,j) = CoefBtype(j,i)
       CoefCtype(i,j) = CoefCtype(j,i)
       IDtable(i,j)   = IDtable(j,i)
       Rcut2(i,j)     = Rcut2(j,i)
       Rmin2(i,j)     = Rmin2(j,i)
       Rbk2 (i,j)     = Rbk2 (j,i)
       Rsw2 (i,j)     = Rsw2 (j,i)
       Swch (i,j)     = Swch (j,i)
       Rc13 (i,j)     = Rc13 (j,i)
       Eps13(i,j)     = Eps13(j,i)
     end do
   end do

! ## Check for Rcutoff (real space cutoff distance for Coulomb interaction)
! ## Rcutoff should be less than that of any vdW interaction 

   if(QPBC) then

   Rcutmin2 = Rcutoff2
   do i = 1, NumAtype - 1
     do j = i, NumAtype

       if(Rcutmin2 > Rcut2(i,j)) then
         Rcutmin2 = Rcut2(i,j)
       end if

     end do
   end do

   if(Rcutmin2 < Rcutoff2) then
     if(QMaster) then
       write(*,*)
       write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       write(*,*) '  WARNING : Rcutoff has been changed !  '
       write(*,*) '  your choice = ', sqrt(Rcutoff2)
       write(*,*) '  new length  = ', sqrt(Rcutmin2)
       write(*,*) '  Make sure if your Ewald parameters make sense with the new cutoff length '
       write(*,*) '  In this software, the cutoff for the real-space Ewald sum is assumed '
       write(*,*) '  shorter than that of vdW interaction.'
       write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       write(*,*)
     end if
     Rcutoff2 = Rcutmin2
     Rcutoff = sqrt( Rcutoff2 )
     Rbook = Rcutoff + Rskin
     Rbook2 = Rbook * Rbook
   end if

   if(Rrespa > sqrt(Rcutoff2)) then
     Rrespa = sqrt(Rcutoff2) - Rheal
     if(QMaster) then
       write(*,'(a)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       write(*,'(a)') 'WARNING: RRESPA is set to RCUTOFF - RHEAL'
       write(*,'(a)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     end if
   end if
   Rcut_short2 = Rrespa ** 2
   Rbk_short2  = (Rrespa + Rheal + Rskin) ** 2

   end if

! ----------------------------------------------
! ## Wall potential

   if(QCGWALL) then

     open(19,file=trim(WallFile),status='old')

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

     if(num<NumAtype) then
       write(*,*) 'error : the number of wall files is fewer than needed'
       call Finalize
     end if

     allocate( invdz_size(NumAtype) )
     allocate( TabWall(2,ngrid_wall,NumAtype) )
     allocate( Rwmin(NumAtype) )
     allocate( Rwmax(NumAtype) )

     do i = 1, NumAtype
iwn:   do j = 1, num
         if(Tmpname(j)==AtomtypeList(i)) then
           open(99,file=trim(Tmpfile(j)),status='old')
           Rwmin(i) = TmpRmin(j)
           Rwmax(i) = TmpRmax(j)
           l = 0
           do k = 1, ngrid_wall
             read(99,*) x, vp, fp
             if(x>(Rwmin(i)-1.d-10)) then
               l = l + 1
               TabWall(1,l,i) = vp * kb ! U
               TabWall(2,l,i) = fp * kb ! dU/dZ
               if(l==1) then
                 Rwmin(i) = x
               end if
               Rmx = x
             end if
           end do
           if(Rmx < Rwmax(i)) then
             Rwmax(i) = Rmx
             write(*,*) 'WARNING: Rwmax for Atomtype',i,' has been changed to ',Rwmax(i)
           end if
           gridsize = (Rmx - Rwmin(i)) / (dble(l) - 1.d0)
           invdz_size(i) = 1.d0 / gridsize
           exit iwn
         end if
         if(j==num) then
           write(*,*) 'error : cannot find a wall file'
           call Finalize
         end if
       end do iwn
     end do

     deallocate( Tmpname, Tmpfile, TmpRmin, TmpRmax )

     do i = 1, NumAtype
       Rwmax(i) = Rwmax(i) - Rwmin(i)
     end do

   end if

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

     if(num<NumAtype) then
       write(*,*) 'error : the number of wall files is fewer than needed'
       call Finalize
     end if

     allocate( invdr_size(NumAtype) )
     allocate( TabCyl(ngrid_cyl,3,NumAtype) )
     allocate( Rcylmin2(NumAtype) )
     allocate( Rcylmax2(NumAtype) )

     do i = 1, NumAtype
irn:   do j = 1, num
         if(Tmpname(j)==AtomtypeList(i)) then
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

! ----------------------------------------------
! ## Cylinder potential (Finite Length)
   if(QFSCyl) then

     open(19,file=trim(CylFile),status='old')

     allocate( Tmpname(100) )
     allocate( Tmpfile(100) )
     allocate( TmpRmax(100) )
     allocate( TmpZmax(100) )

     num = 0
     do
       read(19,'(a80)',iostat=eofile) String1
       if(eofile == -1) exit
       String = trim(adjustl(String1))
       if(String(1:1)=='#'.or.String(1:1)=='!') cycle
       num = num + 1
       read(String,*) Tmpname(num), TmpRmax(num), TmpZmax(num), Tmpfile(num)
     end do

     close(19)

     if(num<NumAtype) then
       write(*,*) 'error : the number of wall files is fewer than needed'
       call Finalize
     end if

     allocate( TabCylZ(4,ngrid_cyl,ngrid_cylZ,5,NumAtype) )
     allocate( FscylR(ngrid_cyl ,NumAtype) )
     allocate( FscylZ(ngrid_cylZ,NumAtype) )
     allocate( Rcylmax2(NumAtype) )
     allocate( Zcylmax(NumAtype) )

     do i = 1, NumAtype

irn2:  do j = 1, num

         if(Tmpname(j)==AtomtypeList(i)) then
           open(99,file=trim(Tmpfile(j)),status='old')
           Rcylmax2(i) = TmpRmax(j)*TmpRmax(j)
           Zcylmax(i) = TmpZmax(j)
           l = 0

           do k = 1, ngrid_cylZ
             do kk = 1, ngrid_cyl
               read(99,*) x, z, Upot, dUdr, dUdz, d2Udr2, d2Udz2, d2Udrdz,  &
               &          d3Ud2rdz, d3Udrd2z
               !fr, fh

               FscylR(kk,i) = x!x
               FscylZ(k ,i) = z
               TabCylZ(1,kk,k,1,i) = Upot * ExParam     ! U
               TabCylZ(1,kk,k,2,i) = dUdr * ExParam     ! dU/dR=Fr
               TabCylZ(1,kk,k,3,i) = dUdz * ExParam     ! dU/dZ=Fz
               TabCylZ(2,kk,k,1,i) = TabCylZ(1,kk,k,2,i)
               TabCylZ(3,kk,k,1,i) = TabCylZ(1,kk,k,3,i)
               TabCylZ(4,kk,k,1,i) = d2Udrdz * ExParam
               TabCylZ(2,kk,k,2,i) = d2Udr2 * ExParam
               TabCylZ(3,kk,k,2,i) = TabCylZ(4,kk,k,1,i)
               TabCylZ(4,kk,k,2,i) = d3Ud2rdz * ExParam
               TabCylZ(2,kk,k,3,i) = TabCylZ(4,kk,k,1,i)
               TabCylZ(3,kk,k,3,i) = d2Udz2 * ExParam
               TabCylZ(4,kk,k,3,i) = d3Udrd2z * ExParam
!               TabCylZ(1,kk,k,4,i) = fr * kb      ! dU/dRadius
!               TabCylZ(1,kk,k,5,i) = fh * kb      ! dU/dH

             end do
           end do

           Rmx = FscylR(ngrid_cyl,i)**2
           if(Rmx < Rcylmax2(i)) then
             Rcylmax2(i) = Rmx
             write(*,*) 'WARNING: Rcylmax for Atomtype',i,' has been changed to ',&
             &          sqrt(Rcylmax2(i))
           end if
           exit irn2
         end if

         if(j==num) then
           write(*,*) 'error : cannot find a cylinder file'
           call Finalize
         end if

       end do irn2

     end do

!     do i = 1, NumAtype

!       do k = 1, ngrid_cylZ
!         do kk = 1, ngrid_cyl
!           do j = 1, 5
!             ixp = kk + 1
!             ixm = kk - 1
!             iyp = k  + 1
!             iym = k  - 1
!             if(kk==1) then
!               ixm = 1
!             else if(kk==ngrid_cyl) then
!               ixp = ngrid_cyl
!             end if
!             if(k==1) then
!               iym = 1
!             else if(k==ngrid_cylZ) then
!               iyp = ngrid_cylZ
!             end if
!             dx12 = FscylR(kk ,i) - FscylR(ixm,i)
!             dx23 = FscylR(ixp,i) - FscylR(kk ,i)
!             dy12 = FscylZ(k  ,i) - FscylZ(iym,i)
!             dy23 = FscylZ(iyp,i) - FscylZ(k  ,i)
!             if(dx12==0.) dx12=1.d-99
!             if(dx23==0.) dx23=1.d-99
!             if(dy12==0.) dy12=1.d-99
!             if(dy23==0.) dy23=1.d-99
!             invx = 1.d0/(FscylR(ixp,i)-FscylR(ixm,i))
!             invy = 1.d0/(FscylZ(iyp,i)-FscylZ(iym,i))

!             TabCylZ(2,kk,k,j,i) = &
!             &     ((TabCylZ(1,ixp,k  ,j,i)-TabCylZ(1,kk ,k  ,j,i))/dx23*dx12 + &
!             &      (TabCylZ(1,kk ,k  ,j,i)-TabCylZ(1,ixm,k  ,j,i))/dx12*dx23 ) * invx

!             TabCylZ(3,kk,k,j,i) = &
!             &     ((TabCylZ(1,kk ,iyp,j,i)-TabCylZ(1,kk ,k  ,j,i))/dy23*dy12 + &
!             &      (TabCylZ(1,kk ,k  ,j,i)-TabCylZ(1,kk ,iym,j,i))/dy12*dy23 ) * invy

!             f11 = (TabCylZ(1,kk ,k  ,j,i) - TabCylZ(1,kk ,iym,j,i) &
!             &     -TabCylZ(1,ixm,k  ,j,i) + TabCylZ(1,ixm,iym,j,i))/(dx12*dy12)
!             f12 = (TabCylZ(1,kk ,iyp,j,i) - TabCylZ(1,kk ,k  ,j,i) &
!             &     -TabCylZ(1,ixm,iyp,j,i) + TabCylZ(1,ixm,k  ,j,i))/(dx12*dy23)
!             f21 = (TabCylZ(1,ixp,k  ,j,i) - TabCylZ(1,ixp,iym,j,i) &
!             &     -TabCylZ(1,kk ,k  ,j,i) + TabCylZ(1,kk ,iym,j,i))/(dx23*dy12)
!             f22 = (TabCylZ(1,ixp,iyp,j,i) - TabCylZ(1,ixp,k  ,j,i) &
!             &     -TabCylZ(1,kk ,iyp,j,i) + TabCylZ(1,kk ,k  ,j,i))/(dx23*dy23)
!             TabCylZ(4,kk,k,j,i) = (dx23*dy23*f11 + dx23*dy12*f12 + &
!             &                      dx12*dy23*f21 + dx12*dy12*f22) * invx*invy

!           end do
!         end do
!       end do

!     end do

     deallocate( Tmpname, Tmpfile, TmpZmax, TmpRmax )

     FCylRs = 0.d0
     FCylH  = 0.d0

   end if

!-------------------------------------------------------
! ## Allocate intramolecular parameters

   if(NumBond /= 0) then

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

     allocate( FTypeBond( NumBond ) )
     allocate( kBond( NumBond ) )
     allocate( rBond( NumBond ) )

     do k = 1, NumBond

       i = BondI(k)
       j = BondJ(k)

       ii = NBAtomType(i)
       jj = NBAtomType(j)

inn4:  do kk = 1, NumBondParam
         if(((AtomtypeList(ii)==BondPairAtoms(1,kk)).and.  &
         &   (AtomtypeList(jj)==BondPairAtoms(2,kk)))  .or.&
         &  ((AtomtypeList(ii)==BondPairAtoms(2,kk)).and.  &
         &   (AtomtypeList(jj)==BondPairAtoms(1,kk)))) then
           kBond(k) = kBondParam(kk)
           rBond(k) = rBondParam(kk)
           FTypeBond(k) = BondType(kk)
           exit inn4
         end if
         if(kk == NumBondParam) then
           if(QMaster) then
             write(*,*) 'error : assign of bond parameter'
             write(*,*) 'Bond : No.',k,'{',i,j,'}'
           end if
           call Finalize
         end if
       end do inn4

     end do

     call NoLJList(1)

! ----------------------------------------
     if(QSHAKE) then

       call MakeSHAKEList

     else

       NSHAKEGroup = 0

     end if
! ----------------------------------------

   else

     call NoLJList(0)

   end if

   if(NumAngle /= 0) then

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

     allocate( FTypeAngle(NumAngle) )
     allocate( kTheta(NumAngle) )
     allocate( Theta0(NumAngle) )
     allocate( vdwSubtAng(NumAngle) )

     vdwSubtAng = .False.

     do l = 1, NumAngle

       i = AngleI(l)
       j = AngleJ(l)
       k = AngleK(l)

       ii = NBAtomType(i)
       jj = NBAtomType(j)
       kk = NBAtomType(k)

inn5:  do ll = 1, NumAngleParam
         if(((AtomtypeList(ii)==AnglePairAtoms(1,ll)).and.  &
         &   (AtomtypeList(jj)==AnglePairAtoms(2,ll)).and.  &
         &   (AtomtypeList(kk)==AnglePairAtoms(3,ll)))  .or.&
         &  ((AtomtypeList(ii)==AnglePairAtoms(3,ll)).and.  &
         &   (AtomtypeList(jj)==AnglePairAtoms(2,ll)).and.  &
         &   (AtomtypeList(kk)==AnglePairAtoms(1,ll)))) then

           FTypeAngle(l) = AngleType(ll)
           kTheta(l) = kThetaParam(ll)
           Theta0(l) = Theta0Param(ll)

           exit inn5

         end if
         if(ll == NumAngleParam) then
           if(QMaster) then
             write(*,*) 'error : assign of angle parameter'
             write(*,*) 'Angle : No.',l,'{',i,j,k,'}'
           end if
           call Finalize
         end if
       end do inn5

     end do

     call NoLJList(2)

     do i = 1, NumAngle

       i1 = AngleI(i)
       i2 = AngleK(i)

       do j = 1, NumBond

         j1 = BondI(j)
         j2 = BondJ(j)

         if( ( ( i1 == j1 ) .and. (i2 == j2 ) ) .or. &
         &   ( ( i1 == j2 ) .and. (i2 == j1 ) ) ) then

           vdwSubtAng(i) = .True.

           exit

         end if

       end do

       do j = i+1 , NumAngle

         j1 = AngleI(j)
         j2 = AngleK(j)

         if( ( ( i1 == j1 ) .and. (i2 == j2 ) ) .or. &
         &   ( ( i1 == j2 ) .and. (i2 == j1 ) ) ) then

           vdwSubtAng(i) = .True.

           exit

         end if

       end do

     end do

   end if

   ntmp = NumDihedral

   if(ntmp/=0.and.QRigidBody) then

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

     NumDihedral = num

   end if

   ntmp = NumDihedral

   if(ntmp/=0) then

     allocate( iFlagDihed(NumDihedral) )
     iFlagDihed(:) = 0

     num = 0
     do ii = 1, NumDihedral

       i = TmpDihedI(ii)
       j = TmpDihedJ(ii)
       k = TmpDihedK(ii)
       l = TmpDihedL(ii)

       it = NBAtomType(i)
       jt = NBAtomType(j)
       kt = NBAtomType(k)
       lt = NBAtomType(l)

inn7:  do ll = 1, NumDihedParam
         If(((AtomtypeList(it)==DihedPairAtoms(1,ll)).and.  &
         &   (AtomtypeList(jt)==DihedPairAtoms(2,ll)).and.  &
         &   (AtomtypeList(kt)==DihedPairAtoms(3,ll)).and.  &
         &   (AtomtypeList(lt)==DihedPairAtoms(4,ll)))  .or.&
         &  ((AtomtypeList(it)==DihedPairAtoms(4,ll)).and.  &
         &   (AtomtypeList(jt)==DihedPairAtoms(3,ll)).and.  &
         &   (AtomtypeList(kt)==DihedPairAtoms(2,ll)).and.  &
         &   (AtomtypeList(lt)==DihedPairAtoms(1,ll)))) then
           num = num + 1
           iFlagDihed(ii) = ll
           exit inn7
         end if
       end do inn7
     end do

     if(num/=0) then
       allocate( DihedI(num) )
       allocate( DihedJ(num) )
       allocate( DihedK(num) )
       allocate( DihedL(num) )

       allocate( FTypeDihed(num) )
       allocate( kChi(num) )
       allocate( DeltaDih(num) )
       allocate( NDih(num) )
       allocate( DupFlag(num) )
       allocate( DupkChi(6,num) )
       allocate( DupDeltaDih(6,num) )
       allocate( DupNDih(6,num) )

       DupFlag = 0

       jj = 0
       do ii = 1, NumDihedral
         if(iFlagDihed(ii)/=0) then
           jj = jj + 1
           DihedI(jj) = TmpDihedI(ii)
           DihedJ(jj) = TmpDihedJ(ii)
           DihedK(jj) = TmpDihedK(ii)
           DihedL(jj) = TmpDihedL(ii)

           ll = iFlagDihed(ii)
           FTypeDihed(jj) = DihedType(ll)
           kChi(jj)     = kChiParam(ll)
           DeltaDih(jj) = DeltaDihParam(ll)

           it = NBAtomType(DihedI(jj))
           jt = NBAtomType(DihedJ(jj))
           kt = NBAtomType(DihedK(jj))
           lt = NBAtomType(DihedL(jj))

           if(DihedType(ll) == 1) then
             NDih(jj) = NdihParam(ll)
             m = 0
             do kk = ll+1, NumDihedParam
               if(((AtomtypeList(it)==DihedPairAtoms(1,kk)).and.  &
               &   (AtomtypeList(jt)==DihedPairAtoms(2,kk)).and.  &
               &   (AtomtypeList(kt)==DihedPairAtoms(3,kk)).and.  &
               &   (AtomtypeList(lt)==DihedPairAtoms(4,kk)))  .or.&
               &  ((AtomtypeList(it)==DihedPairAtoms(4,kk)).and.  &
               &   (AtomtypeList(jt)==DihedPairAtoms(3,kk)).and.  &
               &   (AtomtypeList(kt)==DihedPairAtoms(2,kk)).and.  &
               &   (AtomtypeList(lt)==DihedPairAtoms(1,kk)))) then
                 m = m + 1
                 DupFlag(jj)                 = m - 1
                 DupkChi(DupFlag(jj),jj)     = kChiParam(kk)
                 DupDeltaDih(DupFlag(jj),jj) = DeltaDihParam(kk)
                 DupNDih(DupFlag(jj),jj)     = NdihParam(kk)
               end if
             end do                 
           end if

         end if
       end do

     end if

     NumDihedral = num

   end if

   ntmp = NumImproper

   if(ntmp/=0.and.QRigidBody) then

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

     if(NumImproper/=0) then

       allocate( ImproI(NumImproper) )
       allocate( ImproJ(NumImproper) )
       allocate( ImproK(NumImproper) )
       allocate( ImproL(NumImproper) )

       do ii = 1, NumImproper
         ImproI(ii) = TmpImproI(ii)
         ImproJ(ii) = TmpImproJ(ii)
         ImproK(ii) = TmpImproK(ii)
         ImproL(ii) = TmpImproL(ii)
       end do

     end if

     deallocate( TmpImproI, TmpImproJ, TmpImproK, TmpImproL )

   end if

   ntmp = NumImproper

   if(ntmp/=0) then

     allocate( FTypeImpro( NumImproper ) )
     allocate( kImp( NumImproper ) )
     allocate( NImp( NumImproper ) )
     allocate( DeltaImp(NumImproper) )

     do i1 = 1, NumImproper

       i = ImproI(i1)
       j = ImproJ(i1)
       k = ImproK(i1)
       l = ImproL(i1)

       ii = NBAtomType(i)
       jj = NBAtomType(j)
       kk = NBAtomType(k)
       ll = NBAtomType(l)

inn8:  do i2 = 1, NumImproParam
         if((AtomtypeList(ii)==ImproPairAtoms(1,i2)).and.  &
         &  (AtomtypeList(jj)==ImproPairAtoms(2,i2)).and.  &
         &  (AtomtypeList(kk)==ImproPairAtoms(3,i2)).and.  &
         &  (AtomtypeList(ll)==ImproPairAtoms(4,i2))) then
           kImp(i1) = kPsiParam(i2)
           DeltaImp(i1) = Psi0Param(i2)
           FTypeImpro(i1) = ImproType(i2)
           if(ImproType(i2)==1) then
             NImp(i1) = NimpParam(i2)
           else
             NImp(i1) = 0
           end if
           exit inn8
         end if
         if(i2 == NumImproParam) then
           if(QMaster) then
             write(*,*) 'error : assign of improper parameter'
             write(*,*) 'Bond : No.',i1,'{',i,j,k,l,'}'
           end if
           call Finalize
         end if
       end do inn8

     end do

   end if

end subroutine AllocateCGdata


!######################################################################
!######################################################################


subroutine Read_TableFunctions(FileHeader,Rmin,ii)

use UnitExParam, only : ExParam
use TableFuncs

character(len=13) :: FileHeader
character(len=17) :: TabFileName
integer :: ii, kk, i, ntable
real(8) ::r, U, F, r2
real(8) :: Rmin,Rmax, GridSize

   write(TabFileName,'(a,a)') trim(adjustl(FileHeader)),'.ftable'

   open(99,file=trim(TabFileName))

   kk = 0
   do i = 1, 2000
     read(99,*) r, U, F
     r2 = r * r
     if(r2>(Rmin-1.d-10)) then
       kk = kk + 1
       TabFunc(kk,1,ii) = U * ExParam
       TabFunc(kk,2,ii) = F * ExParam / r ! (dU/dR)*(1/R) 
       if(kk==1) then
         Rmin = r2
       end if
       Rmax = r2
     end if
   end do
   ntable = kk

   GridSize  = ( Rmax - Rmin ) / (dble(ntable) - 1.d0)
   Invgs(ii) = 1.d0 / GridSize

   close(99)

end subroutine Read_TableFunctions
