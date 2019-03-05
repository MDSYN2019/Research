! ############################
! ## SUBROUTINE LIST 
! ## -- SetupMPI 
! ## -- Synch 
! ## -- FinMPI 
! ## -- AllocPara 
! ## -- SumFrc 
! ## -- SumEnergy 
! ## -- SumVir 
! ## -- BcastE 
! ## -- BcastR 
! ## -- BcastRH 
! ## -- BcastRgQuat 
! ## -- BcastRgQuatH 
! ## -- SumRDF 
! ## -- SumEneInsertion 
! ## -- SumRho 
! ## -- SumDistPI
! ## -- SumFrcPI
! ## -- SumFrcPI_DCV
! ## -- RnmPI
! ## -- VnmPI
! ## -- RVBathPI
! ## -- SumEnePI
! ## -- SumTempPI
! ## -- SumFrcgA
! ## -- Prep_Atom_to_Mesh 
! ## -- SumChargeDens 
! ## -- DistChargeDens 
! ## -- SumEneTICR 
! ## -- SumPMF_CGball 
! ## -- SumHC 
! ## -- FFT_ChAxisF 
! ## -- FFT_ChAxisB 
! ## -- SumER1
! ############################


!######################################################################
!######################################################################


subroutine SetupMPI

use CommonBlocks, only : QMaster
use CommonMPI, only : ierror, NProcs, MyRank

implicit none

include 'mpif.h'

   call Mpi_Init(ierror)
   call Mpi_Comm_Size(Mpi_Comm_World, NProcs, ierror)
   call Mpi_Comm_Rank(Mpi_Comm_World, MyRank, ierror)

#ifdef SCINT
   call m64_init()
#endif

   QMaster = .False.
   if( MyRank == 0 ) QMaster = .True.

end subroutine SetupMPI


!######################################################################
!######################################################################


subroutine Synch(NN)

use Numbers, only : N
use CommonBlocks, only : QPathInt, cCOULOMB
use BathParam, only : NHchain
use CommonPI, only : Nbead
use CommonMPI
use PMEparam, only : NfftDim

implicit none

include 'mpif.h'

integer :: NN,NF,NDIM

   if(NN==0) then
     NDIM = 3 * N + 10

     if(cCOULOMB == 'PME') then
       NF = NfftDim(1) * NfftDim(2) * NfftDim(3)
       if(NF > NDIM) NDIM = NF + 1
     end if

     if(QPathInt) then
       NDIM = 3 * N * Nbead * NHchain
     end if

     allocate( Buff1( NDIM ) )
     allocate( Buff2( NDIM ) )

     Buff1 = 0.
     Buff2 = 0.
   end if

   call Mpi_Barrier(Mpi_Comm_World, ierror)

end subroutine Synch


!######################################################################
!######################################################################


subroutine FinMPI

use CommonMPI

implicit none

include 'mpif.h'

   call Mpi_Barrier(Mpi_Comm_World, ierror)

   call Mpi_Finalize(ierror)

end subroutine FinMPI


!######################################################################
!######################################################################


subroutine AllocPara

use CommonBlocks, only : QOption, ForceField, QPathInt
use CommonMPI
use CommonPI
use OptConstraintParam, only : NumOptC, OptCI, OptCJ, kOptC, rOptC
use BondedParam, only : NumBond, NumAngle, NumUB, NumDihedral, NumImproper, &
&   BondI, BondJ, kBond, rBond, AngleI, AngleJ, AngleK, kTheta, Theta0,     &
&   UB_I, UB_J, Kub, S0, DihedI, DihedJ, DihedK, DihedL, vdWSubtDih,        &
&   kChi, DeltaDih, NDih, DupFlag, DupkChi, DupDeltaDih, DupNDih, Ts,       &
&   CsDelDih, DupCsDelDih, ImproI, ImproJ, ImproK, ImproL, kPsi, PsiImp,    &
&   kImp, DeltaImp, NImp, CsDelImp, FTypeBond, FTypeAngle, FTypeDihed, FTypeImpro

implicit none

include 'mpif.h'

integer :: Count, i, NBD, NAN, NUB, NDH, NIM, NOC
integer :: TmpNBond, TmpNAngle, TmpNUB, TmpNDihedral, TmpNImproper
integer :: NAddition
integer, dimension(:), allocatable :: TmpBondI, TmpBondJ, TmpBondType
real(8), dimension(:), allocatable :: TmpkBond, TmprBond
integer, dimension(:), allocatable :: TmpAngleI, TmpAngleJ, TmpAngleK
integer, dimension(:), allocatable :: TmpAngleType
real(8), dimension(:), allocatable :: TmpkTheta, TmpTheta0
integer, dimension(:), allocatable :: TmpUB_I,TmpUB_J
real(8), dimension(:), allocatable :: TmpKub,TmpS0
integer, dimension(:), allocatable :: TmpDihedI,TmpDihedJ,TmpDihedK,TmpDihedL
integer, dimension(:), allocatable :: TmpDihedType
real(8), dimension(:), allocatable :: TmpkChi,TmpDeltaDih
integer, dimension(:), allocatable :: TmpNDih,TmpDupFlag
real(8), dimension(:,:), allocatable :: TmpDupkChi,TmpDupDeltaDih
integer, dimension(:,:), allocatable :: TmpDupNDih
logical, dimension(:), allocatable :: TmpvdWSubtDih
integer, dimension(:), allocatable :: TmpImproI,TmpImproJ,TmpImproK,TmpImproL
integer, dimension(:), allocatable :: TmpImproType
real(8), dimension(:), allocatable :: TmpkPsi,TmpPsiImp
integer :: TmpNOptC
integer, dimension(:), allocatable :: TmpOptCI,TmpOptCJ
real(8), dimension(:), allocatable :: TmpkOptC,TmprOptC
real(8), dimension(:), allocatable :: TmpCsDelDih
real(8), dimension(:,:), allocatable :: TmpDupCsDelDih
real(8), dimension(:), allocatable :: TmpkImp,TmpDeltaImp,TmpCsDelImp
integer, dimension(:), allocatable :: TmpNImp

integer :: MyRankTemp, NProcsTemp

   NBD = NumBond
   NAN = NumAngle
   NUB = NumUB
   NDH = NumDihedral
   NIM = NumImproper
   NOC = NumOptC

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   call Mpi_Barrier(Mpi_Comm_World,ierror)

! ## Bond Parameters

   if(NumBond /= 0) then

     allocate( TmpBondI(NBD) )
     allocate( TmpBondJ(NBD) )
     allocate( TmpkBond(NBD) )
     allocate( TmprBond(NBD) )
     if(ForceField(1:2)=='CG') then
       allocate( TmpBondType(NBD) )
       TmpBondType = FTypeBond
       deallocate( FTypeBond )
     end if

     TmpBondI = BondI
     TmpBondJ = BondJ
     TmpkBond = kBond
     TmprBond = rBond

     TmpNBond  = NumBond / NProcsTemp
     NAddition = mod(NumBond,NProcsTemp)

     if(MyRankTemp < NAddition) TmpNBond = TmpNBond + 1

     deallocate( BondI, BondJ, kBond, rBond )

     if(TmpNBond /= 0) then

       allocate( BondI(TmpNBond) )
       allocate( BondJ(TmpNBond) )
       allocate( kBond(TmpNBond) )
       allocate( rBond(TmpNBond) )
       if(ForceField(1:2) == 'CG') then
         allocate( FTypeBond(TmpNBond) )
       end if

       Count = 0

       do i = MyRankTemp+1, NumBond, NProcsTemp

         Count = Count + 1
         BondI(Count) = TmpBondI(i)
         BondJ(Count) = TmpBondJ(i)
         kBond(Count) = TmpkBond(i)
         rBond(Count) = TmprBond(i)
         if(ForceField(1:2)=='CG') FTypeBond(Count) = TmpBondType(i)

       end do

       deallocate( TmpBondI, TmpBondJ, TmpkBond, TmprBond )
       if(ForceField(1:2)=='CG') deallocate( TmpBondType )

     end if

     NumBond = TmpNBond

   end if

! ## Angle Parameters

   if(NumAngle /= 0) then

     allocate( TmpAngleI(NAN) )
     allocate( TmpAngleJ(NAN) )
     allocate( TmpAngleK(NAN) )
     allocate( TmpkTheta(NAN) )
     allocate( TmpTheta0(NAN) )

     TmpAngleI = AngleI
     TmpAngleJ = AngleJ
     TmpAngleK = AngleK
     TmpkTheta = kTheta
     TmpTheta0 = Theta0

     if(ForceField(1:2)=='CG') then
       allocate( TmpAngleType(NAN) )
       TmpAngleType = FTypeAngle
       deallocate( FTypeAngle)
     end if

     TmpNAngle = NumAngle / NProcsTemp
     NAddition = mod(NumAngle, NProcsTemp)

     if(MyRankTemp < NAddition) TmpNAngle = TmpNAngle + 1

     deallocate( AngleI, AngleJ, AngleK, kTheta, Theta0 )

     if(TmpNAngle /= 0) then

       allocate( AngleI(TmpNAngle) )
       allocate( AngleJ(TmpNAngle) )
       allocate( AngleK(TmpNAngle) )
       allocate( kTheta(TmpNAngle) )
       allocate( Theta0(TmpNAngle) )
       if(ForceField(1:2)=='CG') allocate( FTypeAngle(TmpNAngle) )

       Count = 0

       do i = MyRankTemp+1, NumAngle, NProcsTemp

         Count = Count + 1
         AngleI(Count) = TmpAngleI(i)
         AngleJ(Count) = TmpAngleJ(i)
         AngleK(Count) = TmpAngleK(i)
         kTheta(Count) = TmpkTheta(i)
         Theta0(Count) = TmpTheta0(i)
         if(ForceField(1:2)=='CG') FTypeAngle(Count) = TmpAngleType(i)

       end do

       deallocate( TmpAngleI, TmpAngleJ, TmpAngleK, TmpkTheta, TmpTheta0 )
       if(ForceField(1:2)=='CG') deallocate( TmpAngleType )

     end if

     NumAngle = TmpNAngle

   end if

! ## UB Parameters

   if(NumUB /= 0) then

     allocate( TmpUB_I(NUB) )
     allocate( TmpUB_J(NUB) )
     allocate( TmpKub(NUB) )
     allocate( TmpS0(NUB) )

     TmpUB_I = UB_I
     TmpUB_J = UB_J
     TmpKub  = Kub
     TmpS0   = S0

     TmpNUB    = NumUB / NProcsTemp
     NAddition = mod(NumUB,NProcsTemp)

     if(MyRankTemp < NAddition) TmpNUB = TmpNUB + 1

     deallocate( UB_I, UB_J, Kub, S0 )

     if(TmpNUB /= 0) then

       allocate( UB_I(TmpNUB) )
       allocate( UB_J(TmpNUB) )
       allocate( Kub (TmpNUB) )
       allocate( S0  (TmpNUB) )

       Count = 0

       do i = MyRankTemp+1, NumUB, NProcsTemp

         Count = Count + 1
         UB_I(Count) = TmpUB_I(i)
         UB_J(Count) = TmpUB_J(i)
         Kub (Count) = TmpKub (i)
         S0  (Count) = TmpS0  (i)

       end do

       deallocate( TmpUB_I, TmpUB_J, TmpKub, TmpS0 )

     end if

     NumUB = TmpNUB

   end if

! ## Dihedral parameters

   if(NumDihedral /= 0) then

     if(ForceField(1:5)=='CHARM') then

       allocate( TmpDihedI(NDH) )
       allocate( TmpDihedJ(NDH) )
       allocate( TmpDihedK(NDH) )
       allocate( TmpDihedL(NDH) )
       allocate( TmpkChi(NDH) )
       allocate( TmpDeltaDih(NDH) )
       allocate( TmpNDih(NDH) )
       allocate( TmpDupFlag(NDH) )
       allocate( TmpDupkChi(3,NDH) )
       allocate( TmpDupDeltaDih(3,NDH) )
       allocate( TmpDupNDih(3,NDH) )
       allocate( TmpvdWSubtDih(NDH) )

       TmpDihedI   = DihedI
       TmpDihedJ   = DihedJ
       TmpDihedK   = DihedK
       TmpDihedL   = DihedL
       TmpkChi     = kChi
       TmpDeltaDih = DeltaDih
       TmpNDih     = NDih

       TmpDupFlag     = DupFlag
       TmpDupkChi     = DupkChi
       TmpDupDeltaDih = DupDeltaDih
       TmpDupNDih     = DupNDih

       TmpvdWSubtDih  = vdWSubtDih

       TmpNDihedral  = NumDihedral / NProcsTemp
       NAddition = mod(NumDihedral,NProcsTemp)

       if(MyRankTemp < NAddition) TmpNDihedral = TmpNDihedral + 1

       deallocate( DihedI, DihedJ, DihedK, DihedL, kChi, DeltaDih, NDih )
       deallocate( DupFlag, DupkChi, DupDeltaDih, DupNDih, vdWSubtDih )

       if(TmpNDihedral /= 0) then

         allocate( DihedI(TmpNDihedral) )
         allocate( DihedJ(TmpNDihedral) )
         allocate( DihedK(TmpNDihedral) )
         allocate( DihedL(TmpNDihedral) )
         allocate( kChi    (TmpNDihedral) )
         allocate( DeltaDih(TmpNDihedral) )
         allocate( NDih    (TmpNDihedral) )

         allocate( DupFlag    (  TmpNDihedral) )
         allocate( DupkChi    (3,TmpNDihedral) )
         allocate( DupDeltaDih(3,TmpNDihedral) )
         allocate( DupNDih    (3,TmpNDihedral) )

         allocate( vdWSubtDih(TmpNDihedral) )

         Count = 0

         do i = MyRankTemp+1, NumDihedral, NProcsTemp

           Count = Count + 1

           DihedI(Count) = TmpDihedI(i)
           DihedJ(Count) = TmpDihedJ(i)
           DihedK(Count) = TmpDihedK(i)
           DihedL(Count) = TmpDihedL(i)

           kChi    (Count) = TmpkChi    (i)
           DeltaDih(Count) = TmpDeltaDih(i)
           NDih    (Count) = TmpNDih    (i)

           DupFlag    (  Count) = TmpDupFlag    (  i)
           DupkChi    (:,Count) = TmpDupkChi    (:,i)
           DupDeltaDih(:,Count) = TmpDupDeltaDih(:,i)
           DupNDih    (:,Count) = TmpDupNDih    (:,i)

           vdWSubtDih(Count) = TmpvdWSubtDih(i)

         end do

         deallocate( TmpDihedI, TmpDihedJ, TmpDihedK, TmpDihedL, TmpkChi, TmpDeltaDih )
         deallocate( TmpNDih, TmpDupFlag, TmpDupkChi, TmpDupDeltaDih, TmpDupNDih )
         deallocate( TmpvdWSubtDih )

       end if

       NumDihedral = TmpNDihedral

     else if(ForceField(1:4) == 'OPLS') then

       allocate( TmpDihedI(NDH) )
       allocate( TmpDihedJ(NDH) )
       allocate( TmpDihedK(NDH) )
       allocate( TmpDihedL(NDH) )
       allocate( TmpkChi(NDH) )
       allocate( TmpDeltaDih(NDH) )
       allocate( TmpNDih(NDH) )
       allocate( TmpDupFlag(NDH) )
       allocate( TmpDupkChi(3,NDH) )
       allocate( TmpDupDeltaDih(3,NDH) )
       allocate( TmpDupNDih(3,NDH) )
       allocate( TmpCsDelDih(NDH) )
       allocate( TmpDupCsDelDih(3,NDH) )
       allocate( TmpvdWSubtDih(NDH) )

       TmpDihedI   = DihedI
       TmpDihedJ   = DihedJ
       TmpDihedK   = DihedK
       TmpDihedL   = DihedL
       TmpkChi     = kChi
       TmpDeltaDih = DeltaDih
       TmpNDih     = NDih
       TmpCsDelDih = CsDelDih

       TmpDupFlag     = DupFlag
       TmpDupkChi     = DupkChi
       TmpDupDeltaDih = DupDeltaDih
       TmpDupNDih     = DupNDih
       TmpDupCsDelDih = DupCsDelDih

       TmpvdWSubtDih  = vdWSubtDih

       TmpNDihedral  = NumDihedral / NProcsTemp
       NAddition = mod(NumDihedral,NProcsTemp)

       if(MyRankTemp < NAddition) TmpNDihedral = TmpNDihedral + 1

       deallocate( DihedI, DihedJ, DihedK, DihedL, kChi, DeltaDih, NDih )
       deallocate( DupFlag, DupkChi, DupDeltaDih, DupNDih, vdWSubtDih )
       deallocate( CsDelDih, DupCsDelDih )

       if(TmpNDihedral /= 0) then

         allocate( DihedI(TmpNDihedral) )
         allocate( DihedJ(TmpNDihedral) )
         allocate( DihedK(TmpNDihedral) )
         allocate( DihedL(TmpNDihedral) )
         allocate( kChi    (TmpNDihedral) )
         allocate( DeltaDih(TmpNDihedral) )
         allocate( NDih    (TmpNDihedral) )
         allocate( CsDelDih(TmpNDihedral) )

         allocate( DupFlag    (  TmpNDihedral) )
         allocate( DupkChi    (3,TmpNDihedral) )
         allocate( DupDeltaDih(3,TmpNDihedral) )
         allocate( DupNDih    (3,TmpNDihedral) )
         allocate( DupCsDelDih(3,TmpNDihedral) )

         allocate( vdWSubtDih(TmpNDihedral) )

         Count = 0

         do i = MyRankTemp+1, NumDihedral, NProcsTemp

           Count = Count + 1

           DihedI(Count) = TmpDihedI(i)
           DihedJ(Count) = TmpDihedJ(i)
           DihedK(Count) = TmpDihedK(i)
           DihedL(Count) = TmpDihedL(i)

           kChi    (Count) = TmpkChi    (i)
           DeltaDih(Count) = TmpDeltaDih(i)
           NDih    (Count) = TmpNDih    (i)
           CsDelDih(Count) = TmpCsDelDih(i)

           DupFlag    (  Count) = TmpDupFlag    (  i)
           DupkChi    (:,Count) = TmpDupkChi    (:,i)
           DupDeltaDih(:,Count) = TmpDupDeltaDih(:,i)
           DupNDih    (:,Count) = TmpDupNDih    (:,i)
           DupCsDelDih(:,Count) = TmpDupCsDelDih(:,i)

           vdWSubtDih(Count) = TmpvdWSubtDih(i)

         end do

         deallocate( TmpDihedI, TmpDihedJ, TmpDihedK, TmpDihedL, TmpkChi, TmpDeltaDih )
         deallocate( TmpNDih, TmpDupFlag, TmpDupkChi, TmpDupDeltaDih, TmpDupNDih )
         deallocate( TmpvdWSubtDih )
         deallocate( TmpCsDelDih, TmpDupCsDelDih )

       end if

       NumDihedral = TmpNDihedral

     else if(ForceField(1:2)=='CG') then

       allocate( TmpDihedI(NDH) )
       allocate( TmpDihedJ(NDH) )
       allocate( TmpDihedK(NDH) )
       allocate( TmpDihedL(NDH) )
       allocate( TmpDihedType(NDH) )
       allocate( TmpkChi(NDH) )
       allocate( TmpDeltaDih(NDH) )
       allocate( TmpNDih(NDH) )
       allocate( TmpDupFlag(NDH) )
       allocate( TmpDupkChi(3,NDH) )
       allocate( TmpDupDeltaDih(3,NDH) )
       allocate( TmpDupNDih(3,NDH) )

       TmpDihedI   = DihedI
       TmpDihedJ   = DihedJ
       TmpDihedK   = DihedK
       TmpDihedL   = DihedL
       TmpDihedType= FTypeDihed
       TmpkChi     = kChi
       TmpDeltaDih = DeltaDih
       TmpNDih     = NDih

       TmpDupFlag     = DupFlag
       TmpDupkChi     = DupkChi
       TmpDupDeltaDih = DupDeltaDih
       TmpDupNDih     = DupNDih

       TmpNDihedral  = NumDihedral / NProcsTemp
       NAddition = mod(NumDihedral,NProcsTemp)

       if(MyRankTemp < NAddition) TmpNDihedral = TmpNDihedral + 1

       deallocate( DihedI, DihedJ, DihedK, DihedL, kChi, DeltaDih, NDih )
       deallocate( DupFlag, DupkChi, DupDeltaDih, DupNDih, FTypeDihed )

       if(TmpNDihedral /= 0) then

         allocate( DihedI(TmpNDihedral) )
         allocate( DihedJ(TmpNDihedral) )
         allocate( DihedK(TmpNDihedral) )
         allocate( DihedL(TmpNDihedral) )
         allocate( FTypeDihed(TmpNDihedral) )
         allocate( kChi    (TmpNDihedral) )
         allocate( DeltaDih(TmpNDihedral) )
         allocate( NDih    (TmpNDihedral) )

         allocate( DupFlag    (  TmpNDihedral) )
         allocate( DupkChi    (3,TmpNDihedral) )
         allocate( DupDeltaDih(3,TmpNDihedral) )
         allocate( DupNDih    (3,TmpNDihedral) )

         Count = 0

         do i = MyRankTemp+1, NumDihedral, NProcsTemp

           Count = Count + 1

           DihedI(Count) = TmpDihedI(i)
           DihedJ(Count) = TmpDihedJ(i)
           DihedK(Count) = TmpDihedK(i)
           DihedL(Count) = TmpDihedL(i)

           FTypeDihed(Count) = TmpDihedType(i)

           kChi    (Count) = TmpkChi    (i)
           DeltaDih(Count) = TmpDeltaDih(i)
           NDih    (Count) = TmpNDih    (i)

           DupFlag    (  Count) = TmpDupFlag    (  i)
           DupkChi    (:,Count) = TmpDupkChi    (:,i)
           DupDeltaDih(:,Count) = TmpDupDeltaDih(:,i)
           DupNDih    (:,Count) = TmpDupNDih    (:,i)

         end do

         deallocate( TmpDihedI, TmpDihedJ, TmpDihedK, TmpDihedL, TmpkChi, TmpDeltaDih )
         deallocate( TmpNDih, TmpDupFlag, TmpDupkChi, TmpDupDeltaDih, TmpDupNDih )
         deallocate( TmpDihedType )

       end if

       NumDihedral = TmpNDihedral

     end if

   end if

! ## Improper Parameters

   if(NumImproper /= 0) then

     if(ForceField(1:5)=='CHARM') then

       allocate( TmpImproI(NIM) )
       allocate( TmpImproJ(NIM) )
       allocate( TmpImproK(NIM) )
       allocate( TmpImproL(NIM) )
       allocate( TmpkPsi(NIM) )
       allocate( TmpPsiImp(NIM) )

       TmpImproI = ImproI
       TmpImproJ = ImproJ
       TmpImproK = ImproK
       TmpImproL = ImproL
       TmpkPsi   = kPsi
       TmpPsiImp = PsiImp

       TmpNImproper = NumImproper / NProcsTemp
       NAddition = mod(NumImproper,NProcsTemp)

       if(MyRankTemp < NAddition) TmpNImproper = TmpNImproper + 1

       deallocate( ImproI, ImproJ, ImproK, ImproL, kPsi, PsiImp )

       if(TmpNImproper /= 0) then

         allocate( ImproI(TmpNImproper) )
         allocate( ImproJ(TmpNImproper) )
         allocate( ImproK(TmpNImproper) )
         allocate( ImproL(TmpNImproper) )
         allocate( kPsi  (TmpNImproper) )
         allocate( PsiImp(TmpNImproper) )

         Count = 0

         do i = MyRankTemp+1, NumImproper, NProcsTemp

           Count = Count + 1
           ImproI(Count) = TmpImproI(i)
           ImproJ(Count) = TmpImproJ(i)
           ImproK(Count) = TmpImproK(i)
           ImproL(Count) = TmpImproL(i)
           kPsi  (Count) = TmpkPsi  (i)
           PsiImp(Count) = TmpPsiImp(i)

         end do

         deallocate( TmpImproI, TmpImproJ, TmpImproK, TmpImproL, TmpkPsi, TmpPsiImp )

       end if

       NumImproper = TmpNImproper

     else if(ForceField(1:4) == 'OPLS') then

       allocate( TmpImproI(NIM) )
       allocate( TmpImproJ(NIM) )
       allocate( TmpImproK(NIM) )
       allocate( TmpImproL(NIM) )
       allocate( TmpkImp(NIM) )
       allocate( TmpDeltaImp(NIM) )
       allocate( TmpNImp(NIM) )
       allocate( TmpCsDelImp(NIM) )

       TmpImproI = ImproI
       TmpImproJ = ImproJ
       TmpImproK = ImproK
       TmpImproL = ImproL

       TmpkImp     = kImp
       TmpDeltaImp = DeltaImp
       TmpNImp     = NImp
       TmpCsDelImp = CsDelImp

       TmpNImproper = NumImproper / NProcsTemp
       NAddition = mod(NumImproper,NProcsTemp)

       if(MyRankTemp < NAddition) TmpNImproper = TmpNImproper + 1

       deallocate( ImproI, ImproJ, ImproK, ImproL, kImp, DeltaImp, NImp, CsDelImp )

       if(TmpNImproper /= 0) then

         allocate( ImproI(TmpNImproper) )
         allocate( ImproJ(TmpNImproper) )
         allocate( ImproK(TmpNImproper) )
         allocate( ImproL(TmpNImproper) )

         allocate( kImp    (TmpNImproper) )
         allocate( DeltaImp(TmpNImproper) )
         allocate( NImp    (TmpNImproper) )
         allocate( CsDelImp(TmpNImproper) )

         Count = 0

         do i = MyRankTemp+1, NumImproper, NProcsTemp

           Count = Count + 1
           ImproI(Count) = TmpImproI(i)
           ImproJ(Count) = TmpImproJ(i)
           ImproK(Count) = TmpImproK(i)
           ImproL(Count) = TmpImproL(i)

           kImp    (Count) = TmpkImp    (i)
           DeltaImp(Count) = TmpDeltaImp(i)
           NImp    (Count) = TmpNImp    (i)
           CsDelImp(Count) = TmpCsDelImp(i)

         end do

         deallocate( TmpImproI, TmpImproJ, TmpImproK, TmpImproL, TmpkImp, TmpDeltaImp )
         deallocate( TmpNImp, TmpCsDelImp )

       end if

       NumImproper = TmpNImproper

     else if(ForceField(1:2) == 'CG') then

       allocate( TmpImproI(NIM) )
       allocate( TmpImproJ(NIM) )
       allocate( TmpImproK(NIM) )
       allocate( TmpImproL(NIM) )
       allocate( TmpImproType(NIM) )
       allocate( TmpkImp(NIM) )
       allocate( TmpDeltaImp(NIM) )
       allocate( TmpNImp(NIM) )

       TmpImproI = ImproI
       TmpImproJ = ImproJ
       TmpImproK = ImproK
       TmpImproL = ImproL

       TmpImproType = FTypeImpro

       TmpkImp     = kImp
       TmpDeltaImp = DeltaImp
       TmpNImp     = NImp

       TmpNImproper = NumImproper / NProcsTemp
       NAddition = mod(NumImproper,NProcsTemp)

       if(MyRankTemp < NAddition) TmpNImproper = TmpNImproper + 1

       deallocate( ImproI, ImproJ, ImproK, ImproL, FtypeImpro, kImp, DeltaImp, NImp )

       if(TmpNImproper /= 0) then

         allocate( ImproI(TmpNImproper) )
         allocate( ImproJ(TmpNImproper) )
         allocate( ImproK(TmpNImproper) )
         allocate( ImproL(TmpNImproper) )

         allocate( FTypeImpro(TmpNImproper) )

         allocate( kImp    (TmpNImproper) )
         allocate( DeltaImp(TmpNImproper) )
         allocate( NImp    (TmpNImproper) )

         Count = 0

         do i = MyRankTemp+1, NumImproper, NProcsTemp

           Count = Count + 1
           ImproI(Count) = TmpImproI(i)
           ImproJ(Count) = TmpImproJ(i)
           ImproK(Count) = TmpImproK(i)
           ImproL(Count) = TmpImproL(i)

           FTypeImpro(Count) = TmpImproType(i)

           kImp    (Count) = TmpkImp    (i)
           DeltaImp(Count) = TmpDeltaImp(i)
           NImp    (Count) = TmpNImp    (i)

         end do

         deallocate( TmpImproI, TmpImproJ, TmpImproK, TmpImproL, TmpkImp, TmpDeltaImp )
         deallocate( TmpNImp, TmpImproType )

       end if

       NumImproper = TmpNImproper

     end if

   end if

! ## Optional Bond Parameters

   if(QOption) then

     allocate( TmpOptCI(NOC) )
     allocate( TmpOptCJ(NOC) )
     allocate( TmpkOptC(NOC) )
     allocate( TmprOptC(NOC) )

     TmpOptCI = OptCI
     TmpOptCJ = OptCJ
     TmpkOptC = kOptC
     TmprOptC = rOptC

     TmpNOptC  = NumOptC / NProcsTemp
     NAddition = mod(NumOptC,NProcsTemp)

     if(MyRankTemp < NAddition) TmpNOptC = TmpNOptC + 1

     deallocate( OptCI, OptCJ, kOptC, rOptC )

     if(TmpNOptC /= 0) then

       allocate( OptCI(TmpNOptC) )
       allocate( OptCJ(TmpNOptC) )
       allocate( kOptC(TmpNOptC) )
       allocate( rOptC(TmpNOptC) )

       Count = 0

       do i = MyRankTemp+1, NumOptC, NProcsTemp

         Count = Count + 1
         OptCI(Count) = TmpOptCI(i)
         OptCJ(Count) = TmpOptCJ(i)
         kOptC(Count) = TmpkOptC(i)
         rOptC(Count) = TmprOptC(i)

       end do

       deallocate( TmpOptCI, TmpOptCJ, TmpkOptC, TmprOptC )

     end if

     NumOptC = TmpNOptC

   end if


end subroutine AllocPara


!######################################################################
!######################################################################


subroutine SumFrc( Acc )

use Numbers, only : N
use CommonBlocks, only : QMaster
use CommonMPI

implicit none

include 'mpif.h'

integer :: i , j
integer :: Nall
real(8), dimension(3,N) :: Acc

   Nall = 3 * N

   do i = 1, N

     j = (i-1) * 3
     Buff1(j+1) = Acc(1,i)
     Buff1(j+2) = Acc(2,i)
     Buff1(j+3) = Acc(3,i)

   end do

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(QMaster) then

     do i = 1, N

       j = (i-1) * 3
       Acc(1,i) = Buff2(j+1)
       Acc(2,i) = Buff2(j+2)
       Acc(3,i) = Buff2(j+3)

     end do

   end if

end subroutine SumFrc


!######################################################################
!######################################################################


subroutine SumEnergy(E_Bond, E_Angle, E_UB, E_Dihed, E_Impro, &
&                    E_LJ, E_Elec, E_Ersp, E_Eksp, E_OptC)

use CommonMPI

implicit none

include 'mpif.h'

integer, parameter :: Nall = 10
real(8) :: E_Bond, E_Angle, E_UB, E_Dihed, E_Impro
real(8) :: E_LJ, E_Elec, E_Ersp, E_Eksp, E_OptC

   Buff1( 1) = E_Bond
   Buff1( 2) = E_Angle
   Buff1( 3) = E_UB
   Buff1( 4) = E_Dihed
   Buff1( 5) = E_Impro
   Buff1( 6) = E_LJ
   Buff1( 7) = E_Elec
   Buff1( 8) = E_Ersp
   Buff1( 9) = E_Eksp
   Buff1(10) = E_OptC

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call mpi_reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     E_Bond  = Buff2( 1)
     E_Angle = Buff2( 2)
     E_UB    = Buff2( 3)
     E_Dihed = Buff2( 4)
     E_Impro = Buff2( 5)
     E_LJ    = Buff2( 6)
     E_Elec  = Buff2( 7)
     E_Ersp  = Buff2( 8)
     E_Eksp  = Buff2( 9)
     E_OptC  = Buff2(10)

   end if


end subroutine SumEnergy


!######################################################################
!######################################################################


subroutine SumVir( Vir )

use CommonMPI

implicit none

include 'mpif.h'

integer :: i , j, k
integer, parameter :: Nall = 9
real(8), dimension(3,3) :: Vir

   k = 0

   do i = 1 , 3

     do j = 1 , 3

       k = k + 1

       Buff1(k) = Vir(i,j)

     end do

   end do

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call mpi_reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     k = 0

     do i = 1 , 3

       do j = 1 , 3

         k = k + 1

         Vir(i,j) = Buff2(k)

       end do

     end do

   end if

end subroutine SumVir


!######################################################################
!######################################################################


subroutine BcastE(deltaE,Energy)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall
real(8) :: deltaE, Energy

   Nall = 2

   Buff1(1) = deltaE
   Buff1(2) = Energy

   call Mpi_Bcast( Buff1, Nall, Mpi_Double_Precision, 0, &
   &               Mpi_Comm_World, ierror )

   if(MyRank /= 0) then

     deltaE = Buff1(1)
     Energy = Buff1(2)

   end if

end subroutine BcastE


!######################################################################
!######################################################################


subroutine BcastR

use Numbers, only : N
use CommonBlocks, only : QMaster
use Configuration, only : R
use CommonMPI

implicit none

include 'mpif.h'

integer :: i, j, Nall

   Nall = 3 * N

   if(QMaster) then

     do i = 1 , N

       j = (i-1) * 3
       Buff1(j+1) = R(1,i)
       Buff1(j+2) = R(2,i)
       Buff1(j+3) = R(3,i)

     end do

   end if

   call Mpi_Bcast( Buff1, Nall, Mpi_Double_Precision, 0, &
   &               Mpi_Comm_World, ierror )

   if(.not.QMaster) then

     do i = 1 , N

       j = (i-1) * 3
       R(1,i) = Buff1(j+1)
       R(2,i) = Buff1(j+2)
       R(3,i) = Buff1(j+3)

     end do

   end if

end subroutine BcastR



!######################################################################
!######################################################################


subroutine BcastRH

use Numbers, only : N
use CommonBlocks, only : QMaster
use Configuration, only : R
use CellParam, only : H
use CommonMPI

implicit none

include 'mpif.h'

integer :: i, j, k, Nall

   Nall = 3 * N + 9

   if(QMaster) then

     do i = 1 , N

       j = (i-1) * 3
       Buff1(j+1) = R(1,i)
       Buff1(j+2) = R(2,i)
       Buff1(j+3) = R(3,i)

     end do

     k = 3 * N

     do i = 1 , 3

       do j = 1 , 3

         k = k + 1
         Buff1(k) = H(i,j)

       end do

     end do

   end if

   call Mpi_Bcast( Buff1, Nall, Mpi_Double_Precision, 0, &
   &               Mpi_Comm_World, ierror )

   if(.not.QMaster) then

     do i = 1 , N

       j = (i-1) * 3
       R(1,i) = Buff1(j+1)
       R(2,i) = Buff1(j+2)
       R(3,i) = Buff1(j+3)

     end do

     k = 3 * N

     do i = 1 , 3

       do j = 1 , 3

         k = k + 1
         H(i,j) = Buff1(k)

       end do

     end do

   end if

end subroutine BcastRH


!######################################################################
!######################################################################


subroutine BcastRGQuat

use CommonMPI
use RBparam, only : NumRB, R_RB, Quaternion

implicit none

include 'mpif.h'

integer :: i, j, Nall

   Nall = 7 * NumRB

   if(MyRank == 0) then

     do i = 1 , NumRB

       j = (i-1) * 7
       Buff1(j+1) = R_RB(1,i)
       Buff1(j+2) = R_RB(2,i)
       Buff1(j+3) = R_RB(3,i)
       Buff1(j+4) = Quaternion(1,i)
       Buff1(j+5) = Quaternion(2,i)
       Buff1(j+6) = Quaternion(3,i)
       Buff1(j+7) = Quaternion(4,i)

     end do

   end if

   call Mpi_Bcast( Buff1, Nall, Mpi_Double_Precision, 0, &
   &               Mpi_Comm_World, ierror )

   if(MyRank /= 0) then

     do i = 1 , NumRB

       j = (i-1) * 7

       R_RB      (1,i) = Buff1(j+1)
       R_RB      (2,i) = Buff1(j+2)
       R_RB      (3,i) = Buff1(j+3)
       Quaternion(1,i) = Buff1(j+4)
       Quaternion(2,i) = Buff1(j+5)
       Quaternion(3,i) = Buff1(j+6)
       Quaternion(4,i) = Buff1(j+7)

     end do

   end if

end subroutine BcastRGQuat


!######################################################################
!######################################################################


subroutine BcastRgQuatH

use CellParam, only : H
use CommonMPI
use RBparam, only : NumRB, R_RB, Quaternion

implicit none

include 'mpif.h'

integer :: i, j, k, Nall

   Nall = 7 * NumRB + 9

   if(MyRank == 0) then

     do i = 1 , NumRB

       j = (i-1) * 7
       Buff1(j+1) = R_RB(1,i)
       Buff1(j+2) = R_RB(2,i)
       Buff1(j+3) = R_RB(3,i)
       Buff1(j+4) = Quaternion(1,i)
       Buff1(j+5) = Quaternion(2,i)
       Buff1(j+6) = Quaternion(3,i)
       Buff1(j+7) = Quaternion(4,i)

     end do

     k = 7 * NumRB

     do i = 1 , 3

       do j = 1 , 3

         k = k + 1
         Buff1(k) = H(i,j)

       end do

     end do

   end if

   call Mpi_Bcast( Buff1, Nall, Mpi_Double_Precision, 0, &
   &               Mpi_Comm_World, ierror )

   if(MyRank /= 0) then

     do i = 1 , NumRB

       j = (i-1) * 7
       R_RB      (1,i) = Buff1(j+1)
       R_RB      (2,i) = Buff1(j+2)
       R_RB      (3,i) = Buff1(j+3)
       Quaternion(1,i) = Buff1(j+4)
       Quaternion(2,i) = Buff1(j+5)
       Quaternion(3,i) = Buff1(j+6)
       Quaternion(4,i) = Buff1(j+7)

     end do

     k = 7 * NumRB

     do i = 1 , 3

       do j = 1 , 3

         k = k + 1
         H(i,j) = Buff1(k)

       end do

     end do

   end if

end subroutine BcastRgQuatH


!######################################################################
!######################################################################


subroutine SumRDF( Acc )

use CommonMPI
use ParamAnalyze, only : NumGR, IRcut

implicit none

include 'mpif.h'

integer :: i, j, k
integer :: Nall
integer, dimension(IRcut,NumGR) :: Acc
integer, dimension(IRcut*NumGR) :: Ibuff1, Ibuff2

   Nall = IRcut * NumGR

   k = 0
   do j = 1 , NumGR
     do i = 1, IRcut

       k = k + 1
       Ibuff1(k) = Acc(i,j)

     end do
   end do

   call Mpi_Reduce(Ibuff1,Ibuff2,Nall,Mpi_Integer,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     k = 0
     do j = 1 , NumGR
       do i = 1, IRcut

         k = k + 1
         Acc(i,j) = Ibuff2(k)

       end do
     end do

   end if

end subroutine SumRDF


!######################################################################
!######################################################################


subroutine SumEneInsertion

use CommonMPI
use ParamAnalyze, only : NumSpecInsert, EneInsert, EneInsertPT

implicit none

include 'mpif.h'

integer :: i, j
integer :: Nall

   Nall = NumSpecInsert*2

   do i = 1, NumSpecInsert

     j = (i-1) * 2

     Buff1(j+1) = EneInsert(i)
     Buff1(j+2) = EneInsertPT(i)
!     Buff1(j+3) = EneInsertLJ(i)
!     Buff1(j+4) = EneInsertEl(i)

   end do

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     do i = 1, NumSpecInsert

       j = (i-1) * 2

       EneInsert(i)   = Buff2(j+1)
       EneInsertPT(i) = Buff2(j+2)
!       EneInsertLJ(i) = Buff2(j+3)
!       EneInsertEl(i) = Buff2(j+4)

     end do

   end if

end subroutine SumEneInsertion


!######################################################################
!######################################################################


#ifdef MEAM

subroutine SumRho( Acc )

use Numbers, only : N
use CommonMPI
use EAM_param, only : NumEAM

implicit none

include 'mpif.h'

integer :: i, k
integer :: Nall
real(8), dimension(N,NumEAM) :: Acc

   Nall = N

   do k = 1, NumEAM

      do i = 1, N

        Buff1(i) = Acc(i,k)

      end do

      call Mpi_AllReduce(Buff1,Buff2,Nall,Mpi_Double_Precision, &
      &                  Mpi_Sum,Mpi_Comm_World,ierror)

      do i = 1, N

        Acc(i,k) = Buff2(i)

      end do

   end do

end subroutine SumRho

#else

subroutine SumRho( Acc )

use Numbers, only : N
use CommonMPI

implicit none

include 'mpif.h'

integer :: i
integer :: Nall
real(8), dimension(N) :: Acc

   Nall = N

   do i = 1, N

     Buff1(i) = Acc(i)

   end do

   call Mpi_AllReduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &                  Mpi_Comm_World,ierror)

   do i = 1, N

     Acc(i) = Buff2(i)

   end do

end subroutine SumRho

#endif


!######################################################################
!######################################################################


subroutine SumDistPI( Nall, Vect )

use CommonMPI

implicit none

include 'mpif.h'

integer :: i
integer :: Nall
real(8), dimension(Nall) :: Vect

   do i = 1, Nall
     Buff1(i) = Vect(i)
   end do

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_AllReduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &                  Mpi_Comm_World,ierror)

   do i = 1, Nall
     Vect(i) = Buff2(i)
   end do

end subroutine SumDistPI


!######################################################################
!######################################################################


subroutine SumFrcPI( Vect )

use Numbers, only : N
use CommonPI
use CommonMPI

implicit none

include 'mpif.h'

integer :: i, j, k
integer :: Nall, Nrec
real(8), dimension(3,N,Nbead) :: Vect
integer, dimension(0:NProcs-1) :: Nsend, Displs

   Nall = 3 * N * Nbead

   do i = 0, NProcs-1
     Nsend(i) = Nassi(i+1) * 3 * N
   end do

   Displs(0) = 0

   do i = 1, NProcs-1
     Displs(i) = Displs(i-1) + Nsend(i-1)
   end do

   k = 0

   do j = 1, Nbead

     do i = 1, N

       Buff1(k+1) = Vect(1,i,j)
       Buff1(k+2) = Vect(2,i,j)
       Buff1(k+3) = Vect(3,i,j)

       k = k + 3

     end do

   end do

!   call Mpi_AllReduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
!   &                  Mpi_Comm_World,ierror)

!   do i = 1, Nall
!     Vect(i) = Buff2(i)
!   end do

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(QMasterPI) then

     Nrec = 3 * N * BeadNum

   else

     Nrec = 0

   end if

   call Mpi_Scatterv(Buff2,Nsend,Displs,Mpi_Double_Precision, &
   &                 Buff1,Nrec,        Mpi_Double_Precision, &
   &                 0,Mpi_Comm_World,ierror)

   if(QMasterPI) then

     k = 0

     do j = IniBead, FinBead

       do i = 1, N

         Vect(1,i,j) = Buff1(k+1)
         Vect(2,i,j) = Buff1(k+2)
         Vect(3,i,j) = Buff1(k+3)

         k = k + 3

       end do

     end do

   end if

end subroutine SumFrcPI


!######################################################################
!######################################################################


subroutine SumFrcPI_DCV( Vect )

use Numbers, only : N
use CommonBlocks, only : QMaster
use CommonMPI
use CommonPI

implicit none

include 'mpif.h'

integer :: i , j, k
integer :: Nall
real(8), dimension(3,N,Nbead) :: Vect

   Nall = 3 * N * Nbead

   k = 0

   do j = 1, Nbead

     do i = 1, N

       Buff1(k+1) = Vect(1,i,j)
       Buff1(k+2) = Vect(2,i,j)
       Buff1(k+3) = Vect(3,i,j)

       k = k + 3

     end do

   end do

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(QMaster) then

     k = 0

     do j = IniBead, FinBead

       do i = 1, N

         Vect(1,i,j) = Buff1(k+1)
         Vect(2,i,j) = Buff1(k+2)
         Vect(3,i,j) = Buff1(k+3)

         k = k + 3

       end do

     end do

   end if

end subroutine SumFrcPI_DCV


!######################################################################
!######################################################################


subroutine RnmPI

use Numbers, only : N
use CommonBlocks, only : QMaster
use CommonMPI
use CommonPI

implicit none

include 'mpif.h'

integer :: i, j, k, Nall
integer, dimension(0:NProcs-1) :: Nreciev, Displs

   do i = 0, NProcs-1
     Nreciev(i) = Nassi(i+1) * 3 * N
   end do

   Displs(0) = 0

   do i = 1, NProcs-1
     Displs(i) = Displs(i-1) + Nreciev(i-1)
   end do

   if(QMasterPI) then

     Nall = 3 * N * BeadNum

     k = 0

     do j = IniBead, FinBead

       do i = 1 , N

         Buff1(k+1) = Rnm(1,i,j)
         Buff1(k+2) = Rnm(2,i,j)
         Buff1(k+3) = Rnm(3,i,j)

         k = k + 3

       end do

     end do

   else

     Nall = 0
     Buff1(1) = 0.d0

   end if

   call Mpi_Gatherv( Buff1, Nall,            Mpi_Double_Precision, &
   &                 Buff2, Nreciev, Displs, Mpi_Double_Precision, &
   &                 0, Mpi_Comm_World, ierror )

   if(QMaster) then

     k = 0

     do j = 1, Nbead

       do i = 1 , N

         Rnm(1,i,j) = Buff2(k+1)
         Rnm(2,i,j) = Buff2(k+2)
         Rnm(3,i,j) = Buff2(k+3)

         k = k + 3

       end do

     end do

   end if

end subroutine RnmPI


!######################################################################
!######################################################################


subroutine VnmPI

use Numbers, only : N
use CommonBlocks, only : QMaster
use CommonMPI
use CommonPI

implicit none

include 'mpif.h'

integer :: i, j, k, Nall
integer, dimension(0:NProcs-1) :: Nreciev, Displs

   do i = 0, NProcs-1
     Nreciev(i) = Nassi(i+1) * 3 * N
   end do

   Displs(0) = 0

   do i = 1, NProcs-1
     Displs(i) = Displs(i-1) + Nreciev(i-1)
   end do

   if(QMasterPI) then

     Nall = 3 * N * BeadNum

     k = 0

     do j = IniBead, FinBead

       do i = 1 , N

         Buff1(k+1) = Vnm(1,i,j)
         Buff1(k+2) = Vnm(2,i,j)
         Buff1(k+3) = Vnm(3,i,j)

         k = k + 3

       end do

     end do

   else

     Nall = 0
     Buff1(1) = 0.d0

   end if

   call Mpi_Gatherv( Buff1, Nall,            Mpi_Double_Precision, &
   &                 Buff2, Nreciev, Displs, Mpi_Double_Precision, &
   &                 0, Mpi_Comm_World, ierror )

   if(QMaster) then

     k = 0

     do j = 1, Nbead

       do i = 1 , N

         Vnm(1,i,j) = Buff2(k+1)
         Vnm(2,i,j) = Buff2(k+2)
         Vnm(3,i,j) = Buff2(k+3)

         k = k + 3

       end do

     end do

   end if

end subroutine VnmPI


!######################################################################
!######################################################################


subroutine RVBathPI

use Numbers, only : N
use CommonBlocks, only : QMaster
use BathParam, only : NHchain
use CommonMPI
use CommonPI

implicit none

include 'mpif.h'

integer :: i, j, k, l, Nall
integer, dimension(0:NProcs-1) :: Nreciev, Displs

   Nreciev(0) = (Nassi(1) - 1) * 3 * N * NHchain

   do i = 1, NProcs-1
     Nreciev(i) = Nassi(i+1) * 3 * N * NHchain
   end do

   Displs(0) = 0

   do i = 1, NProcs-1
     Displs(i) = Displs(i-1) + Nreciev(i-1)
   end do

   if(QMasterPI) then

     Nall = 3 * N * BeadNum * NHchain

     if(QMaster) Nall = 3 * N * (BeadNum-1) * NHchain

     k = 0

     do j = IniBead, FinBead

       if(j==1) cycle

       do l = 1, NHchain

         do i = 1 , N

           Buff1(k+1) = Rbath(1,i,l,j)
           Buff1(k+2) = Rbath(2,i,l,j)
           Buff1(k+3) = Rbath(3,i,l,j)

           k = k + 3

         end do

       end do

     end do

   else

     Nall = 0
     Buff1(1) = 0.d0

   end if

   call Mpi_Gatherv( Buff1, Nall,            Mpi_Double_Precision, &
   &                 Buff2, Nreciev, Displs, Mpi_Double_Precision, &
   &                 0, Mpi_Comm_World, ierror )

   if(QMaster) then

     k = 0

     do j = 2, Nbead

       do l = 1 , NHchain

         do i = 1 , N

           Rbath(1,i,l,j) = Buff2(k+1)
           Rbath(2,i,l,j) = Buff2(k+2)
           Rbath(3,i,l,j) = Buff2(k+3)

           k = k + 3

         end do

       end do

     end do

   end if

   if(QMasterPI) then

     Nall = 3 * N * BeadNum * NHchain

     if(QMaster) Nall = 3 * N * (BeadNum-1) * NHchain

     k = 0

     do j = IniBead, FinBead

       if(j==1) cycle

       do l = 1, NHchain

         do i = 1 , N

           Buff1(k+1) = Vbath(1,i,l,j)
           Buff1(k+2) = Vbath(2,i,l,j)
           Buff1(k+3) = Vbath(3,i,l,j)

           k = k + 3

         end do

       end do

     end do

   else

     Nall = 0
     Buff1(1) = 0.d0

   end if

   call Mpi_Gatherv( Buff1, Nall,            Mpi_Double_Precision, &
   &                 Buff2, Nreciev, Displs, Mpi_Double_Precision, &
   &                 0, Mpi_Comm_World, ierror )

   if(QMaster) then

     k = 0

     do j = 2, Nbead

       do l = 1 , NHchain

         do i = 1 , N

           Vbath(1,i,l,j) = Buff2(k+1)
           Vbath(2,i,l,j) = Buff2(k+2)
           Vbath(3,i,l,j) = Buff2(k+3)

           k = k + 3

         end do

       end do

     end do

   end if

end subroutine RVBathPI


!######################################################################
!######################################################################


subroutine SumEnePI(E,Qkin,Eth)

use CommonMPI

implicit none

include 'mpif.h'

integer :: i
integer :: Nall
real(8) :: E,Qkin,Eth

   Nall = 3

   Buff1(1) = E
   Buff1(2) = Qkin
   Buff1(3) = Eth

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     E    = Buff2(1)
     Qkin = Buff2(2)
     Eth  = Buff2(3)

   end if

end subroutine SumEnePI


!######################################################################
!######################################################################


subroutine SumTempPI

use CommonBlocks, only : QMaster
use ThermoData, only : Pkinp
use CommonMPI
use CommonPI

implicit none

include 'mpif.h'

integer :: i
integer :: Nall

   Nall = 6

   Buff1(1) = Pkinp(1,1)
   Buff1(2) = Pkinp(1,2)
   Buff1(3) = Pkinp(1,3)
   Buff1(4) = Pkinp(2,2)
   Buff1(5) = Pkinp(2,3)
   Buff1(6) = Pkinp(3,3)

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(QMaster) then

     Pkinp(1,1) = Buff2(1)
     Pkinp(1,2) = Buff2(2)
     Pkinp(1,3) = Buff2(3)
     Pkinp(2,2) = Buff2(4)
     Pkinp(2,3) = Buff2(5)
     Pkinp(3,3) = Buff2(6)

   end if

end subroutine SumTempPI


!######################################################################
!######################################################################


! >> F monitor ##

subroutine SumFrcgA

use CommonMPI
use F_monitor, only : NgA, NiniF, NfinF, Fint

implicit none

include 'mpif.h'

integer :: i, j
integer :: Nall

   Nall = 3 * NgA

   j = 0

   do i = NiniF + 1, NfinF

     Buff1(j+1) = Fint(1,i)
     Buff1(j+2) = Fint(2,i)
     Buff1(j+3) = Fint(3,i)

     j = j + 3

   end do

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     j = 0

     do i = NiniF + 1, NfinF

       Fint(1,i) = Buff2(j+1)
       Fint(2,i) = Buff2(j+2)
       Fint(3,i) = Buff2(j+3)

       j = j + 3

     end do

   end if

end subroutine SumFrcgA
! << F monitor ##


!######################################################################
!######################################################################


subroutine Prep_Atom_to_Mesh

use CommonBlocks, only : QPathInt
use CommonMPI
use PMEparam, only : NfftDim, Nscnt, Ndisp, Nrenum
use CommonPI, only : MyRankPI, NumProcess

implicit none

include 'mpif.h'

integer :: i , j, k, l
integer :: nz, icpu
integer, dimension(NProcs) :: Numnz
integer :: NProcsTemp

   if(QPathInt) then
     NProcsTemp = NumProcess
   else
     NProcsTemp = NProcs
   end if

   allocate( Nscnt(0:NProcsTemp-1) )
   allocate( Ndisp(0:NProcsTemp-1) )

   allocate( Nrenum(NfftDim(3)) )

! ## 

   i = NfftDim(3) / NProcsTemp
   Numnz = i

   j = mod( NfftDim(3), NProcsTemp )
   if(j/=0) then
     do i = NProcsTemp-j+1, NProcsTemp
       NumNz(i) = NumNz(i) + 1
     end do
   end if

! ## 

   l = 0

   do i = 1, NProcsTemp

     do j = 1, NumNz(i)

       l = l + 1
       nz = 0

       do k = 1, NfftDim(3)

         icpu = NProcsTemp - mod(k-1,NProcsTemp)

         if(icpu == i) then
           nz = nz + 1
           if(nz == j) then
             Nrenum(l) = k
             exit
           end if
         end if

       end do

     end do

   end do

   do i = 1, NProcsTemp

     Nscnt(i-1) = Numnz(i) * NfftDim(2) * NfftDim(1)

   end do

   Ndisp(0) = 0

   do i = 1, NProcsTemp-1

     Ndisp(i) = Ndisp(i-1) + Nscnt(i-1)

   end do


end subroutine Prep_Atom_to_Mesh


!######################################################################
!######################################################################


subroutine SumChargeDens( Q )

use CommonBlocks, only : QPathInt
use CommonMPI
use PMEparam, only : NfftDim, Nscnt, Ndisp, Nrenum
use CommonPI, only : MyRankPI, NumProcess

implicit none

include 'mpif.h'

integer :: i , j, k, l, kk
integer :: Nall, Nas
real(8), dimension(NfftDim(1),NfftDim(2),NfftDim(3)) :: Q
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   Nall = NfftDim(1) * NfftDim(2) * NfftDim(3)

   l = 0

   do k = 1, NfftDim(3)

     kk = Nrenum(k)

     do j = 1, NfftDim(2)

       do i = 1, NfftDim(1)

         l = l + 1
         Buff1(l) = Q(i,j,kk)

       end do

     end do

   end do

!   call Mpi_Barrier(Mpi_Comm_World,ierror)

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0, Mpi_Comm_World,ierror)


   call Mpi_ScatterV(Buff2,Nscnt,Ndisp,Mpi_Double_Precision,  &
   &                 Buff1,Nscnt(MyRankTemp),Mpi_Double_Precision, &
   &                 0,Mpi_Comm_World,ierror)

   l = 0

   do k = Nas, NfftDim(3), NProcsTemp

     do j = 1, NfftDim(2)

       do i = 1, NfftDim(1)

         l = l + 1
         Q(i,j,k) = Buff1(l)

       end do

     end do

   end do

end subroutine SumChargeDens


!######################################################################
!######################################################################


subroutine DistChargeDens(Q)

use CommonBlocks, only : QPathInt
use CommonMPI
use PMEparam, only : NfftDim, Nscnt, Ndisp, Nrenum
use CommonPI, only : MyRankPI, NumProcess

implicit none

include 'mpif.h'

integer :: i , j, k, l, kk
integer :: Nas
real(8), dimension(NfftDim(1),NfftDim(2),NfftDim(3)) :: Q
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   l = 0

   do k = Nas, NfftDim(3), NProcsTemp

     do j = 1, NfftDim(2)

       do i = 1, NfftDim(1)

         l = l + 1
         Buff1(l) = Q(i,j,k)

       end do

     end do

   end do

   call Mpi_Allgatherv(Buff1,Nscnt(MyRankTemp),Mpi_Double_Precision, &
   &                   Buff2,Nscnt,Ndisp, Mpi_Double_Precision, &
   &                   Mpi_Comm_World, ierror)

   l = 0

   do k = 1, NfftDim(3)

     kk = Nrenum(k)

     do j = 1, NfftDim(2)

       do i = 1, NfftDim(1)

         l = l + 1
         Q(i,j,kk) = Buff2(l)

       end do

     end do

   end do

end subroutine DistChargeDens


!######################################################################
!######################################################################


subroutine SumEneTICR

use CommonMPI
use FEparam, only : E_TICR

implicit none

include 'mpif.h'

integer :: Nall

   Nall = 1

   Buff1(1) = E_TICR

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     E_TICR = Buff2(1)

   end if

end subroutine SumEneTICR


!######################################################################
!######################################################################


subroutine SumPMF_CGball(PMF,Nall)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall, i
real(8), dimension(Nall) :: PMF

   do i = 1, Nall
     Buff1(i) = PMF(i)
   end do

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then
     do i = 1, Nall
       PMF(i) = Buff2(i)
     end do
   end if

end subroutine SumPMF_CGball


!######################################################################
!######################################################################


subroutine SumPMF_Cyl(PMF)

use CommonMPI

implicit none

include 'mpif.h'

real(8) :: PMF

   Buff1(1) = PMF

   call Mpi_Reduce(Buff1,Buff2,1,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then
     PMF = Buff2(1)
   end if

end subroutine SumPMF_Cyl


!######################################################################
!######################################################################


subroutine SumHC(EnePI, CVtrace, Vtrace, HessCVTrace, HessTrace)

use CommonMPI

implicit none

include 'mpif.h'

integer, parameter :: Nall = 5
real(8) :: EnePI, CVtrace, Vtrace, HessCVTrace, HessTrace

   Buff1(1) = EnePI
   Buff1(2) = CVtrace
   Buff1(3) = Vtrace
   Buff1(4) = HessCVTrace
   Buff1(5) = HessTrace

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then

     EnePI       = Buff2(1)
     CVtrace     = Buff2(2)
     Vtrace      = Buff2(3)
     HessCVTrace = Buff2(4)
     HessTrace   = Buff2(5)

   end if

end subroutine SumHC


!######################################################################
!######################################################################


subroutine FFT_ChAxisF(w,n1,n2,n3,ld1,ld2)

use CommonBlocks, only : QPathInt
use PMEparam, only: MaxGrid
use CommonMPI, only: NProcs, MyRank, ierror
use CommonPI, only : MyRankPI, NumProcess

implicit none

include 'mpif.h'

integer :: n1, n2, n3, ld1, ld2
real(8), dimension(2,ld1,ld2,n3) :: w
integer, dimension(NProcs) :: NSdata, NRdata
real(8), dimension(MaxGrid,NProcs) :: Dsend, Drecv
real(8), dimension(MaxGrid) :: Vsend, Vrecv
integer :: i, j, k, ii, jj, kk, Nas
integer :: ireq
integer, dimension(Mpi_Status_Size) :: istatus
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

! ## Send buffer

   NSdata = 0

   do k = Nas, n3, NProcsTemp

     do i = 1, n1

       ii = NProcsTemp - mod(i-1,NProcsTemp) ! NProcs, NProcs-1, ..., 2, 1

       if(ii==(MyRankTemp+1)) cycle

       do j = 1, n2

         jj = NSdata(ii)
         Dsend( jj+1, ii ) = w(1,i,j,k)
         Dsend( jj+2, ii ) = w(2,i,j,k)
         NSdata(ii) = NSdata(ii) + 2

       end do

     end do

   end do

! ## Recieve Count

   NRdata = 0

   do i = Nas, n1, NProcsTemp

     do k = 1, n3

       kk = NProcsTemp - mod(k-1,NProcsTemp)

       if(kk==(MyRankTemp+1)) cycle

       NRdata(kk) = NRdata(kk) + n2 * 2

     end do

   end do

! ## 

   do i = 1, NProcsTemp

     ii = i - 1

     if(NSdata(i) /= 0) then

       do j = 1, NSdata(i)
         Vsend(j) = Dsend(j,i)
       end do

       call Mpi_Isend(Vsend,NSdata(i),Mpi_Double_Precision,ii,i, &
       &              Mpi_Comm_World,ireq,ierror)

     end if

     if(NRdata(i) /= 0) then

       call Mpi_Recv(Vrecv,NRdata(i),Mpi_Double_Precision,ii,Mpi_Any_Tag, &
       &              Mpi_Comm_World,istatus,ierror)

       do j = 1, NRdata(i)
         Drecv(j,i) = Vrecv(j)
       end do

     end if

     if( NSdata(i) /= 0) call Mpi_Wait(ireq,istatus,ierror)

   end do

   NRdata = 0

   do k = 1, n3

     kk = NProcsTemp - mod(k-1, NProcsTemp)

     if( kk == (MyRankTemp + 1) ) cycle

     do i = Nas, n1, NProcsTemp

       do j = 1, n2

         jj = NRdata(kk)
         w(1,i,j,k) = Drecv( jj+1, kk )
         w(2,i,j,k) = Drecv( jj+2, kk )
         NRdata(kk) = NRdata(kk) + 2

       end do

     end do

   end do


end subroutine FFT_ChAxisF


!######################################################################
!######################################################################


subroutine FFT_ChAxisB(w,n1,n2,n3,ld1,ld2)

use CommonBlocks, only : QPathInt
use PMEparam, only: MaxGrid
use CommonMPI, only: NProcs, MyRank, ierror
use CommonPI, only : MyRankPI, NumProcess

implicit none

include 'mpif.h'

integer :: n1, n2, n3, ld1, ld2
real(8), dimension(2,ld1,ld2,n3) :: w
integer, dimension(NProcs) :: NSdata, NRdata
real(8), dimension(MaxGrid,NProcs) :: Dsend, Drecv
real(8), dimension(MaxGrid) :: Vsend, Vrecv
integer :: i, j, k, ii, jj, kk, Nas
integer :: ireq
integer, dimension(Mpi_Status_Size) :: istatus
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

! ## Send buffer

   NSdata = 0

   do i = Nas, n1, NProcsTemp

     do k = 1, n3

       kk = NProcsTemp - mod(k-1,NProcsTemp) ! NProcs, NProcs-1, ..., 2, 1

       if( kk==(MyRankTemp+1) ) cycle

       do j = 1, n2

         jj = NSdata(kk)
         Dsend( jj+1, kk ) = w(1,i,j,k)
         Dsend( jj+2, kk ) = w(2,i,j,k)
         NSdata(kk) = NSdata(kk) + 2

       end do

     end do

   end do

! ## Recieve Count

   NRdata = 0

   do k = Nas, n3, NProcsTemp

     do i = 1, n1

       ii = NProcsTemp - mod(i-1,NProcsTemp)

       if( ii==(MyRankTemp+1) ) cycle

       NRdata(ii) = NRdata(ii) + n2 * 2

     end do

   end do

! ## 

   do i = 1, NProcsTemp

     ii = i - 1

     if(NSdata(i) /= 0) then

       do j = 1, NSdata(i)
         Vsend(j) = Dsend(j,i)
       end do

       call Mpi_Isend(Vsend,NSdata(i),Mpi_Double_Precision,ii,i, &
       &              Mpi_Comm_World,ireq,ierror)

     end if

     if(NRdata(i) /= 0) then

       call Mpi_Recv(Vrecv,NRdata(i),Mpi_Double_Precision,ii,Mpi_Any_Tag, &
       &              Mpi_Comm_World,istatus,ierror)

       do j = 1, NRdata(i)
         Drecv(j,i) = Vrecv(j)
       end do

     end if

     if( NSdata(i) /= 0) call Mpi_Wait(ireq,istatus,ierror)

   end do

! ## 

   NRdata = 0

   do i = 1, n1

     ii = NProcsTemp - mod(i-1, NProcsTemp)

     if( ii == (MyRankTemp + 1) ) cycle

     do k = Nas, n3, NProcsTemp

       do j = 1, n2

         jj = NRdata(ii)
         w(1,i,j,k) = Drecv( jj+1, ii )
         w(2,i,j,k) = Drecv( jj+2, ii )
         NRdata(ii) = NRdata(ii) + 2

       end do

     end do

   end do

end subroutine FFT_ChAxisB


!######################################################################
!######################################################################


subroutine SumList( nn, iarray )

use CommonMPI

implicit none

include 'mpif.h'

integer :: i
integer :: nn
integer, dimension(nn) :: iarray
integer, dimension(nn) :: Ibuff1, Ibuff2

   do i = 1, nn
     Ibuff1(i) = iarray(i)
   end do

   call Mpi_Reduce(Ibuff1,Ibuff2,nn,Mpi_Integer,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then
     do i = 1 , nn
       iarray(i) = Ibuff2(i)
     end do
   end if

end subroutine SumList


!######################################################################
!######################################################################


subroutine SumConstF( nn, array, array2 )

use CommonMPI

implicit none

include 'mpif.h'

integer :: i
integer :: nn, nn2
real(8), dimension(nn) :: array, array2

   nn2 = nn * 2

   do i = 1, nn
     Buff1(i) = array(i)
   end do
   do i = 1, nn
     Buff1(i+nn) = array2(i)
   end do

   call Mpi_Reduce(Buff1,Buff2,nn2,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then
     do i = 1 , nn
       array(i) = Buff2(i)
     end do
     do i = 1 , nn
       array2(i) = Buff2(i+nn)
     end do
   end if

end subroutine SumConstF


!######################################################################
!######################################################################


subroutine Sum_StressProf

use CommonMPI
use ParamAnalyze, only : Nbin, VirProXY,VirProZZ

implicit none

include 'mpif.h'

integer :: Nall, ii, i, j

   Nall = 2 * Nbin
   ii = 0
   do i = 1, Nbin
     Buff1(ii+1) = VirProXY(i)
     Buff1(ii+2) = VirProZZ(i)
     ii = ii + 2
   end do

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then
     ii = 0
     do i = 1, Nbin
       VirProXY(i) = Buff2(ii+1)
       VirProZZ(i) = Buff2(ii+2)
       ii = ii + 2
     end do
   end if

end subroutine Sum_StressProf


!######################################################################
!######################################################################


subroutine Sum_StressProf2

use CommonMPI
use ParamAnalyze, only : Nbin,VirProN,VirProT

implicit none

include 'mpif.h'

integer :: Nall, ii, i, j

   Nall = 2 * Nbin
   ii = 0
   do i = 1, Nbin
     Buff1(ii+1) = VirProN(i)
     Buff1(ii+2) = VirProT(i)
     ii = ii + 2
   end do

   call Mpi_Reduce(Buff1,Buff2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank == 0) then
     ii = 0
     do i = 1, Nbin
       VirProN(i) = Buff2(ii+1)
       VirProT(i) = Buff2(ii+2)
       ii = ii + 2
     end do
   end if

end subroutine Sum_StressProf2


!######################################################################
!######################################################################


subroutine SumEne(B1,Nall)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall, i
real(8), dimension(Nall) :: B1,B2

   call Mpi_Reduce(B1,B2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank==0) then
     do i = 1, Nall
       B1(i) = B2(i)
     end do
   end if

end subroutine SumEne

!######################################################################
!######################################################################


#ifdef EnergyRep
subroutine Sum_ER1(avslf,engnorm,engsmpl,eself,uvmax)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall, ii, i, j
integer :: uvmax
real(8) :: avslf,engnorm,engsmpl
real(8), dimension(uvmax+2) :: eself
real(8), dimension(:), allocatable :: B1,B2

   Nall = uvmax + 5

   allocate(B1(Nall),B2(Nall))

   B1(1) = avslf
   B1(2) = engnorm
   B1(3) = engsmpl
   do i = 1, uvmax+2
     B1(i+3)=eself(i)
   end do

   call Mpi_Reduce(B1,B2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank==0) then
     avslf   = B2(1)
     engnorm = B2(2)
     engsmpl = B2(3)
     do i = 1, uvmax+2
       eself(i)=B2(i+3)
     end do
   end if

   deallocate(B1,B2)

end subroutine Sum_ER1


!######################################################################
!######################################################################


subroutine Sum_ER2(B1,Nall)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall, i
real(8), dimension(Nall) :: B1,B2

   call Mpi_Reduce(B1,B2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               0,Mpi_Comm_World,ierror)

   if(MyRank==0) then
     do i = 1, Nall
       B1(i) = B2(i)
     end do
   end if

end subroutine Sum_ER2


!######################################################################
!######################################################################


subroutine Sum_ER3(B1,Nall)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall,i
real(8), dimension(Nall) :: B1, B2

   call Mpi_AllReduce(B1,B2,Nall,Mpi_Double_Precision,Mpi_Sum, &
   &               Mpi_Comm_World,ierror)

   do i = 1, Nall
     B1(i) = B2(i)
   end do

end subroutine Sum_ER3


!######################################################################
!######################################################################


subroutine Sum_ER4(IB1,Nall)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall,i
integer, dimension(Nall) :: IB1,IB2

   call Mpi_AllReduce(IB1,IB2,Nall,Mpi_Integer,Mpi_Sum, &
   &               Mpi_Comm_World,ierror)

   do i = 1, Nall
     IB1(i) = IB2(i)
   end do

end subroutine Sum_ER4


!######################################################################
!######################################################################


subroutine Bcast1(a)

use CommonMPI

implicit none

include 'mpif.h'

integer :: Nall
real(8) :: a

   Nall = 1

   Buff1(1) = a

   call Mpi_Bcast( Buff1, Nall, Mpi_Double_Precision, 0, &
   &               Mpi_Comm_World, ierror )

   if(MyRank /= 0) then

     a = Buff1(1)

   end if

end subroutine Bcast1

#endif
