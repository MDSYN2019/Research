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
! ## -- BcastRPI 
! ## -- BcastRH 
! ## -- BcastRgQuat 
! ## -- BcastRgQuatH 
! ## -- SumRDF 
! ## -- SumEneInsertion 
! ## -- SumRho 
! ## -- SumFrcPI
! ## -- SumDistPI
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
! ############################


!######################################################################
!######################################################################


subroutine SetupMPI

use CommonBlocks, only : QMaster
use CommonMPI, only : NProcs, MyRank

   QMaster = .True.
   NProcs  = 1
   MyRank  = 0

#ifdef SCINT
   call m64_init()
#endif

end subroutine SetupMPI


!######################################################################
!######################################################################


subroutine Synch(NN)

implicit none

integer :: NN

end subroutine Synch


!######################################################################
!######################################################################


subroutine FinMPI

implicit none

end subroutine FinMPI


!######################################################################
!######################################################################


subroutine AllocPara

implicit none

end subroutine AllocPara


!######################################################################
!######################################################################


subroutine SumFrc( Acc )

use Numbers, only : N

implicit none

real(8), dimension(3,N) :: Acc

end subroutine SumFrc


!######################################################################
!######################################################################


subroutine SumEnergy(E_Bond, E_Angle, E_UB, E_Dihed, E_Impro, &
&                    E_LJ, E_Elec, E_Ersp, E_Eksp, E_OptC)

implicit none

real(8) :: E_Bond, E_Angle, E_UB, E_Dihed, E_Impro
real(8) :: E_LJ, E_Elec, E_Ersp, E_Eksp, E_OptC

end subroutine SumEnergy


!######################################################################
!######################################################################


subroutine SumVir( Vir )

implicit none

real(8), dimension(9) :: Vir

end subroutine SumVir


!######################################################################
!######################################################################


subroutine BcastE(deltaE,Energy)

implicit none

real(8) :: deltaE, Energy

end subroutine BcastE


!######################################################################
!######################################################################


subroutine BcastR

implicit none

end subroutine BcastR


!######################################################################
!######################################################################


subroutine BcastRPI

implicit none

end subroutine BcastRPI


!######################################################################
!######################################################################


subroutine BcastRH

implicit none

end subroutine BcastRH


!######################################################################
!######################################################################


subroutine BcastRgQuat

implicit none

end subroutine BcastRgQuat


!######################################################################
!######################################################################


subroutine BcastRgQuatH

implicit none

end subroutine BcastRgQuatH


!######################################################################
!######################################################################


subroutine SumRDF( Acc )

use ParamAnalyze

implicit none

integer, dimension(IRcut,NumGR) :: Acc

end subroutine SumRDF


!######################################################################
!######################################################################


subroutine SumEneInsertion

implicit none

end subroutine SumEneInsertion


!######################################################################
!######################################################################


subroutine SumRho( Acc )

use Numbers, only : N

implicit none

real(8), dimension(N) :: Acc

end subroutine SumRho


!######################################################################
!######################################################################


subroutine SumDistPI( Nall, Vect )

implicit none

integer :: Nall
real(8), dimension(Nall) :: Vect

end subroutine SumDistPI


!######################################################################
!######################################################################


subroutine SumFrcPI( Vect )

use Numbers, only : N
use CommonPI, only : Nbead

implicit none

real(8), dimension(3,N,Nbead) :: Vect

end subroutine SumFrcPI


!######################################################################
!######################################################################


subroutine SumFrcPI_DCV( Vect )

use Numbers, only : N
use CommonPI, only : Nbead

implicit none

real(8), dimension(3,N,Nbead) :: Vect

end subroutine SumFrcPI_DCV

!######################################################################
!######################################################################


subroutine RnmPI

end subroutine RnmPI


!######################################################################
!######################################################################


subroutine VnmPI

end subroutine VnmPI


!######################################################################
!######################################################################


subroutine RVBathPI

end subroutine RVBathPI


!######################################################################
!######################################################################


subroutine SumEnePI(dum1,dum2,dum3)

real(8) :: dum1,dum2,dum3

end subroutine SumEnePI


!######################################################################
!######################################################################


subroutine SumTempPI

end subroutine SumTempPI


!######################################################################
!######################################################################


! >> F monitor ##

subroutine SumFrcgA

end subroutine SumFrcgA
! << F monitor ##


!######################################################################
!######################################################################


subroutine Prep_Atom_to_Mesh

end subroutine Prep_Atom_to_Mesh


!######################################################################
!######################################################################


subroutine SumChargeDens( Q )

use PMEparam, only : NfftDim

real(8), dimension(NfftDim(1),NfftDim(2),NfftDim(3)) :: Q

end subroutine SumChargeDens


!######################################################################
!######################################################################


subroutine DistChargeDens( Q )

use PMEparam, only : NfftDim

real(8), dimension(NfftDim(1),NfftDim(2),NfftDim(3)) :: Q

end subroutine DistChargeDens


!######################################################################
!######################################################################


subroutine SumEneTICR

end subroutine SumEneTICR


!######################################################################
!######################################################################


subroutine SumPMF_CGball(PMF,Nall)

implicit none
integer :: Nall
real(8), dimension(Nall) :: PMF

end subroutine SumPMF_CGball

!######################################################################
!######################################################################


subroutine SumPMF_Cyl(PMF)

implicit none
real(8) :: PMF

end subroutine SumPMF_Cyl


!######################################################################
!######################################################################


subroutine SumHC(A, B, C, D, E)

implicit none

real(8) :: A, B, C, D, E

end subroutine SumHC


!######################################################################
!######################################################################


subroutine FFT_ChAxisF(w,n1,n2,n3,ld1,ld2)

implicit none

integer :: n1, n2, n3, ld1, ld2
complex(8), dimension(ld1,ld2,n3) :: w

end subroutine FFT_ChAxisF


!######################################################################
!######################################################################


subroutine FFT_ChAxisB(w,n1,n2,n3,ld1,ld2)

implicit none

integer :: n1, n2, n3, ld1, ld2
complex(8), dimension(ld1,ld2,n3) :: w

end subroutine FFT_ChAxisB


!######################################################################
!######################################################################


subroutine SumList( nn, iarray )

integer :: nn
integer, dimension(nn) :: iarray

end subroutine SumList


!######################################################################
!######################################################################


subroutine SumConstF( nn, array, array2 )

implicit none

integer :: nn
real(8), dimension(nn) :: array, array2

end subroutine SumConstF


!######################################################################
!######################################################################


subroutine Sum_StressProf

end subroutine Sum_StressProf


!######################################################################
!######################################################################


subroutine Sum_StressProf2

end subroutine Sum_StressProf2


!######################################################################
!######################################################################


subroutine SumEne(B1,Nall)

implicit none

integer :: Nall
real(8), dimension(Nall) :: B1

end subroutine SumEne


!######################################################################
!######################################################################


#ifdef EnergyRep
subroutine Sum_ER1(avslf,engnorm,engsmpl,eself,uvmax)

implicit none

integer :: uvmax
real(8) :: avslf,engnorm,engsmpl
real(8), dimension(uvmax+2) :: eself

end subroutine Sum_ER1


!######################################################################
!######################################################################


subroutine Sum_ER2(B1,Nall)

implicit none

integer :: Nall
real(8), dimension(Nall) :: B1

end subroutine Sum_ER2


!######################################################################
!######################################################################


subroutine Sum_ER3(B1,Nall)

implicit none

integer :: Nall
real(8), dimension(Nall) :: B1

end subroutine Sum_ER3


!######################################################################
!######################################################################


subroutine Sum_ER4(IB1,Nall)

implicit none

integer :: Nall
integer, dimension(numslv*(uvmax+2)) :: IB1

end subroutine Sum_ER4

!######################################################################
!######################################################################


subroutine Bcast1(a)

implicit none

real(8) :: a

end subroutine Bcast1

#endif
