! ############################
! ## SUBROUTINE LIST 
! ## -- SetupPI 
! ## -- NormalMode 
! ## -- SetMass
! ## -- PIVelocity
! ## -- GetNormalModeForce
! ## -- Assign_Process
! ## -- NormalModeTransform
! ## -- GetForce_Ref
! ############################


!######################################################################
!######################################################################


subroutine SetupPI

implicit none

! ## prepare transformation matrix
   call NormalMode

! ## prepare masses
   call SetMass

end subroutine SetupPI


!#####################################################################
!#####################################################################


! ***********************************************************
! **  Making the matrix of the normal mode transformation  **
! **  tnm: r -> u                                          **
! ***********************************************************

subroutine NormalMode

use CommonPI
use UnitExParam, only : pi2

implicit none

real(8), dimension(Nbead,Nbead) :: U, InvU
real(8) :: Psq, InvPsq, Norm, dummy
integer :: i, j

   Psq    = sqrt(Pbead)
   InvPsq = 1.d0 / sqrt(Pbead)
   Norm   = sqrt(2.d0*InvP)

!     /*  making unitary matrix for diagnalizing the spring matrix *
!      *  using analytical expressions                             */

   dummy = -1.d0

   do i = 1, Nbead

     U(i,1)     = InvPsq
     U(i,Nbead) = dummy * InvPsq
     dummy      = dummy * (-1.d0)

   end do

   do i = 1, (Nbead-2)/2

     do j = 1, Nbead

       U(j,2*i)   = Norm * cos( pi2 * i * j * InvP )
       U(j,2*i+1) = Norm * sin( pi2 * i * j * InvP )

     end do

   end do

   do i = 1, Nbead

     do j = 1, Nbead

       InvU(j,i) = U(i,j)

     end do

   end do

   Tnm    = Psq    * U
   InvTnm = InvPsq * InvU

end subroutine NormalMode


!#####################################################################
!#####################################################################


! *************************
! **  Initialize masses  **
! *************************

subroutine SetMass

use Numbers, only : N
use CommonPI
use UnitExParam, only : pi2
use BathParam, only : kT
use AtomParam, only : Mass

implicit none

integer :: i, j

   allocate( Qmass(Nbead) )
   allocate( InvQmass(Nbead) )

! bath parameters for path integral MD
   do i = 2, Nbead
     Qmass(i) = kT / OmegaP2
   end do
! For centroid MD, Qmass should be scaled by Gamma2, since
! natural frequencies of modes are Omega_p**2/Gamma**2
   do i = 2, Nbead
     Qmass(i) = GammaC2 * Qmass(i)
   end do

   do i = 2, Nbead
     InvQmass(i) = 1.d0 / Qmass(i)
   end do

   do j = 1, N

     NmMass(j,1)     = 0.d0
     NmMass(j,Nbead) = 4.d0 * Pbead * Mass(j)

     do i = 1, (Nbead-2)/2

       NmMass(j,2*i)   = 2.d0 * (1.d0 - cos( pi2 * i * InvP )) * Pbead * Mass(j)
       NmMass(j,2*i+1) = NmMass(j,2*i)

     end do

   end do

! fictitious mass for centroid MD

   do j = 1, N

     FictMass(j,1) = Mass(j)

     do i = 2, Nbead

       FictMass(j,i) = GammaC2 * NmMass(j,i)

     end do

   end do

   InvFictMass = 1.d0 / FictMass

end subroutine SetMass


!#####################################################################
!#####################################################################


! *******************************************************************
! ** Translational velocities from maxwell-boltzmann distribution. **
! ** The distribution is determined by Temperature and mass.       **
! *******************************************************************


subroutine PIVelocity

use Numbers, only : N
use CommonPI
use BathParam, only : kT

implicit none

integer :: i, j
real(8), dimension(3) :: SumV
real(8) :: vv, Gauss
external Gauss

   do j = 1, Nbead

!     /*  vsigma: standard devation of Maxwell distribution  */

     do i = 1, N

       vv = sqrt( kT / FictMass(i,j) )

       Vnm(1,i,j) = vv * Gauss()
       Vnm(2,i,j) = vv * Gauss()
       Vnm(3,i,j) = vv * Gauss()

     end do

   end do

!     /*  remove net momentum  */

   SumV = 0.d0

   do i = 1, N

     SumV = SumV + Vnm(:,i,1)

   end do

   SumV = SumV / dble(N)

   do i = 1, N

     Vnm(:,i,1) = Vnm(:,i,1) - SumV

   end do

end subroutine PIVelocity


!#####################################################################
!#####################################################################


! *******************************************************************
! ** transform force in real space into that in normal mode space  **
! *******************************************************************

subroutine GetNormalModeForce(FF,FFnm)

use Numbers, only : N
use CommonPI

implicit none

integer :: atom, i, j
real(8), dimension(3,N,Nbead) :: FF, FFnm
real(8) :: xx

!     /*  initialize array  */
   FFnm = 0.d0

!     /*  transformation                      *
!      *  fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)  */

   do i = 1, Nbead

     do j = IniBead, FinBead

       xx = Tnm(j,i)

       do atom = 1, N

         FFnm(:,atom,i) = FFnm(:,atom,i) + FF(:,atom,j) * xx

       end do

     end do

   end do

end subroutine GetNormalModeForce


!#####################################################################
!#####################################################################


! *********************************************************
! * For parallel path integral MD calculation             *
! * 1. when the number of processors is less than Nbead,  *
! * interactions associated with more than two bead       *
! * elements should be calculated in a single node        *
! * In this case, all nodes behave as if they are master  *
! * 2. when the number of processors is larger than Nbead,*
! * parallel force calculation between bead elements is   *
! * carried out.                                          *
! *********************************************************


subroutine Assign_Process

use CommonPI
use CommonMPI

implicit none

integer :: i, ii, jj, kk, ia, ib
integer, dimension(NProcs) :: Numb

   allocate( Nassi(NProcs) )

   if(NProcs <= Nbead) then

     Numb(:) = Nbead / NProcs  !  always >= 1 
     ii = mod(Nbead,NProcs)

     allocate(BeadOrder(Nbead))
     QBead = .False.

     if(ii /= 0) then
       do i = 1, ii
         Numb(i) = Numb(i) + 1
       end do
     end if

     ia = 0
     do i = 1, NProcs
       ib = ia + Numb(i)
       if(i == (MyRank+1)) then
         IniBead = ia + 1
         FinBead = ib
         BeadNum = Numb(i)
       end if
       ia = ib
     end do

     ii = 0
     do i = IniBead, FinBead
       ii = ii + 1
       BeadOrder(ii) = i
     end do

     NumProcess = 1
     MyRankPI   = 0
     QMasterPI  = .True.

     Nassi(:) = Numb(:)

   else

     allocate(BeadOrder(1))
     QBead = .True.
     BeadOrder(1) = mod( MyRank, Nbead ) + 1
     BeadNum = 1

     ii = NProcs / Nbead
     jj = ii * Nbead
     kk = NProcs - jj
     if( mod(MyRank,Nbead) < kk) ii = ii + 1
     NumProcess = ii
     MyRankPI   = MyRank/Nbead

     if(MyRankPI == 0) then
       QMasterPI = .True.
     else
       QMasterPI = .False.
     end if

     IniBead = BeadOrder(1)
     FinBead = BeadOrder(1)

     do i = 1, Nbead
       Nassi(i) = 1
     end do
     do i = Nbead + 1, NProcs
       Nassi(i) = 0
     end do

   end if

end subroutine Assign_Process


!#####################################################################
!#####################################################################


subroutine NormalModeTransform

use Numbers, only : N
use CommonPI

implicit none

integer :: atom, i, j, Nall !, ii

   Nall = 3 * N * Nbead

!   if (ii == 0) then
!     /*  from normal mode variables to real variables  *
!      *  x(i) = x(i) + sum_j tnm(i,j)*u(j)              */
!     /*  initialize array  */
   Rpi = 0.d0

   if(QMasterPI) then

     do atom = 1, N

       do i = 1, Nbead

         do j = IniBead, FinBead

           Rpi(:,atom,i) = Rpi(:,atom,i) + Tnm(i,j) * Rnm(:,atom,j)

         end do

       end do

     end do

   end if

   call SumDistPI(Nall, Rpi)

!   else if (ii == 1) then
!     /*  from real variables to normal mode variables  *
!      *  u(i) = u(i) + sum_j tnminv(i,j)*x(j)          */
!     /*  initialize array  */

!     Rnm = 0.d0

!     do atom = 1, N

!       do i = 1, Nbead
!         do j = 1, Nbead

!           Rnm(:,atom,i) = Rnm(:,atom,i) + InvTnm(i,j) * Rpi(:,atom,j)

!         end do
!       end do

!     end do

!   endif

end subroutine NormalModeTransform


!#####################################################################
!#####################################################################


! ************************************
! **  get reference harmonic force  **
! ************************************

subroutine GetForce_Ref

use Numbers, only : N
use CommonPI, only : Fnm_ref, Nbead, NmMass, OmegaP2, Rnm, IniBead, FinBead

implicit none

integer :: atom, mode

!     /*  centroid force in reference system is zero  */

   Fnm_ref(:,:,1) = 0.d0

   do mode = IniBead, FinBead

     if(mode == 1) cycle

     do atom = 1, N

       Fnm_ref(:,atom,mode) = - NmMass(atom,mode) * &
       &                        OmegaP2 * Rnm(:,atom,mode)

     end do

   end do

end subroutine GetForce_Ref
