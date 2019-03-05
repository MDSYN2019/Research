! ############################
! ## SUBROUTINE LIST 
! ## -- HeatCapacity 
! ## -- Calc_Hessian 
! ############################


!######################################################################
!######################################################################


! ********************************************
! * Heat capacity estimator                  *
! * Ref. K. R. Glaesemann and L. E. Fried,   *
! *      J. Chem. Phys. 117, 3020 (2002).    *
! ********************************************

subroutine HeatCapacity(EnePI, CVtrace, Vtrace, HessCVTrace, HessTrace, istep, ioutput)

use CommonBlocks, only: SimMethod, QPBC, ForceField, &
&                       QInitial, QCorrectCutoff
use ThermoData, only : Ene_kin, Virial
use Numbers, only : N
use TimeParam, only : Nstep, isampleHC
use AtomParam, only : Mass
use TailCorrect, only: Ene_LJ_co
use EwaldParam, only : Ene_Eslf
use BathParam, only : Beta, kT
use CommonPI, only: Nbead, Rpi, OmegaP2
use HeatCap
use UnitExParam, only : cvol

implicit none

! input
real(8) :: EnePI, CVtrace, Vtrace, HessCVTrace, HessTrace
integer :: istep, ioutput

! output

real(8) :: const
real(8), dimension(3) :: dR
real(8) :: riacc
integer :: iacc, igrid, i, j, k
integer, save :: irestart
real(8) :: K_T, K_V, K_CV, Vbar, E_T, E_V, E_CV, dE_T_dBeta
real(8) :: dE_V_dBeta, dE_CV_dBeta, B_V_E_V, B_CV_E_CV
real(8) :: ROH, RHH
real(8) :: E_T_E_T, E_T_E_V, E_T_E_CV, E_V_E_V, E_CV_E_CV
real(8) :: Cv_T, Cv_V, Cv_CV, Cv_DV, Cv_DCV

character(len=80) :: Filename_E, Filename_C, Filename_CDiv
character(len=10) :: RestartNum

   if(istep == 1) then

     if(QInitial) then

       irestart = - 5 * ioutput

       acc_K_T         = 0.d0
       acc_K_V         = 0.d0
       acc_K_CV        = 0.d0
       acc_Vbar        = 0.d0
       acc_E_T         = 0.d0
       acc_E_V         = 0.d0
       acc_E_CV        = 0.d0
       acc_E_T_E_T     = 0.d0
       acc_E_T_E_V     = 0.d0
       acc_E_T_E_CV    = 0.d0
       acc_E_V_E_V     = 0.d0
       acc_E_CV_E_CV   = 0.d0
       acc_dE_T_dBeta  = 0.d0
       acc_B_V_E_V     = 0.d0
       acc_B_CV_E_CV   = 0.d0

       if(.not.QPBC) then
         do igrid = 1, 1000
           acc_ROH(igrid) = 0.d0
           acc_RHH(igrid) = 0.d0
         end do
       end if

       open ( 81, file = 'Ene0.data', status = 'unknown' )
       open ( 82, file = 'Cv0.data', status = 'unknown' )
       open ( 85, file = 'CvDiv0.data', status = 'unknown' )

     else

       open(84,file='HCestimater.restart',form='unformatted',status='unknown')

       read(84) irestart
       read(84) acc_K_T
       read(84) acc_K_V
       read(84) acc_K_CV
       read(84) acc_Vbar
       read(84) acc_E_T
       read(84) acc_E_V
       read(84) acc_E_CV
       read(84) acc_E_T_E_T
       read(84) acc_E_T_E_V
       read(84) acc_E_T_E_CV
       read(84) acc_E_V_E_V
       read(84) acc_E_CV_E_CV
       read(84) acc_dE_T_dbeta
       read(84) acc_B_V_E_V
       read(84) acc_B_CV_E_CV
       if(.not.QPBC) then
         read(84) acc_ROH
         read(84) acc_RHH
       end if

       close(84)

       write(RestartNum,'(i10)') irestart

       write(Filename_E,'(a,a,a)') 'Ene',trim(adjustl(RestartNum)),'.data'
       write(Filename_C,'(a,a,a)') 'Cv',trim(adjustl(RestartNum)),'.data'
       write(Filename_CDiv,'(a,a,a)') 'CvDiv',trim(adjustl(RestartNum)),'.data'

       open ( 81, file = trim(adjustl(Filename_E)), status = 'unknown' )
       open ( 82, file = trim(adjustl(Filename_C)), status = 'unknown' )
       open ( 85, file = trim(adjustl(Filename_CDiv)), status = 'unknown' )

     end if

   end if

   iacc = istep + irestart

   if (iacc <= 0 ) return

   if(SimMethod == 'PIMD') then

! ## Estimate K_T 

     K_T = 0.5d0 * kT * dble(3 * N * Nbead)

     do i = 1, N

       const = - 0.5d0 * OmegaP2 * Mass(i)

       do j = 1, Nbead-1
         dR(:) = Rpi(:,i,j) - Rpi(:,i,j+1)
#ifdef PCC
         K_T = K_T + const * ( dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) )
#else
         K_T = K_T + const * dot_product( dR, dR )
#endif
       end do

       dR(:) = Rpi(:,i,Nbead) - Rpi(:,i,1)
#ifdef PCC
       K_T = K_T + const * ( dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) )
#else
       K_T = K_T + const * dot_product( dR, dR )
#endif

     end do

! ## Estimate K_V 

     if(QPBC) then

       K_V = - 0.5d0 * ( Virial(1,1) + Virial(2,2) + Virial(3,3) )

     else

       K_V = 0.5d0 * Vtrace

     end if

! ## Estimate K_CV 

     K_CV = 3.d0 * N * kT * 0.5d0 + CVtrace * 0.5d0

! ## Potential energy 

     if(QPBC) then

       if(QCorrectCutoff) then

         Vbar = EnePI + Ene_Eslf + Ene_LJ_co

#ifdef MEAM
       else if((ForceField(1:3) == 'EAM').or.(ForceField(1:4) == 'MEAM')) then
#else
       else if(ForceField(1:3) == 'EAM') then
#endif

         Vbar = EnePI

       else

         Vbar = EnePI + Ene_Eslf

       end if

     else

       Vbar = EnePI

     end if

! ## Energy

     E_T  = K_T  + Vbar
     E_V  = K_V  + Vbar
     E_CV = K_CV + Vbar

! ## d(E_T)/d(beta) 

     dE_T_dBeta = - 3.d0 * N * Nbead * kT * kT * 0.5d0

     do i = 1, N

       const = Mass(i) * kT * OmegaP2

       do j = 1, Nbead-1

         dR(:) = Rpi(:,i,j) - Rpi(:,i,j+1)
#ifdef PCC
         dE_T_dBeta = dE_T_dBeta + const * ( dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) )
#else
         dE_T_dBeta = dE_T_dBeta + const * dot_product( dR, dR )
#endif

       end do

       dR(:) = Rpi(:,i,Nbead) - Rpi(:,i,1)
#ifdef PCC
       dE_T_dBeta = dE_T_dBeta + const * ( dR(1)*dR(1) + dR(2)*dR(2) + dR(3)*dR(3) )
#else
       dE_T_dBeta = dE_T_dBeta + const * dot_product( dR, dR )
#endif

     end do

! ## d(E_V)/d(beta) 

     dE_V_dBeta = 0.d0

! ## d(E_CV)/d(beta) 

     dE_CV_dBeta = - 3.d0 * N * kT * kT * 0.5d0

! ## B_V*E_V 

     B_V_E_V = 0.d0

     if(.not.QPBC) then

       B_V_E_V = B_V_E_V + 1.5d0 * Vtrace + HessTrace

     end if

! ## B_CV *E_CV

     B_CV_E_CV = 1.5d0 * CVtrace + HessCVtrace

! ## accumulative
     acc_K_T  = acc_K_T  + K_T
     acc_K_V  = acc_K_V  + K_V
     acc_K_CV = acc_K_CV + K_CV

     acc_Vbar  = acc_Vbar  + Vbar

     acc_E_T  = acc_E_T  + E_T
     acc_E_V  = acc_E_V  + E_V
     acc_E_CV = acc_E_CV + E_CV

     acc_E_T_E_T   = acc_E_T_E_T   + E_T*E_T
     acc_E_T_E_V   = acc_E_T_E_V   + E_T*E_V
     acc_E_T_E_CV  = acc_E_T_E_CV  + E_T*E_CV
     acc_E_V_E_V   = acc_E_V_E_V   + E_V*E_V
     acc_E_CV_E_CV = acc_E_CV_E_CV + E_CV*E_CV

     acc_dE_T_dBeta  = acc_dE_T_dBeta  + dE_T_dBeta

     acc_B_V_E_V   = acc_B_V_E_V   + B_V_E_V
     acc_B_CV_E_CV = acc_B_CV_E_CV + B_CV_E_CV

     if(.not.QPBC) then

       do k = 1, Nbead

         dR(:) = Rpi(:,1,k) - Rpi(:,2,k)
         ROH = sqrt( dot_product( dR, dR ) )
         igrid = ROH * 100.d0
         acc_ROH(igrid) = acc_ROH(igrid) + 1.d0

         dR(:) = Rpi(:,1,k) - Rpi(:,3,k)
         ROH = sqrt( dot_product( dR, dR ) )
         igrid = ROH * 100.d0
         acc_ROH(igrid) = acc_ROH(igrid) + 1.d0

         dR(:) = Rpi(:,2,k) - Rpi(:,3,k)
         RHH = sqrt( dot_product( dR, dR ) )
         igrid = RHH * 100.d0
         acc_RHH(igrid) = acc_RHH(igrid) + 1.d0

       end do

     end if

     if ( mod(iacc,ioutput) == 0 ) then

       riacc = 1.d0 / dble(iacc)

       K_T  = acc_K_T *riacc
       K_V  = acc_K_V *riacc
       K_CV = acc_K_CV*riacc

       Vbar = acc_Vbar*riacc

       E_T  = acc_E_T *riacc
       E_V  = acc_E_V *riacc
       E_CV = acc_E_CV*riacc

       E_T_E_T   = acc_E_T_E_T  *riacc
       E_T_E_V   = acc_E_T_E_V  *riacc
       E_T_E_CV  = acc_E_T_E_CV *riacc
       E_V_E_V   = acc_E_V_E_V  *riacc
       E_CV_E_CV = acc_E_CV_E_CV*riacc

       dE_T_dBeta  = acc_dE_T_dBeta *riacc

       B_V_E_V   = acc_B_V_E_V  *riacc
       B_CV_E_CV = acc_B_CV_E_CV*riacc

       Cv_T =  (E_T_E_T  - E_T*E_T   - dE_T_dBeta )*Beta*Beta
       Cv_V =  (E_T_E_V  - E_V*E_V   - dE_V_dBeta )*Beta*Beta
       Cv_CV = (E_T_E_CV - E_CV*E_CV - dE_CV_dBeta)*Beta*Beta

       Cv_DV  = ( E_V_E_V   - E_V*E_V   - dE_V_dBeta &
       &        - 0.5d0 * kT * B_V_E_V   )*Beta*Beta
       Cv_DCV = ( E_CV_E_CV - E_CV*E_CV - dE_CV_dBeta &
       &        - 0.5d0 * kT * B_CV_E_CV )*Beta*Beta

       write (81,'(i8,4d13.5)') iacc, E_T*cvol, E_V*cvol, E_CV*cvol, Vbar*cvol
       if(QPBC) then
         write (82,'(i8,5d13.5)') iacc, Cv_T, Cv_V, Cv_CV, Cv_DCV
       else
         write (82,'(i8,5d13.5)') iacc, Cv_T, Cv_V, Cv_CV, Cv_DV, Cv_DCV
       end if

       write (85,'(i8,10d13.5)') iacc, E_T_E_T*Beta*Beta, -E_T*E_T*Beta*Beta, &
       &   -dE_T_dBeta*Beta*Beta, E_T_E_V*Beta*Beta, -E_V*E_V*Beta*Beta,      &
       &   E_T_E_CV*Beta*Beta, -E_CV*E_CV*Beta*Beta, E_V_E_V*Beta*Beta,       &
       &   - 0.5d0 * kT * B_V_E_V *Beta*Beta, - 0.5d0 * kT * B_CV_E_CV *Beta*Beta

     end if

   else if(SimMethod == 'MD') then

     if(QPBC) then
       Vbar = EnePI + Ene_Eslf + Ene_LJ_co
     else
       Vbar = EnePI
     end if

     call CalcTemp

     K_T = Ene_kin * 0.5d0

     E_T = K_T + Vbar

     acc_K_T  = acc_K_T  + K_T
     acc_Vbar = acc_Vbar + Vbar

     acc_E_T     = acc_E_T     + E_T
     acc_E_T_E_T = acc_E_T_E_T + E_T*E_T

     if ( mod(iacc,ioutput) == 0 ) then

       riacc = 1.d0 / dble(iacc)

       K_T  = acc_K_T  * riacc
       Vbar = acc_Vbar * riacc

       E_T     = acc_E_T     * riacc
       E_T_E_T = acc_E_T_E_T * riacc

       Cv_T =  (E_T_E_T  - E_T*E_T ) * Beta * Beta

       write (81,'(i8,3d13.5)') iacc, E_T*cvol, K_T*cvol, Vbar*cvol
       write (82,'(i8,d13.5)') iacc, Cv_T

     end if

   end if

   if( istep == (Nstep/isampleHC) ) then

     if(.not.QPBC) then
       open (83, file = 'Gr.out')
       do igrid = 1, 1000
         ROH = acc_ROH(igrid)*riacc
         RHH = acc_RHH(igrid)*riacc
         write (83,'(3d23.15)') dble(igrid)*0.01d0, ROH, RHH
       end do
       close (83)
     end if

     close (81)
     close (82)
     close (85)

     irestart = iacc

     open(84,file='HCestimater.restart',form='unformatted',status='unknown')

     write(84) irestart
     write(84) acc_K_T
     write(84) acc_K_V
     write(84) acc_K_CV
     write(84) acc_Vbar
     write(84) acc_E_T
     write(84) acc_E_V
     write(84) acc_E_CV
     write(84) acc_E_T_E_T
     write(84) acc_E_T_E_V
     write(84) acc_E_T_E_CV
     write(84) acc_E_V_E_V
     write(84) acc_E_CV_E_CV
     write(84) acc_dE_T_dbeta
     write(84) acc_B_V_E_V
     write(84) acc_B_CV_E_CV
     if(.not.QPBC) then
       write(84) acc_ROH
       write(84) acc_RHH
     end if

     close(84)

   end if


end subroutine HeatCapacity


!######################################################################
!######################################################################


subroutine Calc_Hessian(HessCVTrace,HessTrace)

use Numbers, only : N, NumSpec, NumMol, NumAtm
use CommonBlocks, only : QPBC, QPathInt
use Configuration, only : R
use CommonMPI
use CommonPI
use HeatCap
use NoLJparam
use BookParam, only : Npair, ListIJ
use UnitExParam, only : reng, pi
use NonbondParam, only : Charge, SgmLJ, EpsLJ
use CellParam, only : H, InvH
use AtomParam, only : ResidName
use CutoffParam, only : Rcutoff2

implicit none

integer :: i, j, Nas, ii, k, l, jj, isp
integer :: iOO, iH1, iH2
real(8), parameter :: ScaleParam = 1.d-18 * reng
real(8), parameter :: a_const =  2.361d0*2.361d0*0.708d0*ScaleParam !3.947d0 * ScaleParam
real(8), parameter :: b_const =  1.803d0 * ScaleParam * 0.5d0
real(8), parameter :: c_const = -1.469d0 * ScaleParam
real(8), parameter :: d_const =  0.776d0 * ScaleParam
real(8), parameter :: b_OHeq = 1.d0
real(8) :: b_HHeq
real(8), dimension(3) :: dOH1, dOH2, dHH0
real(8), dimension(3) :: sOH1, sOH2, sHH0
real(8) :: R2OH1, R2OH2, R2HH0, R1OH1, R1OH2, R1HH0
real(8) :: InvOH1, InvOH2, InvHH0, dR_OH1, dR_OH2, dR_HH0
integer :: NProcsTemp, MyRankTemp
integer :: Na, Na2
real(8) :: aconst2, bconst2, HES, ha, hb, hc, hd
real(8) :: X_OH1, Y_OH1, X_OH2, Y_OH2, X_HH0, Y_HH0
real(8) :: Z_O12, Z_O21, Z_1H0, Z_2H0, Z_HH0, xx1, xx2

real(8) :: Sgm2, Eps

real(8) :: R2, InvR2
real(8) :: SR2, SR6, SR12
real(8) :: cf
real(8), dimension(3) :: Rij, Sij, Dij
real(8), dimension(3,N) :: ScR

integer :: N2, non, noa
real(8) :: InvR1, dfkLJ, fkcLJ, dfk, fkc

real(8) :: vxx, vxy, vxz, vyy, vyx, vyz, vzx, vzy, vzz
real(8), dimension(3) :: Ri, Rk

real(8) :: HessCVTrace, HessTrace

   b_HHeq = b_OHeq * sin(pi*54.d0/180.d0) * 2.d0

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

! ## Zero clear

   Hessian = 0.d0

! ## Intramolecular SPCF

   jj = 0

   do i = 1, NumSpec

     if(ResidName(jj+1) == 'SPCF') then

       isp = i
       exit

     end if

     jj = jj + NumMol(i) * NumAtm(i)

   end do

   Na  = NumMol(isp) * NumAtm(isp)
   Na2 = Na * 2

   aconst2 = 2.d0 * a_const
   bconst2 = 2.d0 * b_const

   do j = Nas, NumMol(isp), NProcsTemp

     ii = (j-1) * 3 + jj

     iOO = ii + 1
     iH1 = ii + 2
     iH2 = ii + 3

     dOH1(:) = R(:,iOO) - R(:,iH1)
     dOH2(:) = R(:,iOO) - R(:,iH2)
     dHH0(:) = R(:,iH1) - R(:,iH2)

#ifdef PCC
     R2OH1 = dOH1(1) * dOH1(1) + dOH1(2) * dOH1(2) + dOH1(3) * dOH1(3)
     R2OH2 = dOH2(1) * dOH2(1) + dOH2(2) * dOH2(2) + dOH2(3) * dOH2(3)
     R2HH0 = dHH0(1) * dHH0(1) + dHH0(2) * dHH0(2) + dHH0(3) * dHH0(3)
#else
     R2OH1 = dot_product(dOH1, dOH1)
     R2OH2 = dot_product(dOH2, dOH2)
     R2HH0 = dot_product(dHH0, dHH0)
#endif
     R1OH1 = sqrt( R2OH1 )
     R1OH2 = sqrt( R2OH2 )
     R1HH0 = sqrt( R2HH0 )

     InvOH1 = 1.d0 / R1OH1
     InvOH2 = 1.d0 / R1OH2
     InvHH0 = 1.d0 / R1HH0

     dR_OH1 = R1OH1 - b_OHeq
     dR_OH2 = R1OH2 - b_OHeq
     dR_HH0 = R1HH0 - b_HHeq

     sOH1(:) = dOH1(:) * InvOH1
     sOH2(:) = dOH2(:) * InvOH2
     sHH0(:) = dHH0(:) * InvHH0

     X_OH1 = InvOH1*dR_OH1
     Y_OH1 = 1.d0 - X_OH1
     X_OH2 = InvOH2*dR_OH2
     Y_OH2 = 1.d0 - X_OH2
     X_HH0 = InvHH0*dR_HH0
     Y_HH0 = 1.d0 - X_HH0

     Z_O12 = InvOH1*dR_OH2
     Z_O21 = InvOH2*dR_OH1

     Z_1H0 = InvOH1*dR_HH0
     Z_2H0 = InvOH2*dR_HH0
     Z_HH0 = InvHH0*(dR_OH1+dR_OH2)

! ## 
     k = iOO
     l = iOO

     xx1 = sOH1(1) * sOH1(1)
     xx2 = sOH2(1) * sOH2(1)

     ha =   aconst2 * ( X_OH1 + xx1*Y_OH1 &
     &                + X_OH2 + xx2*Y_OH2 )
     hc =   c_const * ( Z_1H0*(1.d0 - xx1) &
     &                + Z_2H0*(1.d0 - xx2) )
     hd =   d_const * ( Z_O12*(1.d0 - xx1) &
     &                + Z_O21*(1.d0 - xx2) &
     &                + 2.d0*sOH2(1)*sOH1(1) )

     Hessian(k,l) = Hessian(k,l) + ha + hc + hd

! ## 
     k = iOO
     l = iOO + Na

     xx1 = sOH1(1) * sOH1(2)
     xx2 = sOH2(1) * sOH2(2)

     ha =   aconst2 * ( xx1*Y_OH1 + xx2*Y_OH2 )
     hc = - c_const * ( xx1*Z_1H0 + xx2*Z_2H0 )
     hd =   d_const * ( sOH1(2)*sOH2(1)       &
     &                + sOH1(1)*sOH2(2)       &
     &                - xx1*Z_O12 - xx2*Z_O21 )

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO
     l = iOO + Na2

     xx1 = sOH1(1) * sOH1(3)
     xx2 = sOH2(1) * sOH2(3)

     ha =   aconst2 * ( xx1*Y_OH1 + xx2*Y_OH2 )
     hc = - c_const * ( xx1*Z_1H0 + xx2*Z_2H0 )
     hd =   d_const * ( sOH1(1)*sOH2(3) &
     &                + sOH1(3)*sOH2(1) &
     &                - xx1*Z_O12 - xx2*Z_O21 )

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na
     l = iOO + Na

     xx1 = sOH1(2) * sOH1(2)
     xx2 = sOH2(2) * sOH2(2)

     ha =   aconst2 * ( X_OH1 + xx1*Y_OH1 &
     &                + X_OH2 + xx2*Y_OH2 )
     hc =   c_const * ( (1.d0 - xx1)*Z_1H0 &
     &                + (1.d0 - xx2)*Z_2H0 )
     hd =   d_const * ( Z_O12 * (1.d0-xx1)   &
     &                + Z_O21 * (1.d0-xx2)   &
     &                + sOH1(2)*sOH2(2)*2.d0 )

     Hessian(k,l) = Hessian(k,l) + ha + hc + hd

! ## 
     k = iOO + Na
     l = iOO + Na2

     xx1 = sOH1(2) * sOH1(3)
     xx2 = sOH2(2) * sOH2(3)

     ha =   aconst2 * ( xx1*Y_OH1 + xx2*Y_OH2 )
     hc = - c_const * ( xx1*Z_1H0 + xx2*Z_2H0 )
     hd =   d_const * ( sOH1(2)*sOH2(3) &
     &                + sOH2(2)*sOH1(3) &
     &                - xx1*Z_O12 - xx2*Z_O21 )

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na2
     l = iOO + Na2

     xx1 = sOH1(3) * sOH1(3)
     xx2 = sOH2(3) * sOH2(3)

     ha =   aconst2 * ( X_OH1 + xx1*Y_OH1 &
     &                + X_OH2 + xx2*Y_OH2 )
     hc =   c_const * ( (1.d0 - xx1)*Z_1H0 &
     &                + (1.d0 - xx2)*Z_2H0 )
     hd =   d_const * ( sOH1(3)*sOH2(3)*2.d0 &
     &                + (1.d0 - xx1)*Z_O12 &
     &                + (1.d0 - xx2)*Z_O21 )

     Hessian(k,l) = Hessian(k,l) + ha + hc + hd

!***********************************************************************

! ## 
     k = iOO
     l = iH1

     xx1 = sOH1(1)*sOH1(1)

     ha = - aconst2 * (  xx1*Y_OH1 + X_OH1 )
     hc =   c_const * ( (xx1 - 1.d0)*Z_1H0 &
     &                + (sOH1(1)+sOH2(1))*sHH0(1) )
     hd =   d_const * ( Z_O12*(xx1 - 1.d0) &
     &                - sOH2(1)*sOH1(1) )

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO
     l = iH1 + Na

     xx1 = sOH1(1)*sOH1(2)

     ha = - aconst2 * xx1 * Y_OH1
     hc =   c_const * ( xx1*Z_1H0     &
     &                + sHH0(2)*(sOH1(1)+sOH2(1)) )
     hd =   d_const * (xx1*Z_O12 - sOH1(2)*sOH2(1))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO
     l = iH1 + Na2

     xx1 = sOH1(1)*sOH1(3)

     ha = - aconst2 * xx1 * Y_OH1
     hc =   c_const * ( xx1*Z_1H0    &
     &                + sHH0(3)*(sOH1(1)+sOH2(1)) )
     hd =   d_const * (xx1*Z_O12 - sOH1(3)*sOH2(1))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iOO

     hc =   c_const * ( xx1*Z_1H0    &
     &                + sHH0(3)*(sOH1(1)+sOH2(1)) )
     hd =   d_const * (xx1*Z_O12 - sOH1(3)*sOH2(1))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iOO + Na
     l = iH1

     xx1 = sOH1(2)*sOH1(1)

     ha = - aconst2 * xx1 * Y_OH1
     hc =   c_const * ( xx1*Z_1H0   &
     &                + sHH0(1)*(sOH1(2)+sOH2(2)) )
     hd =   d_const * (xx1*Z_O12-sOH1(1)*sOH2(2))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1
     l = iOO + Na

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iOO + Na
     l = iH1 + Na

     xx1 = sOH1(2)*sOH1(2)

     ha = - aconst2 * ( xx1*Y_OH1 + X_OH1)
     hc =   c_const * ( (xx1 - 1.d0)*Z_1H0 &
    &                 + sHH0(2)*(sOH1(2) + sOH2(2)))
     hd =   d_const * ( (xx1 - 1.d0)*Z_O12 &
     &                - sOH2(2)*sOH1(2) )

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na
     l = iH1 + Na2

     xx1 = sOH1(2)*sOH1(3)

     ha = - aconst2 * xx1 * Y_OH1
     hc =   c_const * ( xx1*Z_1H0       &
     &                + sHH0(3)*(sOH1(2)+sOH2(2)) )
     hd =   d_const * (xx1*Z_O12 - sOH1(3)*sOH2(2))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iOO + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iOO + Na2
     l = iH1

     xx1 = sOH1(1)*sOH1(3)

     ha = - aconst2 * xx1 * Y_OH1
     hc =   c_const * ( xx1 * Z_1H0   &
     &                + (sOH1(3)+sOH2(3))* sHH0(1) )
     hd =   d_const * (xx1*Z_O12-sOH1(1)*sOH2(3))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na2
     l = iH1 + Na

     xx1 = sOH1(2)*sOH1(3)

     ha = - aconst2 * xx1 * Y_OH1
     hc =   c_const * ( xx1 * Z_1H0   &
     &                + sHH0(2)*(sOH1(3) + sOH2(3)) )
     hd =   d_const * (xx1*Z_O12 - sOH1(2)*sOH2(3))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na2
     l = iH1 + Na2

     xx1 = sOH1(3)*sOH1(3)

     ha = - aconst2 * (xx1 * Y_OH1 + X_OH1)
     hc =   c_const * ( Z_1H0*(xx1-1.d0) &
     &                + (sOH1(3)+sOH2(3))*sHH0(3))
     hd =   d_const * ((xx1-1.d0)*Z_O12-sOH1(3)*sOH2(3))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

!***********************************************************************

! ## 
     k = iOO
     l = iH2

     xx1 = sOH2(1)*sOH2(1)

     ha = - aconst2 * ( xx1*Y_OH2+X_OH2 )
     hc =   c_const * ( Z_2H0*(xx1-1.d0) &
     &                - sHH0(1)*(sOH2(1)+sOH1(1)) )
     hd =   d_const * ((xx1-1.d0)*Z_O21-sOH2(1)*sOH1(1))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO
     l = iH2 + Na

     xx1 = sOH2(1)*sOH2(2)

     ha = - aconst2 * xx1 * Y_OH2
     hc =   c_const * ( xx1*Z_2H0   &
     &                - (sOH2(1)+sOH1(1)) * sHH0(2) )
     hd =   d_const * (xx1*Z_O21-sOH2(2)*sOH1(1))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO
     l = iH2 + Na2

     xx1 = sOH2(1)*sOH2(3)

     ha = - aconst2 * xx1 * Y_OH2
     hc =   c_const * ( xx1 * Z_2H0   &
     &                - (sOH2(1)+sOH1(1))*sHH0(3) )
     hd =   d_const * (xx1*Z_O21-sOH2(3)*sOH1(1))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na
     l = iH2

     xx1 = sOH2(2)*sOH2(1)

     ha = - aconst2 * xx1 * Y_OH2
     hc =   c_const * ( xx1*Z_2H0   &
     &                - (sOH2(2)+sOH1(2))*sHH0(1) )
     hd =   d_const * (xx1*Z_O21-sOH2(1)*sOH1(2))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2
     l = iOO + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iOO + Na
     l = iH2 + Na

     xx1 = sOH2(2)*sOH2(2)

     ha = - aconst2 * (xx1*Y_OH2 + X_OH2)
     hc =   c_const * ( Z_2H0*(xx1-1.d0) &
     &                - (sOH2(2)+sOH1(2))*sHH0(2) )
     hd =   d_const * ((xx1-1.d0)*Z_O21-sOH2(2)*sOH1(2))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na
     l = iH2 + Na2

     xx1 = sOH2(2)*sOH2(3)

     ha = - aconst2 * xx1*Y_OH2
     hc =   c_const * ( xx1*Z_2H0 &
     &                - (sOH2(2)+sOH1(2))*sHH0(3) )
     hd =   d_const * (xx1*Z_O21-sOH2(3)*sOH1(2))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iOO + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iOO + Na2
     l = iH2

     xx1 = sOH2(1)*sOH2(3)

     ha = - aconst2 * xx1 * Y_OH2
     hc =   c_const * ( xx1 * Z_2H0 &
     &                - (sOH2(3)+sOH1(3))*sHH0(1) )
     hd =   d_const * (xx1*Z_O21-sOH2(1)*sOH1(3))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na2
     l = iH2 + Na

     xx1 = sOH2(2)*sOH2(3)

     ha = - aconst2 * xx1 * Y_OH2
     hc =   c_const * ( xx1 * Z_2H0 &
     &                - (sOH2(3)+sOH1(3)) *sHH0(2) )
     hd =   d_const * (xx1*Z_O21-sOH2(2)*sOH1(3))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

! ## 
     k = iOO + Na2
     l = iH2 + Na2

     xx1 = sOH2(3)*sOH2(3)

     ha = - aconst2 * (xx1*Y_OH2 + X_OH2)
     hc =   c_const * ( Z_2H0*(xx1-1.d0) &
     &                - (sOH2(3)+sOH1(3))*sHH0(3) )
     hd =   d_const * ((xx1-1.d0)*Z_O21-sOH2(3)*sOH1(3))

     HES = ha + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     Hessian(l,k) = Hessian(l,k) + HES

!***********************************************************************

! ## 
     k = iH1
     l = iH1

     xx1 = sOH1(1)*sOH1(1)
     xx2 = sHH0(1)*sHH0(1)

     ha =   aconst2 * (xx1*Y_OH1+X_OH1)
     hb =   bconst2 * (xx2*Y_HH0+X_HH0 )
     hc =   c_const * ( (1.d0 - xx1)*Z_1H0 &
     &                - sOH1(1)*sHH0(1)*2.d0           &
     &                + (1.d0 - xx2)*Z_HH0 )
     hd =   d_const * (1.d0 - xx1)*Z_O12

     Hessian(k,l) = Hessian(k,l) + ha + hb + hc + hd

! ## 
     k = iH1
     l = iH1 + Na

     xx1 = sOH1(1)*sOH1(2)
     xx2 = sHH0(1)*sHH0(2)

     ha =   aconst2 * xx1*Y_OH1
     hb =   bconst2 * xx2*Y_HH0
     hc = - c_const * ( sOH1(1)*sHH0(2)+xx1*Z_1H0 &
     &                + sHH0(1)*sOH1(2)+xx2*Z_HH0 )
     hd = - d_const * xx1*Z_O12

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na
     l = iH1

     hc = - c_const * ( sOH1(2)*sHH0(1)+xx1*Z_1H0 &
     &                + sHH0(2)*sOH1(1)+xx2*Z_HH0 )

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1
     l = iH1 + Na2

     ha =   aconst2 * sOH1(1)*sOH1(3)*Y_OH1
     hb =   bconst2 * sHH0(1)*sHH0(3)*Y_HH0
     hc = - c_const * ( sOH1(1)*(sHH0(3)+sOH1(3)*Z_1H0) &
     &                + sHH0(1)*(sOH1(3)+sHH0(3)*Z_HH0) )
     hd = - d_const * sOH1(1)*sOH1(3)*Z_O12

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iH1

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na
     l = iH1 + Na

     ha =   aconst2 * (sOH1(2)*sOH1(2)*Y_OH1+X_OH1)
     hb =   bconst2 * (sHH0(2)*sHH0(2)*Y_HH0+X_HH0)
     hc =   c_const * ( (1.d0-sOH1(2)*sOH1(2))*Z_1H0 &
     &                - sOH1(2)*sHH0(2)*2.d0         &
     &                + (1.d0-sHH0(2)*sHH0(2))*Z_HH0 )
     hd =   d_const * (1.d0 - sOH1(2)*sOH1(2))*Z_O12

     Hessian(k,l) = Hessian(k,l) + ha + hb + hc + hd

! ## 
     k = iH1 + Na
     l = iH1 + Na2

     ha =   aconst2 * sOH1(2)*sOH1(3)*Y_OH1
     hb =   bconst2 * sHH0(2)*sHH0(3)*Y_HH0
     hc = - c_const * ( sOH1(2)*(sHH0(3)+sOH1(3)*Z_1H0) &
     &                + sHH0(2)*(sOH1(3)+sHH0(3)*Z_HH0) )
     hd = - d_const * sOH1(2)*sOH1(3)*Z_O12

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iH1 + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iH1 + Na2

     ha =   aconst2 * (sOH1(3)*sOH1(3)*Y_OH1+X_OH1)
     hb =   bconst2 * (sHH0(3)*sHH0(3)*Y_HH0+X_HH0)
     hc =   c_const * ( Z_1H0*(1.d0-SOH1(3)*SOH1(3))            &
     &                - sOH1(3)*sHH0(3)*2.d0          &
     &                + (1.d0 - sHH0(3)*sHH0(3))*Z_HH0 )
     hd =   d_const * (1.d0 - sOH1(3)*sOH1(3))*Z_O12

     Hessian(k,l) = Hessian(k,l) + ha + hb + hc + hd

!***********************************************************************

! ## 
     k = iH2
     l = iH2

     ha =   aconst2 * (sOH2(1)*sOH2(1)*Y_OH2+X_OH2)
     hb =   bconst2 * (sHH0(1)*sHH0(1)*Y_HH0+X_HH0)
     hc =   c_const * ( sOH2(1)*sHH0(1)*2.d0           &
     &                + (1.d0 - sOH2(1)*sOH2(1))*Z_2H0 &
     &                + (1.d0 - sHH0(1)*sHH0(1))*Z_HH0)
     hd =   d_const * (1.d0-sOH2(1)*sOH2(1))*Z_O21

     Hessian(k,l) = Hessian(k,l) + ha + hb + hc + hd

! ## 
     k = iH2
     l = iH2 + Na

     ha =   aconst2 * sOH2(1)*sOH2(2)*Y_OH2
     hb =   bconst2 * sHH0(1)*sHH0(2)*Y_HH0
     hc =   c_const * ( sOH2(1)*sHH0(2)       &
     &                + sOH2(2)*sHH0(1)       &
     &                - sOH2(1)*sOH2(2)*Z_2H0 &
     &                - sHH0(1)*sHH0(2)*Z_HH0 )
     hd = - d_const * sOH2(1)*sOH2(2)*Z_O21

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na
     l = iH2

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2
     l = iH2 + Na2

     ha =   aconst2 * sOH2(1)*sOH2(3)*Y_OH2
     hb =   bconst2 * sHH0(1)*sHH0(3)*Y_HH0
     hc =   c_const * ( sOH2(1)*sHH0(3) &
     &                + sOH2(3)*sHH0(1) &
     &                - sOH2(1)*sOH2(3)*Z_2H0    &
     &                - sHH0(1)*sHH0(3)*Z_HH0 )
     hd = - d_const * sOH2(1)*sOH2(3)*Z_O21

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iH2

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na
     l = iH2 + Na

     ha =   aconst2 * (sOH2(2)*sOH2(2)*Y_OH2+X_OH2)
     hb =   bconst2 * (sHH0(2)*sHH0(2)*Y_HH0+X_HH0)
     hc =   c_const * ( sOH2(2)*sHH0(2)*2.d0 &
     &                + (1.d0-sOH2(2)*sOH2(2))*Z_2H0 &
     &                + (1.d0-sHH0(2)*sHH0(2))*Z_HH0 )
     hd =   d_const * (1.d0-sOH2(2)*sOH2(2))*Z_O21

     Hessian(k,l) = Hessian(k,l) + ha + hb + hc + hd

! ## 
     k = iH2 + Na
     l = iH2 + Na2

     ha =   aconst2 * sOH2(2)*sOH2(3)*Y_OH2
     hb =   bconst2 * sHH0(2)*sHH0(3)*Y_HH0
     hc =   c_const * ( sOH2(2)*(sHH0(3)-sOH2(3)*Z_2H0) &
     &                + sHH0(2)*(sOH2(3)-sHH0(3)*Z_HH0) )
     hd = - d_const * sOH2(2)*sOH2(3)*Z_O21

     HES = ha + hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iH2 + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iH2 + Na2

     ha =   aconst2 * (sOH2(3)*sOH2(3)*Y_OH2+X_OH2)
     hb =   bconst2 * (sHH0(3)*sHH0(3)*Y_HH0+X_HH0)
     hc =   c_const * ( sOH2(3)*sHH0(3)*2.d0         &
     &                + (1.d0-sOH2(3)*sOH2(3))*Z_2H0 &
     &                + (1.d0-sHH0(3)*sHH0(3))*Z_HH0 )
     hd =   d_const * (1.d0-sOH2(3)*sOH2(3))*Z_O21

     Hessian(k,l) = Hessian(k,l) + ha + hb + hc + hd

!***********************************************************************

! ## 
     k = iH1
     l = iH2

     hb = - bconst2 * ( sHH0(1)*sHH0(1)*Y_HH0+X_HH0)
     hc =   c_const * ( (sOH1(1)-sOH2(1))*sHH0(1) &
     &                + (sHH0(1)*sHH0(1)-1.d0)*Z_HH0 )
     hd =   d_const * sOH1(1)*sOH2(1)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2
     l = iH1

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1
     l = iH2 + Na

     hb = - bconst2 * sHH0(1)*sHH0(2)*Y_HH0
     hc =   c_const * ( sOH1(1)*sHH0(2) &
     &                - sOH2(2)*sHH0(1) &
     &                + sHH0(1)*sHH0(2)*Z_HH0 )
     hd =   d_const * sOH1(1)*sOH2(2)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na
     l = iH1

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1
     l = iH2 + Na2

     hb = - bconst2 * sHH0(1)*sHH0(3)*Y_HH0
     hc =   c_const * ( sOH1(1)*sHH0(3)       &
     &                - sOH2(3)*sHH0(1)       &
     &                + sHH0(1)*sHH0(3)*Z_HH0 )
     hd =   d_const * sOH1(1)*sOH2(3)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iH1

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na
     l = iH2

     hb = - bconst2 * sHH0(2)*sHH0(1)*Y_HH0
     hc =   c_const * ( sOH1(2)*sHH0(1)       &
     &                - sOH2(1)*sHH0(2)       &
     &                + sHH0(2)*sHH0(1)*Z_HH0 )
     hd =   d_const * sOH1(2)*sOH2(1)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2
     l = iH1 + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na
     l = iH2 + Na

     hb = - bconst2 * (sHH0(2)*sHH0(2)*Y_HH0+X_HH0)
     hc =   c_const * ( (sOH1(2)-sOH2(2))*sHH0(2) &
     &                + (sHH0(2)*sHH0(2)-1.d0)*Z_HH0 )
     hd =   d_const * sOH1(2)*sOH2(2)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na
     l = iH1 + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na
     l = iH2 + Na2

     hb = - bconst2 * sHH0(2)*sHH0(3)*Y_HH0
     hc =   c_const * ( sOH1(2)*sHH0(3)       &
     &                - sOH2(3)*sHH0(2)       &
     &                + sHH0(2)*sHH0(3)*Z_HH0 )
     hd =   d_const * sOH1(2)*sOH2(3)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iH1 + Na

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iH2

     hb = - bconst2 * sHH0(3)*sHH0(1)*Y_HH0
     hc =   c_const * ( sOH1(3)*sHH0(1)       &
     &                - sOH2(1)*sHH0(3)       &
     &                + sHH0(3)*sHH0(1)*Z_HH0 )
     hd =   d_const * sOH1(3)*sOH2(1)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2
     l = iH1 + Na2

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iH2 + Na

     hb = - bconst2 * sHH0(3)*sHH0(2)*Y_HH0
     hc =   c_const * ( sOH1(3)*sHH0(2)       &
     &                - sOH2(2)*sHH0(3)       &
     &                + sHH0(3)*sHH0(2)*Z_HH0 )
     hd =   d_const * sOH1(3)*sOH2(2)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na
     l = iH1 + Na2

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH1 + Na2
     l = iH2 + Na2

     hb = - bconst2 * (sHH0(3)*sHH0(3)*Y_HH0+X_HH0)
     hc =   c_const * ( (sOH1(3)-sOH2(3))*sHH0(3) &
     &                + (sHH0(3)*sHH0(3)-1.d0)*Z_HH0 )
     hd =   d_const * sOH1(3)*sOH2(3)

     HES = hb + hc + hd

     Hessian(k,l) = Hessian(k,l) + HES

! ## 
     k = iH2 + Na2
     l = iH1 + Na2

     Hessian(k,l) = Hessian(k,l) + HES

   end do

! ## Intermolecular

   N2 = N * 2

   if(QPBC) then

     do i = 1 , N

#ifdef PCC
       ScR(1,i) = InvH(1,1) * R(1,i) + InvH(1,2) * R(2,i) + InvH(1,3) * R(3,i)
       ScR(2,i) = InvH(2,1) * R(1,i) + InvH(2,2) * R(2,i) + InvH(2,3) * R(3,i)
       ScR(3,i) = InvH(3,1) * R(1,i) + InvH(3,2) * R(2,i) + InvH(3,3) * R(3,i)
#else
       ScR(:,i) = matmul( InvH , R(:,i) )
#endif

     end do

     do l = 1 , Npair

       i = ListIJ(1,l)
       j = ListIJ(2,l)

       Sij = ScR(:,i) - ScR(:,j)
#ifdef PCC
       if(Sij(1) >  0.5d0) Sij(1) = Sij(1) - 1.d0
       if(Sij(1) < -0.5d0) Sij(1) = Sij(1) + 1.d0
       if(Sij(2) >  0.5d0) Sij(2) = Sij(2) - 1.d0
       if(Sij(2) < -0.5d0) Sij(2) = Sij(2) + 1.d0
       if(Sij(3) >  0.5d0) Sij(3) = Sij(3) - 1.d0
       if(Sij(3) < -0.5d0) Sij(3) = Sij(3) + 1.d0
       Rij(1) = H(1,1) * Sij(1) + H(1,2) * Sij(2) + H(1,3) * Sij(3)
       Rij(2) = H(2,1) * Sij(1) + H(2,2) * Sij(2) + H(2,3) * Sij(3)
       Rij(3) = H(3,1) * Sij(1) + H(3,2) * Sij(2) + H(3,3) * Sij(3)
       R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
       Sij = Sij - nint( Sij )
       Rij = matmul( H , Sij )
       R2  = dot_product( Rij, Rij )
#endif

       if(R2 <= Rcutoff2) then

         InvR2 = 1.d0 / R2
         Eps   = EpsLJ(i) * EpsLJ(j)
         InvR1 = sqrt( InvR2 )
         Dij(:) = Rij(:) * InvR1

         if(Eps /= 0.) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkcLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           dfkLJ = Eps * 24.d0 * ( 26.d0 * SR12 - 7.d0 * SR6 ) * InvR2

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(1) - fkcLJ
           Hessian(i,i) = Hessian(i,i) + HES
           Hessian(i,j) = Hessian(i,j) - HES
           Hessian(j,i) = Hessian(j,i) - HES
           Hessian(j,j) = Hessian(j,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(2)
           Hessian(i,i+N) = Hessian(i,i+N) + HES
           Hessian(i,j+N) = Hessian(i,j+N) - HES
           Hessian(j,i+N) = Hessian(j,i+N) - HES
           Hessian(j,j+N) = Hessian(j,j+N) + HES
           Hessian(i+N,i) = Hessian(i+N,i) + HES
           Hessian(i+N,j) = Hessian(i+N,j) - HES
           Hessian(j+N,i) = Hessian(j+N,i) - HES
           Hessian(j+N,j) = Hessian(j+N,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(3)
           Hessian(i,i+N2) = Hessian(i,i+N2) + HES
           Hessian(i,j+N2) = Hessian(i,j+N2) - HES
           Hessian(j,i+N2) = Hessian(j,i+N2) - HES
           Hessian(j,j+N2) = Hessian(j,j+N2) + HES
           Hessian(i+N2,i) = Hessian(i+N2,i) + HES
           Hessian(i+N2,j) = Hessian(i+N2,j) - HES
           Hessian(j+N2,i) = Hessian(j+N2,i) - HES
           Hessian(j+N2,j) = Hessian(j+N2,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(2) * Dij(2) - fkcLJ
           Hessian(i+N,i+N) = Hessian(i+N,i+N) + HES
           Hessian(i+N,j+N) = Hessian(i+N,j+N) - HES
           Hessian(j+N,i+N) = Hessian(j+N,i+N) - HES
           Hessian(j+N,j+N) = Hessian(j+N,j+N) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(2) * Dij(3)
           Hessian(i+N,i+N2) = Hessian(i+N,i+N2) + HES
           Hessian(i+N,j+N2) = Hessian(i+N,j+N2) - HES
           Hessian(j+N,i+N2) = Hessian(j+N,i+N2) - HES
           Hessian(j+N,j+N2) = Hessian(j+N,j+N2) + HES
           Hessian(i+N2,i+N) = Hessian(i+N2,i+N) + HES
           Hessian(i+N2,j+N) = Hessian(i+N2,j+N) - HES
           Hessian(j+N2,i+N) = Hessian(j+N2,i+N) - HES
           Hessian(j+N2,j+N) = Hessian(j+N2,j+N) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(3) * Dij(3) - fkcLJ
           Hessian(i+N2,i+N2) = Hessian(i+N2,i+N2) + HES
           Hessian(i+N2,j+N2) = Hessian(i+N2,j+N2) - HES
           Hessian(j+N2,i+N2) = Hessian(j+N2,i+N2) - HES
           Hessian(j+N2,j+N2) = Hessian(j+N2,j+N2) + HES

         end if

         cf = Charge(i) * Charge(j)

         fkc = cf * InvR2 * InvR1
         dfk = 3.d0 * fkc

         HES = dfk * Dij(1) * Dij(1) - fkc
         Hessian(i,i) = Hessian(i,i) + HES
         Hessian(i,j) = Hessian(i,j) - HES
         Hessian(j,i) = Hessian(j,i) - HES
         Hessian(j,j) = Hessian(j,j) + HES

         HES = dfk * Dij(1) * Dij(2)
         Hessian(i,i+N) = Hessian(i,i+N) + HES
         Hessian(i,j+N) = Hessian(i,j+N) - HES
         Hessian(j,i+N) = Hessian(j,i+N) - HES
         Hessian(j,j+N) = Hessian(j,j+N) + HES
         Hessian(i+N,i) = Hessian(i+N,i) + HES
         Hessian(i+N,j) = Hessian(i+N,j) - HES
         Hessian(j+N,i) = Hessian(j+N,i) - HES
         Hessian(j+N,j) = Hessian(j+N,j) + HES

         HES = dfk * Dij(1) * Dij(3)
         Hessian(i,i+N2) = Hessian(i,i+N2) + HES
         Hessian(i,j+N2) = Hessian(i,j+N2) - HES
         Hessian(j,i+N2) = Hessian(j,i+N2) - HES
         Hessian(j,j+N2) = Hessian(j,j+N2) + HES
         Hessian(i+N2,i) = Hessian(i+N2,i) + HES
         Hessian(i+N2,j) = Hessian(i+N2,j) - HES
         Hessian(j+N2,i) = Hessian(j+N2,i) - HES
         Hessian(j+N2,j) = Hessian(j+N2,j) + HES

         HES = dfk * Dij(2) * Dij(2) - fkc
         Hessian(i+N,i+N) = Hessian(i+N,i+N) + HES
         Hessian(i+N,j+N) = Hessian(i+N,j+N) - HES
         Hessian(j+N,i+N) = Hessian(j+N,i+N) - HES
         Hessian(j+N,j+N) = Hessian(j+N,j+N) + HES

         HES = dfk * Dij(2) * Dij(3)
         Hessian(i+N,i+N2) = Hessian(i+N,i+N2) + HES
         Hessian(i+N,j+N2) = Hessian(i+N,j+N2) - HES
         Hessian(j+N,i+N2) = Hessian(j+N,i+N2) - HES
         Hessian(j+N,j+N2) = Hessian(j+N,j+N2) + HES
         Hessian(i+N2,i+N) = Hessian(i+N2,i+N) + HES
         Hessian(i+N2,j+N) = Hessian(i+N2,j+N) - HES
         Hessian(j+N2,i+N) = Hessian(j+N2,i+N) - HES
         Hessian(j+N2,j+N) = Hessian(j+N2,j+N) + HES

         HES = dfk * Dij(3) * Dij(3) - fkc
         Hessian(i+N2,i+N2) = Hessian(i+N2,i+N2) + HES
         Hessian(i+N2,j+N2) = Hessian(i+N2,j+N2) - HES
         Hessian(j+N2,i+N2) = Hessian(j+N2,i+N2) - HES
         Hessian(j+N2,j+N2) = Hessian(j+N2,j+N2) + HES

       end if

     end do

   else  ! for isolated system

     do i = Nas, N, NProcsTemp

       non = NumNoLJ(i)

       do j = i-2 , 1, -2

         noa = 0
hep6:    do k = 1, non
           if(j == NoLJ(i,k)) then
             noa=1
             exit hep6
           end if
         end do hep6

         if(noa == 1) cycle

         Rij = R(:,i) - R(:,j)
#ifdef PCC
         R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
         R2  = dot_product( Rij, Rij )
#endif

         InvR2 = 1.d0 / R2
         InvR1 = sqrt( InvR2 )
         Dij(:) = Rij(:) * InvR1

         Eps   = EpsLJ(i) * EpsLJ(j)

         if(Eps /= 0.) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkcLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           dfkLJ = Eps * 24.d0 * ( 26.d0 * SR12 - 7.d0 * SR6 ) * InvR2

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(1) - fkcLJ
           Hessian(i,i) = Hessian(i,i) + HES
           Hessian(i,j) = Hessian(i,j) - HES
           Hessian(j,i) = Hessian(j,i) - HES
           Hessian(j,j) = Hessian(j,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(2)
           Hessian(i,i+N) = Hessian(i,i+N) + HES
           Hessian(i,j+N) = Hessian(i,j+N) - HES
           Hessian(j,i+N) = Hessian(j,i+N) - HES
           Hessian(j,j+N) = Hessian(j,j+N) + HES
           Hessian(i+N,i) = Hessian(i+N,i) + HES
           Hessian(i+N,j) = Hessian(i+N,j) - HES
           Hessian(j+N,i) = Hessian(j+N,i) - HES
           Hessian(j+N,j) = Hessian(j+N,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(3)
           Hessian(i,i+N2) = Hessian(i,i+N2) + HES
           Hessian(i,j+N2) = Hessian(i,j+N2) - HES
           Hessian(j,i+N2) = Hessian(j,i+N2) - HES
           Hessian(j,j+N2) = Hessian(j,j+N2) + HES
           Hessian(i+N2,i) = Hessian(i+N2,i) + HES
           Hessian(i+N2,j) = Hessian(i+N2,j) - HES
           Hessian(j+N2,i) = Hessian(j+N2,i) - HES
           Hessian(j+N2,j) = Hessian(j+N2,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(2) * Dij(2) - fkcLJ
           Hessian(i+N,i+N) = Hessian(i+N,i+N) + HES
           Hessian(i+N,j+N) = Hessian(i+N,j+N) - HES
           Hessian(j+N,i+N) = Hessian(j+N,i+N) - HES
           Hessian(j+N,j+N) = Hessian(j+N,j+N) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(2) * Dij(3)
           Hessian(i+N,i+N2) = Hessian(i+N,i+N2) + HES
           Hessian(i+N,j+N2) = Hessian(i+N,j+N2) - HES
           Hessian(j+N,i+N2) = Hessian(j+N,i+N2) - HES
           Hessian(j+N,j+N2) = Hessian(j+N,j+N2) + HES
           Hessian(i+N2,i+N) = Hessian(i+N2,i+N) + HES
           Hessian(i+N2,j+N) = Hessian(i+N2,j+N) - HES
           Hessian(j+N2,i+N) = Hessian(j+N2,i+N) - HES
           Hessian(j+N2,j+N) = Hessian(j+N2,j+N) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(3) * Dij(3) - fkcLJ
           Hessian(i+N2,i+N2) = Hessian(i+N2,i+N2) + HES
           Hessian(i+N2,j+N2) = Hessian(i+N2,j+N2) - HES
           Hessian(j+N2,i+N2) = Hessian(j+N2,i+N2) - HES
           Hessian(j+N2,j+N2) = Hessian(j+N2,j+N2) + HES

         end if

         cf = Charge(i) * Charge(j)

         fkc = cf * InvR2 * InvR1
         dfk = 3.d0 * fkc

         HES = dfk * Dij(1) * Dij(1) - fkc
         Hessian(i,i) = Hessian(i,i) + HES
         Hessian(i,j) = Hessian(i,j) - HES
         Hessian(j,i) = Hessian(j,i) - HES
         Hessian(j,j) = Hessian(j,j) + HES

         HES = dfk * Dij(1) * Dij(2)
         Hessian(i,i+N) = Hessian(i,i+N) + HES
         Hessian(i,j+N) = Hessian(i,j+N) - HES
         Hessian(j,i+N) = Hessian(j,i+N) - HES
         Hessian(j,j+N) = Hessian(j,j+N) + HES
         Hessian(i+N,i) = Hessian(i+N,i) + HES
         Hessian(i+N,j) = Hessian(i+N,j) - HES
         Hessian(j+N,i) = Hessian(j+N,i) - HES
         Hessian(j+N,j) = Hessian(j+N,j) + HES

         HES = dfk * Dij(1) * Dij(3)
         Hessian(i,i+N2) = Hessian(i,i+N2) + HES
         Hessian(i,j+N2) = Hessian(i,j+N2) - HES
         Hessian(j,i+N2) = Hessian(j,i+N2) - HES
         Hessian(j,j+N2) = Hessian(j,j+N2) + HES
         Hessian(i+N2,i) = Hessian(i+N2,i) + HES
         Hessian(i+N2,j) = Hessian(i+N2,j) - HES
         Hessian(j+N2,i) = Hessian(j+N2,i) - HES
         Hessian(j+N2,j) = Hessian(j+N2,j) + HES

         HES = dfk * Dij(2) * Dij(2) - fkc
         Hessian(i+N,i+N) = Hessian(i+N,i+N) + HES
         Hessian(i+N,j+N) = Hessian(i+N,j+N) - HES
         Hessian(j+N,i+N) = Hessian(j+N,i+N) - HES
         Hessian(j+N,j+N) = Hessian(j+N,j+N) + HES

         HES = dfk * Dij(2) * Dij(3)
         Hessian(i+N,i+N2) = Hessian(i+N,i+N2) + HES
         Hessian(i+N,j+N2) = Hessian(i+N,j+N2) - HES
         Hessian(j+N,i+N2) = Hessian(j+N,i+N2) - HES
         Hessian(j+N,j+N2) = Hessian(j+N,j+N2) + HES
         Hessian(i+N2,i+N) = Hessian(i+N2,i+N) + HES
         Hessian(i+N2,j+N) = Hessian(i+N2,j+N) - HES
         Hessian(j+N2,i+N) = Hessian(j+N2,i+N) - HES
         Hessian(j+N2,j+N) = Hessian(j+N2,j+N) + HES

         HES = dfk * Dij(3) * Dij(3) - fkc
         Hessian(i+N2,i+N2) = Hessian(i+N2,i+N2) + HES
         Hessian(i+N2,j+N2) = Hessian(i+N2,j+N2) - HES
         Hessian(j+N2,i+N2) = Hessian(j+N2,i+N2) - HES
         Hessian(j+N2,j+N2) = Hessian(j+N2,j+N2) + HES

       end do

       do j = i+1 , N, 2

         noa = 0
hep7:    do k = 1, non
           if(j == NoLJ(i,k)) then
             noa=1
             exit hep7
           end if
         end do hep7

         if(noa == 1) cycle

         Rij = R(:,i) - R(:,j)
#ifdef PCC
         R2  = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
#else
         R2  = dot_product( Rij, Rij )
#endif

         InvR2 = 1.d0 / R2
         InvR1 = sqrt( InvR2 )
         Dij(:) = Rij(:) * InvR1

         Eps   = EpsLJ(i) * EpsLJ(j)

         if(Eps /= 0.) then

           Sgm2  = SgmLJ(i) * SgmLJ(j)

           SR2  = Sgm2 * InvR2                    !(sigma/r)^2
           SR6  = SR2 * SR2 * SR2                 !         ^6
           SR12 = SR6 * SR6                       !         ^12
           fkcLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
           dfkLJ = Eps * 24.d0 * ( 26.d0 * SR12 - 7.d0 * SR6 ) * InvR2

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(1) - fkcLJ
           Hessian(i,i) = Hessian(i,i) + HES
           Hessian(i,j) = Hessian(i,j) - HES
           Hessian(j,i) = Hessian(j,i) - HES
           Hessian(j,j) = Hessian(j,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(2)
           Hessian(i,i+N) = Hessian(i,i+N) + HES
           Hessian(i,j+N) = Hessian(i,j+N) - HES
           Hessian(j,i+N) = Hessian(j,i+N) - HES
           Hessian(j,j+N) = Hessian(j,j+N) + HES
           Hessian(i+N,i) = Hessian(i+N,i) + HES
           Hessian(i+N,j) = Hessian(i+N,j) - HES
           Hessian(j+N,i) = Hessian(j+N,i) - HES
           Hessian(j+N,j) = Hessian(j+N,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(1) * Dij(3)
           Hessian(i,i+N2) = Hessian(i,i+N2) + HES
           Hessian(i,j+N2) = Hessian(i,j+N2) - HES
           Hessian(j,i+N2) = Hessian(j,i+N2) - HES
           Hessian(j,j+N2) = Hessian(j,j+N2) + HES
           Hessian(i+N2,i) = Hessian(i+N2,i) + HES
           Hessian(i+N2,j) = Hessian(i+N2,j) - HES
           Hessian(j+N2,i) = Hessian(j+N2,i) - HES
           Hessian(j+N2,j) = Hessian(j+N2,j) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(2) * Dij(2) - fkcLJ
           Hessian(i+N,i+N) = Hessian(i+N,i+N) + HES
           Hessian(i+N,j+N) = Hessian(i+N,j+N) - HES
           Hessian(j+N,i+N) = Hessian(j+N,i+N) - HES
           Hessian(j+N,j+N) = Hessian(j+N,j+N) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(2) * Dij(3)
           Hessian(i+N,i+N2) = Hessian(i+N,i+N2) + HES
           Hessian(i+N,j+N2) = Hessian(i+N,j+N2) - HES
           Hessian(j+N,i+N2) = Hessian(j+N,i+N2) - HES
           Hessian(j+N,j+N2) = Hessian(j+N,j+N2) + HES
           Hessian(i+N2,i+N) = Hessian(i+N2,i+N) + HES
           Hessian(i+N2,j+N) = Hessian(i+N2,j+N) - HES
           Hessian(j+N2,i+N) = Hessian(j+N2,i+N) - HES
           Hessian(j+N2,j+N) = Hessian(j+N2,j+N) + HES

           HES = ( dfkLJ + fkcLJ ) * Dij(3) * Dij(3) - fkcLJ
           Hessian(i+N2,i+N2) = Hessian(i+N2,i+N2) + HES
           Hessian(i+N2,j+N2) = Hessian(i+N2,j+N2) - HES
           Hessian(j+N2,i+N2) = Hessian(j+N2,i+N2) - HES
           Hessian(j+N2,j+N2) = Hessian(j+N2,j+N2) + HES

         end if

         cf = Charge(i) * Charge(j)

         fkc = cf * InvR2 * InvR1
         dfk = 3.d0 * fkc

         HES = dfk * Dij(1) * Dij(1) - fkc
         Hessian(i,i) = Hessian(i,i) + HES
         Hessian(i,j) = Hessian(i,j) - HES
         Hessian(j,i) = Hessian(j,i) - HES
         Hessian(j,j) = Hessian(j,j) + HES

         HES = dfk * Dij(1) * Dij(2)
         Hessian(i,i+N) = Hessian(i,i+N) + HES
         Hessian(i,j+N) = Hessian(i,j+N) - HES
         Hessian(j,i+N) = Hessian(j,i+N) - HES
         Hessian(j,j+N) = Hessian(j,j+N) + HES
         Hessian(i+N,i) = Hessian(i+N,i) + HES
         Hessian(i+N,j) = Hessian(i+N,j) - HES
         Hessian(j+N,i) = Hessian(j+N,i) - HES
         Hessian(j+N,j) = Hessian(j+N,j) + HES

         HES = dfk * Dij(1) * Dij(3)
         Hessian(i,i+N2) = Hessian(i,i+N2) + HES
         Hessian(i,j+N2) = Hessian(i,j+N2) - HES
         Hessian(j,i+N2) = Hessian(j,i+N2) - HES
         Hessian(j,j+N2) = Hessian(j,j+N2) + HES
         Hessian(i+N2,i) = Hessian(i+N2,i) + HES
         Hessian(i+N2,j) = Hessian(i+N2,j) - HES
         Hessian(j+N2,i) = Hessian(j+N2,i) - HES
         Hessian(j+N2,j) = Hessian(j+N2,j) + HES

         HES = dfk * Dij(2) * Dij(2) - fkc
         Hessian(i+N,i+N) = Hessian(i+N,i+N) + HES
         Hessian(i+N,j+N) = Hessian(i+N,j+N) - HES
         Hessian(j+N,i+N) = Hessian(j+N,i+N) - HES
         Hessian(j+N,j+N) = Hessian(j+N,j+N) + HES

         HES = dfk * Dij(2) * Dij(3)
         Hessian(i+N,i+N2) = Hessian(i+N,i+N2) + HES
         Hessian(i+N,j+N2) = Hessian(i+N,j+N2) - HES
         Hessian(j+N,i+N2) = Hessian(j+N,i+N2) - HES
         Hessian(j+N,j+N2) = Hessian(j+N,j+N2) + HES
         Hessian(i+N2,i+N) = Hessian(i+N2,i+N) + HES
         Hessian(i+N2,j+N) = Hessian(i+N2,j+N) - HES
         Hessian(j+N2,i+N) = Hessian(j+N2,i+N) - HES
         Hessian(j+N2,j+N) = Hessian(j+N2,j+N) + HES

         HES = dfk * Dij(3) * Dij(3) - fkc
         Hessian(i+N2,i+N2) = Hessian(i+N2,i+N2) + HES
         Hessian(i+N2,j+N2) = Hessian(i+N2,j+N2) - HES
         Hessian(j+N2,i+N2) = Hessian(j+N2,i+N2) - HES
         Hessian(j+N2,j+N2) = Hessian(j+N2,j+N2) + HES

       end do

     end do

     do i = 1, N

       Ri(:) = R(:,i)

       do j = 1, N

         vxx = Hessian(i   ,j   )
         vxy = Hessian(i   ,j+N )
         vxz = Hessian(i   ,j+N2)
         vyx = Hessian(i+N ,j   )
         vyy = Hessian(i+N ,j+N )
         vyz = Hessian(i+N ,j+N2)
         vzx = Hessian(i+N2,j   )
         vzy = Hessian(i+N2,j+N )
         vzz = Hessian(i+N2,j+N2)

         Rk(:) = R(:,j)
         HessTrace = HessTrace + Ri(1) * Rk(1) * vxx &
         &                     + Ri(1) * Rk(2) * vxy &
         &                     + Ri(1) * Rk(3) * vxz &
         &                     + Ri(2) * Rk(1) * vyx &
         &                     + Ri(2) * Rk(2) * vyy &
         &                     + Ri(2) * Rk(3) * vyz &
         &                     + Ri(3) * Rk(1) * vzx &
         &                     + Ri(3) * Rk(2) * vzy &
         &                     + Ri(3) * Rk(3) * vzz

       end do

     end do

   end if

   do i = 1, N

     Ri(:) = R(:,i) - Rnm(:,i,1)

     do j = 1, N

       vxx = Hessian(i   ,j   )
       vxy = Hessian(i   ,j+N )
       vxz = Hessian(i   ,j+N2)
       vyx = Hessian(i+N ,j   )
       vyy = Hessian(i+N ,j+N )
       vyz = Hessian(i+N ,j+N2)
       vzx = Hessian(i+N2,j   )
       vzy = Hessian(i+N2,j+N )
       vzz = Hessian(i+N2,j+N2)

       Rk(:) = R(:,j) - Rnm(:,j,1)
       HessCVTrace = HessCVTrace + Ri(1) * Rk(1) * vxx &
       &                         + Ri(1) * Rk(2) * vxy &
       &                         + Ri(1) * Rk(3) * vxz &
       &                         + Ri(2) * Rk(1) * vyx &
       &                         + Ri(2) * Rk(2) * vyy &
       &                         + Ri(2) * Rk(3) * vyz &
       &                         + Ri(3) * Rk(1) * vzx &
       &                         + Ri(3) * Rk(2) * vzy &
       &                         + Ri(3) * Rk(3) * vzz

     end do

   end do

end subroutine Calc_Hessian
