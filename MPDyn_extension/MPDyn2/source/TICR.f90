

!######################################################################
!######################################################################


subroutine SetTICR(ii,jj)

use CommonBlocks, only : cCOULOMB
use FEparam
use NonbondParam, only : Charge

implicit none

integer :: ii,i,j,jj
real(8) :: Lmd

   if(jj>iTIstep) then
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     write(*,*) 'WARNING : TI calculation has been already finished !'
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     Return
   end if

   if(ii == 0) then

     open(15,file='TICRconv.dat',status='unknown')
     open(16,file='TICRvalue.dat',status='unknown')

     if(QCharge) then

       allocate( ChargeOrg(NumTIpart) )

       j = 0
       do i = NiniTI+1, NfinTI
         j = j + 1
         ChargeOrg(j) = Charge(i)
       end do

     else ! ## Charge OFF

       do i = NiniTI+1, NfinTI
         Charge(i) = 0.d0
       end do

       if((trim(cCOULOMB) == 'EWALD').or.&
       &  (trim(cCOULOMB) == 'PME')  ) then
         call Ewald_SelfTerm
       end if

     end if

   end if

   AcAvE_TICR = 0.d0

   if(QCharge) then

     if(QCreation) then
       Lmd = TIlambda(jj)
     else
       Lmd = ( 1.d0 - TIlambda(jj) )
     end if

     j = 0
     do i = NiniTI+1, NfinTI
       j = j + 1
       Charge(i) = Lmd * ChargeOrg(j)
     end do

     if((trim(cCOULOMB) == 'EWALD').or.&
     &  (trim(cCOULOMB) == 'PME')  ) then
       call Ewald_SelfTerm
     end if

   end if

end subroutine SetTICR


!######################################################################
!######################################################################


subroutine TICRSampling(niter,npoint)

use CommonBlocks, only : QMaster
use FEparam
use UnitExParam, only : cvol
use TimeParam, only : Timeps

implicit none

integer :: niter, npoint
real(8) :: CurrEne
integer, save :: samplenumber

   if(niter==0) then
     samplenumber = 0
   end if

   niter = niter + 1

   if(niter > iequil_TICR(npoint)) then

     samplenumber = samplenumber + 1

     call EneTICR(npoint)
     call SumEneTICR

     if(QMaster) then

       AcAvE_TICR = AcAvE_TICR + E_TICR

       CurrEne = AcAvE_TICR / dble(samplenumber)

       write(15,'(f9.4,e16.8)') Timeps, CurrEne*cvol

     end if

     if(samplenumber == isample_TICR(npoint)) then
       write(16,'(f8.5,e16.8)') TIlambda(npoint), &
       &                        AcAvE_TICR / dble(isample_TICR(npoint)) * cvol
       npoint = npoint + 1
       niter  = 0
       call SetTICR(1,npoint)
     end if

   end if

end subroutine TICRSampling


!######################################################################
!######################################################################


subroutine EneTICR(npoint)

use FEparam

implicit none

real(8) :: Lmd, ShftR
integer :: npoint

! ##

   if(QCreation) then
     Lmd = TIlambda(npoint)
   else
     Lmd = 1.d0 - TIlambda(npoint)
   end if

   ShftR = Shift_Param * (1.d0 - Lmd)

! ##

   E_TICR = 0.d0


   if(QLJ) then

     if(QCharge) then

!     --------------------------------
       call Calc_LJ_Charge(Lmd,ShftR)
!     --------------------------------

     else

!     --------------------------
       call Calc_LJ(Lmd,ShftR)
!     --------------------------

     end if

   else if(QCharge) then

!   -----------------------------
     call Calc_Charge
!   -----------------------------

   end if

   if(.not.QCreation) then
     E_TICR = - E_TICR
   end if

end subroutine EneTICR


!######################################################################
!######################################################################


subroutine Calc_LJ_Charge(Lmd,ShftR)

use Numbers, only : N
use CommonBlocks, only : QPathInt, QSwitch, ForceField
use Configuration, only : R
use CommonMPI
use FEparam
use CommonPI
use NonbondParam, only : Charge, Rminh, SgmLJ, EpsLJ
use CellParam, only : H, InvH
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

integer :: i, j, k, ii
integer :: Nas
real(8), dimension(3,N) :: ScR
real(8), dimension(3) :: Rij, Sij
real(8) :: R2, Sgm, Sgm2, Sgm6, Sgm12, Eps
real(8) :: Sc, Sc3, Sc4, Sc6, Sc7, fk1, ek, cf
real(8) :: AA, BB
real(8) :: Dummy
real(8) :: Qfrac
real(8) :: Lmd, ShftR
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   do i = 1 , N

#ifdef PCC
     ScR(1,i) = InvH(1,1) * R(1,i) + InvH(1,2) * R(2,i) + InvH(1,3) * R(3,i)
     ScR(2,i) = InvH(2,1) * R(1,i) + InvH(2,2) * R(2,i) + InvH(2,3) * R(3,i)
     ScR(3,i) = InvH(3,1) * R(1,i) + InvH(3,2) * R(2,i) + InvH(3,3) * R(3,i)
#else
     ScR(:,i) = matmul( InvH , R(:,i) )
#endif

   end do

   if(ForceField(1:5) == 'CHARM') then

     ii = 0

     do i = NiniTI+1, NfinTI

       ii = ii + 1

       Qfrac = ChargeOrg(ii)

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij(:) = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then

             Sgm   = Rminh(i) + Rminh(j)
             Sgm2  = Sgm * Sgm
             Eps   = EpsLJ(i) * EpsLJ(j)
             Sc    = 1.d0 / ( R2 + ShftR )

             Sc3  = Sc  * Sc  * Sc
             Sc6  = Sc3 * Sc3
             Sc4  = Sc3 * Sc
             Sc7  = Sc6 * Sc

             Sgm6 = Sgm2 * Sgm2 * Sgm2
             Sgm12 = Sgm6 * Sgm6

             AA =          Eps * Sgm12
             BB = - 2.d0 * Eps * Sgm6

             ek = 3.d0 * Lmd * Shift_Param * ( 2.d0 * AA * Sc7 + BB * Sc4 ) &
             &  + ( AA * Sc6 + BB * Sc3 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,Dummy,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             cf = Qfrac * Charge(j)

             fk1 = cf / sqrt( R2 )

             E_TICR = E_TICR + ek + fk1

           end if

         end if

       end do

     end do

   else if(ForceField(1:4) == 'OPLS') then

     ii = 0

     do i = NiniTI+1, NfinTI

       ii = ii + 1

       Qfrac = ChargeOrg(ii)

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then

             Eps   = EpsLJ(i) * EpsLJ(j)

             Sgm2  = SgmLJ(i) * SgmLJ(j)
             Sgm6  = Sgm2 * Sgm2 * Sgm2
             Sgm12 = Sgm6 * Sgm6

             Sc  = 1.d0 / ( R2 + ShftR )           !(sigma/r)^2
             Sc3 = Sc * Sc * Sc                 !         ^6
             Sc6 = Sc3 * Sc3
             Sc4 = Sc3 * Sc
             Sc7 = Sc6 * Sc

             AA =   4.d0 * Eps * Sgm12
             BB = - 4.d0 * Eps * Sgm6

             ek = 3.d0 * Lmd * Shift_Param * ( 2.d0 * AA * Sc7 + BB * Sc4 ) &
             &  + ( AA * Sc6 + BB * Sc3 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,Dummy,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             cf = Qfrac * Charge(j)

             fk1 = cf / sqrt( R2 )

             E_TICR = E_TICR + ek + fk1

           end if

         end if

       end do

     end do

   else

     write(*,*) 'ERROR : "TICR" option is valid only for CHARMM or OPLS forcefield'
     call Finalize

   end if

end subroutine Calc_LJ_Charge


! #################################################################
! #################################################################


subroutine Calc_LJ(Lmd,ShftR)

use Numbers, only : N
use CommonBlocks, only : QPathInt, QSwitch, ForceField
use Configuration, only : R
use CommonMPI
use FEparam
use CommonPI
use NonbondParam, only : Rminh, EpsLJ, SgmLJ
use CellParam, only : H, InvH
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

integer :: i, j, k
integer :: Nas
real(8), dimension(3,N) :: ScR
real(8), dimension(3) :: Rij, Sij
real(8) :: R2, Sgm, Sgm2, Sgm6, Sgm12, Eps
real(8) :: Sc, Sc3, Sc4, Sc6, Sc7, ek
real(8) :: AA, BB
real(8) :: Dummy
real(8) :: Lmd, ShftR
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   do i = 1 , N

#ifdef PCC
     ScR(1,i) = InvH(1,1) * R(1,i) + InvH(1,2) * R(2,i) + InvH(1,3) * R(3,i)
     ScR(2,i) = InvH(2,1) * R(1,i) + InvH(2,2) * R(2,i) + InvH(2,3) * R(3,i)
     ScR(3,i) = InvH(3,1) * R(1,i) + InvH(3,2) * R(2,i) + InvH(3,3) * R(3,i)
#else
     ScR(:,i) = matmul( InvH , R(:,i) )
#endif

   end do

   if(ForceField(1:5) == 'CHARM') then

     do i = NiniTI+1, NfinTI

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij(:) = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then

             Sgm   = Rminh(i) + Rminh(j)
             Sgm2  = Sgm * Sgm
             Eps   = EpsLJ(i) * EpsLJ(j)
             Sc    = 1.d0 / ( R2 + ShftR )

             Sc3  = Sc  * Sc  * Sc
             Sc6  = Sc3 * Sc3
             Sc4  = Sc3 * Sc
             Sc7  = Sc6 * Sc

             Sgm6 = Sgm2 * Sgm2 * Sgm2
             Sgm12 = Sgm6 * Sgm6

             AA =          Eps * Sgm12
             BB = - 2.d0 * Eps * Sgm6

             ek = 3.d0 * Lmd * Shift_Param * ( 2.d0 * AA * Sc7 + BB * Sc4 ) &
             &  + ( AA * Sc6 + BB * Sc3 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,Dummy,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             E_TICR = E_TICR + ek

           end if

         end if

       end do

     end do

   else if(ForceField(1:4) == 'OPLS') then

     do i = NiniTI+1, NfinTI

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then

             Eps   = EpsLJ(i) * EpsLJ(j)

             Sgm2  = SgmLJ(i) * SgmLJ(j)
             Sgm6  = Sgm2 * Sgm2 * Sgm2
             Sgm12 = Sgm6 * Sgm6

             Sc  = 1.d0 / ( R2 + ShftR )           !(sigma/r)^2
             Sc3 = Sc * Sc * Sc                 !         ^6
             Sc6 = Sc3 * Sc3
             Sc4 = Sc3 * Sc
             Sc7 = Sc6 * Sc

             AA =   4.d0 * Eps * Sgm12
             BB = - 4.d0 * Eps * Sgm6

             ek = 3.d0 * Lmd * Shift_Param * ( 2.d0 * AA * Sc7 + BB * Sc4 ) &
             &  + ( AA * Sc6 + BB * Sc3 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,Dummy,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             E_TICR = E_TICR + ek

           end if

         end if

       end do

     end do

   else

     write(*,*) 'ERROR : "TICR" option is valid only for CHARMM or OPLS forcefield'
     call Finalize

   end if

end subroutine Calc_LJ


! #################################################################
! #################################################################


subroutine Calc_Charge

use Numbers, only : N
use CommonBlocks, only : QPathInt
use Configuration, only : R
use CommonMPI
use CommonPI
use FEparam
use NonbondParam, only : Charge
use CellParam, only : H, InvH
use CutoffParam, only : Rcutoff2

implicit none

integer :: i, j, k, ii
integer :: Nas
real(8), dimension(3,N) :: ScR
real(8), dimension(3) :: Rij, Sij
real(8) :: R2
real(8) :: fk1, cf
real(8) :: Qfrac
integer :: NProcsTemp, MyRankTemp

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

   Nas = NProcsTemp - MyRankTemp

   do i = 1 , N

#ifdef PCC
     ScR(1,i) = InvH(1,1) * R(1,i) + InvH(1,2) * R(2,i) + InvH(1,3) * R(3,i)
     ScR(2,i) = InvH(2,1) * R(1,i) + InvH(2,2) * R(2,i) + InvH(2,3) * R(3,i)
     ScR(3,i) = InvH(3,1) * R(1,i) + InvH(3,2) * R(2,i) + InvH(3,3) * R(3,i)
#else
     ScR(:,i) = matmul( InvH , R(:,i) )
#endif

   end do

   ii = 0

   do i = NiniTI+1, NfinTI

     ii = ii + 1

     Qfrac = ChargeOrg(ii)

     do j = Nas, N, NProcsTemp

       if(QTICRcalc(j)) then

         Sij(:) = ScR(:,i) - ScR(:,j)
         Rij(:) = R(:,i) - R(:,j)
         call MImageF(Sij,Rij,R2,k)

         if(R2 <= Rcutoff2) then

           cf = Qfrac * Charge(j)

           fk1 = cf / sqrt( R2 )

           E_TICR = E_TICR + fk1

         end if

       end if

     end do

   end do

end subroutine Calc_Charge


!######################################################################
!######################################################################


subroutine Force_TICR(npoint)

use Numbers, only : N
use CommonBlocks, only : QRigidBody, QPathInt, QSwitch, ForceField
use Configuration, only : R
use FEparam
use CommonMPI
use CommonPI
use RBparam, only : NumRB, R_RB, AtomUnitNum, ScG
use EwaldParam, only : Alpha, ar2
use NonbondParam, only : Charge, Rminh, EpsLJ, SgmLJ, Frc_Ersp, Ene_Ersp, &
&  Ene_LJ, Vir_Ersp
use CellParam, only : H, InvH
use CutoffParam, only : Ron2, Rcutoff2, swf1

implicit none

integer :: i, j, k, npoint
integer :: IRB, JRB
integer :: Nas
integer, dimension(3) :: Nij
real(8), dimension(3,N) :: ScR
real(8), dimension(3) :: Rij, Sij
real(8), dimension(3) :: RGij, SGij
real(8), dimension(3) :: FijLJ, FijEL
real(8) :: R2, R1, InvR2, Sgm, Sgm2, Sgm6, Sgm12, Eps
real(8) :: Sc, Sc3, Sc4, Sc6, Sc7, fk1, ek, cf
real(8) :: SR2, SR6, SR12
real(8) :: AA, BB
real(8) :: fk, x, xtm, fk2, fkLJ, ErrorFunc
real(8) :: Lmd, ShftR
integer :: NProcsTemp, MyRankTemp
real(8) :: Error_Function
external Error_Function

   if(QPathInt) then
     MyRankTemp = MyRankPI
     NProcsTemp = NumProcess
   else
     MyRankTemp = MyRank
     NProcsTemp = NProcs
   end if

! ##

   if(QCreation) then
     Lmd = TIlambda(npoint)
   else
     Lmd = 1.d0 - TIlambda(npoint)
   end if

   ShftR = Shift_Param * (1.d0 - Lmd)

! ##

   do i = 1 , N

#ifdef PCC
     ScR(1,i) = InvH(1,1) * R(1,i) + InvH(1,2) * R(2,i) + InvH(1,3) * R(3,i)
     ScR(2,i) = InvH(2,1) * R(1,i) + InvH(2,2) * R(2,i) + InvH(2,3) * R(3,i)
     ScR(3,i) = InvH(3,1) * R(1,i) + InvH(3,2) * R(2,i) + InvH(3,3) * R(3,i)
#else
     ScR(:,i) = matmul( InvH , R(:,i) )
#endif

   end do

   if(QRigidBody) then

     do i = 1 , NumRB

#ifdef PCC
       ScG(1,i) = InvH(1,1) * R_RB(1,i) + InvH(1,2) * R_RB(2,i) + InvH(1,3) * R_RB(3,i)
       ScG(2,i) = InvH(2,1) * R_RB(1,i) + InvH(2,2) * R_RB(2,i) + InvH(2,3) * R_RB(3,i)
       ScG(3,i) = InvH(3,1) * R_RB(1,i) + InvH(3,2) * R_RB(2,i) + InvH(3,3) * R_RB(3,i)
#else
       ScG(:,i) = matmul( InvH , R_RB(:,i) )
#endif

     end do

   end if

   Nas = NProcsTemp - MyRankTemp

if(QLJ) then

   if(Lmd == 0.) Return

   if(ForceField(1:5) == 'CHARM') then

     do i = NiniTI+1, NfinTI

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij(:) = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then

             IRB = AtomUnitNum(i)
             JRB = AtomUnitNum(j)

             SGij(:) = ScG(:,IRB) - ScG(:,JRB) + Nij(:)
#ifdef PCC
             RGij(1) = H(1,1) * SGij(1) + H(1,2) * SGij(2) + H(1,3) * SGij(3)
             RGij(2) = H(2,1) * SGij(1) + H(2,2) * SGij(2) + H(2,3) * SGij(3)
             RGij(3) = H(3,1) * SGij(1) + H(3,2) * SGij(2) + H(3,3) * SGij(3)
#else
             RGij = matmul( H , SGij )
#endif
             Sgm   = Rminh(i) + Rminh(j)
             Sgm2  = Sgm * Sgm
             Eps   = EpsLJ(i) * EpsLJ(j)
             Sc    = 1.d0 / ( R2 + ShftR )

             Sc3  = Sc  * Sc  * Sc
             Sc6  = Sc3 * Sc3
             Sc4  = Sc3 * Sc
             Sc7  = Sc6 * Sc

             Sgm6 = Sgm2 * Sgm2 * Sgm2
             Sgm12 = Sgm6 * Sgm6

             AA =          Eps * Sgm12
             BB = - 2.d0 * Eps * Sgm6

             ek = Lmd * ( AA * Sc6 + BB * Sc3 )

             fkLJ = Lmd * ( 12.d0 * AA * Sc7 + 6.d0 * BB * Sc4 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             FijLJ = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

             Ene_LJ = Ene_LJ + ek

             cf = Charge(i) * Charge(j)

             if(cf /= 0.) then

               InvR2 = 1.d0 / R2

               R1  = sqrt( R2 )
               x   = Alpha * R1

               ErrorFunc = Error_Function(x)

               xtm = -x * x
               fk1 = cf * ErrorFunc / R1
               fk2 = cf * ar2 * exp(xtm)
               fk  = ( fk1 + fk2 ) * InvR2

               FijEL = fk * Rij

               Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
               Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

               Ene_Ersp = Ene_Ersp + fk1

             else

               FijEL = 0.d0

             end if

             Vir_Ersp(1,1) = Vir_Ersp(1,1) + ( FijLJ(1) + FijEL(1) ) * RGij(1)
             Vir_Ersp(1,2) = Vir_Ersp(1,2) + ( FijLJ(1) + FijEL(1) ) * RGij(2)
             Vir_Ersp(1,3) = Vir_Ersp(1,3) + ( FijLJ(1) + FijEL(1) ) * RGij(3)
             Vir_Ersp(2,1) = Vir_Ersp(2,1) + ( FijLJ(2) + FijEL(2) ) * RGij(1)
             Vir_Ersp(2,2) = Vir_Ersp(2,2) + ( FijLJ(2) + FijEL(2) ) * RGij(2)
             Vir_Ersp(2,3) = Vir_Ersp(2,3) + ( FijLJ(2) + FijEL(2) ) * RGij(3)
             Vir_Ersp(3,1) = Vir_Ersp(3,1) + ( FijLJ(3) + FijEL(3) ) * RGij(1)
             Vir_Ersp(3,2) = Vir_Ersp(3,2) + ( FijLJ(3) + FijEL(3) ) * RGij(2)
             Vir_Ersp(3,3) = Vir_Ersp(3,3) + ( FijLJ(3) + FijEL(3) ) * RGij(3)

           end if

         end if

       end do

     end do

   else if(ForceField(1:4) == 'OPLS') then

     do i = NiniTI+1, NfinTI

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then


             IRB = AtomUnitNum(i)
             JRB = AtomUnitNum(j)

             SGij(:) = ScG(:,IRB) - ScG(:,JRB) + Nij(:)
#ifdef PCC
             RGij(1) = H(1,1) * SGij(1) + H(1,2) * SGij(2) + H(1,3) * SGij(3)
             RGij(2) = H(2,1) * SGij(1) + H(2,2) * SGij(2) + H(2,3) * SGij(3)
             RGij(3) = H(3,1) * SGij(1) + H(3,2) * SGij(2) + H(3,3) * SGij(3)
#else
             RGij = matmul( H , SGij )
#endif

             Eps   = EpsLJ(i) * EpsLJ(j)

             Sgm2  = SgmLJ(i) * SgmLJ(j)
             Sgm6  = Sgm2 * Sgm2 * Sgm2
             Sgm12 = Sgm6 * Sgm6

             Sc  = 1.d0 / ( R2 + ShftR )           !(sigma/r)^2
             Sc3 = Sc * Sc * Sc                 !         ^6
             Sc6 = Sc3 * Sc3
             Sc4 = Sc3 * Sc
             Sc7 = Sc6 * Sc

             AA =   4.d0 * Eps * Sgm12
             BB = - 4.d0 * Eps * Sgm6

             ek   = Lmd * ( AA * Sc6 + BB * Sc3 )

             fkLJ = Lmd * ( 12.d0 * AA * Sc7 + 6.d0 * BB * Sc4 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             Ene_LJ = Ene_LJ + ek

             FijLJ = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

             cf = Charge(i) * Charge(j)

             if(cf /= 0.) then

               InvR2 = 1.d0 / R2

               R1  = sqrt( R2 )
               x   = Alpha * R1
               ErrorFunc = Error_Function(x)

               xtm = -x * x
               fk1 = cf * ErrorFunc / R1
               fk2 = cf * ar2 * exp(xtm)
               fk  = ( fk1 + fk2 ) * InvR2

               FijEL = fk * Rij

               Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
               Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

               Ene_Ersp = Ene_Ersp + fk1

             else

               FijEL = 0.d0

             end if

             Vir_Ersp(1,1) = Vir_Ersp(1,1) + ( FijLJ(1) + FijEL(1) ) * RGij(1)
             Vir_Ersp(1,2) = Vir_Ersp(1,2) + ( FijLJ(1) + FijEL(1) ) * RGij(2)
             Vir_Ersp(1,3) = Vir_Ersp(1,3) + ( FijLJ(1) + FijEL(1) ) * RGij(3)
             Vir_Ersp(2,1) = Vir_Ersp(2,1) + ( FijLJ(2) + FijEL(2) ) * RGij(1)
             Vir_Ersp(2,2) = Vir_Ersp(2,2) + ( FijLJ(2) + FijEL(2) ) * RGij(2)
             Vir_Ersp(2,3) = Vir_Ersp(2,3) + ( FijLJ(2) + FijEL(2) ) * RGij(3)
             Vir_Ersp(3,1) = Vir_Ersp(3,1) + ( FijLJ(3) + FijEL(3) ) * RGij(1)
             Vir_Ersp(3,2) = Vir_Ersp(3,2) + ( FijLJ(3) + FijEL(3) ) * RGij(2)
             Vir_Ersp(3,3) = Vir_Ersp(3,3) + ( FijLJ(3) + FijEL(3) ) * RGij(3)

           end if

         end if

       end do

     end do

   else

     write(*,*) 'ERROR : "TICR" option is valid only for CHARMM or OPLS forcefield'
     call Finalize

   end if

else

   if(ForceField(1:5) == 'CHARM') then

     do i = NiniTI+1, NfinTI

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij(:) = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then

             IRB = AtomUnitNum(i)
             JRB = AtomUnitNum(j)

             SGij(:) = ScG(:,IRB) - ScG(:,JRB) + Nij(:)
#ifdef PCC
             RGij(1) = H(1,1) * SGij(1) + H(1,2) * SGij(2) + H(1,3) * SGij(3)
             RGij(2) = H(2,1) * SGij(1) + H(2,2) * SGij(2) + H(2,3) * SGij(3)
             RGij(3) = H(3,1) * SGij(1) + H(3,2) * SGij(2) + H(3,3) * SGij(3)
#else
             RGij = matmul( H , SGij )
#endif
             InvR2 = 1.d0 / R2

             Sgm   = Rminh(i) + Rminh(j)
             Sgm2  = Sgm * Sgm
             Eps   = EpsLJ(i) * EpsLJ(j)

             SR2  = Sgm2 * InvR2                    !(sigma/r)^2
             SR6  = SR2 * SR2 * SR2                 !         ^6
             SR12 = SR6 * SR6                       !         ^12
             fkLJ = Eps * 12.d0 * ( SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
             ek   = Eps * ( SR12 - 2.d0 * SR6 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

             FijLJ = fkLJ * Rij

             Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
             Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

             Ene_LJ = Ene_LJ + ek

             cf = Charge(i) * Charge(j)

             if(cf /= 0.) then

               R1  = sqrt( R2 )
               x   = Alpha * R1

               ErrorFunc = Error_Function(x)

               xtm = -x * x
               fk1 = cf * ErrorFunc / R1
               fk2 = cf * ar2 * exp(xtm)
               fk  = ( fk1 + fk2 ) * InvR2

               FijEL = fk * Rij

               Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
               Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

               Ene_Ersp = Ene_Ersp + fk1

             else

               FijEL = 0.d0

             end if

             Vir_Ersp(1,1) = Vir_Ersp(1,1) + ( FijLJ(1) + FijEL(1) ) * RGij(1)
             Vir_Ersp(1,2) = Vir_Ersp(1,2) + ( FijLJ(1) + FijEL(1) ) * RGij(2)
             Vir_Ersp(1,3) = Vir_Ersp(1,3) + ( FijLJ(1) + FijEL(1) ) * RGij(3)
             Vir_Ersp(2,1) = Vir_Ersp(2,1) + ( FijLJ(2) + FijEL(2) ) * RGij(1)
             Vir_Ersp(2,2) = Vir_Ersp(2,2) + ( FijLJ(2) + FijEL(2) ) * RGij(2)
             Vir_Ersp(2,3) = Vir_Ersp(2,3) + ( FijLJ(2) + FijEL(2) ) * RGij(3)
             Vir_Ersp(3,1) = Vir_Ersp(3,1) + ( FijLJ(3) + FijEL(3) ) * RGij(1)
             Vir_Ersp(3,2) = Vir_Ersp(3,2) + ( FijLJ(3) + FijEL(3) ) * RGij(2)
             Vir_Ersp(3,3) = Vir_Ersp(3,3) + ( FijLJ(3) + FijEL(3) ) * RGij(3)

           end if

         end if

       end do

     end do

   else if(ForceField(1:4) == 'OPLS') then

     do i = NiniTI+1, NfinTI

       do j = Nas, N, NProcsTemp

         if(QTICRcalc(j)) then

           Sij = ScR(:,i) - ScR(:,j)
           Rij(:) = R(:,i) - R(:,j)
           call MImageF(Sij,Rij,R2,k)

           if(R2 <= Rcutoff2) then


             IRB = AtomUnitNum(i)
             JRB = AtomUnitNum(j)

             SGij(:) = ScG(:,IRB) - ScG(:,JRB) + Nij(:)
#ifdef PCC
             RGij(1) = H(1,1) * SGij(1) + H(1,2) * SGij(2) + H(1,3) * SGij(3)
             RGij(2) = H(2,1) * SGij(1) + H(2,2) * SGij(2) + H(2,3) * SGij(3)
             RGij(3) = H(3,1) * SGij(1) + H(3,2) * SGij(2) + H(3,3) * SGij(3)
#else
             RGij = matmul( H , SGij )
#endif

             InvR2 = 1.d0 / R2

             Eps   = EpsLJ(i) * EpsLJ(j)

             if(Eps /= 0.) then

               Sgm2 = SgmLJ(i) * SgmLJ(j)
               SR2  = Sgm2 * InvR2                    !(sigma/r)^2
               SR6  = SR2 * SR2 * SR2                 !         ^6
               SR12 = SR6 * SR6                       !         ^12
               fkLJ = Eps * 24.d0 * ( 2.d0 * SR12 - SR6 ) * InvR2  !4e(12()^12-6()^6)
               ek   = Eps * 4.d0 * ( SR12 - SR6 )

! ----------------------------------------------------------
 if(QSwitch) call SwitchFunc(R2,fkLJ,ek,Ron2,Rcutoff2,swf1)
! ----------------------------------------------------------

               FijLJ = fkLJ * Rij

               Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijLJ
               Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijLJ

               Ene_LJ = Ene_LJ + ek

             else

               FijLJ = 0.d0

             end if

             cf = Charge(i) * Charge(j)

             if(cf /= 0.) then

               InvR2 = 1.d0 / R2

               R1  = sqrt( R2 )
               x   = Alpha * R1
               ErrorFunc = Error_Function(x)

               xtm = -x * x
               fk1 = cf * ErrorFunc / R1
               fk2 = cf * ar2 * exp(xtm)
               fk  = ( fk1 + fk2 ) * InvR2

               FijEL = fk * Rij

               Frc_Ersp(:,i) = Frc_Ersp(:,i) + FijEL
               Frc_Ersp(:,j) = Frc_Ersp(:,j) - FijEL

               Ene_Ersp = Ene_Ersp + fk1

             else

               FijEL = 0.d0

             end if

             Vir_Ersp(1,1) = Vir_Ersp(1,1) + ( FijLJ(1) + FijEL(1) ) * RGij(1)
             Vir_Ersp(1,2) = Vir_Ersp(1,2) + ( FijLJ(1) + FijEL(1) ) * RGij(2)
             Vir_Ersp(1,3) = Vir_Ersp(1,3) + ( FijLJ(1) + FijEL(1) ) * RGij(3)
             Vir_Ersp(2,1) = Vir_Ersp(2,1) + ( FijLJ(2) + FijEL(2) ) * RGij(1)
             Vir_Ersp(2,2) = Vir_Ersp(2,2) + ( FijLJ(2) + FijEL(2) ) * RGij(2)
             Vir_Ersp(2,3) = Vir_Ersp(2,3) + ( FijLJ(2) + FijEL(2) ) * RGij(3)
             Vir_Ersp(3,1) = Vir_Ersp(3,1) + ( FijLJ(3) + FijEL(3) ) * RGij(1)
             Vir_Ersp(3,2) = Vir_Ersp(3,2) + ( FijLJ(3) + FijEL(3) ) * RGij(2)
             Vir_Ersp(3,3) = Vir_Ersp(3,3) + ( FijLJ(3) + FijEL(3) ) * RGij(3)

           end if

         end if

       end do

     end do

   else

     write(*,*) 'ERROR : "TICR" option is valid only for CHARMM or OPLS forcefield'
     call Finalize

   end if

end if

end subroutine Force_TICR
