

!######################################################################
!######################################################################


subroutine MMcalc

use Numbers, only : N
use CommonBlocks, only : QPBC, QCorrectCutoff, QMacro
use EAM_param, only : Vir_EAM
use UnitExParam, only : kb, cvol
use EwaldParam, only : Vir_Eksp
use OptConstraintParam, only : Vir_OptC
use NonbondParam, only : Vir_Ersp, Vir_NBshrt, Vir_NBlong
use BondedParam, only : Vir_Bond, Vir_Angle, Vir_UB, Vir_Dihed, Vir_Impro
use CellParam, only : Volume
use TailCorrect
use AtomParam, only : Mass
use ThermoData, only : Virial
use CGball, only : NumSphere, FSphRs

implicit none

integer :: i, Nall
real(8), dimension(3,N) :: F1, F2, F3, FF

   Nall = 3 * N

! -------------------------------
   if(QPBC) then
     call GetForce(0,0)
   else
     call GetForceIso(0)
   end if
! -------------------------------

   call GetAcc(F1,1)
   call GetAcc(F2,2)
   call GetAcc(F3,3)

   FF = F1 + F2 + F3

   do i = 1 , N
     FF(:,i) = FF(:,i) * Mass(i)
   end do

   call SumFrc( FF )

   if(QPBC) then

     Virial = Vir_Bond  + Vir_Angle + Vir_UB                &
     &      + Vir_Dihed + Vir_Impro + Vir_Ersp + Vir_NBshrt &
     &      + Vir_Eksp  + Vir_OptC  + Vir_EAM  + Vir_NBlong

     if(QCorrectCutoff) then

       Virial_co = 0.d0

       Virial_co(1,1) = CorrectV / (3.d0*Volume)
       Virial_co(2,2) = Virial_co(1,1)
       Virial_co(3,3) = Virial_co(1,1)

       Ene_LJ_co = CorrectE / Volume

     end if

     call SumVir( Virial )

     if(QCorrectCutoff) then

       Virial = Virial - Virial_co

     end if

     call Print_Energy_MM(FF)

   else

     call Print_Energy_MM_iso(FF)

   end if

!  #ifdef CHECKF
   do i = 1, N
     write(100,*) FF(:,i)/kb
   end do
!  #endif

   if(QMacro.and.NumSphere==1) then
     write( 6,'(/a)') '  --------- dU/dRadius ---------  '
     write(11,'(/a)') '  --------- dU/dRadius ---------  '
     do i = 1, NumSphere
       write( 6,'(5x,i2,a,e16.8,a)') i,' dU/dRad = ',FSphRs(i)*cvol,' [kcal/mol/A] '
       write( 6,'(5x,i2,a,e16.8,a)') i,' dU/dRad = ',FSphRs(i)/kb,' [K/A] '
       write(11,'(5x,i2,a,e16.8,a)') i,' dU/dRad = ',FSphRs(i)*cvol,' [kcal/mol/A] '
       write(11,'(5x,i2,a,e16.8,a)') i,' dU/dRad = ',FSphRs(i)/kb,' [K/A] '
     end do
     write( 6,'(/)')
     write(11,'(/)')
   end if

end subroutine MMcalc
