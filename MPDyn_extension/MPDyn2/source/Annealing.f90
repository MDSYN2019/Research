! ############################
! ## SUBROUTINE LIST 
! ## -- PreAnneal 
! ## -- Annealing 
! ############################


!######################################################################
!######################################################################


subroutine PreAnneal

use SimAnneal
use CommonBlocks, only : QMaster
use BathParam, only : kT, gkT

implicit none

   Nswitch = nint(switch_time)
   if(QMaster) open(50,file='AnnealTemp.data',status='unknown')
   kT_0 = kT
   gkT_0 = gkT

end subroutine PreAnneal


!######################################################################
!######################################################################


subroutine Annealing(istep)

use SimAnneal
use CommonBlocks, only : QMaster
use TimeParam, only : Timeps
use BathParam, only : kT, gkT, Temp_o
use UnitExParam, only : kb

implicit none

integer :: istep

   if(cAnnealMethod=='LINEAR') then
     TempSet = factA * istep + Temp_o
   else if(cAnnealMethod=='HYBRID') then
     if(istep <= Nswitch) then
       TempSet = factA * istep + Temp_o
     else
       TempSet = factB / dble(istep) + factC
     end if
   end if

   kT  =  kT_0 * TempSet / Temp_o
   gkT = gkT_0 * TempSet / Temp_o

   if(QMaster.and.mod(istep,50)==0) write(50,'(f9.4,f7.2)') Timeps, kT/kb

end subroutine Annealing
