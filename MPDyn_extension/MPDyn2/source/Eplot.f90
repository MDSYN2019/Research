Program Eplot

use Forplot

implicit none

   call Datainputs

   if(Method(1:2)=='MD') then

     if((Ensemble=='NPT').or. &
     &  (Ensemble=='NtT').or. &
     &  (Ensemble=='NPH').or. &
     &  (Ensemble=='NtH')) then
       call MD_NPT
     else if((Ensemble=='NVT').or. &
     &       (Ensemble=='NVE')) then
       call MD_NVT
     else if((Ensemble(1:2)=='NT').or. &
     &       (Ensemble(1:2)=='NE')) then
       call MD_iso
     end if

   else if(Method(1:3)=='HMC') then

     if(Ensemble(1:3)=='NPT') then
       call HMC_NPT
     else if(Ensemble(1:3)=='NVT') then
       call HMC_NVT
     end if

   else if((Method(1:4)=='PIMD').or.(Method(1:3)=='CMD')) then

     if(Ensemble(1:3)=='NPT') then
       call PIMD_NPT
     else if(Ensemble(1:3)=='NVT') then
       call PIMD_NVT
     else if(Ensemble(1:2)=='NT') then
       call PIMD_iso
     end if

   end if

end Program Eplot
