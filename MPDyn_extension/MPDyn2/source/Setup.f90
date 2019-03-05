! ############################
! ## SUBROUTINE LIST 
! ## -- Setup 
! ## -- ParamFixed 
! ## -- ParamFix_OptC 
! ## -- ParamFix_SHAKE 
! ## -- PrepDyna
! ## -- InfoConf
! ## -- InfoConc
! ## -- SetFmoni
! ############################


!######################################################################
!######################################################################


! *************
! **  Setup  **
! *************

subroutine Setup

use CommonBlocks, only : QMaster, Qstdout, SimMethod, ForceField, &
&   QFmonitor, QPathInt, QRigidBody, QSTWALL, QMP, QInitial, QPBC,&
&   QThermostat, QGeneconf, QPDB
use ParamAnalyze
use CommonHMC
use CellParam, only : CellShape
use TimeParam, only : Timeps

implicit none

! ## random number
   call PreRandomNumber

! ## Simulation Condition
   call Read_Condition(1)
   if(QMaster.and.Qstdout) print *, 'OK:Condition1'

! ## DPD >>
   if(SimMethod == 'DPD') then
     call SetupDP
   end if

! ## File Names
   call Set_IO
   if(QMaster.and.Qstdout) print *, 'OK:Set_IO'

! ## for DynaLib
   if( SimMethod == 'DynaLib' ) then

     call SetDynaLib(1)

! ## 
   else

! ## Potential Parameters
     if( ForceField(1:5) == 'CHARM' ) then
! Charmm Parameters
       call Read_Charmm_Parameter
     else if(ForceField(1:4) == 'OPLS') then
! OPLS Parameters
       call Read_OPLS_Parameter !<< new >> 
     else if(ForceField(1:3) == 'BKS') then
! BKS Parameters
       call Read_BKS_Parameter !<< new >> 
! ForCG
     else if(ForceField(1:2) == 'CG') then
       QGeneconf = .False.
       call Read_CG_Parameter
     end if
     if(QMaster.and.Qstdout) print *, 'OK:Parameter file'

     if( ForceField(1:5) == 'CHARM' ) then
! Charmm Topology File ## not necessary anymore
!##     call Read_Charmm_Topology
     else if((ForceField(1:4) == 'OPLS').or.(ForceField(1:3) == 'BKS')) then
! OPLS Topology File
       call Read_OPLS_Topology !<< new >> 
! ForCG
     else if(ForceField(1:2) == 'CG') then
       call Read_CG_Topology
     end if
     if(QMaster.and.Qstdout) print *, 'OK:Topology file'

! ## DPD 
     if(SimMethod == 'DPD') then
       call SetupDP_IO
       Return
     end if

     call MemoryAllocation
     if(QMaster.and.Qstdout) print *, 'OK:Memory allocation'

! ## Configurations
     call Read_Config(1)
     if(QMaster.and.Qstdout) print *, 'OK:Configuration'

! ## Charmm PSF File
     if( ForceField(1:5) == 'CHARM' ) then
       call Read_PSF
     else if(ForceField(1:4) == 'OPLS') then !<< new >> 
       call MakeTopologicalData
       call AllocateOPLS
     else if(ForceField(1:3) == 'EAM') then
       call Set_EAM
     else if(ForceField(1:3) == 'BKS') then
       call AllocateBKS
     else if(ForceField(1:2) == 'CG') then
       call AllocateCGdata
     end if
     if(QMaster.and.Qstdout) print *, 'OK:Parameter Assign'
     if(QPDB) call Write_CONECT

! ## Path Integral 
     if(QPathInt) then
       call SetupPI
       call Read_Conf_PI
       if(QMaster.and.Qstdout) print *, 'OK:Path Integral'
     end if

     if(SimMethod /= 'Analysis') then

       if(QRigidBody) then
         call RB_data
         call Read_Config(2)
         if(QMaster.and.Qstdout) print *, 'OK:Rigid-Body'
       end if

     end if

   end if

! ## Simulation Condition
   call Read_Condition(2)
   if(QMaster.and.Qstdout) print *, 'OK:Condition2'

   if(QFmonitor) call SetFmoni

   if(QSTWALL) call Set_Steele

   if(QMP) call Set_MembrPot

! >> Analyses
   if(SimMethod == 'Analysis') then
! ## Simulation Condition
     call Read_Condition(4)
     if(QMaster.and.Qstdout) print *, 'OK:Condition4'
     QConfig = .True.
     QVeloc  = .False.
! ForCG 
     CellShape = 3
     close(7)
     Return
   end if
! << Analyses

! ## for DynaLib
   if( SimMethod == 'DynaLib' ) then
     call SetDynaLib(2)
     Return
   end if

! ## Velocities
   if(QInitial) then
     call Gene_Velocity
   else
     call Read_Velocity
   end if
   if(QMaster.and.Qstdout) print *, 'OK:Velocity'

! ## Bath Parameters
   if((QPBC).or.(QThermostat)) then
     call Read_Bath
     if(QMaster.and.Qstdout) print *, 'OK:Bath'
   end if

   if(QInitial) then
     Timeps = 0.d0
     if(SimMethod == 'HMC') TimeMC = 0
   else
     call Read_Time
   end if

   close(7)

   if(ForceField(1:2) == 'CG') call Set_CGcond

! ## Remove the net momentum of the total system
   if(QMaster.and.(.not.QPathInt)) call Elim_CellMove

   if(QMaster.and.Qstdout) print *, '## Finalize Initial Setup! ##'

end subroutine Setup


!#####################################################################
!#####################################################################


! ************************************
! ** Calculation of the initial     **
! ** temperature and the parameters **
! ** used for Ewald sum or cut-off  **
! ** correction term                **
! ************************************

subroutine PrepDyna

use CommonBlocks, only : QMaster, SimMethod, QPBC, ForceField, QCorrectCutoff, &
&   cCOULOMB, QPathInt

implicit NONE

   if(SimMethod == 'MD') then

     if(QMaster) call CalcTemp

     if(QPBC.and.(ForceField(1:3)/='EAM')) then

       if(QCorrectCutoff) then
         call CorrectPotential        !  correction for the potential cutoff
       end if

       if( trim(cCOULOMB) == 'EWALD' ) then
         call ErrorFuncList           !  make a table of error-function
         call RecLatticeList          !  define reciprocal lattice
         call Ewald_SelfTerm
       else if( trim(cCOULOMB) == 'PME' ) then
         call ErrorFuncList           !  make a table of error-function
         call PMEwald_pre
         call Ewald_SelfTerm
       end if

     end if

   else if(SimMethod == 'HMC') then

     if(QPBC.and.(ForceField(1:3)/='EAM')) then

       if(QCorrectCutoff) then
         call CorrectPotential        !  correction for the potential cutoff
       end if

       if( trim(cCOULOMB) == 'EWALD' ) then
         call ErrorFuncList           !  make a table of error-function
         call RecLatticeList          !  define reciprocal lattice
         call Ewald_SelfTerm
       else if( trim(cCOULOMB) == 'PME' ) then
         call ErrorFuncList           !  make a table of error-function
         call PMEwald_pre
         call Ewald_SelfTerm
       end if

     end if

   else if(QPathInt) then

     call Synch(0)

! ## for parallel
     call Assign_Process

     call CalcTempPI

     if(QPBC.and.(ForceField(1:3)/='EAM')) then

       if(QCorrectCutoff) then
         call CorrectPotential        !  correction for the potential cutoff
       end if

       if( trim(cCOULOMB) == 'EWALD' ) then
         call ErrorFuncList           !  make a table of error-function
         call RecLatticeList          !  define reciprocal lattice
         call Ewald_SelfTerm
       else if( trim(cCOULOMB) == 'PME' ) then
         call ErrorFuncList           !  make a table of error-function
         call PMEwald_pre
         call Ewald_SelfTerm
       end if

     end if

   else if(SimMethod == 'DPD' ) then

     if(QMaster) call CalcTempDP

   else if((SimMethod == 'MM' ).or.(SimMethod == 'NMA')) then

     if(QPBC.and.(ForceField(1:3)/='EAM')) then

       if(QCorrectCutoff) then
         call CorrectPotential        !  correction for the potential cutoff
       end if

       if( trim(cCOULOMB) == 'EWALD' ) then
         call ErrorFuncList           !  make a table of error-function
         call RecLatticeList          !  define reciprocal lattice
         call Ewald_SelfTerm
       else if( trim(cCOULOMB) == 'PME' ) then
         call ErrorFuncList           !  make a table of error-function
         call PMEwald_pre
         call Ewald_SelfTerm
       end if

     end if

   end if

   call TransCellList

end subroutine PrepDyna


!#####################################################################
!#####################################################################


subroutine InfoConf

use Numbers, only : N
use CommonBlocks, only : SimMethod, QSHAKE, QMinimization, QMaster, &
&   Qstdout, QAveTh, QInitial
use CommonDPD
use SHAKEparam, only : R_o
use CommonPI, only : QBead
use CommonMPI

implicit NONE

   if( SimMethod == 'MD') then

     call Synch(0)

! ## for parallel
     call AllocPara

     if(QSHAKE) then

       allocate( R_o(3,N) )

       if((QMinimization).and.(QMaster)) then

         if(Qstdout) write( 6,'(2x,a,a/)') '## Energy Minimization is skipped, ', &
                     &                     'when the SHAKE method is selected. ##'
#ifndef BMONI
         write(11,'(2x,a,a/)') '## Energy Minimization is skipped, ', &
         &                     'when the SHAKE method is selected. ##'
#endif
       end if

     else

       if(QMinimization) then

         allocate( R_o(3,N) )
         call EnergyMinimize

       end if

     end if

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start MD time evolution >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start MD time evolution >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start MD time evolution >##############   '
       end if

     end if

     call Synch(1)

   else if(SimMethod == 'HMC') then

     call Synch(0)

! ## for parallel
     call AllocPara

     if(QMinimization) then

       allocate( R_o(3,N) )
       call EnergyMinimize

     end if

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start Hybrid MC iteration >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start Hybrid MC iteration >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start Hybrid MC iteration >##############   '
       end if

     end if

     call Synch(1)

   else if(SimMethod == 'PIMD') then

     call Synch(1)

! ## for parallel
     if(QBead) call AllocPara

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start Path Integral MD >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start Path Integral MD >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start Path Integral MD >##############   '
       end if

     end if

     call Synch(1)

   else if(SimMethod == 'CMD') then

     call Synch(1)

! ## for parallel
     if(QBead) call AllocPara

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start Centroid MD >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start Centroid MD >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start Centroid MD >##############   '
       end if

     end if

     call Synch(1)

   else if(SimMethod == 'Analysis') then

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start Analysis >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start Analysis >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start Analysis >##############   '
       end if

     end if

     call Synch(0)

   else if(SimMethod == 'DPD' ) then

     if(QInitial.and.(TypeNum(N)==1).and.(a(1,1)/=0.)) then
       call Energy_Minimization_DPD
     end if

     if(NProcs /= 1) then
       write( 6,*) 'ERROR : no parallel calculation can be done for DPD in this version!'
#ifndef BMONI
       write(11,*) 'ERROR : no parallel calculation can be done for DPD in this version!'
#endif
       call Finalize
     end if

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start DPD time evolution >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start DPD time evolution >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start DPD time evolution >##############   '
       end if

     end if

   else if(SimMethod == 'MM') then

     call Synch(0)

! ## for parallel
     call AllocPara

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start MM Energy calculation >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start MM Energy calculation >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start MM Energy calculation >##############   '
       end if

     end if

     call Synch(1)

   else if(SimMethod == 'NMA') then

     call Synch(0)

! ## for parallel
     call AllocPara

     if(QMaster) then

#ifndef BMONI
       write(11,'(a/)')  '  #############< Start Normal Mode Analysis >##############   '
       if(QAveTh) then
         write(12,'(a/)')  '  #############< Start Normal Mode Analysis >##############   '
       end if
#endif
       if(Qstdout) then
         write( 6,'(a/)')  '  #############< Start Normal Mode Analysis >##############   '
       end if

     end if

     call Synch(1)

   end if

end subroutine InfoConf


!#####################################################################
!#####################################################################


subroutine InfoConc(start_time,end_time)

use CommonBlocks, only : QMaster, Qstdout, QAveTh, SimMethod

implicit none

character(len=17) :: start_time, end_time
integer, dimension(3), parameter :: IO = (/11,12,6/)
integer :: i, ii

   if(QMaster) then

     ii = 2
     if(Qstdout) ii = 3

     do i = 1, ii

#ifdef BMONI
       if(i==1.or.i==2) cycle
#endif

       if((.not.QAveTh).and.(i==2)) cycle

       if(SimMethod == 'MD') then

         write(IO(i),'(a)') '-------------------------------'
         write(IO(i),'(a)') ' MD simulation ended normally! '
         write(IO(i),'(a)') '-------------------------------'

       else if(SimMethod == 'HMC') then

         write(IO(i),'(a)') '--------------------------------'
         write(IO(i),'(a)') ' HMC simulation ended normally! '
         write(IO(i),'(a)') '--------------------------------'

       else if(SimMethod == 'PIMD') then

         write(IO(i),'(a)') '---------------------------------'
         write(IO(i),'(a)') ' PIMD simulation ended normally! '
         write(IO(i),'(a)') '---------------------------------'

       else if(SimMethod == 'CMD') then

         write(IO(i),'(a)') '--------------------------------'
         write(IO(i),'(a)') ' CMD simulation ended normally! '
         write(IO(i),'(a)') '--------------------------------'

       else if(SimMethod == 'Analysis') then

         write(IO(i),'(a)') '--------------------------'
         write(IO(i),'(a)') ' Analysis ended normally! '
         write(IO(i),'(a)') '--------------------------'

       else if(SimMethod == 'DPD') then

         write(IO(i),'(a)') '--------------------------------'
         write(IO(i),'(a)') ' DPD simulation ended normally! '
         write(IO(i),'(a)') '--------------------------------'

       else if(SimMethod == 'MM') then

         write(IO(i),'(a)') '--------------------------------'
         write(IO(i),'(a)') ' MM calculation ended normally! '
         write(IO(i),'(a)') '--------------------------------'

       end if

       write(IO(i),'(a23,a17)') 'Program was started at ',start_time
       write(IO(i),'(a23,a17)') '            ended   at ',end_time

     end do

     close(11)
     if(QAveTh) close(12)

   end if

end subroutine InfoConc


! >> F monitor ##

!#####################################################################
!#####################################################################


subroutine SetFmoni

use Numbers, only : N
use CommonBlocks, only : Job_name
use F_monitor, only : NfrcF, IfrcF, force_file, confg_file, Fint, Fext

implicit NONE

   NfrcF = 1
   IfrcF = 0
   write(force_file,'(a,a,i1)') trim(adjustl(Job_name)),'.FRC000',NfrcF
   write(confg_file,'(a,a,i1)') trim(adjustl(Job_name)),'.CRD000',NfrcF

   allocate( Fint(3,N) )
   allocate( Fext(3,N) )

   call ListFmoni

end subroutine SetFmoni
! << F monitor ##
