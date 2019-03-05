! ############################
! ## SUBROUTINE LIST 
! ## -- Read_Condition 
! ## -- Read_Script 
! ## -- Read_Jobname 
! ## -- Read_Simulation_Method 
! ## -- Read_FF 
! ## -- Read_Conf_Flag 
! ## -- Read_System_Cond 
! ## -- Read_PBC_Cond 
! ## -- DATA_Cond 
! ## -- Read_Emin_Cond 
! ## -- Read_Coulomb_Cond 
! ## -- Read_Thermostat_Cond 
! ## -- Read_Annealing 
! ## -- Read_Barostat_Cond 
! ## -- Read_Constraint_Cond 
! ## -- Read_Temperature 
! ## -- Read_StepNumber 
! ## -- Read_Time_Step 
! ## -- Read_Bath_Step 
! ## -- Read_Bead_Length 
! ## -- Read_CMD_param 
! ## -- Mic1 
! ## -- Read_FixedAtom_Cond 
! ## -- Read_Constraint_Mol 
! ## -- Read_OptionalBond_Cond 
! ## -- Read_SPCF_Cond 
! ## -- Read_OptionalConstraint 
! ## -- Read_Jarzynski 
! ## -- Read_wall 
! ## -- Read_Hill 
! ## -- Read_F_monitor 
! ## -- Read_TICR 
! ## -- Read_HeatCapCond 
! ## -- Read_VolScale 
! ## -- Calc_Number 
! ## -- Calc_BathMass 
! ## -- DPD_Condition 
! ## -- ErrorMessage 
! ############################


!######################################################################
!######################################################################


! **********************************************
! **  Conditions : set up simulation details  **
! **********************************************

subroutine Read_Condition(NFL)

use CommonBlocks, only : SimMethod, QPathInt

implicit none

integer :: NFL
logical :: Flag

! ======================================================================
! ======================================================================

   if( NFL == 1 ) then

! ## Save Job_script on memory
     call Read_Script

! ## Job_name --------------------------------------------------------
     call Read_Jobname

! ## SimMethod -------------------------------------------------------
     call Read_Simulation_Method

! ## PBC for Molecular Simulation, not for DPD------------------------
     call Read_PBC_Cond

! ## Initial conf ----------------------------------------------------
     call Read_Conf_Flag

! ## SYSTEM ----------------------------------------------------------
     call Read_System_Cond

! ## 
     if(SimMethod /= 'DynaLib') then

! ## Force Field -----------------------------------------------------
     call Read_FF

! ## ENERGY MINIMIZATION ---------------------------------------------
     call Read_Emin_Cond

! ## COULOMB ---------------------------------------------------------
     call Read_Coulomb_Cond

     end if

! ## DATA_PER_FILE----------------------------------------------------
! ## DATA_FORM_TRJ ---------------------------------------------------
     call DATA_Cond

! ## CALCULATE AVERAGE OF THERMODYNAMIC QUANTITIES -------------------
     call Read_AvThermo

! ## MD or PIMD or CMD -----------------------------------------------
! ## TEMPERATURE -----------------------------------------------------
     if((SimMethod == 'MD').or.(SimMethod == 'MM').or.(QPathInt).or.&
     &  (SimMethod == 'DynaLib')) then
       call Read_Thermostat_Cond(Flag)
     end if

! ## MD or PIMD or CMD or HMC------------------------------------------
! ## PRESSURE --------------------------------------------------------
     if((SimMethod == 'MD').or.(QPathInt).or.(SimMethod == 'HMC').or. &
     &  (SimMethod == 'DynaLib').or.(SimMethod == 'MM')) then
       call Read_Barostat_Cond(Flag)
     end if

! ## MD & HMC---------------------------------------------------------
! ## CONSTRAINT ------------------------------------------------------
     if((SimMethod == 'MD').or.(SimMethod == 'HMC').or. &
     &  (SimMethod == 'DynaLib')) then
       call Read_Constraint_Cond
     end if

! ## DPD -------------------------------------------------------------
! ## HMC -------------------------------------------------------------
! ## TEMPERATURE -----------------------------------------------------
     call Read_Temperature

! ## STEPS -----------------------------------------------------------
     call Read_StepNumber

! ## MD --------------------------------------------------------------
! ## TIME_STEP -------------------------------------------------------
! ## HMC -------------------------------------------------------------
! ## STEP and Random Velocity ----------------------------------------
     call Read_Time_Step

! ## MD or CMD or PIMD ----------------------------------------------
! ## BATH_INT -------------------------------------------------------
     call Read_Bath_Step

! ## PIMD or CMD ----------------------------------------------------
! ## Bead -----------------------------------------------------------
     call Read_Bead_Length
     call Read_CMD_param

! ## MD & HMC ----------------------------------------------------------
     if(SimMethod == 'MD') then
       call Read_Annealing
     end if

! ## Macroparticle
     call Read_Macroparticle

! ## Wall
     call Read_wall

! ## Cylinder
     call Read_Cylinder

     call Mic1

#ifdef EnergyRep
     call Read_EnergyRep
#endif

! ======================================================================
! ======================================================================

   else if( NFL == 2 ) then

! ## FIXED_ATOM ------------------------------------------------------
     call Read_FixedAtom_Cond

     call Read_FixedCOM_Cond

! ## OPTIONALBOND ----------------------------------------------------
     call Read_OptionalBond_Cond

! ## for SPCF 
     call Read_SPCF_Cond

! ## OPTIONALCONSTRAINT
     call Read_OptionalConstraint

! ## CONSTRAINT_MOLECULE
     call Read_Constraint_Mol

! ## Jarzynski
     call Read_Jarzynski

! ## Implicit membrane potential
     call Read_Hill

! ## Controling Moment of Inertia of molecular assembly
     call Read_ContInertia

! ## F_MONITOR
     call Read_F_monitor

! ## E_flux
     call Read_Eflux_mon

! ## ElecticField
     call Read_Efield

! ## TICR
     call Read_TICR

! ## Heat Capacity Estimator
     call Read_HeatCapCond

! ## Volume Scaling
     call Read_VolScale

! ----------------------
! ## degrees of freedom
! ----------------------
     call Calc_Number

! ------------------------
! ## mass of external bath
! ------------------------
     if(SimMethod /= 'Analysis') then
       call Calc_BathMass
       call ErrorMessage
     end if

! ======================================================================
! ======================================================================

   else if( NFL == 3 ) then

     call DPD_Condition
     call ErrorMessage

! ## ANALYSIS --------------------------------------------------------

! ======================================================================
! ======================================================================

   else if( NFL == 4 ) then

! ## ANALYZ_FILE -----------------------------------------------------
     call Read_Analyz_filename

! ## CWHAT -----------------------------------------------------------
     call Read_Analyz_Cond

     call ErrorMessage

   end if

end subroutine Read_Condition


! #####################################################################
! #####################################################################


subroutine Read_Script

use ScriptData
use CommonBlocks, only : QMaster

implicit none


integer :: i
integer :: MaxOption
integer :: MaxScript
character(len=80) :: String, String1
integer :: NumCh, Count
character(len=1) :: Ch
character(len=20), dimension(:), allocatable :: TmpOPTIONS
character(len=80), dimension(:), allocatable :: TmpScript
integer, dimension(:), allocatable :: TmpIst
logical, dimension(:), allocatable :: TmpReadFlag
integer :: eofile

   MaxOption = 200
   MaxScript = 2000

   allocate(OPTIONS(MaxOption))
   allocate(Script(MaxScript))
   allocate(Ist(MaxOption))
   allocate(ReadFlag(MaxOption))

   open(2,file='input.data',status='old')

   Count = 0
   NumOption = 0

   mfile = 6
   Qcheck = .False.

   do

     read(2,'(a)',iostat=eofile) String1

     if(eofile == -1) exit

     String = adjustl(String1)

     if(String(1:5)=='<end>') exit

     if((String(1:1) == '!').or. &
     &  (String(1:1) == '#').or. &
     &  (String(1:1) == ' ')) cycle

     if((String(1:5)=='CHECK').and.(Count==0)) then
       mfile = 99
       Qcheck = .True.
       cycle
     end if

     Count = Count + 1

     if(Count > MaxScript) then
       allocate(TmpScript(MaxScript))
       TmpScript(:) = Script(:)
       deallocate(Script)
       allocate(Script(MaxScript+1000))
       Script(1:MaxScript) = TmpScript(1:MaxScript)
       deallocate(TmpScript)
       MaxScript = MaxScript + 1000
     end if

     NumCh = 0

id0: do i = 1, 80

       Ch = String(i:i)

       if(Ch == '!') exit id0

       NumCh = NumCh + 1

     end do id0

     write(Script(Count),'(a)') adjustl(String(1:NumCh))

     if(String(1:2)=='>>') then

       NumOption = NumOption + 1

       if(NumOption > MaxOption) then
         allocate(TmpOPTIONS(MaxOption))
         allocate(TmpIst(MaxOption))
         allocate(TmpReadFlag(MaxOption))
         TmpOPTIONS(:) = OPTIONS(:)
         TmpIst(:) = Ist(:)
         TmpReadFlag(:) = ReadFlag(:)
         deallocate(OPTIONS,Ist,ReadFlag)
         allocate(OPTIONS(MaxOption+100))
         allocate(Ist(MaxOption+100))
         allocate(ReadFlag(MaxOption+100))
         OPTIONS(1:MaxOption) = TmpOPTIONS(1:MaxOption)
         Ist(1:MaxOption) = TmpIst(1:MaxOption)
         ReadFlag(1:MaxOption) = TmpReadFlag(1:MaxOption)
         deallocate(TmpOPTIONS,TmpIst,TmpReadFlag)
         MaxOption = MaxOption + 100
       end if

       Ist(NumOption) = Count

       String1 = adjustl(String(3:NumCh))

       write(OPTIONS(NumOption),'(a)') trim(String1)

     end if

   end do

   close(2)

   NumScript = Count

   ReadFlag = .False.

end subroutine Read_Script


! #####################################################################
! #####################################################################


subroutine Read_Jobname

use ScriptData
use CommonBlocks, only : QMaster, Job_name
use IOparam, only : DirectoryName

implicit none

integer :: i
character(len=80) :: Check_file
character(len=17) :: ptime
logical :: JobFlag
integer :: icharc, IJob

   JobFlag = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='JOB_NAME') then

       JobFlag = .True.
       IJob = Ist(i) + 1

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(JobFlag) then
     read(Script(IJob),*) Job_name
   else
     call JNAME(ptime)
     write(Job_name,'(a)') trim(ptime)
   end if

! ## for check !

   Check_file = 'CHECK_INPUTS.'//trim(Job_name)

   if(QMaster) then

     if(Qcheck) then

       open(mfile,file=trim(Check_file))

       write(mfile,'(a/)') ' ************** JOB SCRIPTS ************** '
       do i = 1, NumScript
         write(mfile,'(i4,a,2x,a)') i,':',trim(Script(i))
       end do
       write(mfile,'(/a/)') ' ***************** OPTIONS ***************** '
       do i = 1, NumOption
         write(mfile,'(i4,a,2x,a)') i,'.',trim(OPTIONS(i))
       end do
       write(mfile,'(/a)') ' ******************** CONDITIONS ********************** '
       write(mfile,'(a )') ' ***  MPDyn has interpreted the script as follows: **** '
       write(mfile,'(a/)') ' ****************************************************** '
       write(mfile,'(a,a)') 'JOB NAME    = ',trim(Job_name)

     end if

   end if

   DirectoryName = './'

   do i = 1, NumOption

     if(OPTIONS(i)=='DIR_NAME') then

       write(DirectoryName,*) trim(Script(Ist(i)+1))

       ReadFlag(i) = .True.

       exit

     end if

   end do

   icharc = len_trim(DirectoryName)
   if(DirectoryName(icharc:icharc)/='/') then
     DirectoryName(icharc+1:icharc+1)='/'
   end if

   if(QMaster.and.Qcheck) then
     write(mfile,'(a,a)') 'Directory   = ',trim(DirectoryName)
   end if

end subroutine Read_Jobname


! #####################################################################
! #####################################################################


subroutine Read_Simulation_Method

use ScriptData
use CommonBlocks, only : QMaster, SimMethod, QPathInt, QCellList

implicit none

integer :: i
logical :: Flag
integer :: IJob

   Flag = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='METHOD') then

       IJob = Ist(i) + 1
       read(Script(IJob),*) SimMethod

       ReadFlag(i) = .True.

       Flag = .True.

       exit

     end if

   end do

   if(.not.Flag) then
     if(QMaster) write(*,*) 'ERROR : METHOD must be defined in (input.data)'
     call Finalize
   end if

   if(QMaster.and.Qcheck) then
     write(mfile,'(a,a)') 'METHOD      = ',trim(SimMethod)
   end if

   if((SimMethod == 'PIMD').or.(SimMethod == 'CMD')) then
     QPathInt = .True.
   else
     QPathInt = .False.
   end if

   QCellList = .False.

end subroutine Read_Simulation_Method


! #####################################################################
! #####################################################################


subroutine Read_FF

use ScriptData
use CommonBlocks, only : QMaster, ForceField, QSwitch, QCorrectCutoff
use IOparam, only : PSF_file, Topology_file, Parameter_file
#ifdef MEAM
use EAM_param
#endif

implicit none

integer :: i, ii
logical :: psf_on, top_on, par_on
logical :: Flag

   Flag = .False.
   psf_on = .False.
   top_on = .False.
   par_on = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='FORCEFIELD') then

       ii = Ist(i)

lig1:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lig1

         if(Script(ii)(1:3) == 'FF=') then
           read(Script(ii)(4:),*) ForceField
           Flag = .True.
         end if

       end do lig1

       ii = Ist(i)

lig2:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lig2

         if(Script(ii)(1:4) == 'PSF=') then
           read(Script(ii)(5:),*) PSF_file
           psf_on = .True.
         end if

         if(Script(ii)(1:8) == 'TOPFILE=') then
           read(Script(ii)(9:),*) Topology_file
           top_on = .True.
         end if

         if(Script(ii)(1:8) == 'PARFILE=') then
           read(Script(ii)(9:),*) Parameter_file
           par_on = .True.
         end if

#ifdef MEAM
         if(Script(ii)(1:5) == 'MEAM=') then
           read(Script(ii)(6:),*) NumEAM
         end if
#endif

       end do lig2

       ReadFlag(i) = .True.

       if(ForceField(1:5) == 'CHARM') then
         if(.not.par_on) write(Parameter_file,*) 'param/par_charmm27.prm'
         if(.not.top_on) write(Topology_file,*) 'param/top_charmm27.prm'
       else if(ForceField(1:4) == 'OPLS') then
         if(.not.par_on) write(Parameter_file,*) 'param/par_oplsaa.prm'
         if(.not.top_on) write(Topology_file,*) 'param/top_oplsaa.prm'
       else if(ForceField == 'BKS') then
         if(.not.par_on) write(Parameter_file,*) 'param/par_BKS.prm'
         if(.not.top_on) write(Topology_file,*) 'param/top_BKS.prm'
       else if(ForceField == 'CG') then
         if(.not.par_on) write(Parameter_file,*) '~/CGparam/par_CG.prm'
         if(.not.top_on) write(Topology_file,*) '~/CGparam/top_CG.prm'
       else if(ForceField(1:1) == 'F') then
         if(.not.par_on) write(Parameter_file,*) 'param/interaction.param'
       end if

       exit

     end if

   end do

   if(ForceField(1:5) == 'CHARM') then
     if(.not.QdescSW) then
       QSwitch = .True.
       QCorrectCutoff = .False.
     end if
   end if

   if(.not.Flag) then
     if(QMaster) write(*,'(a)') 'ERROR : FORCEFIELD must be defined in (input.data)'
     call Finalize
   end if

   if((ForceField(1:5)=='CHARM').and.(.not.psf_on)) then
     if(QMaster) write(*,'(a)') 'ERROR : the PSF file is needed for the CHARMM force field'
     call Finalize
   end if

   if(QMaster.and.Qcheck) then
     write(mfile,'(a,a)') 'FORCE FIELD = ',trim(ForceField)
   end if

end subroutine Read_FF


! #####################################################################
! #####################################################################


subroutine Read_Conf_Flag

use ScriptData
use CommonBlocks, only : QMaster, QInitial, QResForm, QPDB
use TimeParam, only : irs

implicit none

integer :: i, ii
character(len=11) :: cResForm
character(len=3) :: cPDB

   QInitial=.False.
   write(cResForm,'(a)') 'unformatted'
   irs = 1000
   QPDB = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='STATUS') then

       ii = Ist(i)

ppp1:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit ppp1

         if(Script(ii)(1:7) =='Initial') then
           QInitial=.True.
         else if(Script(ii)(1:7) =='Restart') then
           QInitial=.False.
         else if(Script(ii)(1:5) == 'FORM=') then
           read(Script(ii)(6:),*) cResForm
         else if(Script(ii)(1:7) == 'BACKUP=') then
           read(Script(ii)(8:),*) irs
         else if(Script(ii)(1:4) == 'PDB=') then
           read(Script(ii)(5:),*) cPDB
           if(cPDB(1:2)=='ON') then
             QPDB=.True.
           end if
         else
           if(QMaster) write(*,*) 'ERROR : Wrong flag on STATUS'
           if(QMaster) write(*,*) 'reconsider the line : '
           if(QMaster) write(*,*) Script(ii)
           call Finalize
         end if

       end do ppp1

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(trim(adjustl(cResForm))=='formatted') then
     QResForm = .True.
   else if(trim(adjustl(cResForm))=='unformatted') then
     QResForm = .False.
   else
     if(QMaster) write(*,*) 'ERROR in the defined form of the restart file'
     call Finalize
   end if

   if(QMaster.and.Qcheck) then

     if(QInitial) then
       write(mfile,'(a  )') 'STATUS      = Initial'
       write(mfile,'(a,a)') 'FILE FORM   =', cResForm
     else
       write(mfile,'(a  )') 'STATUS      = Restart'
       write(mfile,'(a,a)') 'FILE FORM   =', cResForm
     end if
   end if

end subroutine Read_Conf_Flag


! #####################################################################
! #####################################################################


subroutine Read_System_Cond

use ScriptData
use CommonBlocks, only : QMaster, PolyFlag, SimMethod, Job_name
use Numbers, only : NumSpec, NumMol, NumAtm, NumMer
use AtomParam, only : MolName

implicit none

integer :: i, ii, j, jj, k
character(len=80) :: String1
logical, dimension(10) :: CheckF
logical :: Flag

   Flag = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='SYSTEM') then

       ii = Ist(i) + 1
       jj = Ist(i+1) - 1

       Flag = .True.

       CheckF = .False.

       do j = ii, jj

         if(Script(j)(1:7)=='SYSMOL=') then

           read(Script(j)(8:),*) NumSpec

           CheckF(1) = .True.

           allocate( MolName(NumSpec) )
           allocate( NumMol(NumSpec) )
           allocate( NumAtm(NumSpec) )
           allocate( NumMer(NumSpec) )
           allocate( PolyFlag(NumSpec) )

           if((SimMethod=='MC').or.(SimMethod=='gcMC')) then

             if((NumSpec/=1)) then

               write(String1,*) 'ERROR_MPDyn_',Job_name
               open(61, file = trim(adjustl(String1)), status = 'unknown')

               write( 6,*) ' ERROR : In this version, MC method can be used only '
               write( 6,*) '         for one component (rigid body) system.      '
               write(61,*) ' ERROR : In this version, MC method can be used only '
               write(61,*) '         for one component (rigid body) system.      '

               call Finalize

             end if

           end if

         end if

       end do

       if(QMaster) then
         if(.not.CheckF(1)) then
           write(*,*) 'ERROR : "SYSMOL" must be defined in SYSTEM'
           call Finalize
         end if
       end if

       do j = ii, jj

         if(Script(j)(1:8)=='MOLNAME=') then

           read(Script(j)(9:),*) (MolName(k),k=1,NumSpec)

           CheckF(2) = .True.

           do k = 1, NumSpec
             if(MolName(k)(1:4) == 'Poly') then
               PolyFlag(k) = .True.
             else
               PolyFlag(k) = .False.
             end if
           end do

         end if

         if(Script(j)(1:7)=='MOLNUM=') then

           read(Script(j)(8:),*) (NumMol(k),k=1,NumSpec)

           CheckF(3) = .True.

         end if

         if(Script(j)(1:7)=='ATMNUM=') then

           read(Script(j)(8:),*) (NumAtm(k),k=1,NumSpec)

           CheckF(4) = .True.

         end if

       end do

       do j = 1, 4

         if(.not.CheckF(j)) then

           if(j==1) then
             if(QMaster) write(*,*) 'ERROR : "SYSMOL" must be defined in SYSTEM'
           else if(j==2) then
             if(QMaster) write(*,*) 'ERROR : "MOLNAME" must be defined in SYSTEM'
           else if(j==3) then
             if(QMaster) write(*,*) 'ERROR : "MOLNUM" must be defined in SYSTEM'
           else if(j==4) then
             if(QMaster) write(*,*) 'ERROR : "ATMNUM" must be defined in SYSTEM'
           end if

           call Finalize

         end if

       end do

       do j = 1, NumSpec

         if(PolyFlag(j)) then
           NumMer(j) = NumAtm(j)
         end if

       end do

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(.not.Flag) then

     if(QMaster) write(*,*) 'ERROR : the flag "SYSTEM" is necessary in (input.data)'
     call Finalize

   end if

   if(QMaster.and.Qcheck) then

     write(mfile,'(/a)') 'SYSTEM:'
     write(mfile,'(2x,a)') '------------------------------'
     write(mfile,'(2x,a)') 'MolName     NumMol    NumAtm  '

     do i = 1, NumSpec
       write(mfile,'(2x,a8,2(i8,2x))') MolName(i), NumMol(i), NumAtm(i)
     end do

     write(mfile,'(2x,a/)') '------------------------------'

   end if

end subroutine Read_System_Cond


! #####################################################################
! #####################################################################


subroutine Read_PBC_Cond

use ScriptData
use CommonBlocks, only : QMaster, QPBC, QSwitch, QCorrectCutoff, SimMethod, &
&   ForceField
use CGdata, only : Rbk_short2, Rcut_short2, Rrespa, Rheal, Rheal2, Fsw1
use BookParam, only : MaxPair
use TimeParam, only : BookFreq
use CutoffParam

implicit none

integer :: i, ii, k
logical :: Qskin
logical :: Qbook
character(len=3) :: cONOFF, cCORRECT
character(len=70) :: Strcp
real(8) :: Ron, Rcutoff, Rbook

   QPBC=.False.
   Qskin=.False.
   Qbook=.False.
   QdescSW=.False.   
   QSwitch=.False.   ! switching function near the cut-off distance
   QCorrectCutoff=.True.  ! corrections for vdW cut-off, useful only 
!                             when switching func is off

! # these two quantities are still used for isolated systems
   Rrespa = 6.d0
   Rheal  = 4.d0

   MaxPair  = 1200000

   do i = 1, NumOption

     if(OPTIONS(i)=='PBC') then

       if(SimMethod == 'DPD') then
         if(QMaster) then
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           write(*,*) 'WARNING : PBC option is not usefule for DPD'
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         end if
         ReadFlag(i) = .True.
         exit
       end if

       QPBC=.True.

       Rcutoff  = 12.d0
       Rbook    = 14.d0
       BookFreq = 1

       ii = Ist(i)

lia0:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia0

         if(Script(ii)(1:3) == 'RC=') then
           read(Script(ii)(4:),*) Rcutoff
         end if
         if(Script(ii)(1:6) == 'RBOOK=') then
           read(Script(ii)(7:),*) Rbook
           Qbook = .True.
         end if
         if(Script(ii)(1:7) == 'RRESPA=') then
           read(Script(ii)(8:),*) Rrespa
         end if
         if(Script(ii)(1:6) == 'RHEAL=') then
           read(Script(ii)(7:),*) Rheal
         end if
         if(Script(ii)(1:6) == 'RSKIN=') then
           read(Script(ii)(7:),*) Rskin
           Qskin = .True.
         end if
         if(Script(ii)(1:11) == 'VDWCORRECT=') then
           read(Script(ii)(12:),*) cCORRECT
           if(trim(cCORRECT)=='OFF') then
             QCorrectCutoff = .False.
           end if
         end if
         if(Script(ii)(1:8) == 'MAXPAIR=') then
           read(Script(ii)(9:),*) MaxPair
         end if
         if(Script(ii)(1:6) == 'TBOOK=') then
           read(Script(ii)(7:),*) BookFreq
         end if

       end do lia0

       ii = Ist(i)

lia1:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia1

         if(Script(ii)(1:10) == 'VDWSWITCH=') then
           read(Script(ii)(11:),*) cONOFF
           if(trim(cONOFF)/='OFF') then
             Ron = Rcutoff - 2.d0
             Strcp = adjustl(Script(ii)(11:80))
lp0:         do k = 1,70
               if(Strcp(k:k)=='!') exit lp0
               if(Strcp(k:k)/=' '.and.Strcp(k-1:k-1)==' ') then
                 read(Script(ii)(11:),*) cONOFF, Ron
                 exit lp0
               end if
             end do lp0
             QdescSW = .True.
             if(trim(cONOFF)=='ON') then
               QSwitch = .True.
             end if
           end if
         end if

       end do lia1

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QPBC) then

     Ron2 = Ron * Ron

     Rcutoff2 = Rcutoff * Rcutoff

     if(.not.Qbook) then
       if(.not.Qskin) Rskin = 2.d0
       Rbook  = Rcutoff + Rskin
     end if
     if(Qbook) Rskin = Rbook - Rcutoff
     Rbook2 = Rbook * Rbook

   end if

! ## currently just for CG models
   if(QPBC) then

     if(Rrespa+Rheal > Rcutoff) then
       Rrespa = Rcutoff - Rheal
       if(QMaster.and.ForceField(1:2)=='CG') write(*,*) 'RRESPA is set to RCUTOFF - RHEAL'
     end if

   end if

   Rcut_short2 = Rrespa ** 2
   Rheal2      = (Rrespa + Rheal) ** 2
   Rbk_short2  = (Rrespa + Rheal + Rskin) ** 2

! ## for force switching between F_fast and F_slow 
   Fsw1 = 1.d0 / (Rheal2 - Rcut_short2)

! ## just for AA models
   if(QPBC) then
     if(Ron2/=Rcutoff2) then 
       swf1 = 1.d0 / (Rcutoff2 - Ron2) ** 3
     else
       swf1 = 0.d0
     end if
   end if

   if(QSwitch) QCorrectCutoff = .False.

   if(QMaster.and.Qcheck) then

     if(QPBC) then

       write(mfile,'(a)') 'Periodic boundary condition'
       write(mfile,'(a,f10.2,a)') '  Rcutoff  = ',Rcutoff,' [A]'
       write(mfile,'(a,f10.2,a)') '  Rbook    = ',Rbook,  ' [A]'
       write(mfile,'(a,f10.2,a)') '  Rskin    = ',Rskin,  ' [A]'
       if(ForceField(1:2) == 'CG') then
       write(mfile,'(a,f10.2,a)') '  Rrespa(CG)= ',Rrespa, ' [A]'
       write(mfile,'(a,f10.2,a)') '  Rheal(CG) = ',Rheal,  ' [A]'
       end if
       write(mfile,'(a,i10    )') '  MaxPair  = ',MaxPair
       write(mfile,'(a,i10/   )') '  BookFreq = ',BookFreq

     else

       write(mfile,'(a/)') 'Isolated system '
       if(ForceField(1:2) == 'CG') then
       write(mfile,'(a,f10.2,a)')  '  Rrespa(CG)= ',Rrespa, ' [A]'
       write(mfile,'(a,f10.2,a/)') '  Rheal(CG) = ',Rheal,  ' [A]'
       end if

     end if

   end if

end subroutine Read_PBC_Cond


! #####################################################################
! #####################################################################


subroutine DATA_Cond

use ScriptData
use CommonBlocks, only : QMaster, iTrjForm, Qstdout, Qdebug, SimMethod
use IOparam, only : StoreTrj, StoreVel

implicit none

integer :: i, ii
character(len=11) :: cTrjForm
character(len=6) :: cINFO

   StoreTrj = 1000
   StoreVel = 1000

   do i = 1, NumOption

     if(OPTIONS(i)=='DATA_PER_FILE') then

       ii = Ist(i)

lia2:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia2

         if(Script(ii)(1:6) == 'COORD=') then
           read(Script(ii)(7:),*) StoreTrj
         end if
         if(Script(ii)(1:5) == 'VELO=') then
           read(Script(ii)(6:),*) StoreVel
         end if

       end do lia2

       ReadFlag(i) = .True.

       exit

     end if

   end do

   write(cTrjForm,'(a)') 'unformatted'

   do i = 1, NumOption

     if(OPTIONS(i)=='DATA_FORM_TRJ') then

       ii = Ist(i) + 1
       if(Script(ii)(1:5) == 'FORM=') then
         read(Script(ii)(6:),*) cTrjForm
       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(cTrjForm=='formatted') then

     iTrjForm = 1

   else if(cTrjForm=='unformatted') then

     iTrjForm = 2

   else if(cTrjForm=='DCD'.or.cTrjForm=='dcd') then

     iTrjForm = 3

#ifdef MOLFILE
   else if(cTrjForm=='TRR'.or.cTrjForm=='trr') then

     iTrjForm = 4

   else if(cTrjForm=='XTC'.or.cTrjForm=='xtc') then

     iTrjForm = 5
#endif

   else

     if(QMaster) write(*,*) 'ERROR : the defined form of trajectory file'
     call Finalize

   end if

   Qstdout = .False.
   Qdebug  = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='STDOUT') then

       read(Script(Ist(i)+1),*) cINFO

       if(cINFO(1:2)=='ON') then
         Qstdout = .True.
       else if(cINFO(1:6)=='DETAIL') then
         Qstdout = .True.
         Qdebug  = .True.
       else if(cINFO(1:3)/='OFF') then
         if(QMaster) write(*,*) 'ERROR : The keyword, STDOUT, has an illegal option'
         call Finalize
       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

#ifdef MOLFILE
   if(SimMethod /= 'Analysis') then
     if(iTrjForm==4) then
       if(QMaster) write(*,*) 'ERROR: TRR file format is not useful except for Analysis'
       call Finalize
     else if(iTrjForm==5) then
       if(QMaster) write(*,*) 'ERROR: XTC file format is not useful except for Analysis'
       call Finalize
     end if
   end if
#endif

end subroutine DATA_Cond


! #####################################################################
! #####################################################################


subroutine Read_Emin_Cond

use ScriptData
use CommonBlocks, only : QMaster, QMinimization
use Eminparam

implicit none

integer :: i, ii

   QMinimization=.False.

   do i = 1, NumOption

     if(OPTIONS(i)=='ENERGYMIN') then

       QMinimization=.True.

       MinTry = 300
       dRmax  = 0.1d0
       dev_relative = 1.d-5

       ii = Ist(i)

lia3:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia3

         if(Script(ii)(1:10) == 'ITERATION=') then
           read(Script(ii)(11:),*) MinTry
         end if
         if(Script(ii)(1:3) == 'DR=') then
           read(Script(ii)(4:),*) dRmax
         end if
         if(Script(ii)(1:6) == 'ERROR=') then
           read(Script(ii)(7:),*) dev_relative
         end if

       end do lia3

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QMinimization) then

       write(mfile,'(a       )') 'Energy Minimization:'
       write(mfile,'(a,i10   )') '  ITERATION = ', MinTry
       write(mfile,'(a,f10.2 )') '  DR        = ', dRmax
       write(mfile,'(a,e10.2/)') '  ERROR     = ', dev_relative

     end if

   end if

end subroutine Read_Emin_Cond


! #####################################################################
! #####################################################################


subroutine Read_Coulomb_Cond

use ScriptData
use CommonBlocks, only : QMaster, cCOULOMB, QPBC, ForceField
use PMEparam, only : Nfft, Bsp_order
use EwaldParam, only : Alpha, ih2mx, alp2, msh, ar2, kmaxx, kmaxy, kmaxz
use UnitExParam, only : pi
use CGdata, only : Kappa, Eps_relative, Ush0, Ush3, Ush4, Fsh2, Fsh3
use CutoffParam, only : Rcutoff2

implicit none

integer :: i, ii
logical :: FlagCoulomb
real(8) :: Rc1, Rc4, Rc5

   FlagCoulomb = .False.

   Alpha = 0.238
   ih2mx = 81
! # for ER analysis
   kmaxx = 9
   kmaxy = 9
   kmaxz = 9
! #

   Kappa = 0.5d0
   Eps_relative = 1.d0

   do i = 1, NumOption

     if(OPTIONS(i)=='COULOMB') then

!       if(.not.QPBC) then
!         if(QMaster) then
!           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!           write(*,*) 'WARNING : COULOMB option is not useful for isolate system !'
!           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!         end if
!       end if

       Nfft(:) = 32
       Bsp_order = 4

       ii = Ist(i)

lia4:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia4

         if(Script(ii)(1:7) == 'METHOD=') then

           read(Script(ii)(8:),*) cCOULOMB

           FlagCoulomb = .True.

         else if(Script(ii)(1:6) == 'ALPHA=') then

           read(Script(ii)(7:),*) Alpha

         else if(Script(ii)(1:6) == 'K2MAX=') then

           read(Script(ii)(7:),*) ih2mx

         else if(Script(ii)(1:6) == 'KMAXX=') then

           read(Script(ii)(7:),*) kmaxx

         else if(Script(ii)(1:6) == 'KMAXY=') then

           read(Script(ii)(7:),*) kmaxy

         else if(Script(ii)(1:6) == 'KMAXZ=') then

           read(Script(ii)(7:),*) kmaxz

         else if(Script(ii)(1:7) == 'NGRIDX=') then

           read(Script(ii)(8:),*) Nfft(1)

         else if(Script(ii)(1:7) == 'NGRIDY=') then

           read(Script(ii)(8:),*) Nfft(2)

         else if(Script(ii)(1:7) == 'NGRIDZ=') then

           read(Script(ii)(8:),*) Nfft(3)

         else if(Script(ii)(1:6) == 'ORDER=') then

           read(Script(ii)(7:),*) Bsp_order

         else if(Script(ii)(1:6) == 'KAPPA=') then

           read(Script(ii)(7:),*) Kappa

         else if(Script(ii)(1:6) == 'EPS_R=') then

           read(Script(ii)(7:),*) Eps_relative

         end if

       end do lia4

       ReadFlag(i) = .True.

       exit

     end if

   end do

   alp2 = 1.d0 / ( Alpha * Alpha )
   msh  = 1000
   ar2 = 2.d0 * Alpha / sqrt( pi )

   if(cCOULOMB=='SHIFT') then
     Rc1  = sqrt(Rcutoff2)
     Rc4  = Rcutoff2 * Rcutoff2
     Rc5  = Rc4 * Rc1
     Ush0 = - 5.d0 / (3.d0 * Rc1)
     Ush3 =   5.d0 / (3.d0 * Rc4)
     Ush4 = - 1.d0 / Rc5
     Fsh2 = - 5.d0 / Rc4
     Fsh3 =   4.d0 / Rc5
   end if

   if(FlagCoulomb) then
     if((cCOULOMB /= 'EWALD').and.(cCOULOMB /= 'PME').and.(cCOULOMB /= 'Cutoff')&
     &  .and.(cCOULOMB /= 'SCREEN').and.(cCOULOMB /= 'SHIFT')) then
       if(QMaster) then
         write(*,*) 'ERROR : METHOD (COULOMB) should be EWALD or PME or Cutoff'
         write(*,*) 'For CG, the options, SCREEN and SHIFT, are available'
       end if
       call Finalize
     end if
   else
     cCOULOMB = 'EWALD'
   end if

   if(QMaster.and.Qcheck) then

     write(mfile,'(a,a)') 'Coulomb interaction = ', trim(cCOULOMB)

     if(trim(cCOULOMB) == 'EWALD') then

       write(mfile,'(2x,a,f8.3)') ' Alpha = ',Alpha
       write(mfile,'(2x,a,i8/ )') ' K2MAX = ',ih2mx

     else if(trim(cCOULOMB) == 'PME') then

       write(mfile,'(2x,a,f8.3)') ' Alpha = ',Alpha
       write(mfile,'(2x,a,i8  )') ' NFFTx = ',Nfft(1)
       write(mfile,'(2x,a,i8  )') ' NFFTy = ',Nfft(2)
       write(mfile,'(2x,a,i8  )') ' NFFTz = ',Nfft(3)
       write(mfile,'(2x,a,i8/ )') ' order = ',Bsp_order

     else if(trim(cCOULOMB) == 'SCREEN') then

       write(mfile,'(2x,a,f8.3)') ' Kappa = ',Kappa
       write(mfile,'(2x,a,f8.3)') ' Eps_r = ',Eps_relative

     else if(trim(cCOULOMB) == 'SHIFT') then

       write(mfile,'(2x,a,f8.3)') ' Eps_r = ',Eps_relative

     else

       write(mfile,*)

     end if

   end if

end subroutine Read_Coulomb_Cond


! #####################################################################
! #####################################################################


subroutine Read_Thermostat_Cond(Flag)

use ScriptData
use CommonBlocks, only : QMaster, QThermostat, cThermostatMethod, QPathInt
use BathParam, only : Temp_o, Tau_s0, Tau_s1, NHchain

implicit none

integer :: i, ii
real(8) :: xi
logical :: Flag

   QThermostat=.False.
   Flag = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='TEMPERATURE') then

       Temp_o = 298.d0
       Tau_s0 = 0.5d0
       NHchain = 1
       xi = 1.d0

       ii = Ist(i)

lia8:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia8

         if(Script(ii)(1:7) == 'METHOD=') then

           read(Script(ii)(8:),*) cThermostatMethod

           if(cThermostatMethod /= 'OFF') then
             QThermostat=.True.
           end if

         else if(Script(ii)(1:5) == 'TEXT=') then

           read(Script(ii)(6:),*) Temp_o

           Flag = .True.

         else if(Script(ii)(1:7) == 'LENGTH=') then

           read(Script(ii)(8:),*) NHchain

         else if(Script(ii)(1:4) == 'TAU=') then

           read(Script(ii)(5:),*) Tau_s0

         else if(Script(ii)(1:6) == 'EXTAU=') then

           read(Script(ii)(5:),*) xi

         end if

       end do lia8

       Tau_s1 = Tau_s0 * xi

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(.not.Flag) then

     if(QMaster) write(*,*) &
     & 'ERROR : the flag "TEMPERATURE(TEXT)" is necessary in (input.data)'
     call Finalize

   end if

   if(QMaster.and.Qcheck) then

     if(QThermostat) then

       write(mfile,'(a)') 'Temperature :'
       write(mfile,'(a,a      )') '  Method = ',trim(cThermostatMethod)
       write(mfile,'(a,f10.1,a)') '  T_ext  = ',Temp_o,' [K] '

       if((cThermostatMethod=='NHC').or.(cThermostatMethod=='MNHC')) then

         write(mfile,'(a,e10.1,a)') '  Tau    = ',Tau_s0,' [ps] '
         write(mfile,'(a,e10.1,a)') '  Tau_1  = ',Tau_s1,' [ps] '
         write(mfile,'(a,i10)')   '  LENGTH = ', NHchain

       else if(cThermostatMethod=='NH') then

         write(mfile,'(a,e10.1,a)') '  Tau    = ',Tau_s0,' [ps] '

       end if

       write(mfile,*)

     else

       write(mfile,'(a)') 'Temperature :'
       write(mfile,'(a        )') '  Method = OFF '
       write(mfile,'(a,f10.1,a)') '  T_ext  = ',Temp_o,' [K] '

     end if

   end if

   if(QThermostat) then

     if(QPathInt) then
       if((cThermostatMethod /= 'MNHC').and.(cThermostatMethod /= 'NHC')) then
         if(QMaster) write(*,*) &
         & 'ERROR   : (massive) Nose-Hoover chain method should be', &
         & 'selected when PIMD or CMD simulation is carried out.'
         call Finalize
       end if
     end if

   end if


end subroutine Read_Thermostat_Cond


! #####################################################################
! #####################################################################


subroutine Read_Annealing

use ScriptData
use CommonBlocks, only : QMaster, QThermostat, cThermostatMethod, QPBC
use SimAnneal, only : QSimAnneal, cAnnealMethod, T_target, switch_time, &
&   factA, factB, factC
use TimeParam, only : Nstep, deltat
use BathParam, only : Temp_o

implicit none

integer :: i, ii

   QSimAnneal=.False.

   do i = 1, NumOption

     if(OPTIONS(i)=='ANNEAL') then

       QSimAnneal = .True.

       T_target = 1.d0
       cAnnealMethod = 'LINEAR'
       switch_time = deltat * Nstep * 0.75

       ii = Ist(i)

lid8:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lid8

         if(Script(ii)(1:7) == 'METHOD=') then

           read(Script(ii)(8:),*) cAnnealMethod

         else if(Script(ii)(1:9) == 'T_TARGET=') then

           read(Script(ii)(10:),*) T_target

         else if(Script(ii)(1:9) == 'SWITCH_T=') then

           read(Script(ii)(10:),*) switch_time

         end if

       end do lid8

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QSimAnneal) then

       write(mfile,'(a)') 'Simulated Annealing :'
       write(mfile,'(a,a      )') '  Method    = ',trim(cAnnealMethod)
       write(mfile,'(a,f10.1,a)') '  T_target  = ',T_target,' [K] '

       if(cAnnealMethod=='HYBRID') then

         write(mfile,'(a,e10.1,a)') ' Switch t = ',switch_time,' [ps] '

       end if

       write(mfile,*)

     end if

   end if

   if(QSimAnneal) then

     if((.not.QThermostat).or.(cThermostatMethod/='NHC')) then
       if(QMaster) write(*,*) &
       & 'ERROR : whenever you use "ANNEAL", you must choose "NHC" thermostat.'
       call Finalize
     end if

     if(QPBC) then
       if(QMaster) write(*,*) 'ERROR : ANNEAL cannot be used with PBC'
       call Finalize
     end if

     if(cAnnealMethod=='LINEAR') then
       factA = (T_target - Temp_o) / dble(Nstep)
     else if(cAnnealMethod=='HYBRID') then
       switch_time = switch_time / deltat
       factB = (Temp_o - T_target) / ( 2.d0 / switch_time - 1.d0 / dble(Nstep) )
       factA = - factB / (switch_time * switch_time)
       factC = T_target - factB / dble(Nstep)
     else
       if(QMaster) write(*,*) 'ERROR : you put a wrong name for METHOD= in ANNEAL'
       call Finalize
     end if

   end if

end subroutine Read_Annealing


! #####################################################################
! #####################################################################


subroutine Read_Barostat_Cond(Flag)

use ScriptData
use CommonBlocks, only : QMaster, QBarostat, QPINPT, QATaup, SimMethod, &
&   cBarostatMethod, QPathInt
use BathParam, only : Pressure_o, Tau_p, prefTaup, Stress_ext, SigmaS, &
&   CoupleEdge
use UnitExParam, only : rprs

implicit none

logical :: Flag
integer :: i, ii, j, k
character(len=2) :: Isodir
! ## Stress >>
real(8), dimension(3,3) :: H_ref, InvH_ref, InvH_ref_T
real(8), dimension(3,3) :: Stress, tempM
real(8) :: Volume_ref, det
logical, dimension(3,3) :: QHref
external det
! ## << Stress

   QBarostat=.False.
   QPINPT=.False.
   QHref =.False.

   do i = 1, NumOption

     if(OPTIONS(i)=='PRESSURE') then

       QBarostat=.True.

       Pressure_o = 0.1d0
       Tau_p = 2.d0

       prefTaup = 1.d0
       QATaup   = .False.

       if((SimMethod=='MD').and.(.not.Flag)) then
         if(QMaster) write(*,*) &
         & 'ERROR : The reference TEMPERATURE must be define even in case of NPH-MD'
         call Finalize
       end if

       ii = Ist(i)

lia5:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia5

         if(Script(ii)(1:7) == 'METHOD=') then

           read(Script(ii)(8:),*) cBarostatMethod

         else if(Script(ii)(1:5) == 'PEXT=') then

           read(Script(ii)(6:),*) Pressure_o

         else if(Script(ii)(1:4) == 'TAU=') then

           read(Script(ii)(5:),*) Tau_p

         end if

       end do lia5

       Pressure_o = 1.d+6 * Pressure_o * rprs

       if(cBarostatMethod == 'ST') then

         ii = Ist(i)

         Stress_ext = 0.d0

lia6:    do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia6

           if(Script(ii)(1:4) == 'SXX=') then

             read(Script(ii)(5:),*) Stress_ext(1,1)

           else if(Script(ii)(1:4) == 'SXY=') then

             read(Script(ii)(5:),*) Stress_ext(1,2)

           else if(Script(ii)(1:4) == 'SXZ=') then

             read(Script(ii)(5:),*) Stress_ext(1,3)

           else if(Script(ii)(1:4) == 'SYX=') then

             read(Script(ii)(5:),*) Stress_ext(2,1)

           else if(Script(ii)(1:4) == 'SYY=') then

             read(Script(ii)(5:),*) Stress_ext(2,2)

           else if(Script(ii)(1:4) == 'SYZ=') then

             read(Script(ii)(5:),*) Stress_ext(2,3)

           else if(Script(ii)(1:4) == 'SZX=') then

             read(Script(ii)(5:),*) Stress_ext(3,1)

           else if(Script(ii)(1:4) == 'SZY=') then

             read(Script(ii)(5:),*) Stress_ext(3,2)

           else if(Script(ii)(1:4) == 'SZZ=') then

             read(Script(ii)(5:),*) Stress_ext(3,3)

           else if(Script(ii)(1:4) == 'HXX=') then

             read(Script(ii)(5:),*) H_ref(1,1)
                 QHref(1,1) = .True.
            else if(Script(ii)(1:4) == 'HXY=') then

             read(Script(ii)(5:),*) H_ref(1,2)
                 QHref(1,2) = .True.
            else if(Script(ii)(1:4) == 'HXZ=') then

             read(Script(ii)(5:),*) H_ref(1,3)
             QHref(1,3) = .True.

           else if(Script(ii)(1:4) == 'HYX=') then

             read(Script(ii)(5:),*) H_ref(2,1)
             QHref(2,1) = .True.

           else if(Script(ii)(1:4) == 'HYY=') then

             read(Script(ii)(5:),*) H_ref(2,2)
             QHref(2,2) = .True.

           else if(Script(ii)(1:4) == 'HYZ=') then

             read(Script(ii)(5:),*) H_ref(2,3)
             QHref(2,3) = .True.

           else if(Script(ii)(1:4) == 'HZX=') then

             read(Script(ii)(5:),*) H_ref(3,1)
             QHref(3,1) = .True.

           else if(Script(ii)(1:4) == 'HZY=') then

             read(Script(ii)(5:),*) H_ref(3,2)
             QHref(3,2) = .True.

           else if(Script(ii)(1:4) == 'HZZ=') then

             read(Script(ii)(5:),*) H_ref(3,3)
             QHref(3,3) = .True.

           end if

         end do lia6

         do j = 1, 3
           do k = 1, 3
             if(.not.QHref(j,k)) then
               if(QMaster) then
                 write(*,*) 'ERRORs in the keyword PRESSURE'
                 write(*,*) 'The reference cell matrix should be defined'
                 write(*,*) 'whenever ST is selected ! '
               end if
               call Finalize
             end if
           end do
         end do

         call InversMatrix(H_ref,InvH_ref)
         tempM = transpose( H_ref )
         call InversMatrix(tempM,InvH_ref_T)

         Stress_ext = 1.d+6 * Stress_ext * rprs

         do j = 1, 3
           do k = 1, 3

             if(j==k) then
               Stress(j,j) = Stress_ext(j,j) - Pressure_o
             else
               Stress(j,k) = Stress_ext(j,k)
             end if

           end do
         end do

         tempM  = matmul(InvH_ref, Stress)
         SigmaS = matmul(tempM, InvH_ref_T)

         Volume_ref = det(H_ref)
         SigmaS = SigmaS * Volume_ref

       end if

       if((cBarostatMethod == 'ST').or.(cBarostatMethod == 'PR').or. &
       &  (cBarostatMethod == 'A3').or.(cBarostatMethod == 'A2')) then

         ii = Ist(i)

lia7:    do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia7

           if(Script(ii)(1:6) == 'EXTAU=') then

             QATaup   = .True.
             read(Script(ii)(7:),*) prefTaup(1,1),prefTaup(2,2),prefTaup(3,3),&
             &                      prefTaup(1,2),prefTaup(2,3),prefTaup(1,3)

             prefTaup(2,1) = prefTaup(1,2)
             prefTaup(3,2) = prefTaup(2,3)
             prefTaup(3,1) = prefTaup(1,3)

           end if

         end do lia7

       end if

       if(cBarostatMethod == 'A2') then

         Isodir = 'XY'

         CoupleEdge(1) = 1
         CoupleEdge(2) = 2

         ii = Ist(i)

loop:    do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit loop

           if(Script(ii)(1:8) == 'ISOPAIR=') then

             QATaup   = .True.
             read(Script(ii)(9:),*) Isodir

             if(Isodir=='XY'.or.Isodir=='YX'.or.Isodir=='xy'.or.Isodir=='yx') then
               CoupleEdge(1) = 1
               CoupleEdge(2) = 2
             else if(Isodir=='YZ'.or.Isodir=='ZY'.or.Isodir=='yz'.or.Isodir=='zy') then
               CoupleEdge(1) = 2
               CoupleEdge(2) = 3
             else if(Isodir=='XZ'.or.Isodir=='ZX'.or.Isodir=='xz'.or.Isodir=='zx') then
               CoupleEdge(1) = 1
               CoupleEdge(2) = 3
             else
               if(QMaster) write(*,*) 'ERROR : ISOPAIR in PRESSURE'
               call Finalize
             end if

           end if

         end do loop

       end if

       if((cBarostatMethod=='A2').and.(prefTaup(2,2)/=prefTaup(1,1))) then
         if(QMaster) then
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           write(*,*) '  WARNING: check EXTAU in PRESSURE keyword       '
           write(*,*) '  Tau_yy has been changed to be equal to Tau_xx  '
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         end if
         prefTaup(2,2) = prefTaup(1,1)
       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QPathInt.and.QBarostat) QPINPT = .True.

   if(QMaster.and.Qcheck) then

     if(QBarostat) then

       write(mfile,'(a)') 'Pressure Control:'
       write(mfile,'(a,a      )') '  Method = ',trim(cBarostatMethod)
       write(mfile,'(a,f10.1,a)') '  P_ext  = ',1.d-6 * Pressure_o / rprs,' [MPa] '
       write(mfile,'(a,e10.1,a)') '  Tau    = ',Tau_p,' [ps] '
       if(QATaup) then
         write(mfile,'(a,6e8.1)') '  Extau  = ', &
         &                        prefTaup(1,1),prefTaup(2,2),prefTaup(3,3),&
         &                        prefTaup(1,2),prefTaup(2,3),prefTaup(1,3)
       end if

       if(cBarostatMethod == 'ST') then

         write(mfile,'(a)') '  Stress_ext = '
         write(mfile,'(2x,3e10.2)') 1.d-6 * Stress_ext(1,:) / rprs
         write(mfile,'(2x,3e10.2)') 1.d-6 * Stress_ext(2,:) / rprs
         write(mfile,'(2x,3e10.2)') 1.d-6 * Stress_ext(3,:) / rprs

       end if

       if(cBarostatMethod == 'A2') then
         write(mfile,'(2a)') '  isotropic in plane ',Isodir
       end if

       write(mfile,*)

     end if

   end if

end subroutine Read_Barostat_Cond


! #####################################################################
! #####################################################################


subroutine Read_Constraint_Cond

use ScriptData
use CommonBlocks, only : QMaster, QRigidBody, QSHAKE, cThermostatMethod, &
&   SimMethod, QThermostat
use RBparam, only : NumRB
use SHAKEparam, only : TolA, TolB, MaxIteration

implicit none

integer :: i, ii
character(len=10) :: cConstraint

   QRigidBody = .False.
   QSHAKE     = .False.
   NumRB      = 0

   TolA         = 1.d-6
   TolB         = 1.d-28
   MaxIteration = 50

   do i = 1, NumOption

     if(OPTIONS(i)=='CONSTRAINT') then

       ii = Ist(i)

lia9:  do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia9

         if(Script(ii)(1:7) == 'METHOD=') then

           read(Script(ii)(8:),*) cConstraint

           if(cConstraint=='SHAKE') then

             if(SimMethod == 'HMC') then
               if(QMaster) write(*,*) 'ERROR : SHAKE cannot be used in HMC simulations.'
               call Finalize
             end if

             if(QThermostat .and. cThermostatMethod == 'MNHC') then
               if(QMaster) write(*,*) &
               &  'ERROR : SHAKE cannot be used with the Massive Nose-Hoover Chain.'
               call Finalize
             end if

             QSHAKE     = .True.
             QRigidBody = .False.
             NumRB      = 0

           else if(cConstraint=='RigidBody') then  !<<new>>

             QRigidBody = .True.
             QSHAKE     = .False.

           else

             if(QMaster) write(*,*) 'ERROR : Constraint condition'

           end if

         else if(Script(ii)(1:5) == 'TOLR=') then

           read(Script(ii)(6:),*) TolA

         else if(Script(ii)(1:5) == 'TOLV=') then

           read(Script(ii)(6:),*) TolB

         else if(Script(ii)(1:5) == 'MAXI=') then

           read(Script(ii)(6:),*) MaxIteration

         end if

       end do lia9

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QRigidBody) then

       write(mfile,'(a/)') 'Rigid-Body Model'

     else if(QSHAKE) then

       write(mfile,'(a)') 'Constraint dynamics : SHAKE'
       write(mfile,'(a,e10.2)') '  Tolerance R  = ',TolA
       write(mfile,'(a,e10.2)') '  Tolerance V  = ',TolB
       write(mfile,'(a,i10/)' ) '  MaxIteration = ',MaxIteration

     else

       write(mfile,'(a/)') 'Flexible model'

     end if

   end if

end subroutine Read_Constraint_Cond


! #####################################################################
! #####################################################################


subroutine Read_Temperature

use ScriptData
use CommonBlocks, only : QMaster, SimMethod, QThermostat
use BathParam, only : Temp_o

implicit none

integer :: i, ii

   if(SimMethod == 'DPD') then

     Temp_o = 1.d0

     do i = 1, NumOption

       if(OPTIONS(i)=='TEMPERATURE') then

         ii = Ist(i) + 1

         if(Script(ii)(1:3) == 'kT=') then

           read(Script(ii)(4:),*) Temp_o

         else

           if(QMaster) write(*,*) 'ERROR : TEMPERATURE for DPD'

         end if

         ReadFlag(i) = .True.

         exit

       end if

     end do

     if(QMaster.and.Qcheck) then

       write(mfile,'(a,f10.3/)') 'kT = ', Temp_o

     end if

   else if(SimMethod == 'HMC') then

     QThermostat=.False.
     Temp_o = 298.d0

     do i = 1, NumOption

       if(OPTIONS(i)=='TEMPERATURE') then

         ii = Ist(i) + 1

         if(Script(ii)(1:5) == 'TEXT=') then

           read(Script(ii)(6:),*) Temp_o

         else

           if(QMaster) write(*,*) 'ERROR : TEMPERATURE for HMC'

         end if

         ReadFlag(i) = .True.

         exit

       end if

     end do

     if(QMaster.and.Qcheck) then

       write(mfile,'(a,f10.3/)') 'TEXT = ', Temp_o

     end if

   end if

end subroutine Read_Temperature


! #####################################################################
! #####################################################################


subroutine Read_StepNumber

use ScriptData
use CommonBlocks, only : QMaster, SimMethod, QAveTh
use QMDynamics, only : Icurrent
use TimeParam, only : Nstep, itgn, ixc, ixv, irs

implicit none

integer :: i, ii

   Nstep = 10000
   itgn  = 1
   ixc   = 10
   ixv   = 10
   Icurrent = -1

   do i = 1, NumOption

     if(OPTIONS(i)=='STEP_NUMBER') then

       ii = Ist(i)

lia10: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia10

         if(Script(ii)(1:6) == 'TOTAL=') then

           read(Script(ii)(7:),*) Nstep

         else if(Script(ii)(1:7) == 'ENERGY=') then

           read(Script(ii)(8:),*) itgn

         else if(Script(ii)(1:6) == 'COORD=') then

           read(Script(ii)(7:),*) ixc

         else if(Script(ii)(1:5) == 'VELO=') then

           read(Script(ii)(6:),*) ixv

         else if(Script(ii)(1:8) == 'CURRENT=') then

           read(Script(ii)(9:),*) Icurrent

         end if

       end do lia10

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(SimMethod == 'DynaLib') then
     if(Icurrent == -1) then
       if(QMaster) write(*,*) '"Icurrent" option must be set correctly!'
       call Finalize
     end if
   end if

   if(itgn==0) itgn = Nstep + 1
   if(ixc==0) ixc = Nstep + 1
   if(ixv==0) ixv = Nstep + 1
   if(irs==0) irs = Nstep + 1

   if(itgn==1) QAveTh = .False.

   if(QMaster.and.Qcheck) then

     write(mfile,'(a)')     'Step numbers'
     write(mfile,'(a)')     '  --------------------'
     write(mfile,'(a,i10)') '  Total  = ',Nstep
     write(mfile,'(a,i10)') '  Energy = ',itgn 
     write(mfile,'(a,i10)') '  Coord. = ',ixc  
     write(mfile,'(a,i10)') '  Velo.  = ',ixv  
     if(SimMethod == 'DynaLib') then
     write(mfile,'(a,i10)') '  Current= ',Icurrent
     end if
     write(mfile,'(a/)')    '  --------------------'

   end if

end subroutine Read_StepNumber


! #####################################################################
! #####################################################################


subroutine Read_Time_Step

use ScriptData
use CommonBlocks, only : QMaster, SimMethod, QPBC, ForceField, QPathInt
use TimeParam, only : deltat, lp, lk, BookFreq, Nstep, itgn, ixc, ixv, irs
use CommonHMC, only : MDstep, QPartial, ThetaPartial
use UnitExParam, only : pi

implicit none

integer :: i, ii
real(8) :: dt_mode, dt_long
character(len=1) :: Ch

   lp = 1
   lk = 1
   if((SimMethod == 'MD').or.(QPathInt).or.(SimMethod == 'HMC').or. &
   &  (SimMethod == 'DynaLib')) then
     deltat = 0.001d0
   else if(SimMethod == 'DPD') then
     deltat = 0.01d0
   else
     Return
   end if

   do i = 1, NumOption

     if(OPTIONS(i)=='TIME_STEP') then

       ii = Ist(i) + 1

       if(Script(ii)(1:4) == 'SNGL') then

         read(Script(ii)(5:),*) deltat

       else if(Script(ii)(1:4) == 'MULT') then

         read(Script(ii)(5:),*) dt_long, dt_mode, deltat

         lp = nint(dt_mode/deltat)
         lk = nint(dt_long/deltat)

         if(SimMethod == 'DPD') then

           if(QMaster) then
             write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
             write(*,*) 'WARNING : MULT cannot be used for DPD'
             write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           end if
           lp = 1
           lk = 1

         end if

       else

         if(QMaster) write(*,*) 'ERROR : TIME_STEP option!'
         call Finalize

       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QPBC.and.((BookFreq < lp).or.(mod(BookFreq,lp) /= 0))) then

#ifndef BMONI
     if(QMaster) write(11,*) 'ERROR : BookFreq must be chosen to be i*lp (dt_mode/dt_short)'
#endif
     if(QMaster) write( *,*) 'ERROR : BookFreq must be chosen to be i*lp (dt_mode/dt_short)'
     call Finalize

   end if

   if(QPBC.and.(ForceField(1:2)=='CG').and.((BookFreq < lk).or.(mod(BookFreq,lk) /= 0))) then

#ifndef BMONI
     if(QMaster) write(11,*) 'ERROR : BookFreq must be chosen to be i*lk (dt_slow/dt_short)'
#endif
     if(QMaster) write( *,*) 'ERROR : BookFreq must be chosen to be i*lk (dt_slow/dt_short)'
     if(QMaster) write( *,*) '        when you chose the forcefield "CG"'
     call Finalize

   end if

   if(SimMethod == 'HMC') then

     MDstep = 10
     QPartial = .False.
     ThetaPartial = 0.5 * pi

     do i = 1, NumOption

       if(OPTIONS(i)=='HMC_PARAM') then

         ii = Ist(i)

lib0:    do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lib0

           if(Script(ii)(1:7) == 'MDSTEP=') then

             read(Script(ii)(8:),*) MDstep

           else if(Script(ii)(1:8) == 'PARTIAL=') then

             write(Ch,'(a1)') trim(adjustl(Script(ii)(9:)))
             if( Ch == 'Y' ) then
               QPartial = .True.
             else if( Ch == 'N' ) then
               QPartial = .False.
             else
               write(*,*) 'ERROR : "QPARTIAL" of "HMC_PARAM"'
             end if

           else if(Script(ii)(1:6) == 'THETA=') then

             read(Script(ii)(7:),*) ThetaPartial
             ThetaPartial = ThetaPartial / 180. * pi

           end if

         end do lib0

         ReadFlag(i) = .True.

         exit

       end if

     end do

   end if

   if((SimMethod == 'MD').or.(QPathInt)) then

     Nstep = lk * Nstep
     itgn  = lk * itgn
     ixc   = lk * ixc
     ixv   = lk * ixv
     irs   = lk * irs

     if(QMaster.and.Qcheck) then

       if(lp == 1 .and. lk == 1) then

         write(mfile,'(a        )') 'Single time step:'
         write(mfile,'(a,f10.3,a/)') '  dt = ',deltat * 1.d+03,' [fs]'

       else

         write(mfile,'(a        )')  'Multiple time step:'
         write(mfile,'(a,f10.3,a )') '  dt_short = ',deltat    * 1.d+03,' [fs]'
         write(mfile,'(a,f10.3,a )') '  dt_mode  = ',deltat*lp * 1.d+03,' [fs]'
         write(mfile,'(a,f10.3,a/)') '  dt_long  = ',deltat*lk * 1.d+03,' [fs]'

       end if

     end if

   else if(SimMethod == 'HMC') then

     MDstep = MDstep * lk

   end if

end subroutine Read_Time_Step


! #####################################################################
! #####################################################################


subroutine Read_Bath_Step

use ScriptData
use CommonBlocks, only : QMaster, SimMethod, QPathInt
use BathParam, only : Nsc, NYoshid, wdti2, wdti4, wdti8
use TimeParam, only : deltat

implicit none

integer :: i, ii
real(8) :: www

   if((SimMethod == 'MD').or.(QPathInt).or.(SimMethod == 'DynaLib')) then

     Nsc = 1
     NYoshid = 1

     do i = 1, NumOption

       if(OPTIONS(i)=='BATH_INT') then

         ii = Ist(i)

lib1:    do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lib1

           if(Script(ii)(1:6) == 'MULTI=') then

             read(Script(ii)(7:),*) Nsc

           else if(Script(ii)(1:6) == 'ORDER=') then

             read(Script(ii)(7:),*) NYoshid

           else

             if(QMaster) write(*,*) 'ERROR : BATH_INT option!'
             call Finalize

           end if

         end do lib1

         ReadFlag(i) = .True.

         exit

       end if

     end do

     allocate( wdti2(NYoshid) )
     allocate( wdti4(NYoshid) )
     allocate( wdti8(NYoshid) )

! -------------------------------------------
! ## order of expansion [NYoshid=5] for bath
! -------------------------------------------

     if(NYoshid == 1) then

       wdti2(1) = deltat / ( 2.d0 * Nsc )
       wdti4(1) = wdti2(1) * 0.5d0
       wdti8(1) = wdti4(1) * 0.5d0

     else if(NYoshid == 3) then

       www = 1.d0 / ( 2.d0 - 2.d0**(1.d0/3.d0) )

       wdti2(1) = www * deltat / ( 2.d0 * Nsc )
       wdti2(2) = ( 1.d0 - 2.d0 * www ) * deltat / ( 2.d0 * Nsc )
       wdti2(3) = wdti2(1)

       wdti4(1) = wdti2(1) * 0.5d0
       wdti4(2) = wdti2(2) * 0.5d0
       wdti4(3) = wdti4(1)

       wdti8(1) = wdti4(1) * 0.5d0
       wdti8(2) = wdti4(2) * 0.5d0
       wdti8(3) = wdti8(1)

     else if(NYoshid == 5) then

       www = 1.d0 / ( 4.d0 - 4.d0**(1.d0/3.d0) )

       wdti2(1) = www * deltat / ( 2.d0 * Nsc )
       wdti2(2) = wdti2(1)
       wdti2(3) = ( 1.d0 - 4.d0 * www ) * deltat / ( 2.d0 * Nsc )
       wdti2(4) = wdti2(1)
       wdti2(5) = wdti2(1)

       wdti4(1) = wdti2(1) * 0.5d0
       wdti4(2) = wdti4(1)
       wdti4(3) = wdti2(3) * 0.5d0
       wdti4(4) = wdti4(1)
       wdti4(5) = wdti4(1)

       wdti8(1) = wdti4(1) * 0.5d0
       wdti8(2) = wdti8(1)
       wdti8(3) = wdti4(3) * 0.5d0
       wdti8(4) = wdti8(1)
       wdti8(5) = wdti8(1)

     end if

   end if

end subroutine Read_Bath_Step


! #####################################################################
! #####################################################################


subroutine Read_Bead_Length

use ScriptData
use CommonBlocks, only : QMaster, QPathInt
use CommonPI, only : Nbead, Ncref, InvP, deltat_ref, Pbead
use TimeParam, only : deltat
implicit none

integer :: i, ii

   if(QPathInt) then

     Nbead = 16
     Ncref = 10

     do i = 1, NumOption

       if(OPTIONS(i)=='PI_PARAM') then

         ii = Ist(i)

lif5:    do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lif5

           if(Script(ii)(1:12) == 'BEAD_LENGTH=') then

             read(Script(ii)(13:),*) Nbead

           else if(Script(ii)(1:9) == 'REF_TIME=') then

             read(Script(ii)(10:),*) Ncref

           else

             if(QMaster) write(*,*) 'ERROR : PI_PARAM option!'
             call Finalize

           end if

         end do lif5

         ReadFlag(i) = .True.

         exit

       end if

     end do

     Pbead = dble(Nbead)
     InvP  = 1.d0 / Pbead

     deltat_ref = deltat / dble(Ncref)

   end if

end subroutine Read_Bead_Length


! #####################################################################
! #####################################################################


subroutine Read_CMD_param

use ScriptData
use CommonBlocks, only : QMaster, QPathInt
use CommonPI, only : GammaC2

implicit none

integer :: i, ii
real(8) :: TimeScalingFactor, GammaC

   if(QPathInt) then

     TimeScalingFactor = 1.d0

     do i = 1, NumOption

       if(OPTIONS(i)=='CMD_PARAM') then

         ii = Ist(i) + 1

         if(Script(ii)(1:11) == 'TIME_SCALE=') then

             read(Script(ii)(12:),*) TimeScalingFactor

         else

           if(QMaster) write(*,*) 'ERROR : CMD_PARAM option!'
           call Finalize

         end if

         ReadFlag(i) = .True.

         exit

       end if

     end do

     GammaC = 1.d0 / TimeScalingFactor
     GammaC2 = GammaC * GammaC

   end if

end subroutine Read_CMD_param


! #####################################################################
! #####################################################################


subroutine Read_AvThermo

use ScriptData
use CommonBlocks, only : QMaster, QAveTh

implicit none

integer :: i, ii

   QAveTh = .True.

   do i = 1, NumOption

     if(OPTIONS(i)=='AVEMONITOR') then

       ii = Ist(i) + 1

       if(Script(ii)(1:3) == 'OFF') then

         QAveTh = .False.

       else if(Script(ii)(1:2) /= 'ON') then

         if(QMaster) write(*,*) 'ERROR : AVEMONITOR option -- Choose "ON" or "OFF"!'
         call Finalize

       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

end subroutine Read_AvThermo


! #####################################################################
! #####################################################################


subroutine Mic1

use ScriptData
use CommonBlocks, only : QMaster, SimMethod, QPathInt
use TimeParam, only : dt2, deltat
use BathParam, only : kT, Temp_o, Beta
use UnitExParam, only : kb, Plank2pi
use CommonPI, only : OmegaP2, Pbead

implicit none

real(8) :: OmegaP

   if((SimMethod /= 'Analysis').and.(SimMethod /= 'MM').and.&
   &  (SimMethod /= 'NMA')) then

     dt2    = deltat * 0.5d0

! ## parameters

     kT   = kb * Temp_o

     Beta = 1.d0 / kT

   end if

! ## PIMD
   if(QPathInt) then
     OmegaP  = sqrt(Pbead) * kT / Plank2pi
     OmegaP2 = OmegaP * OmegaP
   end if

end subroutine Mic1


! #####################################################################
! #####################################################################


subroutine Read_FixedAtom_Cond

use ScriptData
use CommonBlocks, only : QMaster
use Configuration, only : R
use OptConstraintParam, only : NHam, HamR, kHam, RIni, Rrot
use Numbers, only : N, NumMol, NumAtm
use AtomParam, only : AtomName, Mass, InvMass
use UnitExParam, only : ExParam

implicit none

integer :: i, ii, j, k
character(len= 4) :: TResidName
character(len= 4) :: TAtomName
integer :: TResidNum
character(len= 3) :: Dummy
logical :: QRef
character(len=80) :: RefFile
integer :: Natom, FROM, TO
character(len=80) :: String
integer :: Numa, Numb

   NHam = 0

   HamR = .False.

   kHam = 100.d0 * ExParam
   QRef = .False.
   write(RefFile,'(a)') 'initial.crd'

   do i = 1, NumOption

     if(OPTIONS(i)=='FIXEDATOM') then

       ii = Ist(i)

lia11: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia11

         String = Script(ii)

         if(String(1:5) == 'FILE=') then

           QRef = .True.
           read(String(6:),*) RefFile

         else if(String(1:2) == 'k=') then

           read(String(3:),*) kHam

           kHam = kHam * ExParam

         else if(String(1:3) == 'MAS') then

           if(String(5:14)=='COMPONENT ') then

             read(String(15:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               Mass(k)    = Mass(k)    * 1.d+06
               InvMass(k) = InvMass(k) * 1.d-06

             end do

           else if(String(5:18)=='COMPONENT_MAIN') then

             read(String(19:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               if((AtomName(k)=='N' ).or.(AtomName(k)=='HN').or. &
               &  (AtomName(k)=='CA').or.(AtomName(k)=='HA').or. &
               &  (AtomName(k)=='C' ).or.(AtomName(k)=='O' )) then

                 Mass(k)    = Mass(k)    * 1.d+06
                 InvMass(k) = InvMass(k) * 1.d-06

               end if

             end do

           else if(String(5:18)=='COMPONENT_SIDE') then

             read(String(14:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               if((AtomName(k)/='N' ).and.(AtomName(k)/='HN').and. &
               &  (AtomName(k)/='CA').and.(AtomName(k)/='HA').and. &
               &  (AtomName(k)/='C' ).and.(AtomName(k)/='O' )) then

                 Mass(k)    = Mass(k)    * 1.d+06
                 InvMass(k) = InvMass(k) * 1.d-06

               end if

             end do

           else

             read(String,*) Dummy,FROM, TO

             do k = FROM, TO

               Mass(k)    = Mass(k)    * 1.d+06
               InvMass(k) = InvMass(k) * 1.d-06

             end do

           end if

         else if(String(1:3) == 'HAM') then

           if(String(5:14)=='COMPONENT ') then

             read(String(15:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               if(.not.HamR(k)) NHam = NHam + 1

               HamR(k) = .True.

             end do

           else if(String(5:18)=='COMPONENT_MAIN') then

             read(String(19:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               if((AtomName(k)=='N' ).or.(AtomName(k)=='HN').or. &
               &  (AtomName(k)=='CA').or.(AtomName(k)=='HA').or. &
               &  (AtomName(k)=='C' ).or.(AtomName(k)=='O' )) then

                 if(.not.HamR(k)) NHam = NHam + 1
                 HamR(k) = .True.

               end if

             end do

           else if(String(5:19)=='COMPONENT_CHAIN') then

             read(String(20:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               if((AtomName(k)=='N' ).or.(AtomName(k)=='CA').or. &
               &  (AtomName(k)=='C' )) then

                 if(.not.HamR(k)) NHam = NHam + 1
                 HamR(k) = .True.

               end if

             end do

           else if(String(5:13)=='MAINchain') then

             read(String(14:72),*) FROM,TO

             do k = FROM, TO

               if((AtomName(k)=='N' ).or.(AtomName(k)=='CA').or. &
               &  (AtomName(k)=='C' )) then

                 if(.not.HamR(k)) NHam = NHam + 1
                 HamR(k) = .True.

               end if

             end do

           else if(String(5:18)=='COMPONENT_SIDE') then

             read(String(14:72),*) j

             Numa = 0

             do k = 1, j-1

               Numa = Numa + NumMol(k)*NumAtm(k)

             end do

             Numb = Numa + NumMol(j)*NumAtm(j)

             do k = Numa + 1 , Numb

               if((AtomName(k)/='N' ).and.(AtomName(k)/='HN').and. &
               &  (AtomName(k)/='CA').and.(AtomName(k)/='HA').and. &
               &  (AtomName(k)/='C' ).and.(AtomName(k)/='O' )) then

                 if(.not.HamR(k)) NHam = NHam + 1
                 HamR(k) = .True.

               end if

             end do

           else

             read(String,*) Dummy,FROM, TO

             do k = FROM, TO

               if(.not.HamR(k)) NHam = NHam + 1

               HamR(k) = .True.

             end do

           end if

         else

           if(QMaster) write(*,*) 'OPTION ERROR : check the flag for atom constraints'
           if(QMaster) write(*,*) '               in the FIXEDATOM optional file     '
           call Finalize

         end if

       end do lia11

       ReadFlag(i) = .True.

       exit

     end if

   end do

! ## To Remove Interaction between the Fixed Atoms

   if(NHam /= 0) then

     allocate( RIni(3,N) )
     allocate( Rrot(3,N) )

     if(QRef) then

       open(45,file=trim(adjustl(RefFile)),status='old')

         read(45,'(/)')
         read(45,*) Natom

         if(Natom/=N) then
           if(QMaster) then
             write(*,*) 'error : Fixed Atom reference position!'
             write(*,*) 'CHECK the number of atoms in "', RefFile,'"'
           end if
           call Finalize
         end if

         do k = 1, N

         read(45,'(2i5,2(x,a4),3f10.5)')    &
         & j , TResidNum , TResidName ,    &
         & TAtomName , RIni(:,k)

         end do

       close(45)

     else

       RIni = R

     end if

   end if

end subroutine Read_FixedAtom_Cond


! #####################################################################
! #####################################################################


subroutine Read_FixedCOM_Cond

use ScriptData
use CommonBlocks, only : QMaster, QFixCOM
use Configuration, only : R
use OptConstraintParam, only : Nccom, kccom, Idcomf, Idcomt, Xcon, &
&   Rcon, fscom, fcom, InvMscom
use Numbers, only : N, NumMol, NumAtm, NumSpec
use AtomParam, only : AtomName, Mass, InvMass
use UnitExParam, only : ExParam

implicit none

integer :: i, ii, j, k, jj, kk, l, i1, ini
character(len=80) :: String
integer :: Numbf1, Numbl1, Numhf1

   Numbf1 = 0
   Numbl1 = 0
   Numhf1 = 0

   Nccom = 0

   kccom = 100.d0 * ExParam

   QFixCOM = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='FIXEDCOM') then

       ii = Ist(i)

       QFixCOM = .True.

       jj = 0
liq11: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit liq11

         String = Script(ii)

         if(String(1:2) == 'k=') then

           read(String(3:),*) kccom

           kccom = kccom * ExParam

         else if(String(1:9) == 'COMPONENT') then

           jj = jj + 1

         else if(String(1:1) == '(') then

           jj = jj + 1

         end if

       end do liq11

       Nccom = jj

       allocate(Xcon(3,Nccom))
       allocate(Rcon(3,Nccom))
       allocate(Idcomf(Nccom),Idcomt(Nccom))
       allocate(fcom(3,Nccom))
       allocate(fscom(3,Nccom))
       allocate(InvMscom(Nccom))

       fcom(:,:) = 0.d0

       ii = Ist(i)

       jj = 0
liq12: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit liq12

         String = Script(ii)

         if(String(1:9) == 'COMPONENT') then

           read(String(10:),*) kk

           do k = 10, 80
             if(String(k:k) == '(') then
               Numbf1 = k
             else if(String(k:k) == ')') then
               Numbl1 = k
               exit
             end if
           end do

           jj = jj + 1
           if(Numbf1==0.or.Numbl1==0) then
             write(*,*) "error in FIXEDCOM-COMPONENT; parentheses may be missing"
             call Finalize
           end if
           read(String((Numbf1+1):(Numbl1-1)),*) Xcon(1,jj),Xcon(2,jj),Xcon(3,jj)

           l = 0
           do i1 = 1, NumSpec
             if(i1==kk) then
               Idcomf(jj) = l + 1
               Idcomt(jj) = l + NumAtm(i1)*NumMol(i1)
               exit
             end if
             l = l  + NumAtm(i1)*NumMol(i1)
           end do

         else if(String(1:1) == '(') then

           do k = 2, 80
             if(String(k:k) == '-') then
               Numhf1 = k
             else if(String(k:k) == ')') then
               Numbl1 = k
               exit
             end if
           end do

           jj = jj + 1
           if(Numhf1==0.or.Numbl1==0) then
             write(*,*) "error in FIXEDCOM; parentheses may be missing"
             call Finalize
           end if
           read(String(2:(Numhf1-1)),*) Idcomf(jj)
           read(String((Numhf1+1):(Numbl1-1)),*) Idcomt(jj)

           ini = Numbl1+1

           do k = ini, 80
             if(String(k:k) == '(') then
               Numbf1 = k
             else if(String(k:k) == ')') then
               Numbl1 = k
               exit
             end if
           end do

           if(Numbf1==0.or.Numbl1==0) then
             write(*,*) "error in FIXEDCOM; parentheses may be missing"
             call Finalize
           end if
           read(String((Numbf1+1):(Numbl1-1)),*) Xcon(1,jj),Xcon(2,jj),Xcon(3,jj)

         end if

       end do liq12

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QFixCOM) then
     do i = 1, Nccom
       InvMscom(i) = 0.d0
       do j = Idcomf(i), Idcomt(i)
         InvMscom(i) = InvMscom(i) + Mass(j)
       end do
       InvMscom(i) = 1.d0 / InvMscom(i)
     end do
   end if

end subroutine Read_FixedCOM_Cond


! #####################################################################
! #####################################################################


subroutine Read_Constraint_Mol

use ScriptData
use CommonBlocks, only : QMaster, QOption, Qwcformol, QDelCellMove
use wcparam
use Numbers, only : NumMol, NumAtm, NumSpec
use UnitExParam, only : ExParam
use AtomParam, only : Mass
use TimeParam, only : lk

implicit none

integer :: i, ii, j, jj, kk, k, ks
integer :: mol_init, mol_fin
character(len=1) :: Ch
character(len=10) :: Chap
real(8) :: Mconst_atom, wmx
integer :: CIdirection
external CIdirection

   Qwcformol = .False.
   kc_wall = 100.d0 * ExParam
   Nsample_wc = 1
   corder = 2
   Qwc_edge = .False.
   Nmolwc = 0

   do i = 1, NumOption

     if(OPTIONS(i)=='CONSTWALL') then

       ii = Ist(i)
       Qwcformol = .True.

lia11b:do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia11b

         if(Script(ii)(1:2) == 'k=') then
           read(Script(ii)(3:),*) kc_wall
           kc_wall = kc_wall * ExParam
         else if(Script(ii)(1:8) == 'NSAMPLE=') then
           read(Script(ii)(9:),*) Nsample_wc
         else if(Script(ii)(1:6) == 'ORDER=') then
           read(Script(ii)(7:),*) corder
         else if(Script(ii)(1:5) == 'AXIS=') then
           read(Script(ii)(6:),*) Ch
           Iaxis_const= CIdirection(Ch)
         else if(Script(ii)(1:7) == 'REGION=') then
           read(Script(ii)(8:),*) dw_low, Chap
           if((Chap(1:4)=='EDGE').or.(Chap(1:4)=='edge')) then
             Qwc_edge=.True.
           else
             read(Script(ii)(8:),*) dw_low, dw_upp
             dw_mid = 0.5d0*(dw_upp + dw_low)
             rc_low = dw_low - dw_mid
             rc_upp = dw_upp - dw_mid
           end if
         else if(Script(ii)(1:6)=='MOLS=') then
            read(Script(ii)(7:),*) mol_init, mol_fin
            Nmolwc = mol_fin - mol_init + 1
            allocate(idf(Nmolwc))
            allocate(ide(Nmolwc))
            allocate(invMwc(Nmolwc))
            jj = 0
            kk = 0
            ks = 0
            do j = 1, NumSpec
              do k = 1, NumMol(j)
                kk = kk + 1
                if((kk >= mol_init).and.(kk<=mol_fin)) then
                  ks = ks + 1
                  idf(ks) = jj + 1
                  ide(ks) = jj + NumAtm(j)
                end if
                jj = jj + NumAtm(j)
              end do
            end do
            do j = 1, Nmolwc
              wmx = 0.d0
              do k = idf(j), ide(j)
                wmx = wmx + Mass(k)
              end do
              invMwc(j) = 1.d0 / wmx
            end do
         else if(Script(ii)(1:6)=='ATOMS=') then
           read(Script(ii)(7:),*) atom_init, atom_fin
           Mconst_atom = 0.d0
           do j = atom_init, atom_fin
             Mconst_atom = Mconst_atom + Mass(j)
           end do
           InvMconst_atom = 1.d0 / Mconst_atom
         end if

       end do lia11b

       ReadFlag(i) = .True.

     end if

   end do

   if(Qwcformol) QDelCellMove = .True.

   Nsample_wc = Nsample_wc * lk

end subroutine Read_Constraint_Mol


! #####################################################################
! #####################################################################


subroutine Read_OptionalBond_Cond

use ScriptData
use CommonBlocks, only : QMaster, QOption, QPLC, QClust, QMacro, QJarzynski
use OptConstraintParam, only : NumOptC, OptCI, OptCJ, kOptC, rOptC, &
&   NumPLC, kPLC, rPLC, IaxisPLC, PLCI, &
&   NumClust, NAtomClus, Rsh_clust, k_clust, Ipartner
use UnitExParam, only : ExParam, cvol
use CGball, only : NumTypeSphere, RadiusSphere
use SMDparam, only : NumFreqConst, k_steering, Vshift_steering, &
&   R0_initial, R0_steering, FFsm
use TimeParam, only : lk, deltat

implicit none

integer :: i, ii, j, ic, ic0, jj
character(len=80) :: Script1
character(len=1) :: Ch
character(len=3) :: SwitchCls
integer :: CIdirection
external CIdirection

   QOption =.False.
   QPLC =.False.
   QClust =.False.
   QJarzynski =.False.
   NumOptC = 0
   NumPLC = 0
   NumClust = 0

   do i = 1, NumOption

     if(OPTIONS(i)=='OPTIONALBOND') then

       QOption =.True.

       ii = Ist(i)

lia12: do

         ii = ii + 1
         if(Script(ii)(1:2) == '<<') exit lia12

         NumOptC = NumOptC + 1

       end do lia12

       allocate( OptCI(NumOptC) )
       allocate( OptCJ(NumOptC) )
       allocate( kOptC(NumOptC) )
       allocate( rOptC(NumOptC) )

       do j = 1, NumOptC

         ii = Ist(i) + j

         read(Script(ii),*) OptCI(j),OptCJ(j),kOptC(j),rOptC(j)

       end do

       kOptC = kOptC * ExParam

       ReadFlag(i) = .True.

       exit

     end if

   end do


   do i = 1, NumOption

     if(OPTIONS(i)=='FIXONPLANE') then

       QPLC =.True.

       ii = Ist(i)

lia12b:do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia12b

         if(Script(ii)(1:5) == 'AXIS=') then
           read(Script(ii)(6:),*) Ch
           IaxisPLC= CIdirection(Ch)
         else
           NumPLC = NumPLC + 1
         end if

       end do lia12b

       allocate( PLCI(NumPLC) )
       allocate( kPLC(NumPLC) )
       allocate( rPLC(NumPLC) )

       ii = Ist(i)
       j = 0

lia12c:do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia12c

         if(Script(ii)(1:1) /= 'A') then
           j = j + 1
           read(Script(ii),*) PLCI(j), kPLC(j), rPLC(j)
         end if

       end do lia12c

       kPLC(:) = kPLC(:) * ExParam

       ReadFlag(i) = .True.

       exit

     end if

   end do

   do i = 1, NumOption

     if(OPTIONS(i)=='CLUSTER') then

       QClust =.True.

       ii = Ist(i)

       k_clust = 20.d0
       if(QMacro) then 
         if(NumTypeSphere /= 1) then
           write(*,*) "ERROR: you cannot use CLUSTER option when you have macroparticles"
           write(*,*) "       having different sizes."
           call Finalize
         end if
         Rsh_clust = RadiusSphere(1)*2.d0 + 2.3d0
       end if
       QJarzynski = .False.

lie12: do

         ii = ii + 1
         if(Script(ii)(1:2) == '<<') exit lie12

         if(Script(ii)(1:8) == 'NUMATOM=') then

           read(Script(ii)(9:),*) NumClust

         else if(Script(ii)(1:4) == 'RSH=') then

           read(Script(ii)(5:),*) Rsh_clust

         else if(Script(ii)(1:2) == 'k=') then

           read(Script(ii)(3:),*) k_clust

         else if(Script(ii)(1:3) == 'FE=') then

           read(Script(ii)(4:),*) SwitchCls
           if(SwitchCls(1:2)=='ON'.or.SwitchCls(1:2)=='on') then
             QJarzynski = .True.
           end if

         end if

       end do lie12

       allocate( NAtomClus(NumClust) )

       ii = Ist(i)

       ic0 = 0
lie13: do

         ii = ii + 1
         if(Script(ii)(1:2) == '<<') exit lie13

         if(Script(ii)(1:7) == 'MEMBER=') then

           write(Script1,*) trim(adjustl(Script(ii)(8:80)))
           jj = len(Script1)
           ic = 0
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           read(Script1,*) (NAtomClus(j),j=ic0+1,ic0+ic)
           ic0 = ic0 + ic
         end if

       end do lie13

       if(ic0/=NumClust) then
         if(QMaster) then
           write(*,*) "Number of particles in CLUSTER option may have problems"
           write(*,*) "NUMATOM=",NumClust," COUNTED=",ic0
         end if
         call Finalize
       end if

       if(QJarzynski) then

         ii = Ist(i)

lie14:   do

           ii = ii + 1
           if(Script(ii)(1:2) == '<<') exit lie14

           if(Script(ii)(1:5) == 'k_FE=') then

             read(Script(ii)(6:),*) k_steering

           else if(Script(ii)(1:5) == 'PAIR=') then

             read(Script(ii)(6:),*) Ipartner

           else if(Script(ii)(1:8) == 'NSAMPLE=') then

             read(Script(ii)(9:),*) NumFreqConst

           else if(Script(ii)(1:5) == 'RATE=') then

             read(Script(ii)(6:),*) Vshift_steering

           else if(Script(ii)(1:8) == 'INITIAL=') then

             read(Script(ii)(9:),*) R0_initial
             R0_steering = R0_initial

           end if

         end do lie14

         k_steering = k_steering * ExParam
         Vshift_steering = Vshift_steering * 1.d-03 * deltat ! [A/ns] --> [A/time_step_short]
         NumFreqConst = NumFreqConst * lk
         allocate( FFsm(1) )

       end if

       k_clust = k_clust * ExParam

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.QOption.and.Qcheck) then

     write(mfile,'(a)') 'OPTIONALBOND:'
     write(mfile,'(a,i10)')    '  Number = ',NumOptC
     do j = 1, NumOptC
     write(mfile,'(a,f10.3)') '  k_ham  = ',kOptC(j) * cvol,'[kcal/mol]'
     end do
     write(mfile,*)

   end if

   if(QMaster.and.QPLC.and.Qcheck) then

     write(mfile,'(a)') 'FIXONPLANE:'
     write(mfile,'(a,i10)')    '  Number = ',NumPLC
     do j = 1, NumPLC
     write(mfile,'(a,f10.3/)') '  k_ham  = ',kPLC(j) * cvol,'[kcal/mol]'
     end do

   end if

   if(QMaster.and.QClust.and.Qcheck) then

     write(mfile,'(a)') 'CLUSTER:'
     write(mfile,'(a,i10)')    '  Number = ',NumClust
     write(mfile,*) '  Members= ',NAtomClus(:)
     write(mfile,'(a,f10.3/)') '  k      = ',k_clust * cvol,'[kcal/mol]'
     write(mfile,'(a,f10.3/)') '  Rsh    = ',Rsh_clust,'[A]'

   end if

end subroutine Read_OptionalBond_Cond


! #####################################################################
! #####################################################################


subroutine Read_SPCF_Cond

use ScriptData
use CommonBlocks, only : QMaster, QSPCF
use Numbers, only : NumSpec, NumMol, NumAtm
use AtomParam, only : ResidName

implicit none

integer :: i, ii

   QSPCF = .False.

   ii = 0

   do i = 1, NumSpec

     if(ResidName(ii+1) == 'SPCF') then
       QSPCF = .True.
       exit
     else
       ii = ii + NumMol(i)*NumAtm(i)
     end if

   end do

end subroutine Read_SPCF_Cond


! #####################################################################
! #####################################################################


subroutine Read_OptionalConstraint

use ScriptData
use CommonBlocks, only : QMaster, QOpFix, QJarzynski, QOpJarz, QDelCellMove, &
&   QClust, QCorrCone
use SMDparam, only : NumFreqConst, NumOpFixS, NumOpFixP, NumFixF, NumFixL, &
&   MassRG, FixDir, RGConst, FFsm, NumFixFi, NumFixLi, NumFixFj, NumFixLj, &
&   MassRGi, MassRGj, ConstDis, ConstDis2, Rcomij, k_steering,             &
&   Uvec_steering, Vshift_steering, R0_steering, NSelCyl, NSelCylI, NSelCone, &
&   icaxis, icx1, icx2, IdSelC
use TimeParam, only : Nstep, lk
use UnitExParam, only : ExParam, pi
use CGdata, only : NumAtype, NBAtomType, AtomTypeList
use Numbers, only : N
use OptConstraintParam, only : Nccom, kccom, Xcon, Rcon, fscom, fcom, InvMscom, &
&   Rad_Lipo, Thick_Lipo, Vol_Lipo, facta, factb, factc, Invfacta
use AtomParam, only : Mass

implicit none

integer :: i, ii, j, k, jj, nsel, nnsel, ic
integer :: iiS, iiP
character(len=1) :: DIREC, Ch
character(len=80) :: Script1
character(len=6), dimension(:), allocatable :: SelCType
integer, dimension(:), allocatable :: idsel
integer :: Numbf1, Numbl1, Numhf1
integer :: Numbf2, Numbl2, Numhf2
real(8) :: pref, Rmt
integer :: CIdirection
external CIdirection

   QOpFix   =.False.
   QOpJarz  =.False.
   QDelCellMove = .False.
   QCorrCone = .False.
   NumOpFixS = 0
   NumOpFixP = 0

   Numbf1 = 0
   Numbl1 = 0
   Numhf1 = 0
   Numbf2 = 0
   Numbl2 = 0
   Numhf2 = 0

   if(QClust.and.QJarzynski) then
     Return
   end if

   NumFreqConst = Nstep + 1

   do i = 1, NumOption

     if(OPTIONS(i)=='OPTIONALCONSTRAINT') then

       QOpFix =.True.

       ii = Ist(i)

lib12: do

         ii = ii + 1
         if(Script(ii)(1:2) == '<<') exit lib12

         if(Script(ii)(1:4) == 'SNGL') then

           NumOpFixS = NumOpFixS + 1

         else if(Script(ii)(1:4) == 'PAIR') then

           NumOpFixP = NumOpFixP + 1

         else if(Script(ii)(1:8) == 'CYLINDER') then

           write(Script1,*) trim(adjustl(Script(ii)(9:80)))
           jj = len(Script1)
           ic = -2
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           allocate(SelCType(ic))
           read(Script1,*) Ch, R0_steering, (SelCType(j),j=1,ic)
           icaxis = CIdirection(Ch)
           NSelCyl = ic

           k_steering = 1000.
           QJarzynski = .True.

         else if(Script(ii)(1:8) == 'INNERCYL') then

           write(Script1,*) trim(adjustl(Script(ii)(9:80)))
           jj = len(Script1)
           ic = -2
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           allocate(SelCType(ic))
           read(Script1,*) Ch, R0_steering, (SelCType(j),j=1,ic)
           icaxis = CIdirection(Ch)
           NSelCylI = ic

           k_steering = 1000.
           QJarzynski = .True.

         else if(Script(ii)(1:4) == 'CONE') then

           write(Script1,*) trim(adjustl(Script(ii)(5:80)))
           jj = len(Script1)
           ic = -2
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           allocate(SelCType(ic))
           read(Script1,*) Ch, R0_steering, (SelCType(j),j=1,ic)
           icaxis = CIdirection(Ch)
           NSelCone = ic

           k_steering = 1000.
           QJarzynski = .True.

         else if(Script(ii)(1:11) == 'CORRECTCONE') then

           read(Script(ii)(12:),*) Rad_Lipo, Thick_Lipo
           QCorrCone = .True.
           Rmt = Rad_Lipo - Thick_Lipo
           Vol_Lipo = 4.d0 / 3.d0 * pi * (Rad_Lipo**3 - Rmt**3)
           pref  = 2.d0/3.d0*pi
           facta = pref * 3.d0 * Thick_Lipo
           Invfacta = 1.d0/facta
           factb = - facta * Thick_Lipo
           factc = pref * Thick_Lipo * Thick_Lipo * Thick_Lipo

         else if(Script(ii)(1:8) == 'NSAMPLE=') then

           read(Script(ii)(9:),*) NumFreqConst

         else if(Script(ii)(1:2) == 'k=') then

           read(Script(ii)(3:),*) k_steering
           QJarzynski = .True.

         end if

       end do lib12

       if(NumOpFixS /= 0) then

         allocate( NumFixF(NumOpFixS) )
         allocate( NumFixL(NumOpFixS) )
         allocate( MassRG(NumOpFixS) )
         allocate( FixDir(NumOpFixS) )
         allocate( RGConst(NumOpFixS) )
         allocate( FFsm(NumOpFixS) )

       end if

       if(NumOpFixP /= 0) then

         allocate( NumFixFi(NumOpFixP) )
         allocate( NumFixLi(NumOpFixP) )
         allocate( NumFixFj(NumOpFixP) )
         allocate( NumFixLj(NumOpFixP) )
         allocate( MassRGi(NumOpFixP) )
         allocate( MassRGj(NumOpFixP) )
         allocate( ConstDis(NumOpFixP) )
         allocate( ConstDis2(NumOpFixP) )
         allocate( Rcomij(3,NumOpFixP) )

       end if

       nsel = 0
       if(NSelCyl/=0) then
         nsel = NSelCyl
       else if(NSelCylI/=0) then
         nsel = NSelCylI
       else if(NSelCone/=0) then
         nsel = NSelCone
       end if
       if(nsel/=0) then
         allocate( FFsm(1) )
         jj = 0
         do j = 1, 3
           if(icaxis==j) cycle
           jj = jj + 1
           if(jj==1) icx1 = j
           if(jj==2) icx2 = j
         end do

         nnsel = 0
         if(SelCType(1) == 'all') then

           allocate(IdSelC(N))
           do j = 1, N
             IdSelC(j) = j
           end do
           nnsel = N

         else

           allocate(idsel(nsel))
           do j = 1, nsel
           do k = 1, NumAType
             if(SelCType(j)==AtomTypeList(k)) then
               idsel(j) = k
             end if
           end do
           end do
           jj = 0
           do j = 1, N
             do k = 1, nsel
               if(NBAtomType(j)==idsel(k)) then
                 jj = jj + 1
               end if
             end do
           end do
           nnsel = jj
           allocate(IdSelC(nnsel))
           jj = 0
           do j = 1, N
             do k = 1, nsel
               if(NBAtomType(j)==idsel(k)) then
                 jj = jj + 1
                 IdSelC(jj) = j
               end if
             end do
           end do

         end if

       end if

       if(NSelCyl/=0) then
         NSelCyl = nnsel
       else if(NSelCylI/=0) then
         NSelCylI = nnsel
       else if(NSelCone/=0) then
         NSelCone = nnsel
       end if

       if(NSelCone/=0) then
         Nccom = 1
         allocate(Xcon(3,Nccom))
         allocate(Rcon(3,Nccom))
         allocate(fcom(3,Nccom))
         allocate(fscom(3,Nccom))
         allocate(InvMscom(Nccom))

         fcom(:,:) = 0.d0
         fscom(:,:) = 0.d0

         kccom = 100.d0 * ExParam

         InvMscom(Nccom) = 0.d0
         do j = 1, NSelCone
           k = IdSelC(j)
           InvMscom(Nccom) = InvMscom(Nccom) + Mass(k)
         end do
         InvMscom(Nccom) = 1.d0 / InvMscom(Nccom)

         Xcon(:,Nccom) = 0.d0
       end if

       ii = Ist(i)
       iiS = 0
       iiP = 0

       do j = ii + 1, ii + NumOpFixS + NumOpFixP + 1

         if(Script(j)(1:4) == 'SNGL') then

           iiS = iiS + 1

           do k = 5, 80

             if(Script(j)(k:k) == '(') then
               Numbf1 = k
             else if(Script(j)(k:k) == '-') then
               Numhf1 = k
             else if(Script(j)(k:k) == ')') then
               Numbl1 = k
               exit
             end if

           end do

           if(Numbf1==0.or.Numhf1==0.or.Numbl1==0) then
             write(*,*) "error in Read_OptionalConstraint-SNGL; parentheses may be missing"
             call Finalize
           end if
           read(Script(j)((Numbf1+1):(Numhf1-1)),*) NumFixF(iiS)
           read(Script(j)((Numhf1+1):(Numbl1-1)),*) NumFixL(iiS)

           read(Script(j)((Numbl1+1):),*) DIREC, RGConst(iiS)

           FixDir(iiS) = CIdirection(DIREC)

         else if(Script(j)(1:4) == 'PAIR') then

           iiP = iiP + 1

           do k = 5, 80

             if(Script(j)(k:k) == '(') then
               Numbf1 = k
             else if(Script(j)(k:k) == '-') then
               Numhf1 = k
             else if(Script(j)(k:k) == ')') then
               Numbl1 = k
               exit
             end if

           end do

           do k = Numbl1+1, 80

             if(Script(j)(k:k) == '(') then
               Numbf2 = k
             else if(Script(j)(k:k) == '-') then
               Numhf2 = k
             else if(Script(j)(k:k) == ')') then
               Numbl2 = k
               exit
             end if

           end do

           if(Numbf1==0.or.Numhf1==0.or.Numbl1==0.or. &
           &  Numbf2==0.or.Numhf2==0.or.Numbl2==0) then
             write(*,*) "error in Read_OptionalConstraint-PAIR; parentheses may be missing"
             call Finalize
           end if
           read(Script(j)((Numbf1+1):(Numhf1-1)),*) NumFixFi(iiP)
           read(Script(j)((Numhf1+1):(Numbl1-1)),*) NumFixLi(iiP)
           read(Script(j)((Numbf2+1):(Numhf2-1)),*) NumFixFj(iiP)
           read(Script(j)((Numhf2+1):(Numbl2-1)),*) NumFixLj(iiP)

           read(Script(j)((Numbl2+1):),*) ConstDis(iiP)

           ConstDis2(iiP) = ConstDis(iiP) * ConstDis(iiP)

         end if

       end do

       ReadFlag(i) = .True.

     end if

   end do

!# unit exchange
   k_steering = k_steering * ExParam

! ##

   if(QJarzynski) then
     if(NumOpFixS/=1.and.NumOpFixP/=1.and.NSelCyl==0.and.NSelCone==0.and.&
     &  NSelCylI==0) then
       write(*,*) 'error : harmonic constraint is useful for single molecule or single pair'
       call Finalize
     end if
     if(NumOpFixS==1) then
       Uvec_steering(:) = 0.d0
       Uvec_steering(FixDir(1)) = 1.d0
       R0_steering = RGConst(1)
     else if(NumOpFixP==1) then
       R0_steering = ConstDis(1)
       allocate( FFsm(1) )
     end if
     Vshift_steering = 0.d0
   end if

   if(NumOpFixS/=0) QDelCellMove = .True.

   if(QMaster.and.QOpFix.and.Qcheck) then

     write(mfile,'(a    )') 'OPTIONAL CONSTRAINT:'
     write(mfile,'(a,i10)') 'NSAMPLE= ', NumFreqConst
     write(mfile,'(a,i3,a,i3)') 'N_SNGL= ',NumOpFixS, 'N_PAIR= ',NumOpFixP
 
     do i = 1, NumOpFixS

       j = FixDir(i)

       select case(j)
       case(1)
         DIREC='X'
       case(2)
         DIREC='Y'
       case(3)
         DIREC='Z'
       end select

       write(mfile,'(a,i5,a,i5,a,a,a,f7.1)') &
       & 'SNGL   (',NumFixF(i),'-',NumFixL(i),')  ',DIREC,'=',RGConst(i)

     end do

     do i = 1, NumOpFixP

       write(mfile,'(a,i5,a,i5,a,i5,a,i5,a,f7.2)') &
       & 'PAIR   (',NumFixFi(i),'-',NumFixLi(i),')  (',&
       & NumFixFj(i),'-',NumFixLj(i),')   R0=',ConstDis(i)

     end do

   end if

   NumFreqConst = NumFreqConst * lk

end subroutine Read_OptionalConstraint


! #####################################################################
! #####################################################################
!
! ## Reading Steered MD conditions for free energy calculation
! ## using Jarzynski's equality
!

subroutine Read_Jarzynski

use ScriptData
use CommonBlocks, only : QMaster, QOpFix, QJarzynski, QOpJarz, QDelCellMove, &
&   QClust, QCorrCone
use SMDparam, only : NumOpFixS, NumOpFixP, NumFreqConst, k_steering, &
&   Vshift_steering, R0_initial, R0_steering, Uvec_steering, NumFixF, &
&   NumFixL, MassRG, NumFixFi, NumFixLi, NumFixFj, NumFixLj, MassRGi, &
&   MassRGj, NSelCyl, NSelCylI, NSelCone, icaxis, icx1, icx2, IdSelC
use UnitExParam, only : ExParam, pi
use TimeParam, only : deltat, lk
use CGdata, only : NumAtype, NBAtomType, AtomTypeList
use Numbers, only : N
use OptConstraintParam, only : Nccom, kccom, Xcon, Rcon, fscom, fcom, InvMscom, &
&   Rad_Lipo, Thick_Lipo, Vol_Lipo, facta, factb, factc, Invfacta
use AtomParam, only : Mass

implicit none

integer :: i, ii, j, k, jj, nsel, nnsel, ic
integer :: iiS, iiP
character(len=1) :: Ch
integer :: Numbf1, Numbl1, Numhf1
integer :: Numbf2, Numbl2, Numhf2
real(8) :: R1
character(len=80) :: Script1
character(len=6), dimension(:), allocatable :: SelCType
real(8) :: Rmt, pref
integer :: CIdirection
integer, dimension(:), allocatable :: idsel
external CIdirection

   if(QJarzynski.and.QClust) then
     Return
   end if

   if(QOpFix) then
     if(QJarzynski) then
       QOpFix = .False.
       QOpJarz = .True.
     end if
     Return
   end if

   QJarzynski =.False.
   QCorrCone = .False.
   NumOpFixS  = 0
   NumOpFixP  = 0
   NumFreqConst = 10
   k_steering = 10.
   NSelCyl  = 0
   NSelCylI  = 0
   NSelCone = 0

   Numbf1 = 0
   Numbl1 = 0
   Numhf1 = 0
   Numbf2 = 0
   Numbl2 = 0
   Numhf2 = 0

   do i = 1, NumOption

     if(OPTIONS(i)=='JARZYNSKI') then

       QJarzynski =.True.

       ii = Ist(i)

lib13: do

         ii = ii + 1
         if(Script(ii)(1:2) == '<<') exit lib13

         if(Script(ii)(1:4) == 'SNGL') then

           NumOpFixS = NumOpFixS + 1

         else if(Script(ii)(1:4) == 'PAIR') then

           NumOpFixP = NumOpFixP + 1

         else if(Script(ii)(1:8) == 'CYLINDER') then

           write(Script1,*) trim(adjustl(Script(ii)(9:80)))
           jj = len(Script1)
           ic = -1
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           allocate(SelCType(ic))
           read(Script1,*) Ch, (SelCType(j),j=1,ic)
           icaxis = CIdirection(Ch)
           NSelCyl = ic

         else if(Script(ii)(1:8) == 'INNERCYL') then

           write(Script1,*) trim(adjustl(Script(ii)(9:80)))
           jj = len(Script1)
           ic = -1
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           allocate(SelCType(ic))
           read(Script1,*) Ch, (SelCType(j),j=1,ic)
           icaxis = CIdirection(Ch)
           NSelCylI = ic

         else if(Script(ii)(1:4) == 'CONE') then

           write(Script1,*) trim(adjustl(Script(ii)(5:80)))
           jj = len(Script1)
           ic = -1
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               ic = ic + 1
             end if
           end do
           allocate(SelCType(ic))
           read(Script1,*) Ch, (SelCType(j),j=1,ic)
           icaxis = CIdirection(Ch)
           NSelCone = ic

         else if(Script(ii)(1:11) == 'CORRECTCONE') then

           read(Script(ii)(12:),*) Rad_Lipo, Thick_Lipo
           QCorrCone = .True.
           Rmt = Rad_Lipo - Thick_Lipo
           Vol_Lipo = 4.d0 / 3.d0 * pi * (Rad_Lipo**3 - Rmt**3)
           pref  = 2.d0/3.d0*pi
           facta = pref * 3.d0 * Thick_Lipo
           Invfacta = 1.d0/facta
           factb = - facta * Thick_Lipo
           factc = pref * Thick_Lipo * Thick_Lipo * Thick_Lipo

         else if(Script(ii)(1:8) == 'NSAMPLE=') then

           read(Script(ii)(9:),*) NumFreqConst

         else if(Script(ii)(1:2) == 'k=') then

           read(Script(ii)(3:),*) k_steering

         else if(Script(ii)(1:5) == 'RATE=') then

           read(Script(ii)(6:),*) Vshift_steering

         else if(Script(ii)(1:8) == 'INITIAL=') then

           read(Script(ii)(9:),*) R0_initial
           R0_steering = R0_initial

         else if(Script(ii)(1:10) == 'DIRECTION=') then

           read(Script(ii)(11:),*) Uvec_steering(:)

           R1 = sqrt(dot_product(Uvec_steering,Uvec_steering))
           Uvec_steering(:) = Uvec_steering(:) / R1

         end if

       end do lib13

       if(NumOpFixS /= 0) then

         if(NumOpFixS /= 1) then
           write(*,*) 'ERROR : NumOpFixS /= 1 but',NumOpFixS
           call Finalize
         end if

         allocate( NumFixF(NumOpFixS) )
         allocate( NumFixL(NumOpFixS) )
         allocate( MassRG(NumOpFixS) )

       end if

       if(NumOpFixP /= 0) then

         if(NumOpFixP /= 1) then
           write(*,*) 'ERROR : NumOpFixP /= 1 but',NumOpFixP
           call Finalize
         end if

         allocate( NumFixFi(NumOpFixP) )
         allocate( NumFixLi(NumOpFixP) )
         allocate( NumFixFj(NumOpFixP) )
         allocate( NumFixLj(NumOpFixP) )
         allocate( MassRGi(NumOpFixP) )
         allocate( MassRGj(NumOpFixP) )

       end if

       nsel = 0
       if(NSelCyl/=0) then
         nsel = NSelCyl
       else if(NSelCylI/=0) then
         nsel = NSelCylI
       else if(NSelCone/=0) then
         nsel = NSelCone
       end if
       if(nsel/=0) then
         jj = 0
         do j = 1, 3
           if(icaxis==j) cycle
           jj = jj + 1
           if(jj==1) icx1 = j
           if(jj==2) icx2 = j
         end do

         if(SelCType(1) == 'all') then

           allocate(IdSelC(N))
           do j = 1, N
             IdSelC(j) = j
           end do

         else

           allocate(idsel(nsel))
           do j = 1, nsel
           do k = 1, NumAType
             if(SelCType(j)==AtomTypeList(k)) then
               idsel(j) = k
             end if
           end do
           end do
           nnsel = 0
           do j = 1, N
             do k = 1, nsel
               if(NBAtomType(j)==idsel(k)) then
                 nnsel = nnsel + 1
               end if
             end do
           end do
           allocate(IdSelC(nnsel))
           jj = 0
           do j = 1, N
             do k = 1, nsel
               if(NBAtomType(j)==idsel(k)) then
                 jj = jj + 1
                 IdSelC(jj) = j
               end if
             end do
           end do

         end if

       end if

       if(NSelCyl/=0) then
         NSelCyl = nnsel
       else if(NSelCylI/=0) then
         NSelCylI = nnsel
       else if(NSelCone/=0) then
         NSelCone = nnsel
       end if

       if(NSelCone/=0) then
         Nccom = 1
         allocate(Xcon(3,Nccom))
         allocate(Rcon(3,Nccom))
         allocate(fcom(3,Nccom))
         allocate(fscom(3,Nccom))
         allocate(InvMscom(Nccom))

         fcom(:,:) = 0.d0
         fscom(:,:) = 0.d0

         kccom = 100.d0 * ExParam

         InvMscom(Nccom) = 0.d0
         do j = 1, NSelCone
           k = IdSelC(j)
           InvMscom(Nccom) = InvMscom(Nccom) + Mass(k)
         end do
         InvMscom(Nccom) = 1.d0 / InvMscom(Nccom)

         Xcon(:,Nccom) = 0.d0
       end if

       ii = Ist(i)
       iiS = 0
       iiP = 0

       do j = ii + 1, ii + NumOpFixS + NumOpFixP

         if(Script(j)(1:4) == 'SNGL') then

           iiS = iiS + 1

           do k = 5, 80

             if(Script(j)(k:k) == '(') then
               Numbf1 = k
             else if(Script(j)(k:k) == '-') then
               Numhf1 = k
             else if(Script(j)(k:k) == ')') then
               Numbl1 = k
               exit
             end if

           end do

           if(Numbf1==0.or.Numhf1==0.or.Numbl1==0) then
             write(*,*) "error in JARZYNSKI-SNGL; parentheses may be missing"
             call Finalize
           end if
           read(Script(j)((Numbf1+1):(Numhf1-1)),*) NumFixF(iiS)
           read(Script(j)((Numhf1+1):(Numbl1-1)),*) NumFixL(iiS)

         else if(Script(j)(1:4) == 'PAIR') then

           iiP = iiP + 1

           do k = 5, 80

             if(Script(j)(k:k) == '(') then
               Numbf1 = k
             else if(Script(j)(k:k) == '-') then
               Numhf1 = k
             else if(Script(j)(k:k) == ')') then
               Numbl1 = k
               exit
             end if

           end do

           do k = Numbl1+1, 80

             if(Script(j)(k:k) == '(') then
               Numbf2 = k
             else if(Script(j)(k:k) == '-') then
               Numhf2 = k
             else if(Script(j)(k:k) == ')') then
               Numbl2 = k
               exit
             end if

           end do

           if(Numbf1==0.or.Numhf1==0.or.Numbl1==0.or. &
           &  Numbf2==0.or.Numhf2==0.or.Numbl2==0) then
             write(*,*) "error in JARZYNSKI-PAIR; parentheses may be missing"
             call Finalize
           end if
           read(Script(j)((Numbf1+1):(Numhf1-1)),*) NumFixFi(iiP)
           read(Script(j)((Numhf1+1):(Numbl1-1)),*) NumFixLi(iiP)
           read(Script(j)((Numbf2+1):(Numhf2-1)),*) NumFixFj(iiP)
           read(Script(j)((Numhf2+1):(Numbl2-1)),*) NumFixLj(iiP)

         end if

       end do

       ReadFlag(i) = .True.

     end if

   end do

   if(NumOpFixS/=0) QDelCellMove = .True.

! ## output for check
   if(QMaster.and.QJarzynski.and.Qcheck) then

     write(mfile,'(a    )') 'JARZYNSKI (Steered MD):'
     write(mfile,'(a,f10.3)') 'k=       ', k_steering
     write(mfile,'(a,f10.3)') 'INITIAL= ', R0_steering
     write(mfile,'(a,f10.3)') 'RATE=    ', Vshift_steering
     write(mfile,'(a,i10)') 'NSAMPLE= ', NumFreqConst
     write(mfile,'(a,i3,a,i3)') 'N_SNGL= ',NumOpFixS, 'N_PAIR= ',NumOpFixP
     if(NumOpFixS==1) then
       write(mfile,'(a,3f10.7)') 'vector=',Uvec_steering(:)
     end if

     do i = 1, NumOpFixS

       write(mfile,'(a,i5,a,i5,a)') &
       & 'SNGL   (',NumFixF(i),'-',NumFixL(i),')  '

     end do

     do i = 1, NumOpFixP

       write(mfile,'(a,i5,a,i5,a,i5,a,i5,a)') &
       & 'PAIR   (',NumFixFi(i),'-',NumFixLi(i),')  (',&
       & NumFixFj(i),'-',NumFixLj(i),')'

     end do

   end if

! ## unit exchange
   k_steering = k_steering * ExParam
   Vshift_steering = Vshift_steering * 1.d-03 * deltat ! [A/ns] --> [A/time_step_short]

   NumFreqConst = NumFreqConst * lk

end subroutine Read_Jarzynski


! #####################################################################
! #####################################################################


subroutine Read_ContInertia

use ScriptData
use ParamInertia
use CommonBlocks, only : QMaster, QInert
use Numbers, only : NumSpec, NumMol, NumAtm
use UnitExParam, only : ExParam, cvol

implicit none

integer :: i, ii, jj, j, icomp

   QInert = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='CONT_INERT') then

       QInert =.True.

       ii = Ist(i)

lif13: do

         ii = ii + 1
         if(Script(ii)(1:2) == '<<') exit lif13

         if(Script(ii)(1:7) == 'TARGET=') then

           read(Script(ii)(8:),*) icomp
           jj = 0
           do j = 1, NumSpec
             if(j==icomp) exit
             jj = jj + NumMol(j)*NumAtm(j)
           end do
           NumInert = NumMol(icomp)*NumAtm(icomp)
           allocate( ListInert(NumInert) )
           allocate( Msel(NumInert) )
           allocate( Rsel(3,NumInert) )
           do j = 1, NumInert
             jj = jj + 1
             ListInert(j) = jj
           end do

         else if(Script(ii)(1:2) == 'k=') then

           read(Script(ii)(3:),*) k_D
           k_D = k_D * ExParam

         else if(Script(ii)(1:3) == 'd0=') then

           read(Script(ii)(4:),*) TargetD

         end if

       end do lif13

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.QInert.and.Qcheck) then

     write(mfile,'(a)') 'CONT_INERT:'
     write(mfile,'(a,i10)')     '  Component = ',icomp
     write(mfile,'(a,f10.3,a)') '  k_D       = ',k_D * cvol,'[kcal/mol]'
     write(mfile,'(a,f10.3,a)') '  d0        = ',TargetD,'[-]'
     write(mfile,*)

   end if

end subroutine Read_ContInertia


! #####################################################################
! #####################################################################


subroutine Read_wall

use ScriptData
use CommonBlocks, only : QMaster, QEps_r, QPBC, QREPWALL, QSTWALL, &
&   QCGWALL
use WallParam, only : kcap, Ocap, Rcap, Ccap, Lcap, nwall, Rpcap, &
&   Ifcdir, sg_wall, ep_wall, dlattice_wall, STdir, nlayer, Rpore, &
&   Rstwall, FdirST, nstwall, WallFile, ngrid_wall, NSelR, IdSelR
use UnitExParam, only : ExParam, cvol
use Numbers, only : N
use AtomParam, only : NAAType, AATypeList, NBAAType

implicit none

integer :: i, ii, kk, j, nnsel, k, jj, nsel
integer :: CIdirection
real(8) :: dzz
external CIdirection
character(len=1) :: Ch
character(len=80) :: Script1
character(len=6), dimension(:), allocatable :: SelCType
integer, dimension(:), allocatable :: idsel
logical :: QSelect

   QEps_r = .False.

   if(.not.QPBC) then

     do i = 1, NumOption

       if(OPTIONS(i)=='RDEPENDEPS') then

         QEps_r = .True.

         ReadFlag(i) = .True.

         exit

       end if

     end do

   end if

   QREPWALL = .False.
   QSTWALL = .False.
   QCGWALL = .False.
   QSelect = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='REP_WALL') then

       QREPWALL = .True.

       write(Ccap,'(a)') 'SPHERE'
       kcap = 10.d0 * ExParam
       Ocap = 0.d0
       Rcap = 50.d0
       Rpore = 0.d0
       j = 0

       ii = Ist(i)

lia13: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia13

         if(Script(ii)(1:9) == 'BOUNDARY=') then

           read(Script(ii)(10:),*) Ccap

         else if(Script(ii)(1:2) == 'k=') then

           read(Script(ii)(3:),*) kcap
           kcap = kcap * ExParam

         else if(Script(ii)(1:9) == 'ATOMTYPE=') then

           write(Script1,*) trim(adjustl(Script(ii)(10:80)))
           jj = len(Script1)
           nsel = 0
           do j = 2, jj
             if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
               nsel = nsel + 1
             end if
           end do
           allocate(SelCType(nsel))
           read(Script1,*) (SelCType(j),j=1,nsel)
           QSelect = .True.

         end if

       end do lia13

       ii = Ist(i)

       if(Ccap=='SPHERE') then

lia14:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia14

           if(Script(ii)(1:7) == 'CENTER=') then

             read(Script(ii)(8:),*) Ocap

           else if(Script(ii)(1:7) == 'RADIUS=') then

             read(Script(ii)(8:),*) Rcap

           end if

         end do lia14

       else if(Ccap=='PLANER') then

lia15:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia15

           if(Script(ii)(1:10) == 'DIRECTION=') then

             read(Script(ii)(11:),*) Ch

             Lcap = CIdirection(Ch)

           else if(Script(ii)(1:8) == 'NUMWALL=') then

             read(Script(ii)(9:),*) nwall
             allocate( Rpcap(nwall) )
             allocate( Ifcdir(nwall) )

           else if(Script(ii)(1:5) == 'PORE=') then

             read(Script(ii)(6:),*) Rpore

           end if

         end do lia15

         kk = 0
         ii = Ist(i)

lib16:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lib16

           if(Script(ii)(1:5) == 'WALL=') then

             kk = kk + 1

             read(Script(ii)(6:),*) Rpcap(kk), Ch

             if(Ch=='+') then
               Ifcdir(kk)= 1
             else if(Ch=='-') then
               Ifcdir(kk)=-1
             else
               if(QMaster) write(*,*) 'Error : Force direction of WALL should be "+" or "-"'
               call Finalize
             end if

             if(kk==nwall) exit lib16

           end if

         end do lib16

       end if

       if(.not.QSelect) then

         allocate(IdSelR(N))
         do j = 1, N
           IdSelR(j) = j
         end do
         nnsel = N

       else if(SelCType(1) == 'all') then

         allocate(IdSelR(N))
         do j = 1, N
           IdSelR(j) = j
         end do
         nnsel = N

       else

         allocate(idsel(nsel))
         do j = 1, nsel
         do k = 1, NAAType
           if(SelCType(j)==AATypeList(k)) then
             idsel(j) = k
           end if
         end do
         end do
         jj = 0
         do j = 1, N
           do k = 1, nsel
             if(NBAAType(j)==idsel(k)) then
               jj = jj + 1
             end if
           end do
         end do
         nnsel = jj

         allocate(IdSelR(nnsel))
         jj = 0
         do j = 1, N
           do k = 1, nsel
             if(NBAAType(j)==idsel(k)) then
               jj = jj + 1
               IdSelR(jj) = j
             end if
           end do
         end do

       end if

       NSelR = nnsel

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.QREPWALL.and.Qcheck) then

     write(mfile,'(a  )') 'REPULSIVE WALL:'
     write(mfile,'(a,a)')        '  METHOD = ', trim(Ccap)
     write(mfile,'(a,f10.3,a)')  '  k      = ', kcap * cvol, '[kcal/mol]'
     if(Ccap=='SPHERE') then
     write(mfile,'(a,3f10.2)')     '  Center    = ', Ocap
     write(mfile,'(a,f15.5)')       '  Radius    = ', Rcap
     else if(Ccap=='PLANER') then
     write(mfile,'(a,i5,a)')     '  DIRECTION = ', Lcap,'[x:1,y:2,z:3]'
     write(mfile,'(a,i5)')       '  NUMWALL   = ', nwall
     if(Rpore/=0.) then
     write(mfile,'(a,f5.1,a)')     '  PORE      = ', Rpore,'[A]'
     end if
     do i = 1, nwall
     write(mfile,'(a,i1,a,f10.2,i3)') '  WALL',i,'     = ', Rpcap(i), Ifcdir(i)
     end do
     write(mfile,*)
     end if

   end if

   do i = 1, NumOption

     if(OPTIONS(i)=='ST_WALL') then

       QSTWALL = .True.

       sg_wall = 3.4d0
       ep_wall = 0.0556d0 * ExParam
       dlattice_wall = 2.46d0
       dzz = 3.4d0

       ii = Ist(i)

lig13: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lig13

         if(Script(ii)(1:8) == 'NUMWALL=') then

           read(Script(ii)(9:),*) nstwall
           allocate( Rstwall(nstwall) )
           allocate( FdirST(nstwall) )

         else if(Script(ii)(1:10) == 'DIRECTION=') then

           read(Script(ii)(11:),*) Ch

           STdir = CIdirection(Ch)

         else if(Script(ii)(1:6) == 'LAYER=') then

           read(Script(ii)(7:),*) nlayer

         else if(Script(ii)(1:6) == 'SIGMA=') then

           read(Script(ii)(7:),*) sg_wall

         else if(Script(ii)(1:8) == 'EPSILON=') then

           read(Script(ii)(9:),*) ep_wall
           ep_wall = ep_wall * ExParam

         else if(Script(ii)(1:8) == 'LATTICE=') then

           read(Script(ii)(9:),*) dlattice_wall

         else if(Script(ii)(1:3) == 'DZ=') then

           read(Script(ii)(4:),*) dzz

         end if

       end do lig13

       ii = Ist(i)

       kk = 0

lig16: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lig16

         if(Script(ii)(1:5) == 'WALL=') then

           kk = kk + 1

           read(Script(ii)(6:),*) Rstwall(kk), Ch

           if(Ch=='+') then
             FdirST(kk)= 1
           else if(Ch=='-') then
             FdirST(kk)=-1
           else
             if(QMaster) write(*,*) 'Error : Force direction of WALL should be "+" or "-"'
             call Finalize
           end if

           if(kk==nstwall) exit lig16

         end if

       end do lig16

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.QSTWALL.and.Qcheck) then

     write(mfile,'(a  )') 'Steele WALL:'
     write(mfile,'(a,i5,a)')     '  DIRECTION = ', STdir,'[x:1,y:2,z:3]'
     write(mfile,'(a,i5)')       '  NUMWALL   = ', nstwall
     write(mfile,'(a,i5)')       '  LAYER     = ', nlayer
     write(mfile,'(a,f10.3,a)')  '  SIGMA     = ', sg_wall, '[A]'
     write(mfile,'(a,f10.3,a)')  '  EPSILON   = ', ep_wall * cvol, '[kcal/mol]'
     write(mfile,'(a,f10.3,a)')  '  LATTICE   = ', dlattice_wall, '[A]'
     write(mfile,'(a,f10.3,a)')  '  DZ        = ', dzz, '[A]'
     do i = 1, nstwall
     write(mfile,'(a,i1,a,f10.2,i3)') '  WALL',i,'     = ', Rstwall(i), FdirST(i)
     end do
     write(mfile,*)

   end if

! ## CGwall
   do i = 1, NumOption

     if(OPTIONS(i)=='CG_WALL') then

       QCGWALL = .True.

       WallFile = 'cgwall.prm'

       ii = Ist(i)

lig83: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lig83

         if(Script(ii)(1:8) == 'NUMWALL=') then

           read(Script(ii)(9:),*) nstwall
           allocate( Rstwall(nstwall) )
           allocate( FdirST(nstwall) )

         else if(Script(ii)(1:11) == 'GRIDPOINTS=') then

           read(Script(ii)(12:),*) ngrid_wall

         else if(Script(ii)(1:10) == 'DIRECTION=') then

           read(Script(ii)(11:),*) Ch

           STdir = CIdirection(Ch)

         else if(Script(ii)(1:5) == 'FILE=') then

           read(Script(ii)(6:),*) WallFile

         end if

       end do lig83

       ii = Ist(i)

       kk = 0

lig86: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lig86

         if(Script(ii)(1:5) == 'WALL=') then

           kk = kk + 1

           read(Script(ii)(6:),*) Rstwall(kk), Ch

           if(Ch=='+') then
             FdirST(kk)= 1
           else if(Ch=='-') then
             FdirST(kk)=-1
           else
             if(QMaster) write(*,*) 'Error : Force direction of WALL should be "+" or "-"'
             call Finalize
           end if

           if(kk==nstwall) exit lig86

         end if

       end do lig86

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.QCGWALL.and.Qcheck) then

     write(mfile,'(a  )') 'CG-WALL:'
     write(mfile,'(a,i5,a)')     '  DIRECTION = ', STdir,'[x:1,y:2,z:3]'
     write(mfile,'(a,i5)')       '  NUMWALL   = ', nstwall
     do i = 1, nstwall
     write(mfile,'(a,i1,a,f10.2,i3)') '  WALL',i,'     = ', Rstwall(i), FdirST(i)
     end do
     write(mfile,*)

   end if

end subroutine Read_wall


! #####################################################################
! #####################################################################


subroutine Read_Cylinder

use ScriptData
use CommonBlocks, only : QCyl, QFSCyl, ForceField, QMaster
use CylParam

implicit none

integer :: i, ii, j, jj
character(len=1) :: Ch
integer :: CIdirection
external CIdirection

   QCyl = .False.
   QFSCyl = .False.
   ngrid_cyl = 2000
   ngrid_cylZ = 2000

   do i = 1, NumOption

     if(OPTIONS(i)=='CYLINDER') then

       QCyl = .True.

       ii = Ist(i)

liw24: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit liw24

         if(Script(ii)(1:5) == 'FILE=') then

           read(Script(ii)(6:),*) CylFile

         else if(Script(ii)(1:11) == 'GRIDPOINTS=') then

           read(Script(ii)(12:),*) ngrid_cyl

         else if(Script(ii)(1:8) == 'NGRID_R=') then

           read(Script(ii)(9:),*) ngrid_cyl

         else if(Script(ii)(1:8) == 'NGRID_Z=') then

           read(Script(ii)(9:),*) ngrid_cylZ

         else if(Script(ii)(1:10) == 'DIRECTION=') then

           read(Script(ii)(11:),*) Ch

           icaxis = CIdirection(Ch)
           jj = 0
           do j = 1, 3
             if(icaxis==j) cycle
             jj = jj + 1
             if(jj==1) icx1 = j
             if(jj==2) icx2 = j
           end do

         else if(Script(ii)(1:14) == 'FINITE_HEIGHT=') then

           read(Script(ii)(15:),*) Ch

           if(Ch=='Y'.or.Ch=='y') then
             QFSCyl = .True.
             QCyl = .False.
           end if

         end if

       end do liw24

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.(QCyl.or.QFSCyl).and.Qcheck) then

     write(mfile,'(a  )') 'CYLINDER WALL:'
     write(mfile,'(a,i5,a)')     '  DIRECTION   = ', icaxis,'[x:1,y:2,z:3]'
     write(mfile,'(a,i5)')       '  DATA POINTS = ', ngrid_cyl
     write(mfile,*)

   end if

end subroutine Read_Cylinder


! #####################################################################
! #####################################################################


subroutine Read_Macroparticle

use ScriptData
use CommonBlocks, only : QMacro, ForceField, QMaster
use IOparam, only : Sphere_file
use CGball
use UnitExParam, only : ExParam, kb, cvol

implicit none

integer :: i, ii, j

   QMacro = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='MACROSPHERE') then

       QMacro = .True.

       NumSphere = 1
       NumTypeSphere = 1
       RhoSphere = 0.113d0
       Sphere_file = 'cgsphere.prm'
       sigMM = 3.55d0
       epsMM = 35.23d0 * kb * cvol ! [kcal/mol]

       ii = Ist(i)

lic24: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lic24

         if(Script(ii)(1:7) == 'NUMBER=') then
           read(Script(ii)(8:),*) NumSphere
         end if
         if(Script(ii)(1:9) == 'NDENSITY=') then
           read(Script(ii)(10:),*) RhoSphere
         end if
         if(Script(ii)(1:9) == 'MULTIRAD=') then
           read(Script(ii)(10:),*) NumTypeSphere
         end if
         if(Script(ii)(1:5) == 'FILE=') then
           read(Script(ii)(6:),*) Sphere_file
         end if
         if(Script(ii)(1:8)== 'EPSILON=') then
           read(Script(ii)(9:),*) eps_CO
         end if
         if(Script(ii)(1:11)== 'MACROMACRO=') then
           read(Script(ii)(12:),*) epsMM, sigMM
         end if

       end do lic24

       allocate(RadiusSphere(NumTypeSphere))

       ii = Ist(i)

lie25: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lie25

         if(Script(ii)(1:7) == 'RADIUS=') then
           read(Script(ii)(8:),*) (RadiusSphere(j), j = 1 , NumTypeSphere)
         end if

       end do lie25

       ReadFlag(i) = .True.

     end if

   end do

   if(QMacro) then
     if(NumSphere==0) then
       if(QMaster) write(*,*) "ERROR : NumSphere should not be zero"
       call Finalize
     end if
     if(ForceField/='CG'.and.ForceField/='OPLS') then
       if(QMaster) write(*,*) &
       & 'ERROR: MACROSPHERE option is useful only with CG or OPLS force field'
       call Finalize
     end if
     if(ForceField=='OPLS') then
       eps_CO = eps_CO * ExParam
     end if
     epsMM = epsMM * ExParam
   end if


end subroutine Read_Macroparticle


! #####################################################################
! #####################################################################


subroutine Read_Hill

use ScriptData
use CommonBlocks, only : QMaster, QMP, QPBC
use IOparam, only : Extmem_file

implicit none

integer :: i, ii
character(len=80) :: Chline
logical :: Qfhill

   QMP = .False.
   Qfhill = .False.

   if(.not.QPBC) then

     do i = 1, NumOption

       if(OPTIONS(i)=='EXT_MEMP') then

         ii = Ist(i) + 1
         if(Script(ii)(1:5) == 'FILE=') then
           read(Script(ii)(6:),*) Chline
           Qfhill = .True.
           Extmem_file = trim(adjustl(Chline))
         end if

         QMP = .True.

         ReadFlag(i) = .True.

         exit

       end if

     end do

   end if

   if((QMP).and.(.not.Qfhill)) then
     write(Extmem_file,'(a)') './param/par_externalmem.prm'
   end if

end subroutine Read_Hill


! #####################################################################
! #####################################################################


subroutine Read_Eflux_mon

use ScriptData
use CommonBlocks, only : QMaster, QEflux
use Conduct

   QEflux = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='EFLUX') then

       QEflux = .True.

       id_cation = 0
       id_anion = 0

       ii = Ist(i)

lid20: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lid20

         if(Script(ii)(1:7) == 'CATION=') then
           read(Script(ii)(8:),*) id_cation
         end if

         if(Script(ii)(1:6) == 'ANION=') then
           read(Script(ii)(7:),*) id_anion
         end if

       end do lid20

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QEflux) then

       write(mfile,'(a)'     ) 'electrical current flag is valid'
       write(mfile,'(a,i3  )') '  CATION = ',id_cation
       write(mfile,'(a,i3  )') '  ANION  = ',id_anion

     end if

   end if

end subroutine Read_Eflux_mon


! #####################################################################
! #####################################################################


subroutine Read_Efield

use ScriptData
use CommonBlocks, only : QMaster, QEfield
use Conduct
use UnitExParam, only : reng, e

real(8) :: EV
real(8), dimension(3) :: Evector
integer :: i, ii

   QEfield = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='EFIELD') then

       QEfield = .True.

       Efield = 0.

       ii = Ist(i)


lid21: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lid21

         if(Script(ii)(1:7) == 'EFIELD=') then
           read(Script(ii)(8:),*) Efield  ! [V/A]
         end if

         if(Script(ii)(1:10) == 'DIRECTION=') then
           read(Script(ii)(11:),*) Evector(:)
         end if

       end do lid21

       ReadFlag(i) = .True.

       exit

     end if

   end do

   EV = dot_product(Evector,Evector)
   EV = sqrt(EV)
   Evector(:) = Evector(:)/EV

   Eforce(:) = Efield * reng * e * Evector(:)

   if(QMaster.and.QEfield.and.Qcheck) then

     write(mfile,'(a  )') 'Applied Electric Field:'
     write(mfile,'(a,3f7.4,a)')       '  DIRECTION      = {', Evector(:),'}'
     write(mfile,'(a,e12.4,a)') '  Electric Field = ', Efield,'[V/Angstrom]'
     write(mfile,*)

   end if

end subroutine Read_Efield


! #####################################################################
! #####################################################################


subroutine Read_F_monitor

use ScriptData
use CommonBlocks, only : QMaster, QFmonitor
use F_monitor, only : isampleF, StoreFrc, NiniF, NfinF, NgA
use Numbers, only : NumMol, NumAtm
use TimeParam, only : Nstep, lk

implicit none

integer :: i, ii, j
integer :: compo

   QFmonitor = .False.

   isampleF = Nstep + 1

   do i = 1, NumOption

     if(OPTIONS(i)=='F_MONITOR') then

       QFmonitor = .True.

       compo    = 1
       StoreFrc = 10000

       ii = Ist(i)

lic20: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lic20

         if(Script(ii)(1:10) == 'COMPONENT=') then
           read(Script(ii)(11:),*) compo
         end if
         if(Script(ii)(1:5) == 'FREQ=') then
           read(Script(ii)(6:),*) isampleF
         end if
         if(Script(ii)(1:5) == 'DATA=') then
           read(Script(ii)(6:),*) StoreFrc
         end if

       end do lic20

       ReadFlag(i) = .True.

       NiniF = 0
       NfinF = 0

       if(compo /= 1) then
         do j = 1, compo-1
           NiniF = NiniF + NumMol(j)*NumAtm(j)
         end do
       end if

       NgA = NumMol(compo)*NumAtm(compo)

       NfinF = NiniF + NgA

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QFmonitor) then

       write(mfile,'(a)'     ) 'Force-monitoring flag is valid'
       write(mfile,'(a,i3  )') '  COMPONENT = ',compo
       write(mfile,'(a,i3  )') '  FREQ      = ',isampleF
       write(mfile,'(a,i10/)') '  DATA      = ',StoreFrc

     end if

   end if

   isampleF = isampleF * lk

end subroutine Read_F_monitor


! #####################################################################
! #####################################################################


subroutine Read_TICR

use ScriptData
use CommonBlocks, only : QMaster, QTICR
use FEparam
use TimeParam, only : Nstep, lk
use Numbers, only : N, NumSpec, NumMol, NumAtm

implicit none

integer :: i, ii, j
character(len=12) :: TICRflag, INTERflag

   QTICR = .False.

   isampleTI = Nstep + 1

   do i = 1, NumOption

     if(OPTIONS(i)=='TICR') then

       QTICR = .True.
       QCreation = .True.
       QLJ = .False.
       QCharge = .False.

       NTIspec  = NumSpec
       NumTIMol = NumMol(NTIspec)

       iTIstep = 0
       isampleTI = 1

       Shift_Param = 0.2d0

       ii = Ist(i)

lic21: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lic21

         if(Script(ii)(1:10) == 'COMPONENT=') then
           read(Script(ii)(11:),*) NTIspec
         end if
         if(Script(ii)(1:7) == 'ACTION=') then
           read(Script(ii)(8:),*) TICRflag
           if(trim(TICRflag) == 'ANNIHILATION') then
             QCreation = .False.
           else if(trim(TICRflag) == 'CREATION') then
             QCreation = .True.
           else
             if(QMaster) then
               write(*,*) 'ERROR : ACTION in TICR'
               write(*,*) 'select "CREATION" or "ANNIHILATION"'
             end if
             call Finalize
           end if
         end if
         if(Script(ii)(1:12) == 'INTERACTION=') then
           read(Script(ii)(13:),*) INTERflag
           if(trim(INTERflag)=='LJ') then
             QLJ = .True.
             QCharge = .False.
           else if(trim(INTERflag)=='COULOMB') then
             QLJ = .False.
             QCharge = .True.
           else if(trim(INTERflag)=='ALL') then
             QLJ = .True.
             QCharge = .True.
           else
             if(QMaster) then
               write(*,*) 'ERROR : INTERACTION in TICR'
               write(*,*) 'select "LJ", "COULOMB" or "ALL"'
             end if
             call Finalize
           end if
         end if
         if(Script(ii)(1:6) == 'DELTA=') then
           read(Script(ii)(7:),*) Shift_Param
         end if
         if(Script(ii)(1:8) == 'MOLNUMB=') then
           read(Script(ii)(9:),*) NumTIMol
         end if
         if(Script(ii)(1:5) == 'FREQ=') then
           read(Script(ii)(6:),*) isampleTI
         end if
         if(Script(ii)(1:8) == 'NUMSTEP=') then
           read(Script(ii)(9:),*) iTIstep

           allocate( iequil_TICR(iTIstep) )
           allocate( isample_TICR(iTIstep) )
           allocate( TIlambda(iTIstep) )
         end if

       end do lic21

       if(iTIstep == 0) then
         write(*,*) 'ERROR : "NUMSTEP=" must be described in the option, "TICR"'
         call Finalize
       end if

       ii = Ist(i)

       isample_TICR(:) = 0

lic22: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lic22

         if(Script(ii)(1:5) == 'STEP=') then
           read(Script(ii)(6:),*) j, TIlambda(j), iequil_TICR(j), isample_TICR(j)
         end if

       end do lic22

       do j = 1, iTIstep
         if(isample_TICR(j) == 0) then
           write(*,*) 'ERROR : "STEP=" information is incorrect for ',j,'-th step'
           call Finalize
         end if
       end do

       ii = 0

       do j = 1, iTIstep
         ii = ii + iequil_TICR(j) + isample_TICR(j)
       end do

       if(ii /= Nstep) then
         if(QMaster) write(*,*) &
         & 'ERROR : total step number is not consistent with the TICR sampling number'
         call Finalize
       end if

       ReadFlag(i) = .True.

       NiniTI = 0
       NfinTI = 0

       if(NTIspec /= 1) then
         do j = 1, NTIspec-1
           NiniTI = NiniTI + NumMol(j)*NumAtm(j)
         end do
       end if

       if(NumTIMol /= 1) then
         do j = 1, NumTIMol-1
           NiniTI = NiniTI + NumAtm(NTIspec)
         end do
       end if

       NumTIpart = NumAtm(NTIspec)

       NfinTI = NiniTI + NumTIpart

       allocate( QTICRcalc(N) )
       QTICRcalc = .True.
       do j = NiniTI+1, NfinTI
         QTICRcalc(j) = .False.
       end do

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QTICR) then

       write(mfile,'(a)'     ) 'TICR flag is valid'
       write(mfile,'(a,i3  )') '  COMPONENT = ',NTIspec
       write(mfile,'(a,i3  )') '  MOLNUMB   = ',NumTIMol
       write(mfile,'(a,i3  )') '  FREQ      = ',isampleTI
       write(mfile,'(a,i3  )') '  NUMSTEP   = ',iTIstep
       do j = 1, iTIstep
         write(mfile,'(a,f7.4,2i7)') '  STEP = ', TIlambda(j), iequil_TICR(j), &
         &                                                     isample_TICR(j)
       end do

     end if

   end if

   isampleTI = isampleTI * lk

end subroutine Read_TICR


! #####################################################################
! #####################################################################


subroutine Read_HeatCapCond

use ScriptData
use CommonBlocks, only : QMaster, QHeatCap
use HeatCap, only : Hessian
use Numbers, only : N
use TimeParam, only : lk, isampleHC

implicit none

integer :: i, ii

   QHeatCap = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='HEATCAPACITY') then

       QHeatCap = .True.

       isampleHC = 1

       ii = Ist(i)

lic23: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lic23

         if(Script(ii)(1:5) == 'FREQ=') then
           read(Script(ii)(6:),*) isampleHC
         end if

       end do lic23

       ReadFlag(i) = .True.

       allocate( Hessian(3*N,3*N) )

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     if(QHeatCap) then
       write(mfile,'(a)'     ) 'Heat Capacity flag is valid'
       write(mfile,'(a,i3/ )') '  FREQ      = ',isampleHC
     end if

   end if

   isampleHC = isampleHC * lk

end subroutine Read_HeatCapCond


! #####################################################################
! #####################################################################


subroutine Read_VolScale

use ScriptData
use CommonBlocks, only : QMaster, QVolScale, QPBC
use CellParam, only : Vsc_Rate

implicit none

integer :: i, ii

   QVolScale = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='VOL_SCALE') then

       QVolScale = .True.

       Vsc_Rate = 0.d0

       ii = Ist(i)

lic25: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lic25

         if(Script(ii)(1:7) == 'RATE_X=') then
           read(Script(ii)(8:),*) Vsc_Rate(1)
         else if(Script(ii)(1:7) == 'RATE_Y=') then
           read(Script(ii)(8:),*) Vsc_Rate(2)
         else if(Script(ii)(1:7) == 'RATE_Z=') then
           read(Script(ii)(8:),*) Vsc_Rate(3)
         end if

       end do lic25

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QVolScale.and.(.not.QPBC)) then
     if(QMaster) &
     & write(*,*) 'Error : VOL_SCALE is useful only with the priodic boundary condition'
     call Finalize
   end if

   if(QMaster.and.Qcheck) then

     if(QVolScale) then

       write(mfile,'(a)'     )  'Volume Scaling flag is valid'
       write(mfile,'(a,f7.2)')  '  Rate_X      = ',Vsc_Rate(1)
       write(mfile,'(a,f7.2)')  '  Rate_Y      = ',Vsc_Rate(2)
       write(mfile,'(a,f7.2/)') '  Rate_Z      = ',Vsc_Rate(3)

     end if

   end if

end subroutine Read_VolScale


! #####################################################################
! #####################################################################


subroutine Calc_Number

use ScriptData
use CommonBlocks, only : QMaster, QSHAKE, QRigidBody, cThermostatMethod, &
&   QPBC
use SHAKEparam, only : NSHAKEGroup, NCoupleBond
use RBparam, only : NumRB, QSingle, QLinear
use Numbers, only : NfT, NfR, Nf, Number
use BathParam, only : InvNf, kT, gkT

implicit none

integer :: i
integer :: NSK

   NSK= 0

   if(QSHAKE) then

     do i = 1 , NSHAKEGroup

       NSK = NSK + NCoupleBond(i)

     end do

   end if

   if(QRigidBody) then

     NfT = 3 * NumRB

     NfR = 0

     do i = 1 , NumRB

       if(QSingle(i)) cycle

       if(QLinear(i)) then

         NfR = NfR + 2

       else

         NfR = NfR + 3

       end if

     end do

   else

     NfT = 3 * Number - NSK
     NfR = 0

   end if

   Nf = NfT + NfR

   if(cThermostatMethod /= 'MNHC') then
     Nf = Nf - 3
     if(.not.QPBC) Nf = Nf - 3 ! elimination of angular momentum of the system
   end if

   if(QMaster.and.Qcheck) then

     write(mfile,'(a,i10/)') 'Number of degrees of freedom = ',Nf

   end if

   InvNf = 1.d0 / dble(Nf)

   gkT  = Nf * kT

end subroutine Calc_Number


! #####################################################################
! #####################################################################


subroutine Calc_BathMass

use ScriptData
use CommonBlocks, only : QMaster, QBarostat, cBarostatMethod, QThermostat, &
&   cThermostatMethod
use BathParam, only : kT, gkT, Tau_p, Mp, prefTaup, InvMp, Mts, InvMts, &
&   Rss, Vss, NHchain, NumMNHC, RMNHC, VMNHC, MMNHC, InvMMNHC, Tau_s0, &
&   Tau_s1
use Numbers, only : Nf, NfT
use UnitExParam, only : pi

implicit none

integer :: i, j
real(8) :: Omega_p, Omega_s0, Omega_s1

   if(QBarostat) then

     Omega_p  = 2.d0 * pi / Tau_p

     if( cBarostatMethod == 'AN' ) then

       gkT  = ( Nf + 1.d0 ) * kT
       Mp   = ( NfT + 3.d0 ) * kT / (Omega_p * Omega_p)

     else if( cBarostatMethod == 'A2' ) then

       gkT  = ( Nf + 2.d0 ) * kT
       Mp   = ( NfT + 3.d0 ) * kT / (3.d0 * Omega_p * Omega_p)

     else if( cBarostatMethod == 'A3' ) then

       gkT  = ( Nf + 3.d0 ) * kT
       Mp   = ( NfT + 3.d0 ) * kT / (3.d0 * Omega_p * Omega_p)

     else if( (cBarostatMethod == 'PR').or.( cBarostatMethod == 'ST' ) ) then

       gkT  = ( Nf + 3.d0**2 ) * kT
       Mp   = ( NfT + 3.d0 ) * kT / (3.d0 * Omega_p * Omega_p)

     end if

     Mp = Mp * prefTaup * prefTaup

     InvMp = 1.d0 / Mp

   end if

   if(QThermostat) then

   if((cThermostatMethod == 'NH') .or. &
   &  (cThermostatMethod == 'NHC')) then

     allocate( Mts(NHchain) )
     allocate( InvMts(NHchain) )

     Omega_s0 = 2.d0 * pi / Tau_s0

     if(cThermostatMethod == 'NH') then

       Mts(1) = gkT / ( Omega_s0 * Omega_s0 )

       allocate( Rss(NHchain) )
       allocate( Vss(NHchain) )

     else if(cThermostatMethod == 'NHC') then

       Omega_s1 = 2.d0 * pi / Tau_s1
       do i = 2 , NHchain
         Mts(i) = kT / ( Omega_s1 * Omega_s1 )
       end do

       Mts(1) = gkT / ( Omega_s0 * Omega_s0 )

       allocate( Rss(NHchain) )
       allocate( Vss(NHchain) )

     end if

     InvMts = 1.d0 / Mts

   else if(cThermostatMethod == 'MNHC') then

     if(QBarostat) then

       if( cBarostatMethod == 'AN' ) then
         NumMNHC = Nf + 1
       else if( cBarostatMethod == 'A2' ) then
         NumMNHC = Nf + 2
       else if( cBarostatMethod == 'A3' ) then
         NumMNHC = Nf + 3
       else if( (cBarostatMethod == 'PR').or.( cBarostatMethod == 'ST' ) ) then
         NumMNHC = Nf + 6
       end if

     else

       NumMNHC = Nf

     end if

     allocate( RMNHC(NHchain,NumMNHC) )
     allocate( VMNHC(NHchain,NumMNHC) )
     allocate( MMNHC(NHchain,NumMNHC) )
     allocate( InvMMNHC(NHchain,NumMNHC) )

     Omega_s0 = 2.d0 * pi / Tau_s0
     Omega_s1 = 2.d0 * pi / Tau_s1

     do j = 1, NumMNHC
       MMNHC(1,j) = kT / ( Omega_s0 * Omega_s0 )
     end do

     do i = 2 , NHchain
     do j = 1, NumMNHC
       MMNHC(i,j) = kT / ( Omega_s1 * Omega_s1 )
     end do
     end do

     InvMMNHC = 1.d0 / MMNHC

!     if(QBarostat) then
!       do i = 1, NHchain
!         do j = Nf + 1, NumMNHC
!           MMNHC(i,j) = kT / (Omega_p * Omega_p )
!         end do
!       end do
!     end if

   end if

   end if

end subroutine Calc_BathMass


! #####################################################################
! #####################################################################


subroutine DPD_Condition

use ScriptData
use CommonBlocks, only : QMaster, QInitial, QHfunc, Qmixord
use CommonDPD, only : Iequil, gamma, Lambda, Iterate, IntegrMethod, &
&   sigma, sigmt, gammt, TypeNum, GammDP, dens, irc, irdf, ShearRate
use IOparam, only : Topology_File
use Numbers, only : N
use BathParam, only : Temp_o
use TimeParam, only : deltat
use CellParam, only : Volume, CellL, H, InvCL
use CellListMethod

implicit none

integer :: i, ii
logical :: gammafile
integer :: NumAllType, NumAllType1
logical :: Flag

! ## DPD -------------------------------------------------------------

! ## EQUIL_STEP ------------------------------------------------------

   Iequil = 1000

   do i = 1, NumOption

     if(OPTIONS(i)=='EQUIL_STEP') then

       ii = Ist(i) + 1

       if(Script(ii)(1:5) == 'STEP=') then

         read(Script(ii)(6:),*) Iequil

       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(.not.QInitial) Iequil = 0

! ## INTEGRATION_DPD -------------------------------------------------

   gamma = 4.5d0
   Lambda = 0.65d0
   Iterate = 10
   gammafile = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='INTEGRATION_DPD') then

       ii = Ist(i)

       Flag = .False.

lia16: do

         ii = ii + 1

         if(Script(ii)(1:2) == '<<') exit lia16

         if(Script(ii)(1:7) == 'METHOD=') then

           read(Script(ii)(8:),*) IntegrMethod
           Flag = .True.

         else if(Script(ii)(1:6) == 'GAMMA=') then

           read(Script(ii)(7:),*) gamma

         else if(Script(ii)(1:11) == 'GAMMA_FILE=') then

           read(Script(ii)(12:),*) Topology_File
           gammafile = .True.

         else if(Script(ii)(1:7) == 'LAMBDA=') then

           read(Script(ii)(8:),*) Lambda

         else if(Script(ii)(1:10) == 'ITERATION=') then

           read(Script(ii)(11:),*) Iterate

         end if

       end do lia16

       if(.not.Flag) then
         if(QMaster) write(*,*) 'ERROR : METHOD must be defined for INTEGRATION_DPD'
         call Finalize
       end if

       if((IntegrMethod/='VVerlet').and.(IntegrMethod/='VVerletSC').and. &
       &  (IntegrMethod/='Lowe').and.(IntegrMethod/='Peters')) then
         if(QMaster) write(*,*) 'ERROR : wrong METHOD for INTEGRATION_DPD'
         call Finalize
       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then

     write(mfile,'(a)') 'DPD:'
     write(mfile,'(a,a    )')   '  Integration = ',trim(IntegrMethod)
     write(mfile,'(a,f10.3)')   '  Gamma       = ',gamma

     if(IntegrMethod == 'VVerlet' .or. IntegrMethod == 'VVerletSC') then
       write(mfile,'(a,f10.3)') '  Lambda      = ',Lambda
     end if

     if(IntegrMethod == 'LeapFrogSC' .or. IntegrMethod == 'VVerletSC') then
       write(mfile,'(a,f10.3)') '  Iteration   = ',Iterate
     end if

   end if

   sigma = sqrt( 2.d0 * gamma * Temp_o )
   sigmt = 2.d0 * sqrt(3.d0) * sigma / sqrt(deltat)

   if(IntegrMethod == 'Lowe') then
     gammt = gamma * deltat
     sigmt = sqrt( 2.d0 * Temp_o )
   end if

   if(IntegrMethod == 'Peters') then

     NumAllType = TypeNum(N)
     allocate( GammDP(NumAllType,NumAllType) )

     if(gammafile) then

       open(3,file=trim(Topology_file),status='old')

       read(3,*) NumAllType1
       read(3,*)

       if(NumAllType1 /= NumAllType) then

         write(*,*) 'ERROR : in the file chosen by GAMMA_FILE'
         call Finalize

       end if

       do i = 1, NumAllType
         read(3,*) GammDP(i,:)
       end do

       close(3)

     else

       GammDP(:,:) = gamma

     end if

     GammDP = GammDP * deltat
     sigmt = 2.d0 * sqrt( 3.d0 * Temp_o )
   end if

! ## NDENSITY --------------------------------------------------------

   dens = 3.d0

   do i = 1, NumOption

     if(OPTIONS(i)=='NDENSITY') then

       ii = Ist(i) + 1

       if(Script(ii)(1:6) == 'NDENS=') then

       read(Script(ii)(7:),*) dens

       else

       write(*,*) 'ERROR : "NDENS=" must be is explicitly declared with the flag "NDENSITY"'
       call Finalize

       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

   Volume = dble(N) / dens

! assuming cubic shape
! cell length : wlx, wly, wlz
   CellL(1) = Volume**(1.d0/3.d0)
   CellL(2) = CellL(1)
   CellL(3) = CellL(1)

   irc = nint( CellL(1) * 0.5 * 20. )
   allocate( irdf(0:irc) )
   irdf = 0

   H = 0.d0
   H(1,1) = CellL(1)
   H(2,2) = CellL(2)
   H(3,3) = CellL(3)

! cell list parameter
   InvCL = 1.d0 / CellL
   Ndiv = int(CellL)
   Ncell= Ndiv(1) * Ndiv(2) * Ndiv(3)
   Maps = 16 * Ncell

   allocate( Head(Ncell) )
   allocate( NextP(N) )
   allocate( Map(Maps) )

   if(QMaster.and.Qcheck) then
     write(mfile,'(a,f10.3)') '  Ndensity    = ', dens
     write(mfile,'(a,f10.3)') '  Cell        = ', CellL(1)
   end if

! ## SHEAR_RATE ------------------------------------------------------

   ShearRate = 0.d0

   do i = 1, NumOption

     if(OPTIONS(i)=='SHEAR_RATE') then

       ii = Ist(i) + 1

       read(Script(ii),*) ShearRate

       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then
     if(ShearRate /= 0.) write(mfile,'(a,f10.3/)') '  Shear Rate  = ', ShearRate
   end if

! ## Hfunc ------------------------------------------------------

   QHfunc = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='Hfunction') then

       QHfunc = .True.
       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then
     if(QHfunc) write(mfile,'(a/)') '  H-function will be calculated  '
   end if

! ## mixorder ------------------------------------------------------

   Qmixord = .False.

   do i = 1, NumOption

     if(OPTIONS(i)=='MIXORDER') then

       Qmixord = .True.
       ReadFlag(i) = .True.

       exit

     end if

   end do

   if(QMaster.and.Qcheck) then
     if(Qmixord) write(mfile,'(a/)') '  The order paramter for mixing will be calculated  '
   end if


end subroutine DPD_Condition


! #####################################################################
! #####################################################################


subroutine ErrorMessage

use ScriptData
use CommonBlocks, only : QMaster

implicit none

integer :: i

   do i = 1, NumOption
     if(.not.ReadFlag(i)) then
       if(QMaster) write(*,*) OPTIONS(i),' was not correctly read!'
       call Finalize
     end if
   end do

   if(QMaster.and.Qcheck) then
     close(mfile)
   end if

   deallocate(OPTIONS,Script,Ist,ReadFlag)

end subroutine ErrorMessage
