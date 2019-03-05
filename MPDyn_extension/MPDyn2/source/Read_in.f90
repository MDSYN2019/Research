! ############################
! ## SUBROUTINE LIST 
! ## -- Set_IO 
! ## -- Read_OptionC 
! ## -- Read_Config 
! ## -- Read_Conf_PI
! ## -- Read_Conf_Bead 
! ## -- NormalModePosition
! ## -- Gene_Velocity 
! ## -- Read_Velocity 
! ## -- Read_Bath 
! ## -- Read_Time 
! ############################


!######################################################################
!######################################################################


! *****************
! **  file name  **
! *****************

subroutine Set_IO

use CommonBlocks, only : Job_name, QMaster, Qstdout, QAveTh, iTrjForm, &
&   ForceField
use IOparam

implicit none

character(len=80) :: monitor_file, average_file

   if( ForceField(1:5) == 'CHARM' ) then

     if(QMaster.and.Qstdout) print *, 'PSF_file=',trim(PSF_file)

   end if

   write(monitor_file,'(a,a)') trim(adjustl(Job_name)),'.moni'
   write(average_file,'(a,a)') trim(adjustl(Job_name)),'.moni.ave'

   NtrjF = 1
   ItrjF = 0
   NvelF = 1
   IvelF = 0

   if(iTrjForm==1) then
     write(trajectory_file,'(a,a,i1,a)') trim(adjustl(Job_name)),'.r00',NtrjF,'.xyz'
   else if(iTrjForm==2) then
     write(trajectory_file,'(a,a,i1)') trim(adjustl(Job_name)),'.r00',NtrjF
   else if(iTrjForm==3) then
     write(trajectory_file,'(a,a,i1,a)') trim(adjustl(Job_name)),'.r00',NtrjF,'.dcd'
   end if
   write(velocity_file,'(a,a,i1)')   trim(adjustl(Job_name)),'.v00',NvelF

! ----------------------
! ## open monitor files
! ----------------------

   if(QMaster) then
#ifdef BMONI
     open(11,file=trim(DirectoryName)//trim(monitor_file),status='unknown',form='unformatted')
     if(QAveTh) open(12,file=trim(DirectoryName)//trim(average_file),status='unknown',form='unformatted')
#else
     open(11,file=trim(DirectoryName)//trim(monitor_file),status='unknown')
     if(QAveTh) open(12,file=trim(DirectoryName)//trim(average_file),status='unknown')
#endif
   end if

end subroutine Set_IO


!######################################################################
!######################################################################


! *********************************************
! ** read condition for optional constraints **
! *********************************************

subroutine Read_OptionC

use UnitExParam, only : ExParam
use OptConstraintParam, only : NumOptC, OptCI, OptCJ, kOptC, rOptC

implicit none

character(len=40) :: Ch
character(len=40), dimension(1000) :: Chs

integer :: i

   NumOptC = 0

   do

     read(3,'(a40)') Ch

     if( Ch(1:5) == '<end>' ) exit

     if( Ch(1:1) == '#' ) cycle

     NumOptC = NumOptC + 1

     if( NumOptC > 1000 ) then

       write(11,*) 'Too many optional constraints!'
       write( *,*) 'Too many optional constraints!'

       call Finalize

     end if

     Chs(NumOptC) = Ch

   end do

   allocate( OptCI(NumOptC) )
   allocate( OptCJ(NumOptC) )
   allocate( kOptC(NumOptC) )
   allocate( rOptC(NumOptC) )

   do i = 1, NumOptC

     read(Chs(i),*) OptCI(i),OptCJ(i),kOptC(i),rOptC(i)

   end do

   kOptC = kOptC * ExParam

end subroutine Read_OptionC


!######################################################################
!######################################################################


! **************************
! **  Read Configration   **
! **************************

subroutine Read_Config(iflag)

use Numbers, only : N, Number
use CommonBlocks, only : QMaster, QInitial, QResForm
use Configuration, only : R
use TimeParam, only : Timeps
use RBparam, only : NumRB, R_RB, V_RB, F_RB, Torque, Quaternion, Lmoment, ScG
use AtomParam, only : AtomName, ResidName, ResidNum, DefaultAtomName, DefaultResidName

implicit none

integer :: i , j, ii, iflag
character(len=72) :: String, String1

if(iflag == 1) then

   if(QInitial) then

     open( 7, file='initial.crd',status='old')

   else

     if(QResForm) then
       open( 7, file='restart.dat',form='formatted',status='old')
     else
       open( 7, file='restart.dat',form='unformatted',status='old')
     end if

   end if

   if(QInitial) then
     do
       read(7,*) String1
       String = trim(adjustl(String1))
       if(String(1:1)=='*'.or.String(1:1)==' '.or.&
       &  String(1:1)=='!'.or.String(1:1)=='#') cycle
       read(String,*) N
       exit
     end do
   else if(QResForm) then
     read(7,'(/)')
     read(7,*) N
   else
     read(7) N
   end if

   if(N/=Number) then
     if(QMaster) write(*,*) 'ERROR: Defined number of particle is wrong'
     if(QMaster) write(*,*) 'N = ',N, ':  Number = ',Number
     call Finalize
   end if

   allocate( AtomName (N) )
   allocate( ResidName(N) )
   allocate( ResidNum (N) )

   if(QInitial) then

     if(N>100000) then
       do i = 1, N
         read(7,'(2i9,2(x,a4),3f12.5)')  &
         & ii , ResidNum(i) , ResidName(i) ,  &
         & AtomName(i) , R(:,i)
       end do
     else
       do i = 1, N
         read(7,'(2i5,2(x,a4),3f10.5)')  &
         & ii , ResidNum(i) , ResidName(i) ,  &
         & AtomName(i) , R(:,i)
       end do
     end if

     Timeps = 0.d0

   else

     if(QResForm) then

       if(N>100000) then
         do i = 1, N
           read(7,'(2i9,2(x,a4),/3d24.16)')  &
           & ii , ResidNum(i) , ResidName(i) ,   &
           & AtomName(i) , R(:,i)
         end do
       else
         do i = 1, N
           read(7,'(2i5,2(x,a4),/3d24.16)')  &
           & ii , ResidNum(i) , ResidName(i) ,   &
           & AtomName(i) , R(:,i)
         end do
       end if

     else

       read(7) ResidNum
       read(7) ResidName
       read(7) AtomName
       read(7) R

     end if

   end if

   allocate( DefaultAtomName (N) )
   allocate( DefaultResidName (N) )

   DefaultAtomName = AtomName
   DefaultResidName = ResidName

else if(iflag == 2) then

   close(7)

   if(QInitial) then
     open( 7, file='initialRB.dat', status='old')
   else
     if(QResForm) then
       open( 7, file='restartRB.dat', form='formatted', status='old')
     else
       open( 7, file='restartRB.dat', form='unformatted', status='old')
     end if
   end if

   allocate( R_RB(3,NumRB) )
   allocate( ScG(3,NumRB) )
   allocate( V_RB(3,NumRB) )
   allocate( F_RB(3,NumRB) )
   allocate( Torque(3,NumRB) )
   allocate( Quaternion(4,NumRB) )
   allocate( Lmoment(3,NumRB) )

   if(QInitial.or.QResForm) then

     do i = 1 , NumRB
       read(7,'(3d24.16)') ( R_RB(j,i) , j = 1 , 3 )
     end do

     do i = 1 , NumRB
       read(7,'(4d24.16)') ( Quaternion(j,i) , j = 1 , 4 )
     end do

   else

     read(7) R_RB
     read(7) Quaternion

   end if

end if

end subroutine Read_Config


!######################################################################
!######################################################################


! **************************
! **  Read Configration   **
! **************************

subroutine Read_Conf_PI

use Numbers, only : N
use CommonBlocks, only : QInitial
use Configuration, only : R
use CommonPI

implicit NONE

integer :: i

   do i = 1, N

     Rnm(:,i,1) = R(:,i)

   end do

   if(QInitial) then

     call NormalModePosition

   else

     call Read_Conf_Bead

   end if

end subroutine Read_Conf_PI


!#####################################################################
!#####################################################################


subroutine Read_Conf_Bead

use Numbers, only : N
use CommonBlocks, only : QResForm
use CommonPI

implicit none

integer :: i, j

   if(QResForm) then

     do j = 2, Nbead
       do i = 1, N
         read(7,'(3d24.16)') Rnm(:,i,j)
       end do
     end do

   else

     do j = 2, Nbead
       do i = 1, N
         read(7) Rnm(:,i,j)
       end do
     end do

   end if

end subroutine Read_Conf_Bead


!#####################################################################
!#####################################################################


! *********************************************
! **  gaussian distribution of normal modes  **
! **  corresponding to free particle         **
! *********************************************

subroutine NormalModePosition

use Numbers, only : N
use CommonPI
use BathParam, only : kT

implicit none

integer :: i, j
real(8) :: s, Gauss
external Gauss

   do j = 2, Nbead

     do i = 1, N

       s = sqrt( kT / OmegaP2 / NmMass(i,j) )
       Rnm(1,i,j) = s * Gauss()
       Rnm(2,i,j) = s * Gauss()
       Rnm(3,i,j) = s * Gauss()

     end do

   end do

end subroutine NormalModePosition


!#####################################################################
!#####################################################################


! *********************************************
! ** Generate Maxwell-Boltzmann Distribution **
! *********************************************

subroutine Gene_Velocity

use CommonBlocks, only : QPBC, QPathInt, QRigidBody, QSHAKE
use Numbers, only : N
use Configuration, only : Vel
use RBparam, only : NumRB, QSingle, Lmoment, QLinear, RBType, InvMassRB, &
&   InertiaRB, NumRBAtom, V_RB
use SHAKEparam
use BathParam, only : Temp_o, kT
use AtomParam, only : Mass, InvMass
use ThermoData, only : Temp, Temp_Translation, Temp_Rotation

implicit NONE

!integer :: i
integer :: i, j, k, kk
real(8) :: ScV, Scale, SumM
real(8), dimension(3) :: SumP
! integer, dimension(3) :: hms
! integer :: Lnum
integer :: MyType
real(8) :: Gauss
external Gauss

! assign random velocity

   if(QPathInt) then
     call PIVelocity
     Return
   end if

   if(QRigidBody) then

     i = 0
     do j = 1 , NumRB

       if(QSingle(j)) then

         i = i + 1
         ScV = sqrt( kT * InvMass(i) )
         Lmoment(:,j) = 0.d0

       else

         MyType = RBType(j)

         ScV = sqrt( kT * InvMassRB(MyType) )

         if(QLinear(j)) then
           LMoment(1,j) = 0.d0
         else
           LMoment(1,j) = sqrt( InertiaRB(1,MyType) * kT ) * Gauss()
         end if

         LMoment(2,j) = sqrt( InertiaRB(2,MyType) * kT ) * Gauss()
         LMoment(3,j) = sqrt( InertiaRB(3,MyType) * kT ) * Gauss()

         i = i + NumRBAtom(MyType)

       end if

       V_RB(1,j) = ScV * Gauss()
       V_RB(2,j) = ScV * Gauss()
       V_RB(3,j) = ScV * Gauss()

     end do

   else

     do i = 1 , N

       ScV = sqrt( kT / Mass(i) )
       Vel(1,i) = ScV * Gauss()
       Vel(2,i) = ScV * Gauss()
       Vel(3,i) = ScV * Gauss()

     end do

     if(QSHAKE) then

       do k = 1 , NSHAKEGroup

         SumP = 0.d0
         SumM = 0.d0

         do kk = 1, NCoupleAtom(k)

           i = CoupleAtom(k,kk)
           SumP = SumP + Mass(i) * Vel(:,i)
           SumM = SumM + Mass(i)

         end do

         SumP = SumP / SumM

         do kk = 1, NCoupleAtom(k)

           i = CoupleAtom(k,kk)
           Vel(:,i) = SumP

         end do

       end do

     end if

   end if

!    ** remove net momentum **
! --------------------
   call Elim_CellMove
! --------------------

   if(.not.QPBC) call Elim_CellRot

! ---------------
   call CalcTemp
! ---------------

   if(QRigidBody) then

     Scale = sqrt( Temp_o / Temp_Translation )

     do i = 1 , NumRB

       V_RB(:,i) = V_RB(:,i) * Scale

     end do

     Scale = sqrt( Temp_o / Temp_Rotation )

     do i = 1 , NumRB

       if(QSingle(i)) cycle

       if(QLinear(i)) then
         Lmoment(2,i) = Lmoment(2,i) * Scale
         Lmoment(3,i) = Lmoment(3,i) * Scale
       else
         Lmoment(:,i) = Lmoment(:,i) * Scale
       end if

     end do

   else

     Scale = sqrt( Temp_o / Temp )

     do i = 1, N

       Vel(:,i) = Vel(:,i) * Scale

     end do

   end if

end subroutine Gene_Velocity


!#####################################################################
!#####################################################################


! *********************
! **  Read Velocity  **
! *********************

subroutine Read_Velocity

use Numbers, only : N
use CommonBlocks, only : QRigidBody, QPathInt, QResForm
use Configuration, only : Vel
use CommonPI
use RBparam, only : NumRB, V_RB, Lmoment, QSingle, QLinear

implicit NONE

integer :: i, j

   if(QPathInt) then

     if(QResForm) then

       do j = 1, Nbead
         do i = 1, N
           read(7,'(3d24.16)') Vnm(:,i,j)
         end do
       end do

     else

       read(7) Vnm

     end if

   else if(QRigidBody) then

     if(QResForm) then

       do i = 1 , NumRB
         read(7,'(3d24.16)') ( V_RB(j,i) , j = 1 , 3 )
       end do

       do i = 1 , NumRB
         read(7,'(3d24.16)') ( Lmoment(j,i) , j = 1 , 3 )
       end do

     else

       read(7) V_RB
       read(7) Lmoment

     end if

     do i = 1 , NumRB

       if(QSingle(i)) Lmoment(:,i) = 0.d0
       if(QLinear(i)) Lmoment(1,i) = 0.d0

     end do

   else

     if(QResForm) then

       do i = 1 , N
         read(7,'(3d24.16)') ( Vel(j,i) , j = 1 , 3 )
       end do

     else

       read(7) Vel

     end if

   end if

end subroutine Read_Velocity


!#####################################################################
!#####################################################################


! **************************
! **  Read Bath Parameter **
! **************************

subroutine Read_Bath

use Numbers, only : N
use CommonBlocks, only : QThermostat, QInitial, cThermostatMethod, QResForm, &
&   QPathInt, QPBC, cBarostatMethod, QMaster, QVolScale, QBarostat
use CommonPI
use BathParam
use CellParam, only : H, InvH, Volume, CellShape, Vsc_Rate

implicit NONE

integer :: i , j, k, l
real(8) :: det, Gauss, sv
external det, Gauss

   if(QThermostat) then

     if(QInitial) then

       if((cThermostatMethod == 'NHC').or. &
       &  (cThermostatMethod == 'NH')) then
         Rss = 0.d0

         do i = 1, NHchain
           sv = sqrt( kT * InvMts(i) )
           Vss(i) = sv * Gauss()
         end do

       else if(cThermostatMethod == 'MNHC') then
         RMNHC = 0.d0

         do i = 1, NumMNHC
           do j = 1, NHchain
             sv = sqrt( kT * InvMMNHC(j,i) )
             VMNHC(j,i) = sv * Gauss()
           end do
         end do

       end if

     else

       if((cThermostatMethod == 'NHC').or. &
       &  (cThermostatMethod == 'NH')) then

         if(QResForm) then

           do i = 1 , NHchain
             read(7,'(2d23.16)') Rss(i), Vss(i)
           end do

         else

           do i = 1 , NHchain
             read(7) Rss(i), Vss(i)
           end do

         end if

       else if(cThermostatMethod == 'MNHC') then

         if(QResForm) then

           do i = 1 , NumMNHC
             do j = 1 , NHchain
               read(7,'(2d23.16)') RMNHC(j,i), VMNHC(j,i)
             end do
           end do

         else

           do i = 1 , NumMNHC
             do j = 1 , NHchain
               read(7) RMNHC(j,i), VMNHC(j,i)
             end do
           end do

         end if

       end if

     end if

   end if

   if(QPathInt) then

     allocate( Fbath(3,N,NHchain,2:Nbead) )
     allocate( Vbath(3,N,NHchain,2:Nbead) )
     allocate( Rbath(3,N,NHchain,2:Nbead) )

     if(QInitial) then

       Rbath = 0.d0

       do k = 2, Nbead
         do j = 1, NHchain
           do i = 1, N
             sv = sqrt( kT * InvQmass(k) )
             Vbath(1,i,j,k) = sv * Gauss()
             Vbath(2,i,j,k) = sv * Gauss()
             Vbath(3,i,j,k) = sv * Gauss()
           end do
         end do
       end do

     else

       if(QResForm) then

         do k = 2, Nbead
           do j = 1, NHchain
             do i = 1, N
               do l = 1, 3
                 read(7,'(2d23.16)') Rbath(l,i,j,k), Vbath(l,i,j,k)
               end do
             end do
           end do
         end do

       else

         do k = 2, Nbead
           do j = 1, NHchain
             do i = 1, N
               do l = 1, 3
                 read(7) Rbath(l,i,j,k), Vbath(l,i,j,k)
               end do
             end do
           end do
         end do

       end if

     end if

   end if

   if(QPBC) then

     if(QInitial) then

       read(7,*) (H(1,i), i = 1, 3)
       read(7,*) (H(2,i), i = 1, 3)
       read(7,*) (H(3,i), i = 1, 3)

     else

       if(QResForm) then

         read(7,'(3d23.16)') (H(1,i), i = 1, 3)
         read(7,'(3d23.16)') (H(2,i), i = 1, 3)
         read(7,'(3d23.16)') (H(3,i), i = 1, 3)
         read(7,'(3d23.16)') (Vg(1,i), i = 1, 3)
         read(7,'(3d23.16)') (Vg(2,i), i = 1, 3)
         read(7,'(3d23.16)') (Vg(3,i), i = 1, 3)

       else

         read(7) H
         read(7) Vg

       end if

! ## file CHECK

!## Andersen isotropic cell fluctuation
       if( cBarostatMethod == 'AN' ) then

         do i = 1 , 2
           do j = i+1 , 3
             if(H(j,i) /= 0.) then
               if(QMaster) then
                 write(11,*) 'ERROR : cell matrix in the configuration file'
                 write( *,*) 'ERROR : cell matrix in the configuration file'
               end if
               call Finalize
             end if
           end do
         end do
         do i = 1, 3
           do j = 1, 3
             if(i==1.and.j==1) cycle
             if(Vg(i,j)/=0.) then
               if(QMaster) then
               write(*,'(a,i1,a,i1,a)') 'WARNING: Vg(',i,',',j,') is not 0, then changed to 0'
               end if
               Vg(i,j) = 0.d0
             end if
           end do
         end do

       else if( cBarostatMethod == 'A2' ) then
         if(Vg(1,1)/=Vg(2,2)) then
           if(QMaster) then
             write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
             write(*,*) 'WARNING for initial velocities of barostat'
             write(*,*) 'Vg(1,1) and Vg(2,2) should have a same value'
             write(*,*) 'whenever A2 barostat is used'
             write(*,*) 'Now automatically assign the zero velocities'
             write(*,*) 'i.e. Vg(1,1) = Vg(2,2) = Vg(3,3) = 0.0'
             write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           end if
           Vg = 0.d0
         end if
       end if

       if(.not.QBarostat) Vg = 0.d0

     end if

     Volume = det(H)
     call InversMatrix(H,InvH)

! ForCG
     CellShape = 3
     if((H(1,2) == 0.).and.(H(1,3) == 0.).and.&
     &  (H(2,1) == 0.).and.(H(2,3) == 0.).and.&
     &  (H(3,1) == 0.).and.(H(3,2) == 0.)) then
       CellShape = 2
       if((H(1,1) == H(2,2)).and.(H(1,1) == H(3,3))) then
         CellShape = 1
       end if
     end if

     if(QVolScale) then

       do i = 1, 3
         do j = 1, 3
           if(i==j) cycle
           if(H(i,j)/=0.) then
             if(QMaster) write(*,*) 'VOL_SCALE option can be used only for ortholombic cell'
             call Finalize
           end if
         end do
       end do

       if(QBarostat) then
         if(QMaster) then
         write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write(*,*) 'WARNING: VOL_SCALE option cannot be used in combination with BAROSTAT'
         write(*,*) '         BAROSTAT is now automatically switched off                  '
         write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         end if
         QBarostat = .False.
       end if

       if(CellShape==1) then
         if(Vsc_Rate(1)/=Vsc_Rate(2).or.Vsc_Rate(1)/=Vsc_Rate(3)) then
           CellShape = 2
         end if
       end if

       if(QThermostat) then
         if(cThermostatMethod/='VSCALE') then
           cThermostatMethod = 'VSCALE'
           if(QMaster) then
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           write(*,*) 'WARNING: Whenever VOL_SCALE is selected, VSCALE (velocity scaling)'
           write(*,*) 'should be used for thermostat. Now VSCALE method is automatically '
           write(*,*) 'selected for this computation.                                    '
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           end if
         else
           if(QMaster) then
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           write(*,*) 'WARNING: Whenever VOL_SCALE is selected, VSCALE (velocity scaling)'
           write(*,*) 'thermostat should be used. Now VSCALE thermostat is automatically '
           write(*,*) 'switched on for the current computation.                          '
           write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
           end if
           QThermostat = .True.
           cThermostatMethod = 'VSCALE'
         end if
       end if

     end if

   end if

end subroutine Read_Bath


!#####################################################################
!#####################################################################


! **************************
! **  Read Time **
! **************************

subroutine Read_Time

use CommonBlocks, only : SimMethod, QResForm
use TimeParam, only : Timeps
use CommonHMC, only : TimeMC

implicit NONE

   if(SimMethod == 'HMC') then
     if(QResForm) then
       read(7,*) TimeMC
     else
       read(7) TimeMC
     end if
   else
     if(QResForm) then
       read(7,*) Timeps
     else
       read(7) Timeps
     end if
   end if

   close(7)

end subroutine Read_Time
