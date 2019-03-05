! ############################
! ## SUBROUTINE LIST 
! ## -- Select_Ensemble 
! ## -- Select_Ensemble_PI 
! ## -- Select_Ensemble_HMC 
! ## -- Select_Ensemble_DynaLib 
! ## -- Select_DPDalgorithm 
! ############################


!######################################################################
!######################################################################


subroutine Select_Ensemble

use CommonBlocks, only : QPBC, QBarostat, QRigidBody, QSHAKE, cBarostatMethod
use CellParam, only : CellShape

implicit NONE

   if(QPBC) then    ! Periodic Boundary Condition (Bulk System)

     if(QBarostat) then

! NEWDAT
       if((cBarostatMethod=='PR').or.(cBarostatMethod=='ST')) CellShape = 3
       if(((cBarostatMethod=='A3').or.(cBarostatMethod=='A2')) &
       &   .and.(CellShape/=3)) CellShape = 2
       if((cBarostatMethod=='AN').and.(CellShape==1)) CellShape = 1

       if(QRigidBody) then
!       -------------------------------------
         call IntegrEOM_NP_RIGID
!       --------------------------------------
       else if(QSHAKE) then
!       -------------------------------------
         call IntegrEOM_NP_SHAKE
!       --------------------------------------
       else
!       --------------------
         call IntegrEOM_NP
!       --------------------
       end if

     else

       if(QRigidBody) then
!       --------------------------
         call IntegrEOM_NV_RIGID
!       --------------------------
       else
!       --------------------
         call IntegrEOM_NV
!       --------------------
       end if

     end if

   else                   ! Isolated System

     if(QRigidBody) then
!     -----------------------------
       call IntegrEOM_iso_RIGID
!     -----------------------------
     else
!     -----------------------
       call IntegrEOM_iso
!     -----------------------
     end if

   end if

end subroutine Select_Ensemble


!######################################################################
!######################################################################


subroutine Select_Ensemble_PI

use CommonBlocks, only : QPBC, QBarostat, cBarostatMethod
use CellParam, only : CellShape

implicit NONE

   if(QPBC) then    ! Periodic Boundary Condition (Bulk System)

     if(QBarostat) then

! NEWDAT
       if((cBarostatMethod=='PR').or.(cBarostatMethod=='ST')) CellShape = 3
       if(((cBarostatMethod=='A3').or.(cBarostatMethod=='A2')) &
       &   .and.(CellShape/=3)) CellShape = 2
       if((cBarostatMethod=='AN').and.(CellShape==1)) CellShape = 1

!    -------------------------
       call IntegrEOM_NPT_PI
!    -------------------------

     else

!    -------------------------
       call IntegrEOM_NVT_PI
!    -------------------------

     end if

   else                   ! Isolated System

!   -----------------------
     call IntegrEOM_iso_PI
!   -----------------------

   end if

end subroutine Select_Ensemble_PI


!######################################################################
!######################################################################


subroutine Select_Ensemble_HMC

use CommonBlocks, only : QPBC, QBarostat, cBarostatMethod, QRigidBody
use CellParam, only : CellShape

implicit NONE

   if(QPBC) then    ! Periodic Boundary Condition (Bulk System)

     if(QBarostat) then

! NEWDAT
       if((cBarostatMethod=='PR').or.(cBarostatMethod=='ST')) CellShape = 3
       if(((cBarostatMethod=='A3').or.(cBarostatMethod=='A2')) &
       &   .and.(CellShape/=3)) CellShape = 2
       if((cBarostatMethod=='AN').and.(CellShape==1)) CellShape = 1

       if(QRigidBody) then
!       -------------------------------------
         call HMC_NP_RIGID
!       --------------------------------------
       else
!       --------------------
         call HMC_NP
!       --------------------
       end if

     else

       if(QRigidBody) then
!       --------------------------
         call HMC_NV_RIGID
!       --------------------------
       else
!       --------------------
         call HMC_NV
!       --------------------
       end if

     end if

   else                   ! Isolated System

     if(QRigidBody) then
!     -----------------------------
       call HMC_iso_RIGID
!     -----------------------------
     else
!     -----------------------
       call HMC_iso
!     -----------------------
     end if

   end if

end subroutine Select_Ensemble_HMC


!######################################################################
!######################################################################


subroutine Select_Ensemble_DynaLib

use CommonBlocks, only : QPBC, QBarostat, cBarostatMethod
use CellParam, only : CellShape

implicit NONE

   if(QPBC) then    ! Periodic Boundary Condition (Bulk System)

     if(QBarostat) then

! NEWDAT
       if((cBarostatMethod=='PR').or.(cBarostatMethod=='ST')) CellShape = 3
       if(((cBarostatMethod=='A3').or.(cBarostatMethod=='A2')) &
       &   .and.(CellShape/=3)) CellShape = 2
       if((cBarostatMethod=='AN').and.(CellShape==1)) CellShape = 1

       call DynaLib_NP

     else

       call DynaLib_NV

     end if

   else                   ! Isolated System

     call DynaLib_iso

   end if

end subroutine Select_Ensemble_DynaLib


!######################################################################
!######################################################################


subroutine Select_DPDalgorithm

use Numbers, only : N
use IOparam, only : DirectoryName
use CommonBlocks, only : Job_name, QHfunc, Qmixord
use CommonDPD
use BookParam, only : ListIJ

implicit NONE

character(len=72) :: Temp_file, Press_file, Energy_file
character(len=72) :: MSD_file, Hfunc_file, Mix_file

character(len=72) :: Vir_file

   write(Temp_file  ,'(a,a)') trim(Job_Name) , '.Temp'
   write(Press_file ,'(a,a)') trim(Job_Name) , '.Pres'
   write(Energy_file,'(a,a)') trim(Job_Name) , '.Ener'
   write(MSD_file,   '(a,a)') trim(Job_Name) , '.MSD'
   write(Vir_file,   '(a,a)') trim(Job_Name) , '.virial'

   open(22 , file = trim(DirectoryName)//trim(Energy_file) , status = 'unknown' )
   open(23 , file = trim(DirectoryName)//trim(Press_file ) , status = 'unknown' )
   open(24 , file = trim(DirectoryName)//trim(MSD_file   ) , status = 'unknown' )
   open(25 , file = trim(DirectoryName)//trim(Vir_file   ) , status = 'unknown' )
   open(28 , file = trim(DirectoryName)//trim(Temp_file  ) , status = 'unknown' )

   if(QHfunc) then
     write(Hfunc_file,   '(a,a)') trim(Job_Name) , '.Hfunc'
     open(29 , file = trim(DirectoryName)//trim(Hfunc_file) , status = 'unknown' )
   end if

   if(Qmixord) then
     write(Mix_file, '(a,a)') trim(Job_Name), '.morder'
     open(30 , file = trim(DirectoryName)//trim(Mix_file) , status = 'unknown' )
   end if

   if(IntegrMethod=='VVerlet') then

!   ---------------------------
     call IntegrDPD_VVerlet
!   ---------------------------

   else if(IntegrMethod=='VVerletSC') then

     Ndm = 10 * dens * N

     allocate( ListIJ(2,Ndm) )
     allocate( dRList(3,Ndm) )
     allocate( R1List(Ndm) )
     allocate( pfList(Ndm) )
     allocate( SLList(Ndm) )

! ## for safety
     ListIJ = 0
     dRList = 0.d0
     R1List = 0.d0
     pfList = 0.d0
     SLList = 0.d0

     allocate( FrcDPd(3,N) )

     if(QColloid) then
       allocate( FrcCod(3,NumColloid) )
       allocate( Torqued(3,NumColloid) )
     end if

!   ---------------------------
     call IntegrDPD_VVerletSC
!   ---------------------------

   else if(IntegrMethod=='LeapFrogSC') then

!   ---------------------------
     call IntegrDPD_LeapFrogSC
!   ---------------------------

   else if((IntegrMethod=='Lowe').or.(IntegrMethod=='Peters')) then

     Ndm = 10 * dens * N

     allocate( ListIJ(2,Ndm) )
     allocate( R1List(Ndm) )
     allocate( dRList(3,Ndm) )
     allocate( SLList(Ndm) )

! ## for safety
     ListIJ = 0
     R1List = 0.d0
     dRList = 0.d0
     SLList = 0.d0

!   ---------------------------
     call IntegrDPD_Lowe
!   ---------------------------

   end if

   call Print_Average_DP

   close(22)
   close(23)
   close(24)
   close(25)
   close(28)
   if(QHfunc) close(29)
   if(Qmixord) close(30)

end subroutine Select_DPDalgorithm
