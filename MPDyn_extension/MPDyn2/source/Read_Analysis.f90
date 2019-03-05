! ############################
! ## SUBROUTINE LIST 
! ## -- Read_Analyz_filename 
! ## -- Read_Analyz_Cond 
! ############################


! #####################################################################
! #####################################################################


subroutine Read_Analyz_filename

use ScriptData
use CommonBlocks, only : QMaster, Qstdout
use ParamAnalyze, only : NJobs, NTrjStep, Interval, FileHead
#ifdef MOLFILE
use MOLFILEparam
#endif

implicit none

integer :: i, ii, j

   do i = 1, NumOption

     if(OPTIONS(i)=='ANALYZ_FILE') then

       ii = Ist(i) + 1

       read(Script(ii),*) NJobs

       allocate( NTrjStep(NJobs) )
       allocate( Interval(NJobs) )
       allocate( FileHead(NJobs) )
#ifdef MOLFILE
       allocate( Handle(Njobs) )
       Handle(:) = -1
#endif

       do j = 1, NJobs

         FileHead(j) = trim(adjustl(Script(ii+1)))
         read(Script(ii+2),*) NTrjStep(j), Interval(j)
         ii = ii + 2

         if(QMaster.and.Qstdout) print *, FileHead(j)

       end do

       ReadFlag(i) = .True.

       exit

     end if

   end do

end subroutine Read_Analyz_filename


! #####################################################################
! #####################################################################


subroutine Read_Analyz_Cond

use ScriptData
use CommonBlocks, only : QMaster, Qstdout, QSwitch
use ParamAnalyze
use AtomParam, only : MolName
use Numbers, only : NumSpec, NumMol
use BathParam, only : Temp_o, Pressure_o, kT, Beta
use CutoffParam, only : Rcutoff2, Ron2, swf1
use UnitExParam, only : rprs, kb

implicit none

integer :: i, ii, j, k, jj
integer :: NAnaOpt
logical, dimension(:), allocatable :: AnaFlag
integer :: CIdirection
integer :: ibond, iangle, idihed
external CIdirection
character(len=3) :: cONOFF
real(8) :: Rcutoff, Ron
character(len=1) :: Ch, Ch1, Ch2
integer :: FROM, TO, Shift, Ntot
character(len=10) :: CName
character(len=80) :: Script1

   do i = 1, NumOption

     if(OPTIONS(i)=='ANALYZ_NAME') then

       ii = Ist(i) + 1

       if(Script(ii)(1:7) == 'METHOD=') then

         write(cWhat,'(a)') trim(adjustl(Script(ii)(8:)))

         if(QMaster.and.Qcheck) then
           write(mfile,*) 'Analysis : ', cWhat
           write(*,*)  'Analysis : ', cWhat
         end if

       else

         if(QMaster) then
           write(*,*) 'SCRIPT ERROR : METHOD of ANALYSIS (input.data)'
           write(*,*) 'You have to give METHOD= in the first line of '
           write(*,*) 'the keyword option "ANALYZ_NAME"'
         end if
         call Finalize

       end if

! ## for PDB -----------------------------

       if(cWhat == 'PDB') then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia17:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia17

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster) write(*,*) 'number of snapshot to be made : ', Nsnap

           end if

         end do lia17

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "PDB"'
           call Finalize
         end if

         allocate( ResultFile(Nsnap) )

         Shift = 0  ! default

         do j = 1 , Nsnap

           if((j+Shift)<10) then
             write(ResultFile(j),'(a,i1,a)') './Analy/snap000',(j+Shift),'.pdb'
           else if((j+Shift)<100) then
             write(ResultFile(j),'(a,i2,a)') './Analy/snap00', (j+Shift),'.pdb'
           else if((j+Shift)<1000) then
             write(ResultFile(j),'(a,i3,a)') './Analy/snap0',  (j+Shift),'.pdb'
           else if((j+Shift)<10000) then
             write(ResultFile(j),'(a,i4,a)') './Analy/snap',   (j+Shift),'.pdb'
           else
             if(QMaster) write(*,*) 'ERROR :: too many file'
             call Finalize
           end if

         end do

         Ntot = 0

         do j = 1 , NJobs
             Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

         if(QMaster.and.Qcheck) then
           write(mfile,'(a,i10)') '  STEP = ',Nsnap
         end if

! ## for PDBsp -----------------------------

       else if(cWhat == 'PDBsp') then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         Qregion = .False.

lia18:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia18

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qstdout) write(*,*) 'number of snapshot to be made : ', Nsnap

           else if(Script(ii)(1:8) == 'NUMCOMP=') then

             read(Script(ii)(9:),*) Ncomp
             allocate(NumComp(Ncomp))
             if(QMaster.and.Qstdout) write(*,*) 'number of components to be printed : ', Ncomp
             AnaFlag(2) = .True.

           end if

         end do lia18

         ii = Ist(i) + 1

         if(QMaster.and.Qcheck) then
           write(mfile,'(a,i10)') '  STEP    = ',Nsnap
           write(mfile,'(a,i10)') '  NUMCOMP = ',Ncomp
         end if

lia19:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia19

           if(Script(ii)(1:9) == 'SELECTED=') then

             read(Script(ii)(10:),*) (NumComp(j),j=1,Ncomp)
             AnaFlag(3) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The number of selected components : ', NumComp(:)

           else if(Script(ii)(1:7) == 'REGION=') then

             read(Script(ii)(8:),*) Ch

             if(Ch == 'Y') then

               if(QMaster.and.Qcheck) write(mfile,*) '>> the region was specified'
               Qregion = .True.

               Xmin = -10.d0
               Ymin = -10.d0
               Zmin = -10.d0
               Xmax =  10.d0
               Ymax =  10.d0
               Zmax =  10.d0

             end if

           else if(Script(ii)(1:2) == 'X=') then

             read(Script(ii)(3:),*) Xmin, Xmax

           else if(Script(ii)(1:2) == 'Y=') then

             read(Script(ii)(3:),*) Ymin, Ymax

           else if(Script(ii)(1:2) == 'Z=') then

             read(Script(ii)(3:),*) Zmin, Zmax

           end if

         end do lia19

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "PDBsp"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "NUMCOMP" for "PDBsp"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "SELECTED" for "PDBsp"'
           call Finalize
         end if

         if(Qregion) then

           if(QMaster.and.Qcheck) then
             write(mfile,'(2(a,f7.2))') '[X] from ',Xmin,' to ',Xmax
             write(mfile,'(2(a,f7.2))') '[Y] from ',Ymin,' to ',Ymax
             write(mfile,'(2(a,f7.2))') '[Z] from ',Zmin,' to ',Zmax
           end if

         end if

         allocate( ResultFile(Nsnap) )

         Shift = 0  ! default

         do j = 1 , Nsnap

           if((j+Shift)<10) then

             write(ResultFile(j),'(a,i1,a)') './Analy/snap000',(j+Shift),'.pdb'

           else if((j+Shift)<100) then

             write(ResultFile(j),'(a,i2,a)') './Analy/snap00', (j+Shift),'.pdb'

           else if((j+Shift)<1000) then

             write(ResultFile(j),'(a,i3,a)') './Analy/snap0',  (j+Shift),'.pdb'

           else if((j+Shift)<10000) then

             write(ResultFile(j),'(a,i4,a)') './Analy/snap',   (j+Shift),'.pdb'

           else

             if(QMaster) write(*,*) 'ERROR :: too many file'
             call Finalize

           end if

         end do

         Ntot = 0

         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if


! ## for PDBsp -----------------------------

       else if(cWhat == 'PDBex') then

         NAnaOpt = 4
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         Qregion = .False.

lib18:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lib18

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster) write(*,*) 'number of snapshot to be made : ', Nsnap

           else if(Script(ii)(1:2) == 'X=') then

             read(Script(ii)(3:),*) Xmin, Xmax
             AnaFlag(2) = .True.

           else if(Script(ii)(1:2) == 'Y=') then

             read(Script(ii)(3:),*) Ymin, Ymax
             AnaFlag(3) = .True.

           else if(Script(ii)(1:2) == 'Z=') then

             read(Script(ii)(3:),*) Zmin, Zmax
             AnaFlag(4) = .True.

           end if

         end do lib18

         ii = Ist(i) + 1

         if(QMaster.and.Qcheck) then
           write(mfile,'(a,i10)') '  STEP = ',Nsnap
           write(mfile,'(2(a,f7.2))') '[X] from ',Xmin,' to ',Xmax
           write(mfile,'(2(a,f7.2))') '[Y] from ',Ymin,' to ',Ymax
           write(mfile,'(2(a,f7.2))') '[Z] from ',Zmin,' to ',Zmax
         end if

         allocate( ResultFile(Nsnap) )

         Shift = 0  ! default

         do j = 1 , Nsnap

           if((j+Shift)<10) then

             write(ResultFile(j),'(a,i1,a)') './Analy/snap000',(j+Shift),'.pdb'

           else if((j+Shift)<100) then

             write(ResultFile(j),'(a,i2,a)') './Analy/snap00', (j+Shift),'.pdb'

           else if((j+Shift)<1000) then

             write(ResultFile(j),'(a,i3,a)') './Analy/snap0',  (j+Shift),'.pdb'

           else if((j+Shift)<10000) then

             write(ResultFile(j),'(a,i4,a)') './Analy/snap',   (j+Shift),'.pdb'

           else

             if(QMaster) write(*,*) 'ERROR :: too many file'
             call Finalize

           end if

         end do

         Ntot = 0

         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if


! ## for PDBminR -----------------------------

       else if(cWhat == 'PDBminR') then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia20:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia20

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             &  write(mfile,*) '  number of snapshot to be made : ', Nsnap

           else if(Script(ii)(1:8) == 'PROTEIN=') then

             read(Script(ii)(6:),*) Kcomp
             AnaFlag(2) = .True.
             if(QMaster.and.Qcheck) write(mfile,*) &
             & '  The component number of protein is tagged as ', Kcomp

           end if

         end do lia20

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "PDBsp"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "PROTEIN" for "PDBsp"'
           call Finalize
         end if

         allocate( ResultFile(Nsnap) )

         Shift = 0  ! default

         do j = 1 , Nsnap

           if((j+Shift)<10) then
             write(ResultFile(j),'(a,i1,a)') './Analy/snap000',(j+Shift),'.pdb'
           else if((j+Shift)<100) then
             write(ResultFile(j),'(a,i2,a)') './Analy/snap00', (j+Shift),'.pdb'
           else if((j+Shift)<1000) then
             write(ResultFile(j),'(a,i3,a)') './Analy/snap0',  (j+Shift),'.pdb'
           else if((j+Shift)<10000) then
             write(ResultFile(j),'(a,i4,a)') './Analy/snap',   (j+Shift),'.pdb'
           else
             if(QMaster) write(*,*) 'ERROR :: too many file'
             call Finalize
           end if

         end do

         Ntot = 0
         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

! ## for XYZ -----------------------------

       else if( cWhat == 'XYZsp' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         Qregion = .False.

lia21:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia21

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  number of snapshot to be made : ', Nsnap

           else if(Script(ii)(1:8) == 'NUMCOMP=') then

             read(Script(ii)(9:),*) Ncomp
             allocate(NumComp(Ncomp))
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  number of components to be printed : ', Ncomp
             AnaFlag(2) = .True.
             if(Ncomp > NumSpec) then
               if(QMaster) write(*,*) 'INPUT ERROR : the number must be less than NumSpec'
               call Finalize
             end if

           end if

         end do lia21

         ii = Ist(i) + 1

lia22:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia22

           if(Script(ii)(1:9) == 'SELECTED=') then

             read(Script(ii)(10:),*) (NumComp(j),j=1,Ncomp)
             AnaFlag(3) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The number of selected components : ', NumComp(:)

           else if(Script(ii)(1:7) == 'REGION=') then

             read(Script(ii)(8:),*) Ch

             if(Ch == 'Y') then

               if(QMaster.and.Qcheck) write(mfile,*) '>> the region was specified '
               Qregion = .True.

               Xmin = -10.d0
               Ymin = -10.d0
               Zmin = -10.d0
               Xmax =  10.d0
               Ymax =  10.d0
               Zmax =  10.d0

             end if

           else if(Script(ii)(1:2) == 'X=') then

             read(Script(ii)(3:),*) Xmin, Xmax

           else if(Script(ii)(1:2) == 'Y=') then

             read(Script(ii)(3:),*) Ymin, Ymax

           else if(Script(ii)(1:2) == 'Z=') then

             read(Script(ii)(3:),*) Zmin, Zmax

           end if

         end do lia22

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "PDBsp"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "NUMCOMP" for "PDBsp"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "SELECTED" for "PDBsp"'
           call Finalize
         end if

         if(Qregion) then

           if(QMaster.and.Qcheck) then
             write(mfile,'(2(a,f7.2))') '[X] from ',Xmin,' to ',Xmax
             write(mfile,'(2(a,f7.2))') '[Y] from ',Ymin,' to ',Ymax
             write(mfile,'(2(a,f7.2))') '[Z] from ',Zmin,' to ',Zmax
           end if

         end if

         Ntot = 0
         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

! ## for DCD -----------------------------

       else if( cWhat == 'DCD' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         Qregion = .False.

lya21:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lya21

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  number of snapshot to be made : ', Nsnap

           else if(Script(ii)(1:8) == 'NUMCOMP=') then

             read(Script(ii)(9:),*) Ncomp
             allocate(NumComp(Ncomp))
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  number of components to be printed : ', Ncomp
             AnaFlag(2) = .True.
             if(Ncomp > NumSpec) then
               if(QMaster) write(*,*) 'INPUT ERROR : the number must be less than NumSpec'
               call Finalize
             end if

           end if

         end do lya21

         ii = Ist(i) + 1

lya22:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lya22

           if(Script(ii)(1:9) == 'SELECTED=') then

             read(Script(ii)(10:),*) (NumComp(j),j=1,Ncomp)
             AnaFlag(3) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The number of selected components : ', NumComp(:)

           end if

         end do lya22

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "DCD"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "NUMCOMP" for "DCD"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "SELECTED" for "DCD"'
           call Finalize
         end if

         Ntot = 0
         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

! ## for ARC -----------------------------

       else if(cWhat == 'ARC') then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia23:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia23

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The number of snapshot to be made : ', Nsnap

           end if

         end do lia23

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "ARC"'
           call Finalize
         end if

         Ntot = 0
         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

! ## for ARC -----------------------------

       else if(cWhat == 'ARC_P') then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia24:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia24

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The number of snapshot to be made : ', Nsnap

           else if(Script(ii)(1:8) == 'PROTEIN=') then

             read(Script(ii)(6:),*) Kcomp
             AnaFlag(2) = .True.
             if(QMaster.and.Qcheck) write(mfile,*) &
             & '  The component number of protein is tagged as ', Kcomp

           end if

         end do lia24

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "ARC_P"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "PROTEIN" for "ARC_P"'
           call Finalize
         end if

         Ntot = 0
         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

! ## for CRD -----------------------------

       else if( cWhat == 'CRD' ) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia25:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia25

           if(Script(ii)(1:5) == 'STEP=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The number of snapshot to be made : ', Nsnap

           end if

         end do lia25

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "STEP" for "CRD"'
           call Finalize
         end if

         allocate( ResultFile(Nsnap) )

         do j = 1 , Nsnap

           if(j<10) then
             write(ResultFile(j),'(a,i1,a)') './Analy/snap000',j,'.crd'
           else if(j<100) then
             write(ResultFile(j),'(a,i2,a)') './Analy/snap00', j,'.crd'
           else if(j<1000) then
             write(ResultFile(j),'(a,i3,a)') './Analy/snap0',  j,'.crd'
           else if(j<10000) then
             write(ResultFile(j),'(a,i4,a)') './Analy/snap',   j,'.crd'
           else
             if(QMaster) write(*,*) 'ERROR :: too many file'
             call Finalize
           end if

         end do

         Ntot = 0
         do j = 1 , NJobs
           Ntot = Ntot + ( NTrjStep(j) / Interval(j) )
         end do

         if(Ntot/=Nsnap) then
           if(QMaster) write(*,*) 'ERROR :: File number'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if(cWhat=='CBPI') then

         NAnaOpt = 4
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

         write(CBPI_FILE,'(a)') './param/PImolecule.data'
         NumInsertion = 5
         NumOriTry = 10

lia26:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia26

           if(Script(ii)(1:2) == 'T=') then

             read(Script(ii)(3:),*) Temp_o
             AnaFlag(1) = .True.

           else if(Script(ii)(1:2) == 'P=') then

             read(Script(ii)(3:),*) Pressure_o
             AnaFlag(2) = .True.

           else if(Script(ii)(1:5) == 'RCAV=') then

             read(Script(ii)(6:),*) R_cav
             AnaFlag(3) = .True.

           else if(Script(ii)(1:7) == 'REGION=') then

             read(Script(ii)(8:),*) Rxmax, Rymax, Rzmax
             AnaFlag(4) = .True.

           else if(Script(ii)(1:5) == 'FILE=') then

             read(Script(ii)(6:),*) CBPI_FILE

           else if(Script(ii)(1:6) == 'TRIAL=') then

             read(Script(ii)(7:),*) NumInsertion

           else if(Script(ii)(1:7) == 'TRIORI=') then

             read(Script(ii)(8:),*) NumOriTry

           end if

         end do lia26

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "T" for "CBPI"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "P" for "CBPI"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "RCAV" for "CBPI"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "REGION" for "CBPI"'
           call Finalize
         end if

         if(QMaster.and.Qcheck) then

           write(mfile,'(a)') ' #################################### '
           write(mfile,'(a)') '      cavity insertion parameters     '
           write(mfile,'(a)') ' #################################### '
           write(mfile,'(a,a)') '   FILE = ',CBPI_FILE
           write(mfile,'(a,f7.1)') '   Temperature  = ',Temp_o
           write(mfile,'(a,f7.1)') '   Pressure     = ',Pressure_o
           write(mfile,'(a,i5)')   '   NumInsertion = ',NumInsertion
           write(mfile,'(a,i5)')   '   NumOriTry    = ',NumOriTry
           write(mfile,'(a,f5.1)') '   R_cav        = ',R_cav
           write(mfile,'(a,f5.1)') '   Rxmax        = ',Rxmax
           write(mfile,'(a,f5.1)') '   Rymax        = ',Rymax
           write(mfile,'(a,f5.1)') '   Rzmax        = ',Rzmax
           write(mfile,'(a/)') ' ------------------------------------ '

         end if

         Pressure_o = 1.d+6 * Pressure_o * rprs
         kT   = kb * Temp_o
         Beta = 1.d0 / kT

! ## ---------------------------------------------------------

       else if(cWhat=='Cavity') then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia27:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia27

           if(Script(ii)(1:6) == 'LIPID=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.
            if(QMaster.and.Qcheck) &
            & write(mfile,*) '  The component number of lipid is tagged as ', Kcomp

           else if(Script(ii)(1:5) == 'RCAV=') then

             read(Script(ii)(6:),*) R_cav
             AnaFlag(2) = .True.
            if(QMaster.and.Qcheck) &
            & write(mfile,'(a,f9.3)') '  Diameter of the cavity ', R_cav

           else if(Script(ii)(1:7) == 'REGION=') then

             read(Script(ii)(8:),*) Rxmax, Rymax, Rzmax
             AnaFlag(3) = .True.
            if(QMaster.and.Qcheck) &
            & write(mfile,'(a,f9.2)') '  Region = ', Rxmax, Rymax, Rzmax

           end if

         end do lia27

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "LIPID" for "Cavity"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RCAV" for "Cavity"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "REGION" for "Cavity"'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if(( cWhat == 'OriChain'   ).or. &
       &       ( cWhat == 'LenChain'   ).or. &
       &       ( cWhat == 'CorrRotW'   )) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia28:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia28

           if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(1) = .True.
            if(QMaster.and.Qcheck) write(mfile,*) '  Time interval = ', dtime

           end if

         end do lia28

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "DTIME"'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if( cWhat == 'CoTGLipidC' ) then

         NAnaOpt = 5
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lic28:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lic28

           if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(1) = .True.

           else if(Script(ii)(1:6) == 'TOTAL=') then

             read(Script(ii)(7:),*) Nsnap
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'BLOCK=') then

             read(Script(ii)(7:),*) BlockStep
             AnaFlag(3) = .True.

           else if(Script(ii)(1:7) == 'LENGTH=') then

             read(Script(ii)(8:),*) ilength
             AnaFlag(4) = .True.

           else if(Script(ii)(1:7) == 'INTERV=') then

             read(Script(ii)(8:),*) interv
             AnaFlag(5) = .True.

           end if

         end do lic28

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "DTIME" for "CoTGLipidC"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "TOTAL" for "CoTGLipidC"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "BLOCK" for "CoTGLipidC"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "LENGTH" for "CoTGLipidC"'
           call Finalize
         end if
         if(.not.AnaFlag(5)) then
           if(QMaster) write(*,*) 'Missing "INTERV" for "CoTGLipidC"'
           call Finalize
         end if

         if(QMaster.and.Qcheck) then
           write(mfile,*) ' total step number         = ', Nsnap
           write(mfile,*) ' block step number         = ', BlockStep
           write(mfile,*) ' ########<note>#########################################'
           write(mfile,*) ' #  ******************************                     #'
           write(mfile,*) ' #  <---------ilength------------>                     #'
           write(mfile,*) ' #  <--interval-> ***************************          #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #                          ************************** #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #######################################################'
           write(mfile,*) ' ilength  = ', ilength
           write(mfile,*) ' interval = ', interv
           write(mfile,*) ' dtime    = ', dtime
         end if

! ## ---------------------------------------------------------

       else if( cWhat == 'TGchain' ) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lic128:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lic128

           if(Script(ii)(1:6) == 'TOTAL=') then

             read(Script(ii)(7:),*) Nsnap
             AnaFlag(1) = .True.

           end if

         end do lic128

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "TOTAL" for "TGchain"'
           call Finalize
         end if

! ## diffusion coefficient for water from MSD (Z-dependence)------------

       else if( cWhat == 'DiffW_Z' ) then

!#################<note>#########################################
!            ******************************                     #
!            <---------ilength------------>                     #
!            <--interval-> ***************************          #
!                                                               #
!                                    ************************** #
!                                                               #
!################################################################

         NAnaOpt = 9
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia29:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia29

           if(Script(ii)(1:6) == 'WATER=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.

           else if(Script(ii)(1:6) == 'TOTAL=') then

             read(Script(ii)(7:),*) Nsnap
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'BLOCK=') then

             read(Script(ii)(7:),*) BlockStep
             AnaFlag(3) = .True.

           else if(Script(ii)(1:7) == 'LENGTH=') then

             read(Script(ii)(8:),*) ilength
             AnaFlag(4) = .True.

           else if(Script(ii)(1:7) == 'INTERV=') then

             read(Script(ii)(8:),*) interv
             AnaFlag(5) = .True.

           else if(Script(ii)(1:7) == 'DELTAZ=') then

             read(Script(ii)(8:),*) deltaZ
             AnaFlag(6) = .True.

           else if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(7) = .True.

           else if(Script(ii)(1:5) == 'MAXZ=') then

             read(Script(ii)(6:),*) maxdz
             AnaFlag(8) = .True.

           else if(Script(ii)(1:5) == 'MINZ=') then

             read(Script(ii)(6:),*) mindz
             AnaFlag(9) = .True.

           end if

         end do lia29

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "WATER" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "TOTAL" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "BLOCK" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "LENGTH" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(5)) then
           if(QMaster) write(*,*) 'Missing "INTERV" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(6)) then
           if(QMaster) write(*,*) 'Missing "DELTAZ" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(7)) then
           if(QMaster) write(*,*) 'Missing "DTIME" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(8)) then
           if(QMaster) write(*,*) 'Missing "MAXZ" for "DiffW_Z"'
           call Finalize
         end if
         if(.not.AnaFlag(9)) then
           if(QMaster) write(*,*) 'Missing "MINZ" for "DiffW_Z"'
           call Finalize
         end if

         if(QMaster.and.Qcheck) then
           write(mfile,*) ' component number of water = ', Kcomp
           write(mfile,*) ' total step number         = ', Nsnap
           write(mfile,*) ' block step number         = ', BlockStep
           write(mfile,*) ' ########<note>#########################################'
           write(mfile,*) ' #  ******************************                     #'
           write(mfile,*) ' #  <---------ilength------------>                     #'
           write(mfile,*) ' #  <--interval-> ***************************          #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #                          ************************** #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #######################################################'
           write(mfile,*) ' ilength  = ', ilength
           write(mfile,*) ' interval = ', interv
           write(mfile,*) ' deltaZ   = ', deltaZ
           write(mfile,*) ' dtime    = ', dtime
           write(mfile,*) ' maxdz, mindz = ', maxdz, mindz
         end if

! ## radial distribtuion function between COM of lipids ( 2-dimension ) --------------

       else if( cWhat == 'GR_LipidGG' ) then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia32:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia32

           if(Script(ii)(1:6) == 'LIPID=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.

           else if(Script(ii)(1:3) == 'RC=') then

             read(Script(ii)(4:),*) Rcutoff
             AnaFlag(2) = .True.
             Rcutoff2 = Rcutoff * Rcutoff

           end if

         end do lia32

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "LIPID" for "GR_LipidGG"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RC" for "GR_LipidGG"'
           call Finalize
         end if

! ## diffusion coefficient from MSD for the COM of molecules------------

       else if( cWhat == 'DiffMSDcom' ) then

         NAnaOpt = 6
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         MSDdim = 0

lia30:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia30

           if(Script(ii)(1:7) == 'NUMMOL=') then

             read(Script(ii)(8:),*) NumMSD
             AnaFlag(1) = .True.
             allocate( MSDMol(NumMSD) )
             allocate( MSDfile(NumMSD) )

           else if(Script(ii)(1:6) == 'TOTAL=') then

             read(Script(ii)(7:),*) Nsnap
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'BLOCK=') then

             read(Script(ii)(7:),*) BlockStep
             AnaFlag(3) = .True.

           else if(Script(ii)(1:7) == 'LENGTH=') then

             read(Script(ii)(8:),*) ilength
             AnaFlag(4) = .True.

           else if(Script(ii)(1:7) == 'INTERV=') then

             read(Script(ii)(8:),*) interv
             AnaFlag(5) = .True.

           else if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(6) = .True.

           else if(Script(ii)(1:10) == 'DIMENSION=') then

             read(Script(ii)(11:),*) MSDdim

             allocate( MSDdir(MSDdim) )

             if(MSDdim==1) then
               read(Script(ii)(11:),*) jj, Ch1
               MSDdir(1) = CIdirection(Ch1)
             else if(MSDdim==2) then
               read(Script(ii)(11:),*) jj, Ch1, Ch2
               MSDdir(1) = CIdirection(Ch1)
               MSDdir(2) = CIdirection(Ch2)
             else if(MSDdim==3) then
               MSDdir(1) = 1
               MSDdir(2) = 2
               MSDdir(3) = 3
             end if

           end if

         end do lia30

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMMOL" for "DiffMSDcom"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "TOTAL" for "DiffMSDcom"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "BLOCK" for "DiffMSDcom"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "LENGTH" for "DiffMSDcom"'
           call Finalize
         end if
         if(.not.AnaFlag(5)) then
           if(QMaster) write(*,*) 'Missing "INTERV" for "DiffMSDcom"'
           call Finalize
         end if
         if(.not.AnaFlag(6)) then
           if(QMaster) write(*,*) 'Missing "DTIME" for "DiffMSDcom"'
           call Finalize
         end if
         if(MSDdim==0) then
           MSDdim = 3
           allocate( MSDdir(MSDdim) )
           MSDdir(1) = 1
           MSDdir(2) = 2
           MSDdir(3) = 3
         end if

         ii = Ist(i) + 1
         j  = 0

lia31:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia31

           if(Script(ii)(1:4) == 'MOL=') then

             j = j + 1
             write(MSDMol(j),'(a)') trim(adjustl(Script(ii)(5:)))

             if(j>NumMSD) then
               if(QMaster) write(*,*) 'ERROR : The number of MOL exceeds NumMSD'
               call Finalize
             end if

           end if

         end do lia31

         if(j/=NumMSD) then
           if(QMaster) write(*,*) 'ERROR : The total number of MOL is not NumMSD'
           call Finalize
         end if

         do j = 1 , NumMSD
           write(MSDfile(j),'(a,a,a,a,a)') &
           &  './Analy/DMSDcom_',trim(adjustl(MSDMol(j))),'.dat'
         end do

         if(QMaster.and.Qcheck) then
           write(mfile,*) ' component number   = ', NumMSD
           write(mfile,*) ' total step number  = ', Nsnap
           write(mfile,*) ' block step number  = ', BlockStep
           write(mfile,*) ' ########<note>#########################################'
           write(mfile,*) ' #  ******************************                     #'
           write(mfile,*) ' #  <----------length------------>                     #'
           write(mfile,*) ' #  <--interval-> ***************************          #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #                          ************************** #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #######################################################'
           write(mfile,*) ' ilength  = ', ilength
           write(mfile,*) ' interval = ', interv
           write(mfile,*) ' dtime    = ', dtime
           write(mfile,*) ' dimension= ', MSDdim, ' : ', MSDdir
         end if

! ## diffusion coefficient from MSD ------------------------------------

       else if( cWhat == 'DiffMSD' ) then

         NAnaOpt = 6
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         MSDdim = 0

lia33:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia33

           if(Script(ii)(1:7) == 'NUMATM=') then

             read(Script(ii)(8:),*) NumMSD
             AnaFlag(1) = .True.
             allocate( MSDAtom(NumMSD) )
             allocate( MSDResi(NumMSD) )
             allocate( MSDfile(NumMSD) )

           else if(Script(ii)(1:6) == 'TOTAL=') then

             read(Script(ii)(7:),*) Nsnap
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'BLOCK=') then

             read(Script(ii)(7:),*) BlockStep
             AnaFlag(3) = .True.

           else if(Script(ii)(1:7) == 'LENGTH=') then

             read(Script(ii)(8:),*) ilength
             AnaFlag(4) = .True.

           else if(Script(ii)(1:7) == 'INTERV=') then

             read(Script(ii)(8:),*) interv
             AnaFlag(5) = .True.

           else if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(6) = .True.

           else if(Script(ii)(1:10) == 'DIMENSION=') then

             read(Script(ii)(11:),*) MSDdim

             allocate( MSDdir(MSDdim) )

             if(MSDdim==1) then
               read(Script(ii)(11:),*) jj, Ch1
               MSDdir(1) = CIdirection(Ch1)
             else if(MSDdim==2) then
               read(Script(ii)(11:),*) jj, Ch1, Ch2
               MSDdir(1) = CIdirection(Ch1)
               MSDdir(2) = CIdirection(Ch2)
             else if(MSDdim==3) then
               MSDdir(1) = 1
               MSDdir(2) = 2
               MSDdir(3) = 3
             end if

           end if

         end do lia33

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMATM" for "DiffMSD"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "TOTAL" for "DiffMSD"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "BLOCK" for "DiffMSD"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "LENGTH" for "DiffMSD"'
           call Finalize
         end if
         if(.not.AnaFlag(5)) then
           if(QMaster) write(*,*) 'Missing "INTERV" for "DiffMSD"'
           call Finalize
         end if
         if(.not.AnaFlag(6)) then
           if(QMaster) write(*,*) 'Missing "DTIME" for "DiffMSD"'
           call Finalize
         end if
         if(MSDdim==0) then
           MSDdim = 3
           allocate( MSDdir(MSDdim) )
           MSDdir(1) = 1
           MSDdir(2) = 2
           MSDdir(3) = 3
         end if

         ii = Ist(i) + 1
         j  = 0

lia34:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia34

           if(Script(ii)(1:5) == 'ATOM=') then

             j = j + 1
             read(Script(ii)(6:),*) MSDAtom(j), MSDResi(j)

             if(j>NumMSD) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NumMSD'
               call Finalize
             end if

           end if

         end do lia34

         if(j/=NumMSD) then
           if(QMaster) write(*,*) 'ERROR : The total number of ATOM is not NumMSD'
           call Finalize
         end if

         do j = 1 , NumMSD
           write(MSDfile(j),'(a,a,a,a,a)') &
           &  './Analy/DMSD_',trim(adjustl(MSDAtom(j))),'_',trim(adjustl(MSDResi(j))),'.dat'
         end do

         if(QMaster.and.Qcheck) then
           write(mfile,*) ' component number   = ', NumMSD
           write(mfile,*) ' total step number  = ', Nsnap
           write(mfile,*) ' block step number  = ', BlockStep
           write(mfile,*) ' ########<note>#########################################'
           write(mfile,*) ' #  ******************************                     #'
           write(mfile,*) ' #  <----------length------------>                     #'
           write(mfile,*) ' #  <--interval-> ***************************          #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #                          ************************** #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #######################################################'
           write(mfile,*) ' ilength  = ', ilength
           write(mfile,*) ' interval = ', interv
           write(mfile,*) ' dtime    = ', dtime
           write(mfile,*) ' dimension= ', MSDdim, ' : ', MSDdir
         end if

       else if( cWhat == 'GR_LipidGG' ) then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia35:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia35

           if(Script(ii)(1:6) == 'LIPID=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.

           else if(Script(ii)(1:3) == 'RC=') then

             read(Script(ii)(4:),*) Rcutoff
             AnaFlag(2) = .True.
             Rcutoff2 = Rcutoff * Rcutoff

           end if

         end do lia35

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "LIPID" for "GR_LipidGG"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RC" for "GR_LipidGG"'
           call Finalize
         end if

! ## Life time of coordinating water

       else if( cWhat == 'LifeTimeW' ) then

         NAnaOpt = 8
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia36:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia36

           if(Script(ii)(1:6) == 'WATER=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.

           else if(Script(ii)(1:8) == 'NUMCOOR=') then

             read(Script(ii)(9:),*) NumGR
             AnaFlag(2) = .True.
             allocate( GRAtomI(NumGR) )
             allocate( GRResiI(NumGR) )

           else if(Script(ii)(1:3) == 'RC=') then

             read(Script(ii)(4:),*) Rcutoff
             AnaFlag(3) = .True.
             Rcutoff2 = Rcutoff * Rcutoff

           else if(Script(ii)(1:6) == 'TOTAL=') then

             read(Script(ii)(7:),*) Nsnap
             AnaFlag(4) = .True.

           else if(Script(ii)(1:6) == 'BLOCK=') then

             read(Script(ii)(7:),*) BlockStep
             AnaFlag(5) = .True.

           else if(Script(ii)(1:7) == 'LENGTH=') then

             read(Script(ii)(8:),*) ilength
             AnaFlag(6) = .True.

           else if(Script(ii)(1:7) == 'INTERV=') then

             read(Script(ii)(8:),*) interv
             AnaFlag(7) = .True.

           else if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(8) = .True.

           end if

         end do lia36

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "WATER" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "NUMCOOR" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "RC" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "TOTAL" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(5)) then
           if(QMaster) write(*,*) 'Missing "BLOCK" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(6)) then
           if(QMaster) write(*,*) 'Missing "LENGTH" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(7)) then
           if(QMaster) write(*,*) 'Missing "INTERV" for "LifeTimeW"'
           call Finalize
         end if
         if(.not.AnaFlag(8)) then
           if(QMaster) write(*,*) 'Missing "DTIME" for "LifeTimeW"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lia37:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia37

           if(Script(ii)(1:5) == 'ATOM=') then

             j = j + 1
             read(Script(ii)(6:),*) GRAtomI(j), GRResiI(j)

             if(j>NumGR) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NUMCOOR'
               call Finalize
             end if

           end if

         end do lia37

         if(j/=NumGR) then
           if(QMaster) write(*,*) 'ERROR : The total number of ATOM is not NUMCOOR'
           call Finalize
         end if

         if(QMaster.and.Qcheck) then
           write(mfile,*) ' component number of water = ', Kcomp
           write(mfile,*) ' total step number         = ', Nsnap
           write(mfile,*) ' block step number         = ', BlockStep
           write(mfile,*) ' ########<note>#########################################'
           write(mfile,*) ' #  ******************************                     #'
           write(mfile,*) ' #  <---------ilength------------>                     #'
           write(mfile,*) ' #  <--interval-> ***************************          #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #                          ************************** #'
           write(mfile,*) ' #                                                     #'
           write(mfile,*) ' #######################################################'
           write(mfile,*) ' ilength  = ', ilength
           write(mfile,*) ' interval = ', interv
           write(mfile,*) ' dtime    = ', dtime
           write(mfile,*) ' Pair num = ', NumGR
           write(mfile,*) ' PairAtom = ', ( GRAtomI(j),GRResiI(j),j=1,NumGR )
           write(mfile,*) ' Rcutoff  = ', Rcutoff
         end if

! ## radial distribution function  -----------

       else if( cWhat == 'RDF' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         QIntraSubt = .False.

lia38:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia38

           if(Script(ii)(1:8) == 'NUMPAIR=') then

             read(Script(ii)(9:),*) NumGR
             AnaFlag(1) = .True.
             allocate( GRAtomI(NumGR) )
             allocate( GRAtomJ(NumGR) )
             allocate( GRResiI(NumGR) )
             allocate( GRResiJ(NumGR) )

           else if(Script(ii)(1:3) == 'RC=') then

             read(Script(ii)(4:),*) Rcutoff
             AnaFlag(2) = .True.
             Rcutoff2 = Rcutoff * Rcutoff

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(3) = .True.

           else if(Script(ii)(1:8) == 'NOINTRA=') then

             read(Script(ii)(9:),*) cONOFF
             if(trim(cONOFF)=='ON') then
               QIntraSubt = .True.
             end if

           end if

         end do lia38

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMPAIR" for "RDF"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RC" for "RDF"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RDF"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lia39:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia39

           if(Script(ii)(1:5) == 'PAIR=') then

             j = j + 1

             read(Script(ii)(6:),*) &
             &    GRAtomI(j),GRResiI(j),GRAtomJ(j),GRResiJ(j)

             if(j>NumGR) then
               if(QMaster) write(*,*) 'ERROR : The number of PAIR exceeds NUMPAIR'
               call Finalize
             end if

           end if

         end do lia39

         if(j/=NumGR) then
           if(QMaster) write(*,*) 'ERROR : The total number of PAIR is not NUMPAIR'
           call Finalize
         end if

         allocate( GRfile(NumGR) )

         do j = 1 , NumGR

           write(GRfile(j),'(a,a,a,a,a,a,a,a,a)') &
           &  './Analy/GR_',&
           &   trim(adjustl(GRAtomI(j))),'_',trim(adjustl(GRResiI(j))),'--',&
           &   trim(adjustl(GRAtomJ(j))),'_',trim(adjustl(GRResiJ(j))),'.dat'

         end do


! ## radial distribution function between molecular COMs -----------

       else if( cWhat == 'RDFG' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

wia38:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit wia38

           if(Script(ii)(1:8) == 'NUMPAIR=') then

             read(Script(ii)(9:),*) NumGR
             AnaFlag(1) = .True.
             allocate( GRResiI(NumGR) )
             allocate( GRResiJ(NumGR) )

           else if(Script(ii)(1:3) == 'RC=') then

             read(Script(ii)(4:),*) Rcutoff
             AnaFlag(2) = .True.
             Rcutoff2 = Rcutoff * Rcutoff

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(3) = .True.

           end if

         end do wia38

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMPAIR" for "RDFG"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RC" for "RDFG"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RDFG"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

wia39:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit wia39

           if(Script(ii)(1:5) == 'PAIR=') then

             j = j + 1

             read(Script(ii)(6:),*) &
             &    GRResiI(j),GRResiJ(j)

             if(j>NumGR) then
               if(QMaster) write(*,*) 'ERROR : The number of PAIR exceeds NUMPAIR'
               call Finalize
             end if

           end if

         end do wia39

         if(j/=NumGR) then
           if(QMaster) write(*,*) 'ERROR : The total number of PAIR is not NUMPAIR'
           call Finalize
         end if

         allocate( GRfile(NumGR) )

         do j = 1 , NumGR

           write(GRfile(j),'(a,a,a,a,a,a,a,a,a)') &
           &  './Analy/GR_',&
           &   trim(adjustl(GRResiI(j))),'--',trim(adjustl(GRResiJ(j))),'.dat'

         end do


! ## distribution function along the bilayer normal -----------

       else if( cWhat == 'RZ' ) then

         NAnaOpt = 4
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia40:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia40

           if(Script(ii)(1:8) == 'NUMATOM=') then

             read(Script(ii)(9:),*) NumRZ
             AnaFlag(1) = .True.
             allocate( RZAtom(NumRZ) )
             allocate( RZResi(NumRZ) )

           else if(Script(ii)(1:5) == 'AXIS=') then

             read(Script(ii)(6:),*) Ch1
             RZdir = CIdirection(Ch1)
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'RANGE=') then

             read(Script(ii)(7:),*) Xmin, Xmax
             AnaFlag(3) = .True.

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(4) = .True.

           end if

         end do lia40

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMATOM" for "RZ"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "AXIS" for "RZ"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "RANGE" for "RZ"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RZ"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lia41:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia41

           if(Script(ii)(1:5) == 'ATOM=') then

             j = j + 1
             read(Script(ii)(6:),*) RZAtom(j), RZResi(j)

             if(j>NumRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NUMATOM'
               call Finalize
             end if

           end if

         end do lia41

         if(j/=NumRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of ATOM is not NUMATOM'
           call Finalize
         end if

         allocate( RZfile(NumRZ) )

         do j = 1 , NumRZ

           write(RZfile(j),'(a,a,a,a,a)') &
           &  './Analy/RZ_',trim(adjustl(RZAtom(j))),'_',trim(adjustl(RZResi(j))),'.dat'

         end do


! ## distribution function along the bilayer normal -----------

       else if( cWhat == 'RZG' ) then

         NAnaOpt = 4
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

qia40:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit qia40

           if(Script(ii)(1:8) == 'NUMMOL=') then

             read(Script(ii)(9:),*) NumRZ
             AnaFlag(1) = .True.
             allocate( RZResi(NumRZ) )

           else if(Script(ii)(1:5) == 'AXIS=') then

             read(Script(ii)(6:),*) Ch1
             RZdir = CIdirection(Ch1)
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'RANGE=') then

             read(Script(ii)(7:),*) Xmin, Xmax
             AnaFlag(3) = .True.

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(4) = .True.

           end if

         end do qia40

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMMOL" for "RZG"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "AXIS" for "RZG"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "RANGE" for "RZG"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RZG"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

qia41:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit qia41

           if(Script(ii)(1:5) == 'MOL=') then

             j = j + 1
             read(Script(ii)(6:),*) RZResi(j)

             if(j>NumRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NUMATOM'
               call Finalize
             end if

           end if

         end do qia41

         if(j/=NumRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of MOL is not NUMMOL'
           call Finalize
         end if

         allocate( RZfile(NumRZ) )

         do j = 1 , NumRZ

           write(RZfile(j),'(a,a,a)') &
           &  './Analy/RZG_',trim(adjustl(RZResi(j))),'.dat'

         end do


! ## radial distribution function from the center of the simulation box -----------

       else if( cWhat == 'RD' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         NComp = 0

laa40:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit laa40

           if(Script(ii)(1:8) == 'NUMATOM=') then

             read(Script(ii)(9:),*) NumRZ
             AnaFlag(1) = .True.
             allocate( RZAtom(NumRZ) )
             allocate( RZResi(NumRZ) )

           else if(Script(ii)(1:6) == 'RANGE=') then

             read(Script(ii)(7:),*) Xmin, Xmax
             AnaFlag(2) = .True.

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(3) = .True.

           else if(Script(ii)(1:8) == 'CORRCOM=') then

             read(Script(ii)(9:),*) NComp

           end if

         end do laa40

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMATOM" for "RD"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RANGE" for "RD"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RD"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

laa41:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit laa41

           if(Script(ii)(1:5) == 'ATOM=') then

             j = j + 1
             read(Script(ii)(6:),*) RZAtom(j), RZResi(j)

             if(j>NumRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NUMATOM'
               call Finalize
             end if

           end if

         end do laa41

         if(j/=NumRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of ATOM is not NUMATOM'
           call Finalize
         end if

         allocate( RZfile(NumRZ) )

         do j = 1 , NumRZ

           write(RZfile(j),'(a,a,a,a,a)') &
           &  './Analy/RD_',trim(adjustl(RZAtom(j))),'_',trim(adjustl(RZResi(j))),'.dat'

         end do


! ## radial distribution function from the center of the simulation box -----------

       else if( cWhat == 'RDdisk' ) then

         NAnaOpt = 4
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lsa40:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lsa40

           if(Script(ii)(1:8) == 'NUMATOM=') then

             read(Script(ii)(9:),*) NumRZ
             AnaFlag(1) = .True.
             allocate( RZAtom(NumRZ) )
             allocate( RZResi(NumRZ) )

           else if(Script(ii)(1:6) == 'RANGE=') then

             read(Script(ii)(7:),*) Xmin, Xmax
             AnaFlag(2) = .True.

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(3) = .True.

           else if(Script(ii)(1:5) == 'AXIS=') then

             read(Script(ii)(6:),*) Ch1
             RZdir = CIdirection(Ch1)
             AnaFlag(4) = .True.

           end if

         end do lsa40

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMATOM" for "RDdist"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RANGE" for "RDdist"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RDdist"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "AXIS" for "RDdist"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lsa41:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lsa41

           if(Script(ii)(1:5) == 'ATOM=') then

             j = j + 1
             read(Script(ii)(6:),*) RZAtom(j), RZResi(j)

             if(j>NumRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NUMATOM'
               call Finalize
             end if

           end if

         end do lsa41

         if(j/=NumRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of ATOM is not NUMATOM'
           call Finalize
         end if

         allocate( RZfile(NumRZ) )

         do j = 1 , NumRZ

           write(RZfile(j),'(a,a,a,a,a)') &
           &  './Analy/RDdisk_',trim(adjustl(RZAtom(j))),'_',trim(adjustl(RZResi(j))),'.dat'

         end do


! ## radial distribution function from the center of the simulation box -----------

       else if( cWhat == 'RDdiskCOM' ) then

         NAnaOpt = 4
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lta40:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lta40

           if(Script(ii)(1:7) == 'NUMMOL=') then

             read(Script(ii)(8:),*) NumRZ
             AnaFlag(1) = .True.
             allocate( RZResi(NumRZ) )

           else if(Script(ii)(1:6) == 'RANGE=') then

             read(Script(ii)(7:),*) Xmin, Xmax
             AnaFlag(2) = .True.

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(3) = .True.

           else if(Script(ii)(1:5) == 'AXIS=') then

             read(Script(ii)(6:),*) Ch1
             RZdir = CIdirection(Ch1)
             AnaFlag(4) = .True.

           end if

         end do lta40

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMMOL" for "RDdistCOM"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RANGE" for "RDdistCOM"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "RDdistCOM"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "AXIS" for "RDdistCOM"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lta41:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lta41

           if(Script(ii)(1:5) == 'MOL=') then

             j = j + 1
             read(Script(ii)(6:),*) RZResi(j)

             if(j>NumRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of MOL exceeds NUMMOL'
               call Finalize
             end if

           end if

         end do lta41

         if(j/=NumRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of MOL is not NUMMOL'
           call Finalize
         end if

         allocate( RZfile(NumRZ) )

         do j = 1 , NumRZ

           write(RZfile(j),'(a,a,a)') &
           &  './Analy/RDdiskCOM_',trim(adjustl(RZResi(j))),'.dat'

         end do

! ## distribution function along the bilayer normal -----------

       else if( cWhat == 'RZCyl' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia42:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia42

           if(Script(ii)(1:8) == 'NUMATOM=') then

             read(Script(ii)(9:),*) NumRZ
             AnaFlag(1) = .True.
             allocate( RZAtom(NumRZ) )
             allocate( RZResi(NumRZ) )

           else if(Script(ii)(1:7) == 'REFMOL=') then

             read(Script(ii)(8:),*) Kcomp
             AnaFlag(2) = .True.

           else if(Script(ii)(1:7) == 'RADIUS=') then

             read(Script(ii)(8:),*) RadCyl
             AnaFlag(3) = .True.
             RadCyl2 = RadCyl * RadCyl

           end if

         end do lia42

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMATOM" for "RZCyl"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "REFMOL" for "RZCyl"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "RADIUS" for "RZCyl"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lia43:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia43

           if(Script(ii)(1:5) == 'ATOM=') then

             j = j + 1
             read(Script(ii)(6:),*) RZAtom(j), RZResi(j)

             if(j>NumRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of ATOM exceeds NUMATOM'
               call Finalize
             end if

           end if

         end do lia43

         if(j/=NumRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of ATOM is not NUMATOM'
           call Finalize
         end if

         allocate( RZfileI(NumRZ) )
         allocate( RZfileO(NumRZ) )

         do j = 1 , NumRZ

           write(RZfileI(j),'(a,a,a,a,a)') &
           &  './Analy/RZint_',trim(adjustl(RZAtom(j))),'_',trim(adjustl(RZResi(j))),'.dat'
           write(RZfileO(j),'(a,a,a,a,a)') &
           &  './Analy/RZout_',trim(adjustl(RZAtom(j))),'_',trim(adjustl(RZResi(j))),'.dat'

         end do


! ## for Residue -----------------------------

       else if ( ( cWhat == 'RZRes'    ) .or. &
       &         ( cWhat == 'ResidTrj' ) ) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia44:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia44

           if(Script(ii)(1:7) == 'NUMRES=') then

             read(Script(ii)(8:),*) NumResRZ
             AnaFlag(1) = .True.
             allocate( ResRZ(NumRZ) )

           end if

         end do lia44

         if(.not.AnaFlag(1)) then
           if(cWhat == 'RZRes') then
             if(QMaster) write(*,*) 'Missing "NUMRES" for "RZRes"'
           else if(cWhat == 'ResidTrj') then
             if(QMaster) write(*,*) 'Missing "NUMRES" for "ResidTrj"'
           end if
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lia45:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia45

           if(Script(ii)(1:8) == 'RESIDUE=') then

             j = j + 1
             read(Script(ii)(9:),*) ResRZ(j)

             if(j>NumResRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of RESIDUE exceeds NUMRES'
               call Finalize
             end if

           end if

         end do lia45

         if(j/=NumResRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of RESIDUE is not NUMRES'
           call Finalize
         end if

       else if( cWhat == 'RZResG' ) then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia46:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia46

           if(Script(ii)(1:7) == 'NUMRES=') then

             read(Script(ii)(8:),*) NumResRZ
             AnaFlag(1) = .True.
             allocate( ResRZ(NumRZ) )

           else if(Script(ii)(1:7) == 'REFMOL=') then

             read(Script(ii)(8:),*) Kcomp
             AnaFlag(2) = .True.

           end if

         end do lia46

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMRES" for "RZResG"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "REFMOL" for "RZResG"'
           call Finalize
         end if

         ii = Ist(i) + 1
         j  = 0

lia47:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia47

           if(Script(ii)(1:8) == 'RESIDUE=') then

             j = j + 1
             read(Script(ii)(9:),*) ResRZ(j)

             if(j>NumResRZ) then
               if(QMaster) write(*,*) 'ERROR : The number of RESIDUE exceeds NUMRES'
               call Finalize
             end if

           end if

         end do lia47

         if(j/=NumResRZ) then
           if(QMaster) write(*,*) 'ERROR : The total number of RESIDUE is not NUMRES'
           call Finalize
         end if

! ## for water-related -----------------------------

       else if ( ( cWhat == 'WatOriCyl'  ) .or. &
       &         ( cWhat == 'WatOriCylG' ) ) then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia48:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia48

           if(Script(ii)(1:8) == 'PROTEIN=') then

             read(Script(ii)(9:),*) Kcomp
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of protein is tagged as ', Kcomp

           else if(Script(ii)(1:7) == 'RADIUS=') then

            read(Script(ii)(8:),*) RadCyl
            if(QMaster.and.Qcheck) &
            & write(mfile,*) '  radius of the cylinder = ', RadCyl

            RadCyl2 = RadCyl * RadCyl

           end if

         end do lia48

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "PROTEIN" for "WatOriCyl"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RADIUS" for "WatOriCyl"'
           call Finalize
         end if

! ## for water-related -----------------------------

       else if ( ( cWhat == 'PoreWatTrj' ) .or. &
       &         ( cWhat == 'PoreWater'  ) ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia49:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia49

           if(Script(ii)(1:8) == 'PROTEIN=') then

             read(Script(ii)(9:),*) Kcomp
             AnaFlag(1) = .True.

           else if(Script(ii)(1:7) == 'LENGTH=') then

             read(Script(ii)(8:),*) Zsh
             AnaFlag(2) = .True.

           else if(Script(ii)(1:7) == 'RADIUS=') then

             read(Script(ii)(8:),*) RadCyl
             AnaFlag(3) = .True.
             RadCyl2 = RadCyl * RadCyl

           end if

         end do lia49

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "PROTEIN"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "LENGTH"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "RADIUS"'
           call Finalize
         end if

! ## for MSD -----------------------------

       else if ( ( cWhat == 'MSD_P'      ) .or. &
       &         ( cWhat == 'MSD_PD'     ) .or. &
       &         ( cWhat == 'MSD_Pr'     ) .or. &
       &         ( cWhat == 'MSD_PrD'    ) .or. &
       &         ( cWhat == 'MSD_Pm'     ) .or. &
       &         ( cWhat == 'MSD_PmD'    ) .or. &
       &         ( cWhat == 'MSD_PmD_MS' ) ) then

         Nsnap=0

         do j = 1 , NJobs

           Nsnap = Nsnap + NTrjStep(j)

         end do

         ResidFlag = 0

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia50:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia50

           if(Script(ii)(1:8) == 'PROTEIN=') then

             read(Script(ii)(9:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of protein is tagged as ', Kcomp

           else if(Script(ii)(1:3) == 'SP=') then

             read(Script(ii)(4:),*) FROM,TO

             do k = FROM, TO
               ResidFlag(k) = 1
             end do

           end if

         end do lia50

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "PROTEIN" for "MSD"'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if(cWhat=='RofG_Prot') then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia51:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia51

           if(Script(ii)(1:8) == 'PROTEIN=') then

             read(Script(ii)(9:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of protein is tagged as ', Kcomp

           end if

         end do lia51

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "PROTEIN" for "RofG_Prot"'
           call Finalize
         end if

! ## Voronoi tesselation -------------------------------------

       else if((cWhat=='VoronoiL').or.(cWhat=='VoronoiC')) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia52:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia52

           if(Script(ii)(1:6) == 'LIPID=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of LIPID is tagged as ', Kcomp

           end if

         end do lia52

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "LIPID" for "Voronoi"'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if(cWhat=='AreaOccL') then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia53:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia53

           if(Script(ii)(1:6) == 'LIPID=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of LIPID is tagged as ', Kcomp

           else if(Script(ii)(1:5) == 'AREA=') then

             read(Script(ii)(6:),*) AreaL
             AnaFlag(2) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  Averaged molecular area of the LIPID is ', AreaL

           end if

         end do lia53

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "LIPID" for "AreaOccL"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "AREA" for "AreaOccL"'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if(cWhat=='TraceG') then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia54:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia54

           if(Script(ii)(1:7) == 'MOLNUM=') then

             read(Script(ii)(8:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) write(mfile,*) &
             & '  The component number of the molecule to be traced is ', Kcomp

           else if(Script(ii)(1:5) == 'FREQ=') then

             read(Script(ii)(6:),*) Nsnap
             AnaFlag(2) = .True.
             if(QMaster.and.Qcheck) then
               write(mfile,*) '  Trajectory stored at evergy', Nsnap,' steps'
             end if

           end if

         end do lia54

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "MOLNUM" for "TraceG"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "FREQ" for "TraceG"'
           call Finalize
         end if

! ## ---------------------------------------------------------

       else if(cWhat=='LS_anim') then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia55:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia55

           if(Script(ii)(1:9) == 'NUMLIPID=') then

             read(Script(ii)(10:),*) Nlipid
             AnaFlag(1) = .True.

           end if

         end do lia55

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "NUMLIPID" for "LS_anim"'
           call Finalize
         end if

         if(QMaster.and.Qcheck) write(mfile,'(a,i3,a)') &
         & 'The ',Nlipid,'-th lipid molecule was chosen to be visualized'

! ## Coarse Graining -----------------------------------------

       else if((cWhat=='CGdata').or.(cWhat=='CGalkane').or.&
       &       (cWhat=='CGW').or.(cWhat=='CGchainW')) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         Temp_o = 303.15d0

lia56:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia56

           if(Script(ii)(1:3) == 'RC=') then

             read(Script(ii)(4:),*) Rcutoff
             AnaFlag(1) = .True.

           else if(Script(ii)(1:5) == 'TEXT=') then

             read(Script(ii)(6:),*) Temp_o

           end if

         end do lia56

         Rcutoff2 = Rcutoff * Rcutoff

         if(QSwitch) then
           Ron  = Rcutoff - 2.d0
           Ron2 = Ron * Ron
           swf1   = 1.d0 / (Rcutoff2 - Ron2) ** 3
         end if

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "RC" for "CG"'
           call Finalize
         end if


! ## Electron density profile -----------

       else if( cWhat == 'Edens' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia57:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia57

           if(Script(ii)(1:5) == 'AXIS=') then

             read(Script(ii)(6:),*) Ch1
             RZdir = CIdirection(Ch1)
             AnaFlag(1) = .True.

           else if(Script(ii)(1:6) == 'RANGE=') then

             read(Script(ii)(7:),*) Xmin, Xmax
             AnaFlag(2) = .True.

           else if(Script(ii)(1:5) == 'GRID=') then

             read(Script(ii)(6:),*) DRgrid
             AnaFlag(3) = .True.

           end if

         end do lia57

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "AXIS" for "Edens"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "RANGE" for "Edens"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "GRID" for "Edens"'
           call Finalize
         end if

! ## probability distribution of bond, bending, and dihedral pairs -----------

       else if( cWhat == 'DISTinMOL' ) then

         DRgrid = 0.1d0
         DEGgrid = 1.d0
         Rrangemin = 2.d0
         Rrangemax = 8.d0
         ibond = 0
         iangle = 0
         idihed = 0

lia58:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia58

           if(Script(ii)(1:6) == 'RGRID=') then

             read(Script(ii)(7:),*) DRgrid

           else if(Script(ii)(1:6) == 'DGRID=') then

             read(Script(ii)(7:),*) DEGgrid

           else if(Script(ii)(1:7) == 'RRANGE=') then

             read(Script(ii)(8:),*) Rrangemin, Rrangemax

           else if(Script(ii)(1:4) == 'BOND') then

             ibond = ibond + 1

           else if(Script(ii)(1:5) == 'ANGLE') then

             iangle = iangle + 1

           else if(Script(ii)(1:8) == 'DIHEDRAL') then

             idihed = idihed + 1

           end if

         end do lia58

         NBdAna = ibond
         NAnAna = iangle
         NDhAna = idihed
         if(NBdAna/=0) then
           allocate(BdAtomI(NBdAna))
           allocate(BdAtomJ(NBdAna))
           allocate(BdMol(NBdAna))
         end if
         if(NAnAna/=0) then
           allocate(AnAtomI(NAnAna))
           allocate(AnAtomJ(NAnAna))
           allocate(AnAtomK(NAnAna))
           allocate(AnMol(NAnAna))
         end if
         if(NDhAna/=0) then
           allocate(DhAtomI(NDhAna))
           allocate(DhAtomJ(NDhAna))
           allocate(DhAtomK(NDhAna))
           allocate(DhAtomL(NDhAna))
           allocate(DhMol(NDhAna))
         end if

         ibond = 0
         iangle = 0
         idihed = 0

         ii = Ist(i) + 1

lia59:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia59

           if(Script(ii)(1:4) == 'BOND') then

             ibond = ibond + 1
             read(Script(ii)(5:),*) BdMol(ibond), BdAtomI(ibond), BdAtomJ(ibond)

           else if(Script(ii)(1:5) == 'ANGLE') then

             iangle = iangle + 1
             read(Script(ii)(6:),*) AnMol(iangle), AnAtomI(iangle), AnAtomJ(iangle), &
             &                                 AnAtomK(iangle)

           else if(Script(ii)(1:8) == 'DIHEDRAL') then

             idihed = idihed + 1
             read(Script(ii)(9:),*) DhMol(idihed), DhAtomI(idihed), DhAtomJ(idihed), &
             &                                 DhAtomK(idihed), DhAtomL(idihed)

           end if

         end do lia59

       else if( cWhat == 'COORDIS' ) then

         NAnaOpt = 5
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

         mex = 6
         nex = 12
         DRgrid = 0.2d0
         DNgrid = 0.1d0

liw60:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit liw60

           if(Script(ii)(1:9) == 'MOLNAMES=') then

             read(Script(ii)(10:),*) MOLNAME1,MOLNAME2

             imol = 0
             jmol = 0
             do k = 1, NumSpec
               if(MolName(k)==MOLNAME1) then
                 imol = NumMol(k)
               else if(MolName(k)==MOLNAME2) then
                 jmol = NumMol(k)
               end if
             end do
             if(imol==0) then
               write(*,*) "ERROR: wrong name for molecule 1"
               call Finalize
             end if
             if(jmol==0) then
               write(*,*) "ERROR: wrong name for molecule 2"
               call Finalize
             end if

             AnaFlag(1) = .True.

           else if(Script(ii)(1:7) == 'ATOM_I=') then

             write(Script1,*) trim(adjustl(Script(ii)(8:80)))
             jj = len(Script1)
             iatom = 0
lpa:         do j = 2, jj
               if(Script1(j:j)=='!') then
                 if(Script1(j-1:j-1)/=' ') then
                   iatom = iatom + 1
                 end if
                 exit lpa
               end if
               if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
                 iatom = iatom + 1
               end if
             end do lpa
             allocate(NatomI(iatom))
             read(Script1,*) (NAtomI(j),j=1,iatom)

             AnaFlag(2) = .True.

           else if(Script(ii)(1:7) == 'ATOM_J=') then

             write(Script1,*) trim(adjustl(Script(ii)(8:80)))
             jj = len(Script1)
             jatom = 0
lpb:         do j = 2, jj
               if(Script1(j:j)=='!') then
                 if(Script1(j-1:j-1)/=' ') then
                   jatom = jatom + 1
                 end if
                 exit lpb
               end if
               if(Script1(j:j)==' '.and.Script1(j-1:j-1)/=' ') then
                 jatom = jatom + 1
               end if
             end do lpb
             allocate(NatomJ(jatom))
             read(Script1,*) (NatomJ(j),j=1,jatom)

             AnaFlag(3) = .True.

           else if(Script(ii)(1:6) == 'RGCUT=') then

             read(Script(ii)(7:),*) Rcutoff
             Rcutoff2 = Rcutoff*Rcutoff
             AnaFlag(4) = .True.

           else if(Script(ii)(1:2) == 'D=') then

             read(Script(ii)(3:),*) dAB
             AnaFlag(5) = .True.

           else if(Script(ii)(1:2) == 'M=') then

             read(Script(ii)(3:),*) mex

           else if(Script(ii)(1:2) == 'N=') then

             read(Script(ii)(3:),*) nex

           else if(Script(ii)(1:3) == 'DR=') then

             read(Script(ii)(4:),*) DRgrid

           else if(Script(ii)(1:3) == 'DN=') then

             read(Script(ii)(4:),*) DNgrid

           end if

         end do liw60

         ngrid = int(Rcutoff/DRgrid)
         ncgrid = int(dble(iatom*jatom)/DNgrid)

!         print *, "iatom=",iatom
!         print *, "jatom=",jatom
         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "MOLNAMES" for "COORDIS"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "ATOM_I" for "COORDIS"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "ATOM_J" for "COORDIS"'
           call Finalize
         end if
         if(.not.AnaFlag(4)) then
           if(QMaster) write(*,*) 'Missing "RGCUT" for "COORDIS"'
           call Finalize
         end if
         if(.not.AnaFlag(5)) then
           if(QMaster) write(*,*) 'Missing "D" for "COORDIS"'
           call Finalize
         end if

         if(mod(mex,2)/=0) then
           if(QMaster) write(*,*) 'the power "m" should be an even number'
           call Finalize
         end if
         if(mod(nex,2)/=0) then
           if(QMaster) write(*,*) 'the power "n" should be an even number'
           call Finalize
         end if

         mex = mex/2
         nex = nex/2

       else if( cWhat == 'INERTMOM' ) then

         NAnaOpt = 1
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

liw61:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit liw61

           if(Script(ii)(1:7) == 'TARGET=') then

             read(Script(ii)(8:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of the selected molecule is tagged as ', Kcomp

           end if

         end do liw61

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "TARGET" for "INERTMOM"'
           call Finalize
         end if

       else if( cWhat == 'ID_MICELLE' ) then

         NAnaOpt = 2
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

lia60:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia60

           if(Script(ii)(1:6) == 'LIPID=') then

             read(Script(ii)(7:),*) Kcomp
             AnaFlag(1) = .True.
             if(QMaster.and.Qcheck) &
             & write(mfile,*) '  The component number of LIPID is tagged as ', Kcomp

           else if(Script(ii)(1:9) == 'DISTANCE=') then

             read(Script(ii)(10:),*) Rcutoff
             Rcutoff2 = Rcutoff*Rcutoff
             AnaFlag(2) = .True.
             if(QMaster.and.Qcheck) write(mfile,*) &
             & '  The micelle criteria for the distance between carbons is', Rcutoff

           end if

         end do lia60

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "LIPID" for "ID_MICELLE"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "DISTANCE" for "ID_MICELLE"'
           call Finalize
         end if

       else if( cWhat == 'STRESSPROF' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.
         RcutC2 = Rcutoff2
         Qslab = .True.
         NComp = 0

lia62:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia62

           if(Script(ii)(1:8) == 'CONTOUR=') then

             read(Script(ii)(9:),*) CName
             if(CName=='IK') then
               QIK = .True.
             else if(CName(1:1)=='H') then
               QIK = .False.
             else
               if(QMaster) write(*,*) 'error: wrong contour name'
               call Finalize
             end if
             AnaFlag(1) = .True.

           else if(Script(ii)(1:9) == 'GEOMETRY=') then

             read(Script(ii)(10:),*) CName
             if(CName=='SPHERICAL') then
               Qslab = .False.
             else if(CName=='PLANER') then
               Qslab = .True.
             else
               if(QMaster) write(*,*) 'error: wrong geometry name'
               call Finalize
             end if

           else if(Script(ii)(1:8) == 'CORRCOM=') then

             read(Script(ii)(9:),*) NComp

           else if(Script(ii)(1:5) == 'RBIN=') then

             read(Script(ii)(6:),*) dRbin
             InvdRbin = 1.d0/dRbin

           else if(Script(ii)(1:6) == 'NGRID=') then

             read(Script(ii)(7:),*) Nbin
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'EWALD=') then

             read(Script(ii)(7:),*) cONOFF
             if(cONOFF=='ON') then
               QEW = .True.
             else if(cONOFF=='OFF') then
               QEW = .False.
             else
               if(QMaster) write(*,*) 'error : wrong Ewald flag in STRESSPROF'
               call Finalize
             end if
             AnaFlag(3) = .True.

           else if(Script(ii)(1:4) == 'RCC=') then

             read(Script(ii)(5:),*) Rcutoff
             RcutC2 = Rcutoff * Rcutoff

           end if

         end do lia62

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "CONTOUR" for "STRESSPROF"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "NGRID" for "STRESSPROF"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "EWALD" for "STRESSPROF"'
           call Finalize
         end if

       else if( cWhat == 'Rot5Ring' ) then

         NAnaOpt = 3
         allocate( AnaFlag(NAnaOpt) )
         AnaFlag = .False.

         allocate(MSDMol(1))
         allocate(MSDAtom(5))

lia61:   do

           ii = ii + 1

           if(Script(ii)(1:2) == '<<') exit lia61

           if(Script(ii)(1:8) == 'MOLNAME=') then

             read(Script(ii)(9:),*) MSDMol(1)
             AnaFlag(1) = .True.

           else if(Script(ii)(1:6) == 'ATOMS=') then

             read(Script(ii)(7:),*) (MSDAtom(k),k=1,5)
             AnaFlag(2) = .True.

           else if(Script(ii)(1:6) == 'DTIME=') then

             read(Script(ii)(7:),*) dtime
             AnaFlag(3) = .True.

           end if

         end do lia61

         if(.not.AnaFlag(1)) then
           if(QMaster) write(*,*) 'Missing "MOLNAME" for "Rot5Ring"'
           call Finalize
         end if
         if(.not.AnaFlag(2)) then
           if(QMaster) write(*,*) 'Missing "ATOMS" for "Rot5Ring"'
           call Finalize
         end if
         if(.not.AnaFlag(3)) then
           if(QMaster) write(*,*) 'Missing "DTIME" for "Rot5Ring"'
           call Finalize
         end if

#ifdef EnergyRep
       else if( cWhat == 'ENERGYHIST' ) then
         call Read_AnaER(ii)
#endif
       end if

       ReadFlag(i) = .True.

       exit

     end if

   end do

end subroutine Read_Analyz_Cond
