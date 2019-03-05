

!######################################################################
!######################################################################


subroutine Datainputs

use Forplot

integer :: i
character(len=80) :: String
character(len=1) :: Ch

   print *, ' Simulation Method = ( MD / HMC / PIMD / CMD ) ? '
   read(5,'(a)') String
   Method = trim(adjustl(String))

   print *, ' Ensemble = ( NVT / NPT / NtT / NVE / NPH / NE(isolate) / NT(isolate) ) ? '
   read(5,'(a)') String
   Ensemble = trim(adjustl(String))

   print *, ' Calculate accumulative averages ? [ Y or N ]'
   read(5,'(a)') Ch
   if((Ch=='Y').or.(Ch=='y')) then
     QAccAv = .True.
   else
     QAccAv = .False.
   end if

   if(QAccAv) then
   print *, ' Number of steps to be ignored for the average calculation = ? '
   read(5,*) Nskip
   end if

   print *, ' Number of monitor files to be plotted = ? '
   read(5,*) Nfile

   print *, ' The step interval to be plotted = ? '
   read(5,*) Nst

#ifdef SURF
   print *, ' The length of the simulation box in the Z-direction: Lz = ? '
   read(5,*) H(3,3)
#endif

   allocate( Energy_File(Nfile) )

   print *, ' Filenames (  1) = ? '

   do i = 1 , Nfile
     read(5,'(a)') Energy_File(i)
     if(i/=Nfile) print '(a,i3,a)', ' Filenames (',i+1,') = ? '
   end do


   print *, '****************************'
   print *, '**  Method   = ', trim(Method),'      **'
   print *, '**  Ensemble = ', trim(Ensemble),'       **'
   print *, '****************************'

   print *, ''
   print *, '####  Original data files are ...  ###'
   print *, ''

   do i = 1, Nfile
     write(*,'(3x,i3,2x,a)') i,trim(Energy_file(i))
   end do

end subroutine Datainputs


!######################################################################
!######################################################################


