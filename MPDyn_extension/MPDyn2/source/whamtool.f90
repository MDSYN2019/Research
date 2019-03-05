module commonblocks
integer :: nfile
integer :: ngrid
real(8), dimension(:,:), allocatable :: hist
real(8), dimension(:), allocatable :: kconst, R0
integer, dimension(:), allocatable :: nstep
real(8), dimension(:), allocatable :: hpos
end module commonblocks

program WHAM_tool

use commonblocks

implicit none

integer :: i, j, k, n
integer :: eofile
character(len=80), dimension(:), allocatable :: Filename
real(8) :: Rup,Rfm,dr,invdr
real(8) :: Timeps, Fconst, Fav, Rpos
real(8) :: x0, x1
character(len=1) :: Ch
character(len=80) :: String, fname

  print *, "*****************************************"
  print *, "This is a WHAM tool that makes histogram "
  print *, "using the output of MPDyn with the option"
  print *, "'OPTIONALCONSTRAINT'                     "
  print *, "*****************************************"

! ##
  print *, "The range of the histogram ?"
  print *, "Two values are needed; From ?? to ?? in Angstrom"

  read(5,*) Rfm, Rup

! ##
  print *, "The number of grid points ?"

  read(5,*) ngrid

  dr = (Rup-Rfm)/dble(ngrid)

! ##
  print *, "delta R=", dr
  invdr = 1./dr

  if(dr>0.1d0) then
    write(*,*) "WARNING: the grid size is", dr
    write(*,*) "You may increase the number of grid points"
  end if

! ## reading the condition files
  print *, "The number of files?"
  read(*,*) nfile

  allocate(Filename(nfile),nstep(nfile),kconst(nfile),R0(nfile))

! ##
   do i = 1, nfile
     print *, "file name, bias potential parameters, k [kcal/mol/A^2] & R0 [A], for",i
     read(*,*) Filename(i), kconst(i), R0(i)
   end do

! ##
  allocate(hist(ngrid,nfile))
  hist(:,:) = 0.d0
  allocate(hpos(ngrid))

! ##
  print *, "The number of files to be used is",nfile
  print *, "The file names are ..."
   do i = 1, nfile
     print *, Filename(i)
   end do

  print '(/a/)', "Make histogram?[Y] or Read from existing files?[N]"
  read(*,*) Ch

  if(Ch=='Y'.or.Ch=='y') then

    print '(/a/)', "Making histogram ..."

     do i = 1, nfile
       open(1,file=trim(adjustl(Filename(i))))
       print '(a)', "reading ...",Filename(i)
       nstep(i) = 0
       do
         read(1,'(a)',iostat=eofile) String
         if(eofile == -1) exit
         if(String(1:1)=='#') cycle
         nstep(i) = nstep(i) + 1
         read(String,*) Rpos
         if(Rpos<Rfm.or.Rpos>Rup) then
           write(*,*) "WARNING: the distance is out of the range!!"
           nstep(i) = nstep(i) - 1
           cycle
         end if
         x0 = Rpos - Rfm
         x1 = x0*invdr
         k  = int(x1) + 1
         hist(k,i) = hist(k,i) + 1.d0
       end do
       close(1)
       print *, nstep(i)
     end do

     do i = 1, nfile
       do j = 1, ngrid
         hist(j,i) = hist(j,i)/real(nstep(i))
       end do
     end do

     do i = 1, ngrid
       hpos(i) = (i-0.5d0)*dr+Rfm
     end do

     print *, "Histgram is prepared !"

     do i = 1, nfile
       if(i<10) then
       write(fname,'(a,i1,a)') 'Hist_',i,'.dat'
       else if(i<100) then
       write(fname,'(a,i2,a)') 'Hist_',i,'.dat'
       else if(i<1000) then
       write(fname,'(a,i3,a)') 'Hist_',i,'.dat'
       else
       write(*,*) 'ERROR: too many files'
       stop
       end if
       open(1,file=trim(fname))
       write(1,'(a,i12)') '#',nstep(i)
       do j = 1, ngrid
         write(1,'(f9.5,f15.8)') hpos(j),hist(j,i)
       end do
       close(1)
     end do

   else

     do i = 1, nfile
       if(i<10) then
       write(fname,'(a,i1,a)') 'Hist_',i,'.dat'
       else if(i<100) then
       write(fname,'(a,i2,a)') 'Hist_',i,'.dat'
       else if(i<1000) then
       write(fname,'(a,i3,a)') 'Hist_',i,'.dat'
       else
       write(*,*) 'ERROR: too many files'
       stop
       end if
       open(1,file=trim(fname))
       read(1,'(a1,i12)') Ch,nstep(i)
       do j = 1, ngrid
         read(1,'(f9.5,f15.8)') hpos(j),hist(j,i)
       end do
       close(1)
     end do

   end if

   print *, "Histogram analysis will start now"

   call WHAM

end program WHAM_tool

! #######################################################
! #######################################################

subroutine WHAM

use commonblocks

implicit none

integer :: i, j, k, iteration, ii
real(8), dimension(nfile) :: Fi
real(8), dimension(ngrid) :: Rho,PMF
real(8), parameter :: kb = 1.98719137d-03  ! [kcal/mol/K]
real(8) :: Temp, kT, wj
real(8), dimension(ngrid) :: totn
real(8), dimension(ngrid,nfile) :: weigh
logical, dimension(ngrid) :: empty
real(8) :: bnbo, Rho0, PMFmin

   iteration = 100

   Fi(:) = 1.d0  ! Fi = exp(-Fi/kT)

   open(2,file='converge.dat')

   print *, "Temperature = ? [K]"
   read(*,*) Temp

   kT = kb*Temp

   do i = 1, nfile
     kconst(i) = kconst(i) / kT
   end do

   print *, "Histogram analysis ..."

!   open(11,file='allhist.dat')
   do i = 1, ngrid
     totn(i) = 0.d0
     do j = 1, nfile
       totn(i) = totn(i) + nstep(j)*hist(i,j)
     end do
!     write(11,*) hpos(i),totn(i)
   end do
!   close(11)

   do i = 1, ngrid
     do j = 1, nfile
       wj = 0.5d0*kconst(j)*(hpos(i)-R0(j))**2
       if(wj>300.) then
         weigh(i,j) = 0.d0
       else
         weigh(i,j) = exp(-wj)
       end if
     end do
   end do

   do k = 1, iteration

     print *, "iteration cycle =", k

     do i = 1, ngrid
       bnbo = 0.d0
       do j = 1, nfile
         bnbo = bnbo + nstep(j)*weigh(i,j)/Fi(j)
       end do
       Rho(i) = totn(i)/bnbo
     end do

     do j = 1, nfile
       Fi(j) = 0.d0
       do i = 1, ngrid
         Fi(j) = Fi(j) + weigh(i,j)*Rho(i)
       end do
     end do

     write(2,'(i5,100e13.5)') k,Fi(:)

   end do

   close(2)

! ## Free Energy

   open(3,file='FE_wham.dat')

   empty(:) = .False.
   ii = 0
   do i = 1, ngrid
     if(totn(i) < 10.) then
       empty(i) = .True.
       cycle
     end if
     ii = ii + 1
     if(ii==1) then
       Rho0 = Rho(i)
       PMF(i) = 0.d0
     else
       PMF(i) = -kT * log(Rho(i)/Rho0)
     end if
   end do

   PMFmin = 0.d0
   do i = 1, ngrid
     if(empty(i)) cycle
     if(PMF(i)<PMFmin) then
       PMFmin = PMF(i)
     end if
   end do
   do i = 1, ngrid
     if(empty(i)) cycle
     PMF(i) = PMF(i) - PMFmin
     write(3,'(f10.3,e16.8)') hpos(i),PMF(i)
   end do

   close(3)

end subroutine WHAM
