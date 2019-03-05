module checkconfig
logical :: Qbin
integer :: N,iflag,NHchain
real(8), dimension(:,:), allocatable :: R, Vel
real(8), dimension(:), allocatable :: Rss, Vss
character(len=4), dimension(:), allocatable :: AtomName, ResidName
integer, dimension(:), allocatable :: ResidNum
real(8) :: Timeps
real(8), dimension(3,3) :: H, Vg
end module checkconfig

Program CheckX

call Read_Cond
call Read_Conf

end Program CheckX

subroutine Read_Cond

use checkconfig

implicit none

character(len=1) :: Ch

   print '(a)', '##########################################'
   print '(a)', ' Program to check the file information of '
   print '(a)', ' "restart.dat", MPDyn restart file        '
   print '(a)', '##########################################'
   print '(/a/)', ' note: a file, "restart.dat", is needed in the current direcoty!'

   print '(a/)', 'File format?: Ascii or Binary [A or B]  '

   read(*,'(a)') Ch
   if(Ch=='A') then
     Qbin = .False.
     print '(a)', 'Ascii format is selected!'
   else if(Ch=='B') then
     Qbin = .True.
     print '(a)', 'Binary format is selected!'
   else
     print *, 'Incorrect input'
     stop
   end if

   print '(/a)', '##########################################'
   print '(a)', ' Choose one of the following options      '
   print '(a)', ' 1. PDB : generate PDB                    '
   print '(a)', ' 2. Number : Check the number of atoms    '
   print '(a)', ' 3. Cell : Check the simulation box       '
   print '(a/)', '##########################################'
   print '(a)', '[Type the option number]' 

   read(*,*) iflag

   if(iflag/=2) then
     print '(/a)', ' Guess the number of Nose-Hoover chain (NHC)'
     print '(a/)', ' If you used the Massive NHC, give NHC*3*N '
     read(*,*) NHChain
   end if


end subroutine Read_Cond

subroutine Read_Conf

use checkconfig

implicit none

integer :: i, Nm, ii, j
real(8) :: RmaxX, RmaxY, RmaxZ
real(8) :: RminX, RminY, RminZ
real(8) :: Lx_es, Ly_es, Lz_es
real(8), dimension(3) :: Rg

   if(Qbin) then
     open(7,file='restart.dat',status='old',form='unformatted')
   else
     open(7,file='restart.dat',status='old',form='formatted')
   end if

   if(Qbin) then
     read(7) N
   else
     read(7,'(/)')
     read(7,*) N
   end if

   print *, 'N=',N
   if(iflag==2) then
     stop
   end if

   allocate( AtomName (N) )
   allocate( ResidName(N) )
   allocate( ResidNum (N) )
   allocate( R(3,N) )
   allocate( Vel(3,N) )

   if(Qbin) then
     read(7) ResidNum
     read(7) ResidName
     read(7) AtomName
     read(7) R
     read(7) Vel
   else
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
     do i = 1 , N
       read(7,'(3d24.16)') ( Vel(j,i) , j = 1 , 3 )
     end do
   end if

   RmaxX = 0.d0
   RmaxY = 0.d0
   RmaxZ = 0.d0
   RminX = 0.d0
   RminY = 0.d0
   RminZ = 0.d0
   Rg = 0.d0
   Nm = 0
   do i = 1, N
     if((i/=1).and.(ResidNum(i)/=ResidNum(i-1))) then
       Rg = Rg / Nm
       if(Rg(1)>RmaxX) then
         RmaxX=Rg(1)
       else if(Rg(1)<RminX) then
         RminX=Rg(1)
       end if
       if(Rg(2)>RmaxY) then
         RmaxY=Rg(2)
       else if(Rg(2)<RminY) then
         RminY=Rg(2)
       end if
       if(Rg(3)>RmaxZ) then
         RmaxZ=Rg(3)
       else if(Rg(3)<RminZ) then
         RminZ=Rg(3)
       end if
       Rg = 0.d0
       Nm = 0
     end if
     Rg = Rg + R(:,i)
     Nm = Nm + 1
   end do

   Lx_es = RmaxX - RminX
   Ly_es = RmaxY - RminY
   Lz_es = RmaxZ - RminZ

   allocate( Rss(NHChain) )
   allocate( Vss(NHChain) )

   if(Qbin) then
     do i = 1, NHchain
       read(7) Rss(i), Vss(i)
     end do
     read(7) H
     read(7) Vg
     read(7) Timeps
   else
     do i = 1 , NHchain
       read(7,'(2d23.16)') Rss(i), Vss(i)
     end do
     read(7,'(3d23.16)') (H(1,i), i = 1, 3)
     read(7,'(3d23.16)') (H(2,i), i = 1, 3)
     read(7,'(3d23.16)') (H(3,i), i = 1, 3)
     read(7,'(3d23.16)') (Vg(1,i), i = 1, 3)
     read(7,'(3d23.16)') (Vg(2,i), i = 1, 3)
     read(7,'(3d23.16)') (Vg(3,i), i = 1, 3)
     read(7,*) Timeps
   end if

   Lx_es = Lx_es - H(1,1)
   Ly_es = Ly_es - H(2,2)
   Lz_es = Lz_es - H(3,3)

   if((abs(Lx_es)>10.).or.(abs(Ly_es)>10.).or.(abs(Lz_es)>10.)) then
     print '(/a/a/)', ' WARNING: NHC length might be incorrect? ',&
     &                '          Check cell matrix ! '
     write(*,'(a)')            ' #################################################### '
     write(*,'(a,3(f10.3,a))') ' CELL MATRIX = (',H(1,1),',',H(1,2),',',H(1,3),'     '
     write(*,'(a,3(f10.3,a))') '                ',H(2,1),',',H(2,2),',',H(2,3),'     '
     write(*,'(a,3(f10.3,a))') '                ',H(3,1),',',H(3,2),',',H(3,3),')    '
     write(*,'(a)')            ' #################################################### '
   else if(iflag==3) then
     write(*,'(a)')            ' #################################################### '
     write(*,'(a,3(f10.3,a))') ' CELL MATRIX = (',H(1,1),',',H(1,2),',',H(1,3),'     '
     write(*,'(a,3(f10.3,a))') '                ',H(2,1),',',H(2,2),',',H(2,3),'     '
     write(*,'(a,3(f10.3,a))') '                ',H(3,1),',',H(3,2),',',H(3,3),')    '
     write(*,'(a)')            ' #################################################### '
   end if
   if(iflag==1) then
     call GenPDB
     print *, 'Converted.pdb is successfully generated!'
   end if

end subroutine Read_Conf



!#####################################################################
!#####################################################################


! **************************
! ** Write Configration   **
! **************************

subroutine GenPDB

use checkconfig

implicit NONE

integer :: i, ir
character(len=4) :: NameA, RName
character(len=1), dimension(10), parameter :: &
&   Flpr = (/'A','B','C','D','E','F','G','H','I','J'/)
real(8), parameter :: zero=0.d0

real(8), dimension(3) :: va, vb, vc
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca, pi

   pi = atan(1.d0) * 4.d0

open(2,file='Converted.pdb')

   write(2,'(a)') &
   & 'TITLE     A PDB file converted from restart.dat         '

   write(2,'(a)')         'REMARK   2 '
   write(2,'(a,f11.3,a)') 'REMARK   2  Time = ',Timeps,' ps'
   write(2,'(a)')         'REMARK   2 '

   va = H(:,1)
   vb = H(:,2)
   vc = H(:,3)

   LLa = sqrt( dot_product(va,va) )
   LLb = sqrt( dot_product(vb,vb) )
   LLc = sqrt( dot_product(vc,vc) )
   Aab = acos( dot_product(va,vb) / LLa / LLb ) * 180. / pi
   Abc = acos( dot_product(vb,vc) / LLb / LLc ) * 180. / pi
   Aca = acos( dot_product(vc,va) / LLc / LLa ) * 180. / pi

   write(2,'(a,3f9.3,3f7.2)') 'CRYST1',LLa,LLb,LLc,Abc,Aca,Aab

   do i = 1, N

     NameA = AtomName(i)
     RName = ResidName(i)
     if(ResidNum(i)<100000) then
       ir = ResidNum(i)
     else
       ir = mod(ResidNum(i),100000)
     end if

     if(ResidName(i)=='HSD') RName='HIS'
     if(NameA(4:4)/=' ') then

       write(2,'(a4,i7,x,a4,x,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,AtomName(i),RName,ir,  &
       &  R(:,i),zero,zero,NameA(1:1)

     else

       write(2,'(a4,i7,2x,a4,a4,i5,4x,3f8.3,2f6.2,10x,1a)') &
       & 'ATOM',i,AtomName(i),RName,ir, &
       & R(:,i),zero,zero,NameA(1:1)

     end if

   end do

   write(2,'(a)') 'END'

close(2)

end subroutine GenPDB
