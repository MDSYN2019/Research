module convert_init
integer :: N, NumSpec
character(len=50) :: PDB_file
character(len=6), dimension(:), allocatable :: MolName
integer, dimension(:), allocatable :: NumMol, NumAtm
character(len=4), dimension(:), allocatable :: AtomName, ResidName
integer, dimension(:), allocatable :: ResidNum
real(8), dimension(:,:), allocatable :: R
real(8), dimension(3,3) :: H, InvH
end module convert_init

program Gen_init_fm_NAMD

logical :: QCell

call ReadInitCond
call Read_PDB(QCell)
call Write_CRD(QCell)

end program Gen_init_fm_NAMD

subroutine ReadInitCond

use convert_init

implicit none

character(len=80) :: String1, String
integer :: i, ii, j, jj
logical, dimension(80) :: empty
integer,dimension(10) :: first, last

   print *, '**************************************'
   print *, '     Convert : PDB --> CRD file       '
   print *, '**************************************'
   print *, ''

   print *, 'Name of PDB file ?'
   read(5,'(a80)') String1
   String = adjustl(String1)
   read(String,*) PDB_file

   print *, 'number of molecular species ?'
   read(5,*) NumSpec

   allocate( MolName(NumSpec) )
   allocate( NumMol(NumSpec) )
   allocate( NumAtm(NumSpec) )

   print *, 'name of molecules ? write all in one line with space inbetween'
   read(5,'(a80)') String1
   String = trim(adjustl(String1))

   do i = 1, 80
     if(String(i:i) == ' ') then
       empty(i) = .True.
     else
       empty(i) = .False.
     end if
   end do

   ii = 0
   j  = 1
   first(j) = 1
   do i = 2, 80
     if(empty(i).and.(.not.empty(i-1))) then
       ii = ii + 1
       last(ii) = i-1
     end if
     if(empty(i-1).and.(.not.empty(i))) then
       j = j + 1
       first(j) = i
     end if
     if(ii==NumSpec) then
       jj = i
       exit
     end if
   end do

   do i = 1, NumSpec
     MolName(i) = String(first(i):last(i))
   end do

   print *, 'number of molecules ? write all in one line.'
   read(5,'(a80)') String1
   String = trim(adjustl(String1))

   do i = 1, 80
     if(String(i:i) == ' ') then
       empty(i) = .True.
     else
       empty(i) = .False.
     end if
   end do

   ii = 0
   j  = 1
   first(j) = 1
   do i = 2, 80
     if(empty(i).and.(.not.empty(i-1))) then
       ii = ii + 1
       last(ii) = i-1
     end if
     if(empty(i-1).and.(.not.empty(i))) then
       j = j + 1
       first(j) = i
     end if
     if(ii==NumSpec) then
       jj = i
       exit
     end if
   end do

   do i = 1, NumSpec
     read(String(first(i):last(i)),*) NumMol(i)
   end do


   print *, 'number of atoms per molecule ? write all in one line.'
   read(5,'(a80)') String1
   String = trim(adjustl(String1))

   do i = 1, 80
     if(String(i:i) == ' ') then
       empty(i) = .True.
     else
       empty(i) = .False.
     end if
   end do

   ii = 0
   j  = 1
   first(j) = 1
   do i = 2, 80
     if(empty(i).and.(.not.empty(i-1))) then
       ii = ii + 1
       last(ii) = i-1
     end if
     if(empty(i-1).and.(.not.empty(i))) then
       j = j + 1
       first(j) = i
     end if
     if(ii==NumSpec) then
       jj = i
       exit
     end if
   end do

   do i = 1, NumSpec
    read(String(first(i):last(i)),*) NumAtm(i)
   end do

   N = 0
   do i = 1, NumSpec
     N = N + NumMol(i)*NumAtm(i)
   end do

end subroutine ReadInitCond


!######################################################################
!######################################################################


! **************************
! **  Read Configration   **
! **************************

subroutine Read_PDB(QCell)

use convert_init

implicit none

integer :: i, ii
character(len=80) :: String
real(8) :: LLa, LLb, LLc, Aab, Abc, Aca, pi
logical :: QCell
integer :: eofile

   pi = atan(1.d0)*4.d0
   QCell = .False.

   open( 7, file=trim(PDB_file),status='old')

   ii = 0

   allocate( AtomName (N) )
   allocate( ResidName(N) )
   allocate( ResidNum (N) )
   allocate( R(3,N) )

   print *, 'N=', N

   do

     read(7,'(a80)',iostat=eofile) String

     if( eofile == -1 ) exit
     if( String(1:3) == 'END' ) exit

     if((String(1:6) == 'HETATM') .or. (String(1:4) == 'ATOM')) then

       ii = ii + 1

       read(String,'(12x,a4,x,a4,x,i4,4x,3f8.3)') &
       &    AtomName(ii), ResidName(ii), ResidNum(ii), R(:,ii)

     end if

     if(String(1:6)=='CRYST1') then
       read(String,'(6x,3f9.3,3f7.2)') LLa,LLb,LLc,Abc,Aca,Aab
       QCell = .True.
     end if

   end do

   close(7)

   if(QCell) then

     H = 0.

     if((Aab/=90.).or.(Abc/=90.).or.(Aca/=90.)) then

       Aab = pi * Aab / 180.d0
       Abc = pi * Abc / 180.d0
       Aca = pi * Aca / 180.d0

       H(1,1) = LLa
       H(1,2) = LLb * cos(Aab)
       H(2,2) = LLb * sin(Aab)
       H(1,3) = LLc * cos(Aca)
       H(2,3) = LLc * cos(Abc)
       H(3,3) = LLc * sin(Aca) * sin(Abc)

     else

       H(1,1) = LLa
       H(2,2) = LLb
       H(3,3) = LLc

     end if

     call InversMatrix(H,InvH)

   else

     print *, "CRYST line is needed to spcify the cell matrix"

   end if

end subroutine Read_PDB


!######################################################################
!######################################################################


! ***************************
! **  Write Configration   **
! ***************************

subroutine Write_CRD(QCell)

use convert_init, only : N, R, H, InvH, AtomName, ResidName, ResidNum

implicit none

integer :: i
character(len=4) :: ModelName
logical :: QCell

     ModelName = ResidName(1)

     open(7,file='initial.crd')

     write(7,'(a)') '! CRD file'
     write(7,'(a)') '!'
     write(7,'(i10)') N

     if(N>100000) then
       do i = 1, N
         write(7,'(2i9,2(x,a4),3f12.5)')  &
         & i  , ResidNum(i) , adjustl(ResidName(i)) ,  &
         & adjustl(AtomName(i)) , R(:,i)
       end do
     else
       do i = 1, N
         write(7,'(2i5,2(x,a4),3f10.5)')  &
         & i  , ResidNum(i) , adjustl(ResidName(i)) ,  &
         & adjustl(AtomName(i)) , R(:,i)
       end do
     end if

     if(QCell) then

       write(7,'(3f12.4)') (H(1,i), i = 1, 3)
       write(7,'(3f12.4)') (H(2,i), i = 1, 3)
       write(7,'(3f12.4)') (H(3,i), i = 1, 3)

     end if

     close(7)

end subroutine Write_CRD


!######################################################################
!######################################################################


! ************************
! ** A  -->  B = A^(-1) **
! ************************

subroutine InversMatrix(A,B)

implicit none

real(8), dimension(3,3) :: A , B
real(8) :: x, x1, det
external det

   x = det(A)

   x1 = 1.d0 / x

   B(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) * x1
   B(1,2) = ( A(3,2) * A(1,3) - A(3,3) * A(1,2) ) * x1
   B(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * x1
   B(2,1) = ( A(2,3) * A(3,1) - A(2,1) * A(3,3) ) * x1
   B(2,2) = ( A(3,3) * A(1,1) - A(3,1) * A(1,3) ) * x1
   B(2,3) = ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * x1
   B(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) * x1
   B(3,2) = ( A(3,1) * A(1,2) - A(3,2) * A(1,1) ) * x1
   B(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * x1

end subroutine InversMatrix


!######################################################################
!######################################################################


Function det(A)

implicit none

real(8), dimension(3,3) :: A
real(8) :: det

   det = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1)  &
   &   + A(1,3) * A(2,1) * A(3,2) - A(1,3) * A(2,2) * A(3,1)  &
   &   - A(1,1) * A(2,3) * A(3,2) - A(1,2) * A(2,1) * A(3,3)

end Function det
