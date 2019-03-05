! ############################
! ## SUBROUTINE LIST 
! ## -- InversMatrix 
! ## -- det 
! ## -- Jacobi 
! ## -- TransCellList 
! ## -- VecProduct 
! ## -- MatMat 
! ## -- Randp 
! ## -- genran 
! ## -- Gauss 
! ## -- AcTime 
! ## -- JNAME 
! ## -- RanSeed 
! ## -- PreRandomNumber 
! ## -- CIdirection 
! ## -- Error_Function 
! ## -- InterPolate 
! ## -- SwitchFunc 
! ## -- Finalize 
! ############################


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


!######################################################################
!######################################################################


! *********************************************
! **  Eigen Value Solution by Jacobi Method  **
! *********************************************

subroutine Jacobi(aaa,v,d)

implicit none

! a : matrix to be diagonalized
! v : eigen vectors / rotation matrix as a result
! d : eigen values
! irot : number of rotational operations

integer :: ip,iq,irot,i,j
real(8), dimension(3,3) :: aaa,a,v
real(8), dimension(3)   :: d,b,z
real(8) :: sm, tresh, g, h, t, s, c, tau, theta

   a = aaa

   v = 0.d0

   do ip = 1 , 3

     v(ip,ip)=1.d0

   end do

   do ip = 1 , 3

     d(ip) = a(ip,ip)
     b(ip) = d(ip)
     z(ip) = 0.d0

   end do

   irot=0

   do i = 1 , 50

     sm = 0.d0

     do ip = 1 , 2

       do iq = ip + 1 , 3

         sm = sm + abs( a(ip,iq) )

       end do

     end do

     if(sm==0.d0) return

     if(i < 4) then

         tresh = 0.2d0 * ( sm / ( 3 * 3 ) )

     else

         tresh=0.d0

     end if

     do ip = 1 , 2

       do iq = ip + 1 , 3

         g = 1.d2 * abs( a(ip,iq) )

         if((i > 4).and.( (abs(d(ip))+g) == abs(d(ip)) )&
         &         .and.( (abs(d(iq))+g) == abs(d(iq)) )) then

               a(ip,iq) = 0.d0

         else if(abs(a(ip,iq)) > tresh) then

           h = d(iq) - d(ip)

           if( (abs(h)+g) == abs(h) ) then

             t = a(ip,iq) / h

           else

             theta = 0.5d0 * h / a(ip,iq)
             t = 1 / ( abs(theta) + sqrt( 1 + theta * theta ) )

             if(theta < 0.d0) t=-t

           end if

           c = 1.d0 / sqrt( 1.d0 + t * t )

           s     = t * c
           tau   = s / ( 1.d0 + c )
           h     = t * a(ip,iq)
           z(ip) = z(ip) - h
           z(iq) = z(iq) + h
           d(ip) = d(ip) - h
           d(iq) = d(iq) + h
           a(ip,iq) = 0.d0

           do j = 1 , ip - 1

             g       = a(j,ip)
             h       = a(j,iq)
             a(j,ip) = g - s * ( h + g * tau )
             a(j,iq) = h + s * ( g - h * tau )

           end do

           do j = ip+1 , iq-1

             g       = a(ip,j)
             h       = a(j,iq)
             a(ip,j) = g - s * ( h + g * tau )
             a(j,iq) = h + s * ( g - h * tau )

           end do

           do j = iq+1 , 3

             g       = a(ip,j)
             h       = a(iq,j)
             a(ip,j) = g - s * ( h + g * tau )
             a(iq,j) = h + s * ( g - h * tau )

           end do

           do j = 1 , 3

             g       = v(j,ip)
             h       = v(j,iq)
             v(j,ip) = g - s * ( h + g * tau )
             v(j,iq) = h + s * ( g - h * tau )

           end do

           irot=irot+1

         end if

       end do

     end do

     do ip = 1 , 3

       b(ip) = b(ip) + z(ip)
       d(ip) = b(ip)
       z(ip) = 0.0

     end do

   end do

end subroutine Jacobi


!######################################################################
!######################################################################


! ***************************************************
! ** Make the Cell List adjacent to the basic cell **
! ***************************************************

subroutine TransCellList

use CellParam, only : H, CellShft

implicit none

integer :: ix, iy, iz, ii

   do ix = -1, 1
     do iy = -1, 1
       do iz = -1, 1
         ii = 9 * ix + 3 * iy + iz
         CellShft(1,ii) = H(1,1) * ix + H(1,2) * iy + H(1,3) * iz
         CellShft(2,ii) = H(2,1) * ix + H(2,2) * iy + H(2,3) * iz
         CellShft(3,ii) = H(3,1) * ix + H(3,2) * iy + H(3,3) * iz
       end do
     end do
   end do

end subroutine TransCellList


!######################################################################
!######################################################################


subroutine VecProduct(x,y,z)

implicit none

real(8), dimension(3) :: x, y, z
real(8) :: dz

   z(1) = x(2) * y(3) - y(2) * x(3)
   z(2) = x(3) * y(1) - y(3) * x(1)
   z(3) = x(1) * y(2) - y(1) * x(2)

   dz   = sqrt( dot_product( z , z ) )

   z    = z / dz

end subroutine VecProduct


!######################################################################
!######################################################################


subroutine MatMat(X, Y, Z)

implicit none

real(8), dimension(3,3) :: X, Y, Z

   Z(1,1)=X(1,1)*Y(1,1)+X(1,2)*Y(2,1)+X(1,3)*Y(3,1)
   Z(1,2)=X(1,1)*Y(1,2)+X(1,2)*Y(2,2)+X(1,3)*Y(3,2)
   Z(1,3)=X(1,1)*Y(1,3)+X(1,2)*Y(2,3)+X(1,3)*Y(3,3)
   Z(2,1)=X(2,1)*Y(1,1)+X(2,2)*Y(2,1)+X(2,3)*Y(3,1)
   Z(2,2)=X(2,1)*Y(1,2)+X(2,2)*Y(2,2)+X(2,3)*Y(3,2)
   Z(2,3)=X(2,1)*Y(1,3)+X(2,2)*Y(2,3)+X(2,3)*Y(3,3)
   Z(3,1)=X(3,1)*Y(1,1)+X(3,2)*Y(2,1)+X(3,3)*Y(3,1)
   Z(3,2)=X(3,1)*Y(1,2)+X(3,2)*Y(2,2)+X(3,3)*Y(3,2)
   Z(3,3)=X(3,1)*Y(1,3)+X(3,2)*Y(2,3)+X(3,3)*Y(3,3)

end subroutine MatMat


!#####################################################################
!#####################################################################


subroutine Randp(x,y,z)

implicit none

real(8) :: x, y, z, ranf
external ranf

! random number -1 < r < 1

   x = 2.d0 * ranf() - 1.d0
   y = 2.d0 * ranf() - 1.d0
   z = 2.d0 * ranf() - 1.d0

end subroutine Randp


!#####################################################################
!#####################################################################


subroutine genran(seed)
!*
implicit integer(a-z)
!*
!* Period parameters
parameter(N     =  624)
!*
dimension mt(0:N-1)
!*                     the array for the state vector
common /block/mti,mt
save   /block/
!*
!*      setting initial seeds to mt[N] using
!*      the generator Line 25 of Table 1 in
!*      [KNUTH 1981, The Art of Computer Programming
!*         Vol. 2 (2nd Ed.), pp102]
!*
   mt(0)= iand(seed,-1)
   do mti=1,N-1
     mt(mti) = iand(69069 * mt(mti-1),-1)
   end do
!*
end subroutine genran


!#####################################################################
!#####################################################################


function ranf()
!*
implicit integer(a-z)

real(8) :: ranf
!*
!* Period parameters
parameter(N     =  624)
parameter(N1    =  N+1)
parameter(M     =  397)
parameter(MATA  = -1727483681)
!*                                    constant vector a
parameter(UMASK = -2147483648)
!*                                    most significant w-r bits
parameter(LMASK =  2147483647)
!*                                    least significant r bits
!* Tempering parameters
parameter(TMASKB= -1658038656)
parameter(TMASKC= -272236544)
!*
dimension mt(0:N-1)
!*                     the array for the state vector
common /block/mti,mt
save   /block/
data   mti/N1/
!*                     mti==N+1 means mt[N] is not initialized
!*
dimension mag01(0:1)
data mag01/0, MATA/
save mag01
!*                        mag01(x) = x * MATA for x=0,1
!*
   TSHFTU(y)=ishft(y,-11)
   TSHFTS(y)=ishft(y,7)
   TSHFTT(y)=ishft(y,15)
   TSHFTL(y)=ishft(y,-18)
!*
   if(mti.ge.N) then

   do kk = 0 , N-M-1

     y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
     mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))

   end do

   do kk = N-M , N-2

     y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
     mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))

   end do

   y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
   mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
   mti = 0

   endif
!*
   y=mt(mti)
   mti=mti+1
   y=ieor(y,TSHFTU(y))
   y=ieor(y,iand(TSHFTS(y),TMASKB))
   y=ieor(y,iand(TSHFTT(y),TMASKC))
   y=ieor(y,TSHFTL(y))
!*
   if(y.lt.0) then

     ranf=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)

   else

     ranf=dble(y)/(2.0d0**32-1.0d0)

   endif
!*
end function ranf


!********************************************************************
!********************************************************************


!!function Gauss()
!!
!!implicit none
!!
!!real(8) :: Gauss
!!
!!real(8), parameter :: a1=3.949846138d0
!!real(8), parameter :: a3=0.252408784d0
!!real(8), parameter :: a5=0.076542912d0
!!real(8), parameter :: a7=0.008355968d0
!!real(8), parameter :: a9=0.029899776d0
!!real(8)            :: suma,r,r2,ranf
!!integer            :: i
!!external ranf
!!
!!   suma=0.0d0
!!
!!   do i = 1 , 12
!!
!!     suma = suma + ranf()
!!
!!   end do
!!
!!   r  = ( suma - 6.d0 ) / 4.d0
!!   r2 = r * r
!!
!!   Gauss = ( ( ( ( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1  ) * r
!!
!!end function Gauss


!********************************************************************
!********************************************************************


function Gauss()

! *************************************
! *    Generate normal distribution
! *    Numerical recipes.
! *************************************

integer :: Iset
real(8) :: fac, Gset, Rsq, v1, v2, ranf
real(8) :: Gauss
save Iset, Gset
external ranf

data Iset/0/

   if (Iset == 0) then

1    v1 = 2.d0 * ranf() - 1.d0
     v2 = 2.d0 * ranf() - 1.d0
     Rsq = v1*v1 + v2*v2

     if( ( Rsq >= 1.d0 ) .or. ( Rsq == 0.d0 ) ) goto 1

     fac   = sqrt( -2.d0 * log(Rsq) / Rsq )
     Gset  = v1 * fac
     Gauss = v2 * fac
     Iset  = 1

   else

     Gauss = Gset
     Iset=0

   endif

end function Gauss


!#####################################################################
!#####################################################################


subroutine AcTime(ptime)

character(len=17) :: ptime
character(len=2) :: Mon, Date, Hour, Min, Sec
character(len=4) :: Year
integer, dimension(8) :: val

   val = 0

   call date_and_time(VALUES=val)

   write(Year,'(i4)'  ) val(1)
   write(Mon ,'(i2.2)') val(2)
   write(Date,'(i2.2)') val(3)
   write(Hour,'(i2.2)') val(5)
   write(Min ,'(i2.2)') val(6)
   write(Sec ,'(i2.2)') val(7)

   write(ptime,'(5(a2,a1),a2)') Mon,'/',Date,'/',Year(3:4),' ',Hour,':',Min,':',Sec

end subroutine AcTime


!#####################################################################
!#####################################################################


subroutine JNAME(ptime)

character(len=17) :: ptime
character(len=2) :: Mon, Date, Hour, Min
character(len=4) :: Year
integer, dimension(8) :: val

   val = 0

   call date_and_time(VALUES=val)

   write(Year,'(i4)'  ) val(1)
   write(Mon ,'(i2.2)') val(2)
   write(Date,'(i2.2)') val(3)
   write(Hour,'(i2.2)') val(5)
   write(Min ,'(i2.2)') val(6)

   write(ptime,'(2a4,2a2,a1,2a2)') 'JOB_',Year,Mon,Date,'_',Hour,Min

end subroutine JNAME


!#####################################################################
!#####################################################################


subroutine RanSeed(Seed)

use CommonMPI

implicit none

integer :: i, j, Seed
integer, dimension(8) :: val

   val = 0

   call date_and_time(VALUES=val)

   i = val(2)*val(3) + val(5)*val(6) + val(6)*val(7) + val(5)*val(7)
   j = val(1) + val(5) + val(6) + val(7)

   Seed = i + j + MyRank

end subroutine RanSeed


!#####################################################################
!#####################################################################


subroutine PreRandomNumber

use CommonBlocks, only : Qstdout
use CommonMPI

implicit none

!integer :: i, ix
integer :: Seed
!integer, dimension(1000) :: ls
!real(8) :: xi

   call RanSeed(Seed)

   if(Qstdout) write(*,*) 'Seed = ', Seed, ': CPU = ',MyRank

! ----------------------------
! # random number generator
! ----------------------------
   call genran(Seed)

!  open(61,file='ranf.test')
!
!  do i = 1 , 10000000
!
!    xi = ranf()
!    ix = int(xi*1000.)+1
!    if(ix>1000) print *,'error'
!    ls(ix) = ls(ix) + 1
!
!  end do
!
!  do i = 1, 1000
!
!    xi = ls(i) / 10000000. * 1000.
!    write(61,'(2f7.3)') real(i)/1000.,xi
!
!  end do
!
!  close(61)
! -------------------------------------------


end subroutine PreRandomNumber


!######################################################################
!######################################################################


Function CIdirection(Ch)

implicit none

character(len=1) :: Ch
integer :: CIdirection

   if((Ch == 'X').or.(Ch == 'x')) then
     CIdirection = 1
   else if((Ch == 'Y').or.(Ch == 'y')) then
     CIdirection = 2
   else if((Ch == 'Z').or.(Ch == 'z')) then
     CIdirection = 3
   else
     write(*,*) 'ERROR : DIRECTION'
   end if

end Function CIdirection


!######################################################################
!######################################################################


Function Error_Function(x) Result(rs)

use EwaldParam, only: EFList, msh

implicit none

integer :: ii
real(8) :: x, xx, dx, y1, y2, y3, z1, z2
real(8) :: rs

#ifdef CUBIC
real(8) :: y4, z3
real(8), parameter :: sb = 0.3333333333333333d0
real(8), parameter :: sb2 = 0.6666666666666667d0

   xx = x * msh

   ii  = int( xx )

   dx  = xx - dble(ii)

   y1  = EFList( ii-1 )
   y2  = EFList( ii   )
   y3  = EFList( ii+1 )
   y4  = EFList( ii+2 )

   z3  = y2 - y3 + sb * (y4 - y1)
   z2  = y1 - 2.d0 * y2 + y3
   z1  = -sb2 * y1 - y2 + 2.d0 * y3 - sb * y4

   rs = 0.5d0 * dx * ( z1 + dx * ( z2 + dx * z3 ) ) + y2

#else

   xx = x * msh
   ii = nint( xx )

   dx  = xx - dble(ii)

   y1  = EFList( ii-1 )
   y2  = EFList( ii   )
   y3  = EFList( ii+1 )

   z2  = y1 - 2.d0 * y2 + y3
   z1  = - y1 + y3
   rs  = 0.5d0 * (z1 + z2 * dx ) * dx + y2

#endif

end Function Error_Function


!######################################################################
!######################################################################


Function InterPolate(x,drmin,invdr,j,k) Result(rs)

use TableFuncs

implicit none

integer :: ii, j, k
real(8) :: x, xx, dx, y1, y2, y3
real(8) :: z1, z2
real(8) :: rs, drmin, invdr

#ifdef CUBIC
real(8) :: y4, z3
real(8), parameter :: sb = 0.3333333333333333d0
real(8), parameter :: sb2 = 0.6666666666666667d0

   xx  = (x - drmin) * invdr

   ii  = int( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = TabFunc( ii  , j, k )
   y2  = TabFunc( ii+1, j, k )
   y3  = TabFunc( ii+2, j, k )
   y4  = TabFunc( ii+3, j, k )

   z3  = y2 - y3 + sb * (y4 - y1)
   z2  = y1 - 2.d0 * y2 + y3
   z1  = -sb2 * y1 - y2 + 2.d0 * y3 - sb * y4

   rs = 0.5d0 * dx * ( z1 + dx * ( z2 + dx * z3 ) ) + y2

   else

   rs = TabFunc( 1, j, k )
   print *, 'Bad contact'

   end if

#else

   xx  = (x - drmin) * invdr

   ii  = nint( xx )

   if(ii >= 1) then

   dx  = xx - dble(ii)

   y1  = TabFunc( ii  , j, k )
   y2  = TabFunc( ii+1, j, k )
   y3  = TabFunc( ii+2, j, k )

   z2  = y1 - 2.d0 * y2 + y3
   z1  = - y1 + y3
   rs  = 0.5d0 * (z1 + z2 * dx ) * dx + y2

   else

   rs = TabFunc( 1, j, k )
   print *, 'Bad contact'

   end if

#endif

end Function InterPolate


!######################################################################
!######################################################################


subroutine bicubic(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy)!,ansy1,ansy2)

use CommonBlocks, only : QMaster

implicit none

real(8), dimension(4), intent(in) :: y,y1,y2,y12
real(8), intent(in) :: x1l,x1u,x2l,x2u,x1,x2
real(8), intent(out) :: ansy!,ansy1,ansy2
integer :: i
real(8) :: t,u
real(8), dimension(4,4) :: c

   call bicubprep(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

   if (x1u == x1l .or. x2u == x2l) then
     if(QMaster) write(*,*) 'bicubic interpolation is not solved correctly.'
     call Finalize
   end if

   t = (x1-x1l) / (x1u-x1l)
   u = (x2-x2l) / (x2u-x2l)

   ansy  = 0.d0
!   ansy2 = 0.d0
!   ansy1 = 0.d0

   do i = 4, 1, -1
     ansy  = t * ansy  + ((c(i,4) * u + c(i,3)) * u + c(i,2)) * u + c(i,1)
!     ansy2 = t * ansy2 + (3.d0 * c(i,4) * u + 2.d0 * c(i,3)) * u + c(i,2)
!     ansy1 = u * ansy1 + (3.d0 * c(4,i) * t + 2.d0 * c(3,i)) * t + c(2,i)
   end do

!   ansy1 = ansy1 / (x1u-x1l)
!   ansy2 = ansy2 / (x2u-x2l)

end subroutine bicubic


! #######################################################################
! #######################################################################


subroutine bicubprep(y,y1,y2,y12,d1,d2,c)

implicit none

real(8) :: d1,d2
real(8), dimension(4,4) :: c
real(8), dimension(4) :: y, y1, y12, y2

integer :: i,j,l

real(8) :: d12,xx,cl(16),wt(16,16),x(16)

save wt
data wt/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,&
       & 0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,&
       & 0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,&
       & 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,&
       & 0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,&
       & 0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,&
       & 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,&
       & 0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,&
       & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1,&
       & 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1/

   d12 = d1*d2

   do i = 1, 4
     x(i)    = y(i)
     x(i+4)  = y1(i)  * d1
     x(i+8)  = y2(i)  * d2
     x(i+12) = y12(i) * d12
   end do

   do i = 1, 16
     xx = 0.d0
     do j = 1, 16
       xx = xx + wt(i,j) * x(j)
     end do
     cl(i) = xx
   end do

   l=0
   do i = 1, 4
     do j = 1, 4
       l = l + 1
       c(i,j) = cl(l)
     end do
   end do

end subroutine bicubprep


! #######################################################################
! #######################################################################


subroutine locate(xx,n,x,j)

implicit none

integer :: j,n
real(8) :: x
real(8), dimension(n) :: xx
integer :: jl, jm, ju

   jl = 0
   ju = n+1

   do while (ju-jl>1)
     jm=(ju+jl)/2
     if(x>=xx(jm)) then
       jl = jm
     else
       ju = jm
     end if
   end do

   if(x==xx(1)) then
     j = 1
   else if(x==xx(n)) then
     j = n-1
   else
     j = jl
   end if

end subroutine locate


!######################################################################
!######################################################################


subroutine SwitchFunc(R2,fk,ek,Ron2,Rcutoff2,Swpref)

implicit none

real(8) :: R2, fk, ek
real(8) :: Xoff, Xon
real(8) :: Ron2, Rcutoff2, Swpref
real(8) :: Switch, Dswitch

   if(R2 > Ron2) then

     Xoff = R2 - Rcutoff2      ! x-xoff
     Xon  = R2 - Ron2          ! x-xon

     Switch  = ( 3.d0 * Xon - Xoff ) * Xoff * Xoff * Swpref
     Dswitch = 12.d0 * Swpref * Xoff * Xon

     fk = fk * Switch + Dswitch * ek
     ek = ek * Switch

   end if

end subroutine SwitchFunc


!######################################################################
!######################################################################


Function Fswitch(R2)

use CGdata, only: Fsw1, Rcut_short2

implicit none

real(8) :: R2, yr, y3
real(8) :: Fswitch

   yr = (R2 - Rcut_short2) * Fsw1
   y3 = yr * yr * (2.d0 * yr - 3.d0)
   Fswitch = 1.d0 + y3

end Function Fswitch


!######################################################################
!######################################################################


subroutine Finalize

implicit none

   call FinMPI
   stop

end subroutine Finalize
