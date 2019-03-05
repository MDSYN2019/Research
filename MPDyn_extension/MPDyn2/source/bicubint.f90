subroutine bicubprep(y,y1,y2,y12,d1,d2,c)

implicit none

real(8) :: d1,d2
real(8), dimension(4,4) :: c
real(8), dimension(4) :: y, y1, y12, y2

integer :: i,j,k,l

real(8) :: d1d2,xx,cl(16),wt(16,16),x(16)

save wt
data wt/1,0,-3,2,0,0,0,0,-3,0,9,-6,2,0,-6,4,0,0,0,&
& 0,0,0,0,0,3,0,-9,6,-2,0,6,-4,0,0,0,0,0,0,0,0,0,0,&
& 9,-6,0,0,-6,4,0,0,3,-2,0,0,0,0,0,0,-9,6,0,0,6,-4,&
& 0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,0,0,0,0,0,0,&
& 0,0,-1,0,3,-2,1,0,-3,2,0,0,0,0,0,0,0,0,0,0,&
& -3,2,0,0,3,-2,0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2,&
& 0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2,0,0,0,0,0,0,&
& 0,0,0,3,-6,3,0,-2,4,-2,0,0,0,0,0,0,0,0,0,0,&
& -3,3,0,0,2,-2,0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2,&
& 0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1,0,0,0,0,0,0,&
& 0,0,0,-1,2,-1,0,1,-2,1,0,0,0,0,0,0,0,0,0,0,&
& 1,-1,0,0,-1,1,0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1/

   d1d2=d1*d2

   do i = 1, 4
     x(i)    = y(i)
     x(i+4)  = y1(i)  * d1
     x(i+8)  = y2(i)  * d2
     x(i+12) = y12(i) * d1d2
   end do

   do i = 1, 16
     xx = 0.d0
     do k = 1, 16
       xx = xx + wt(i,k) * x(k)
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




subroutine bicubint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)

implicit none

real(8), dimension(4), intent(in) :: y,y1,y2,y12
real(8), intent(in) :: x1l,x1u,x2l,x2u,x1,x2
real(8), intent(out) :: ansy,ansy1,ansy2
integer :: i
real(8) :: t,u
real(8), dimension(4,4) :: c

   call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

   if (x1u == x1l .or. x2u == x2l) then
     write(*,*) 'bicubic interpolation is not solved correctly.'
     stop
   end if

   t = (x1-x1l) / (x1u-x1l)
   u = (x2-x2l) / (x2u-x2l)

   ansy  = 0.d0
   ansy2 = 0.d0
   ansy1 = 0.d0

   do i = 4, 1, -1
     ansy  = t * ansy  + ((c(i,4) * u + c(i,3)) * u + c(i,2)) * u + c(i,1)
     ansy2 = t * ansy2 + (3.d0 * c(i,4) * u + 2.d0 * c(i,3)) * u + c(i,2)
     ansy1 = u * ansy1 + (3.d0 * c(4,i) * t + 2.d0 * c(3,i)) * t + c(2,i)
   end do

   ansy1 = ansy1 / (x1u-x1l)
   ansy2 = ansy2 / (x2u-x2l)

end subroutine bicubint
