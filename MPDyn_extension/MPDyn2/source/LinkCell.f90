! ############################
! ## SUBROUTINE LIST 
! ## -- LinkCell 
! ## -- Mapping 
! ## -- Mapping_TopLine 
! ## -- MappingEQ 
! ############################


!#####################################################################
!#####################################################################


subroutine LinkCell

use Numbers, only : N
use CommonBlocks, only : QMacro
use CellParam, only : CellL, InvCL
use Configuration, only : R
use CellListMethod
use CGball, only : IDsphere

implicit none

integer :: icell, i, Ndd
integer, dimension(3) :: Nij
real(8), dimension(3) :: Sij

! ** zero Head of chain array **

   Head = 0

   Ndd = Ndiv(1) * Ndiv(2)

! ** sort all atoms **

   if(QMacro) then

   do i = 1, N

     if(IDsphere(i)/=0) cycle

     Sij(:) = R(:,i) * InvCL(:)
     Nij(:) = 0
     if(Sij(1) >  0.5d0) then
       Nij(1) = - 1
     else if(Sij(1) < -0.5d0) then
       Nij(1) =   1
     end if
     if(Sij(2) >  0.5d0) then
       Nij(2) = - 1
     else if(Sij(2) < -0.5d0) then
       Nij(2) =   1
     end if
     if(Sij(3) >  0.5d0) then
       Nij(3) = - 1
     else if(Sij(3) < -0.5d0) then
       Nij(3) =   1
     end if
     Sij(:) = Sij(:) + Nij(:)

     SelfShft(i) = 9 * Nij(1) + 3 * Nij(2) + Nij(3)

     icell = 1 + int( (Sij(1)+0.5d0) * Ndiv(1) )                   &
     &         + int( (Sij(2)+0.5d0) * Ndiv(2) ) * Ndiv(1)         &
     &         + int( (Sij(3)+0.5d0) * Ndiv(3) ) * Ndd

     NextP(i) = Head(icell)
     Head(icell) = i

     R(:,i) = Sij(:) * CellL(:)

   end do

   else

   do i = 1, N

     Sij(:) = R(:,i) * InvCL(:)
     Nij(:) = 0
     if(Sij(1) >  0.5d0) then
       Nij(1) = - 1
     else if(Sij(1) < -0.5d0) then
       Nij(1) =   1
     end if
     if(Sij(2) >  0.5d0) then
       Nij(2) = - 1
     else if(Sij(2) < -0.5d0) then
       Nij(2) =   1
     end if
     if(Sij(3) >  0.5d0) then
       Nij(3) = - 1
     else if(Sij(3) < -0.5d0) then
       Nij(3) =   1
     end if
     Sij(:) = Sij(:) + Nij(:)

     SelfShft(i) = 9 * Nij(1) + 3 * Nij(2) + Nij(3)

     icell = 1 + int( (Sij(1)+0.5d0) * Ndiv(1) )                   &
     &         + int( (Sij(2)+0.5d0) * Ndiv(2) ) * Ndiv(1)         &
     &         + int( (Sij(3)+0.5d0) * Ndiv(3) ) * Ndd

     NextP(i) = Head(icell)
     Head(icell) = i

     R(:,i) = Sij(:) * CellL(:)

   end do

   end if

end subroutine LinkCell


!#####################################################################
!#####################################################################


subroutine Mapping

use CellListMethod

implicit NONE

integer :: Ix, Iy, Iz, IMap

!    ** find half the nearest neighbours of each cell **

   do Iz = 1, Ndiv(3)

     do Iy = 1, Ndiv(2)-1

       do Ix = 1, Ndiv(1)

         IMap = ( CellList(Ix,Iy,Iz) - 1 ) * 16

         Map(IMap+ 1) = CellList ( Ix + 1, Iy    , Iz     )
         Map(IMap+ 2) = CellList ( Ix + 1, Iy + 1, Iz     )
         Map(IMap+ 3) = CellList ( Ix    , Iy + 1, Iz     )
         Map(IMap+ 4) = CellList ( Ix - 1, Iy + 1, Iz     )
         Map(IMap+ 5) = CellList ( Ix + 1, Iy    , Iz - 1 )
         Map(IMap+ 6) = CellList ( Ix + 1, Iy + 1, Iz - 1 )
         Map(IMap+ 7) = CellList ( Ix    , Iy + 1, Iz - 1 )
         Map(IMap+ 8) = CellList ( Ix - 1, Iy + 1, Iz - 1 )
         Map(IMap+ 9) = CellList ( Ix + 1, Iy    , Iz + 1 )
         Map(IMap+10) = CellList ( Ix + 1, Iy + 1, Iz + 1 )
         Map(IMap+11) = CellList ( Ix    , Iy + 1, Iz + 1 )
         Map(IMap+12) = CellList ( Ix - 1, Iy + 1, Iz + 1 )
         Map(IMap+13) = CellList ( Ix    , Iy    , Iz + 1 )
         Map(IMap+14) = 0
         Map(IMap+15) = 0
         Map(IMap+16) = 0

       end do

     end do

   end do


 Contains
!    ** statement function to give cell index **

   function CellList(x,y,z) Result(xx)

   integer :: x, y, z, xx

     xx = 1 + mod( x + Ndiv(1)*3 - 1 , Ndiv(1) )                   &
     &      + mod( y + Ndiv(2)*3 - 1 , Ndiv(2) ) * Ndiv(1)         &
     &      + mod( z + Ndiv(3)*3 - 1 , Ndiv(3) ) * Ndiv(1) * Ndiv(2)

   end function CellList

end subroutine Mapping


!#####################################################################
!#####################################################################


subroutine Mapping_TopLine

use CellParam, only : InvCL
use CommonDPD, only : SlideGap
use CellListMethod

implicit none

integer :: iix, ix, iy, iz, IMap

!    ** calculate x offset in cell lengths where SlideGap    **
!    ** is between 0 and +1 and box length = 1.0        **

   iix = int((SlideGap*InvCL(1)+1.d0)*dble(Ndiv(1)))

!    ** find half the nearest neighbours of each cell **

   iy = Ndiv(2)

   do iz = 1 , Ndiv(3)

     do ix = 1 , Ndiv(1)

       IMap = ( CellList(Ix,Iy,Iz) - 1 ) * 16

       Map(IMap+1 ) = CellList(ix+1    ,iy  ,iz  )
       Map(IMap+2 ) = CellList(ix+1-iix,iy+1,iz  )
       Map(IMap+3 ) = CellList(ix  -iix,iy+1,iz  )
       Map(IMap+4 ) = CellList(ix-1-iix,iy+1,iz  )
       Map(IMap+5 ) = CellList(ix+1    ,iy  ,iz-1)
       Map(IMap+6 ) = CellList(ix+1-iix,iy+1,iz-1)
       Map(IMap+7 ) = CellList(ix  -iix,iy+1,iz-1)
       Map(IMap+8 ) = CellList(ix-1-iix,iy+1,iz-1)
       Map(IMap+9 ) = CellList(ix+1    ,iy  ,iz+1)
       Map(IMap+10) = CellList(ix+1-iix,iy+1,iz+1)
       Map(IMap+11) = CellList(ix  -iix,iy+1,iz+1)
       Map(IMap+12) = CellList(ix-1-iix,iy+1,iz+1)
       Map(IMap+13) = CellList(ix      ,iy  ,iz+1)
       Map(IMap+14) = CellList(ix-2-iix,iy+1,iz  )
       Map(IMap+15) = CellList(ix-2-iix,iy+1,iz-1)
       Map(IMap+16) = CellList(ix-2-iix,iy+1,iz+1)

     end do

   end do

 Contains
!    ** statement function to give cell index **

   function CellList(x,y,z) Result(xx)

   integer :: x, y, z, xx

     xx = 1 + mod( x + Ndiv(1)*3 - 1 , Ndiv(1) )                   &
     &      + mod( y + Ndiv(2)*3 - 1 , Ndiv(2) ) * Ndiv(1)         &
     &      + mod( z + Ndiv(3)*3 - 1 , Ndiv(3) ) * Ndiv(1) * Ndiv(2)

   end function CellList

end subroutine Mapping_TopLine


!#####################################################################
!#####################################################################


subroutine MappingEQ

use CellListMethod

implicit NONE

integer :: Ix, Iy, Iz, IMap

!    ** find half the nearest neighbours of each cell **

   MapShft(:) = 0

   do Iz = 1, Ndiv(3)

     do Iy = 1, Ndiv(2)

       do Ix = 1, Ndiv(1)

         IMap = ( CellList(Ix,Iy,Iz) - 1 ) * 13

         Map(IMap+ 1) = CellList ( Ix + 1, Iy    , Iz     )
         Map(IMap+ 2) = CellList ( Ix + 1, Iy + 1, Iz     )
         Map(IMap+ 3) = CellList ( Ix    , Iy + 1, Iz     )
         Map(IMap+ 4) = CellList ( Ix - 1, Iy + 1, Iz     )
         Map(IMap+ 5) = CellList ( Ix + 1, Iy    , Iz - 1 )
         Map(IMap+ 6) = CellList ( Ix + 1, Iy + 1, Iz - 1 )
         Map(IMap+ 7) = CellList ( Ix    , Iy + 1, Iz - 1 )
         Map(IMap+ 8) = CellList ( Ix - 1, Iy + 1, Iz - 1 )
         Map(IMap+ 9) = CellList ( Ix + 1, Iy    , Iz + 1 )
         Map(IMap+10) = CellList ( Ix + 1, Iy + 1, Iz + 1 )
         Map(IMap+11) = CellList ( Ix    , Iy + 1, Iz + 1 )
         Map(IMap+12) = CellList ( Ix - 1, Iy + 1, Iz + 1 )
         Map(IMap+13) = CellList ( Ix    , Iy    , Iz + 1 )

         if(Iz == 1) then
           MapShft(IMap+ 5) = MapShft(IMap+ 5) + 1
           MapShft(IMap+ 6) = MapShft(IMap+ 6) + 1
           MapShft(IMap+ 7) = MapShft(IMap+ 7) + 1
           MapShft(IMap+ 8) = MapShft(IMap+ 8) + 1
         else if(Iz == Ndiv(3)) then
           MapShft(IMap+ 9) = MapShft(IMap+ 9) - 1
           MapShft(IMap+10) = MapShft(IMap+10) - 1
           MapShft(IMap+11) = MapShft(IMap+11) - 1
           MapShft(IMap+12) = MapShft(IMap+12) - 1
           MapShft(IMap+13) = MapShft(IMap+13) - 1
         end if

         if(Iy == Ndiv(2)) then
           MapShft(IMap+ 2) = MapShft(IMap+ 2) - 3
           MapShft(IMap+ 3) = MapShft(IMap+ 3) - 3
           MapShft(IMap+ 4) = MapShft(IMap+ 4) - 3
           MapShft(IMap+ 6) = MapShft(IMap+ 6) - 3
           MapShft(IMap+ 7) = MapShft(IMap+ 7) - 3
           MapShft(IMap+ 8) = MapShft(IMap+ 8) - 3
           MapShft(IMap+10) = MapShft(IMap+10) - 3
           MapShft(IMap+11) = MapShft(IMap+11) - 3
           MapShft(IMap+12) = MapShft(IMap+12) - 3
         end if

         if(Ix == 1) then
           MapShft(IMap+ 4) = MapShft(IMap+ 4) + 9
           MapShft(IMap+ 8) = MapShft(IMap+ 8) + 9
           MapShft(IMap+12) = MapShft(IMap+12) + 9
         else if(Ix == Ndiv(1)) then
           MapShft(IMap+ 1) = MapShft(IMap+ 1) - 9
           MapShft(IMap+ 2) = MapShft(IMap+ 2) - 9
           MapShft(IMap+ 5) = MapShft(IMap+ 5) - 9
           MapShft(IMap+ 6) = MapShft(IMap+ 6) - 9
           MapShft(IMap+ 9) = MapShft(IMap+ 9) - 9
           MapShft(IMap+10) = MapShft(IMap+10) - 9
         end if

       end do

     end do

   end do

 Contains
!    ** statement function to give cell index **

   function CellList(x,y,z) Result(xx)

   integer :: x, y, z, xx

     xx = 1 + mod( x + Ndiv(1)*3 - 1 , Ndiv(1) )                   &
     &      + mod( y + Ndiv(2)*3 - 1 , Ndiv(2) ) * Ndiv(1)         &
     &      + mod( z + Ndiv(3)*3 - 1 , Ndiv(3) ) * Ndiv(1) * Ndiv(2)

   end function CellList

end subroutine MappingEQ


!#####################################################################
!#####################################################################


subroutine MappingAll

use CellListMethod

implicit NONE

integer :: Ix, Iy, Iz, IMap

!    ** find half the nearest neighbours of each cell **

   do Iz = 1, Ndiv(3)

     do Iy = 1, Ndiv(2)

       do Ix = 1, Ndiv(1)

         IMap = ( CellList(Ix,Iy,Iz) - 1 ) * 26

         Map(IMap+ 1) = CellList ( Ix + 1, Iy    , Iz     )
         Map(IMap+ 2) = CellList ( Ix + 1, Iy + 1, Iz     )
         Map(IMap+ 3) = CellList ( Ix    , Iy + 1, Iz     )
         Map(IMap+ 4) = CellList ( Ix - 1, Iy + 1, Iz     )
         Map(IMap+ 5) = CellList ( Ix + 1, Iy    , Iz - 1 )
         Map(IMap+ 6) = CellList ( Ix + 1, Iy + 1, Iz - 1 )
         Map(IMap+ 7) = CellList ( Ix    , Iy + 1, Iz - 1 )
         Map(IMap+ 8) = CellList ( Ix - 1, Iy + 1, Iz - 1 )
         Map(IMap+ 9) = CellList ( Ix + 1, Iy    , Iz + 1 )
         Map(IMap+10) = CellList ( Ix + 1, Iy + 1, Iz + 1 )
         Map(IMap+11) = CellList ( Ix    , Iy + 1, Iz + 1 )
         Map(IMap+12) = CellList ( Ix - 1, Iy + 1, Iz + 1 )
         Map(IMap+13) = CellList ( Ix    , Iy    , Iz + 1 )
         Map(IMap+14) = CellList ( Ix - 1, Iy    , Iz     )
         Map(IMap+15) = CellList ( Ix - 1, Iy - 1, Iz     )
         Map(IMap+16) = CellList ( Ix    , Iy - 1, Iz     )
         Map(IMap+17) = CellList ( Ix + 1, Iy - 1, Iz     )
         Map(IMap+18) = CellList ( Ix - 1, Iy    , Iz + 1 )
         Map(IMap+19) = CellList ( Ix - 1, Iy - 1, Iz + 1 )
         Map(IMap+20) = CellList ( Ix    , Iy - 1, Iz + 1 )
         Map(IMap+21) = CellList ( Ix + 1, Iy - 1, Iz + 1 )
         Map(IMap+22) = CellList ( Ix - 1, Iy    , Iz - 1 )
         Map(IMap+23) = CellList ( Ix - 1, Iy - 1, Iz - 1 )
         Map(IMap+24) = CellList ( Ix    , Iy - 1, Iz - 1 )
         Map(IMap+25) = CellList ( Ix + 1, Iy - 1, Iz - 1 )
         Map(IMap+26) = CellList ( Ix    , Iy    , Iz - 1 )

       end do

     end do

   end do


 Contains
!    ** statement function to give cell index **

   function CellList(x,y,z) Result(xx)

   integer :: x, y, z, xx

     xx = 1 + mod( x + Ndiv(1)*3 - 1 , Ndiv(1) )                   &
     &      + mod( y + Ndiv(2)*3 - 1 , Ndiv(2) ) * Ndiv(1)         &
     &      + mod( z + Ndiv(3)*3 - 1 , Ndiv(3) ) * Ndiv(1) * Ndiv(2)

   end function CellList

end subroutine MappingAll
