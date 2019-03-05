! ############################
! ## SUBROUTINE LIST 
! ## -- NoLJList 
! ## -- check_bonded 
! ############################


!#####################################################################
!#####################################################################


subroutine NoLJList(iFlag)

use Numbers, only : N
use CommonBlocks, only : ForceField
use NoLJparam
use BondedParam, only : NumBond, NumAngle, BondI, BondJ, AngleI, AngleK, &
&   FTypeAngle

implicit none

integer :: i, j, k
integer :: iFlag

if(iFlag == 0) then

   allocate( NumNoLJ(N) )
   allocate( NoLJ(N,1) )

   NumNoLJ = 0
   NoLJ = 0

else if(iFlag == 1) then

   allocate( NumNoLJ(N) )
   allocate( NoLJ(N,16) )

   NumNoLJ = 0
   NoLJ = 0

   do k = 1, NumBond

     i = BondI(k)
     j = BondJ(k)

     NumNoLJ(i) = NumNoLJ(i) + 1
     NumNoLJ(j) = NumNoLJ(j) + 1

     NoLJ(i,NumNoLJ(i)) = j
     NoLJ(j,NumNoLJ(j)) = i

   end do

else if(iFlag == 2) then

   do k = 1, NumAngle

     if(ForceField(1:2) == 'CG') then
       if((FTypeAngle(k)==1).or.(FTypeAngle(k)==2)) then
         cycle
       end if
     end if

     i = AngleI(k)
     j = AngleK(k)

     NumNoLJ(i) = NumNoLJ(i) + 1
     NumNoLJ(j) = NumNoLJ(j) + 1

     NoLJ(i,NumNoLJ(i)) = j
     NoLJ(j,NumNoLJ(j)) = i

   end do

end if

end subroutine NoLJList


!######################################################################
!######################################################################


function check_bonded(i,j,non)

use NoLJparam

integer :: i, j, k, non, check_bonded

  check_bonded = 0
  do k = 1, non
    if(j == NoLJ(i,k)) then
      check_bonded = 1
      exit
    end if
  end do

end function check_bonded
