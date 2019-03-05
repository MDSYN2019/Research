

!#####################################################################
!#####################################################################


subroutine Set_CGcond

use Numbers, only : N
use CommonBlocks, only : QCellList, QMaster, QMacro
use CGdata
use CellListMethod
use CellParam, only : H, CellShape
use CommonMPI, only : NProcs
use BookParam, only : MaxPair, ListIJ, List_shortIJ

implicit none

integer :: i, j, nsmax
real(8) :: Rcut_MAX2, Rbksh

   Rcut_MAX2 = 0.d0
   do i = 1, NumAtype
     do j = i, NumAtype
       if(Rcut2(i,j) > Rcut_MAX2) then
         Rcut_MAX2 = Rcut2(i,j)
       end if
     end do
   end do
   Rcut_MAX = sqrt(Rcut_MAX2)
   InvRcut_MAX = 1.d0 / Rcut_MAX

   if(CellShape /= 3) then

     do i = 1, 3
       Ndiv(i) = int( H(i,i) * InvRcut_MAX )
     end do

     if((Ndiv(1)>=3).and.(Ndiv(2)>=3).and.(Ndiv(3)>=3)) QCellList = .True.

     if(QCellList) then

! ##  adjust the list size
       deallocate( ListIJ, List_shortIJ )
       Rbksh = sqrt(Rbk_short2)
       nsmax = int(Rbksh*Rbk_short2*3.2d-02*N/dble(NProcs))
       if(QMaster.and.(nsmax>MaxPair)) then
         print *, " Max pair number is reset to ", nsmax
         MaxPair = nsmax
       end if
       allocate( List_shortIJ(3,MaxPair) )

       Ncell= Ndiv(1) * Ndiv(2) * Ndiv(3)
       Maps = 13*Ncell

       allocate( Head(Ncell) )
       allocate( NextP(N) )
       allocate( Map(Maps) )
       allocate( MapShft(Maps) )

       call MappingEQ

     end if

   end if

   allocate( SelfShft(N) )
   SelfShft = 0


end subroutine Set_CGcond


!#####################################################################
!#####################################################################


subroutine CheckCellList

use CommonBlocks
use CGdata
use CellListMethod
use CellParam, only : H

implicit none

integer :: i
integer, dimension(3) :: Ich
logical :: QChange

   QChange = .False.

   do i = 1, 3
     Ich(i) = int( H(i,i) * InvRcut_MAX )
   end do

   if((Ich(1)<3).or.(Ich(2)<3).or.(Ich(3)<3)) then
     QCellList = .False.
     Return
   end if

   do i = 1, 3
     if( (Ich(i) < Ndiv(i)).or.(Ich(i) > (Ndiv(i)+1))) then
       QChange = .True.
     end if
   end do

   if(QChange) then
     Ndiv = Ich
     Ncell= Ndiv(1) * Ndiv(2) * Ndiv(3)
     Maps = 13*Ncell

     deallocate( Head )
     deallocate( Map  )
     deallocate( MapShft  )

     allocate( Head(Ncell) )
     allocate( Map(Maps) )
     allocate( MapShft(Maps) )

     call MappingEQ
   end if

end subroutine CheckCellList
