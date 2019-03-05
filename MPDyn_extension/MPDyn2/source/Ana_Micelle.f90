! ############################
! ## SUBROUTINE LIST 
! ## -- IDmicelle
! ############################

subroutine IDmicelle

use Configuration, only : R
use Numbers, only : NumMol, NumAtm
use AtomParam, only : AtomName
use CellParam, only : H, InvH, CellShft
use ParamAnalyze
use CutoffParam, only : Rcutoff2

implicit none

integer :: i, j, k, ii, kk, l, ini, i1, ll, k1
integer :: minrot, jrot, isz
character(len=20) :: Script, Script1
character(len=4), dimension(:), allocatable :: TmpAName
real(8), dimension(:,:,:), allocatable :: ScR, Rlp
integer, dimension(:), allocatable :: Ncatom, Nattrib
integer, dimension(:), allocatable :: Nneighbor
integer, dimension(:,:), allocatable :: Nbooked
integer :: Isolate, Ncarbon
real(8) :: Six, Siy, Siz, Rix, Riy, Riz
real(8) :: Sx, Sy, Sz, Rx, Ry, Rz, R2
integer :: Nx, Ny, Nz, NNN, NAlip, NMlip
integer :: iadj, jatr, minatr
integer :: eofile, iatr, ncls, ncls2, maxsz, minsz, isize
integer, dimension(:), allocatable :: nsize, ndist
real(8) :: sum1, sum2, avsz1, avsz2
integer :: totalstep
character(len=50) :: Filename


   open(12,file='DefineMicelle.data')

   allocate(TmpAName(100))

   ii = 0
   do
     read(12,'(a)',iostat=eofile) Script1
     if(eofile == -1) exit
     Script = trim(adjustl(Script1))

     if(Script(1:1)=='#'.or.Script(1:1)=='!') cycle

     ii = ii + 1
     read(Script,*) TmpAName(ii)

   end do

   close(12)

   Ncarbon = ii
   allocate(Ncatom(Ncarbon))

   NMlip = NumMol(Kcomp)
   NAlip = NumAtm(Kcomp)

   ini = 0
   do i = 1, Kcomp-1
     ini = ini + NumMol(k)*NumAtm(k)
   end do

   do i = 1, Ncarbon
     ii = ini
n00: do j = 1, NAlip
       ii = ii + 1
       if(AtomName(ii)==TmpAName(i)) then
         Ncatom(i) = j
         exit n00
       end if
       if(j==NAlip) then
         write(*,*) 'error: no corresponding atom was found'
         write(*,*) TmpAName(i)
         call Finalize
       end if
     end do n00
   end do

   totalstep = 0


   allocate(Nneighbor(NMlip))
   allocate(Nbooked(NMlip,NMlip))

   allocate(ScR(3,Ncarbon,NMlip))
   allocate(Rlp(3,Ncarbon,NMlip))

   allocate(nsize(NMlip))
   allocate(ndist(NMlip))
   allocate(Nattrib(NMlip))

   open(13,file='./Analy/ClustSize.dat')
   write(13,'(a)') "#-step---*--Ncls---*-Ncls>2--*--AveN---*--Ave^2-*--Nmax---*--Nmin---*"

   do i = 1 , NJobs

     call OpenTraj(i)

     do j = 1 , NTrjStep(i)

!     ---------------------
#ifdef MOLFILE
       call Read_RTraj(i)
#else
       call Read_RTraj
#endif
!     ---------------------
       if(mod(j,Interval(i))/=0) cycle

       call InversMatrix(H,InvH)
       call TransCellList

       totalstep = totalstep + 1

       do k = 1, NMlip
         k1 = (k-1)*NAlip + ini
         do l = 1, Ncarbon
           i1 = k1 + Ncatom(l)
           ScR(:,l,k) = matmul(InvH,R(:,i1))
           Rlp(:,l,k) = R(:,i1)
         end do
       end do

! making neighbor lists

       Nneighbor(:) = 0
       Nbooked(:,:) = 0

       do k = 1, NMlip-1

         do kk = k+1, NMlip

inn:       do l = 1, Ncarbon
             Six = ScR(1,l,k)
             Siy = ScR(2,l,k)
             Siz = ScR(3,l,k)
             Rix = Rlp(1,l,k)
             Riy = Rlp(2,l,k)
             Riz = Rlp(3,l,k)
             do ll = 1, Ncarbon
               Sx = Six - ScR(1,ll,kk)
               Sy = Siy - ScR(2,ll,kk)
               Sz = Siz - ScR(3,ll,kk)
               Rx = Rix - Rlp(1,ll,kk)
               Ry = Riy - Rlp(2,ll,kk)
               Rz = Riz - Rlp(3,ll,kk)
              if(Sx>0.5) then
                Nx = -9
              else if(Sx<-0.5) then
                Nx =  9
              else 
                Nx =  0
              end if
              if(Sy>0.5) then
                Ny = -3
              else if(Sy<-0.5) then
                Ny =  3
              else
                Ny =  0
              end if
              if(Sz>0.5) then
                Nz = -1
              else if(Sz<-0.5) then
                Nz =  1
              else
                Nz =  0
              end if
              NNN = Nx + Ny + Nz
              Rx = Rx + CellShft(1,NNN)
              Ry = Ry + CellShft(2,NNN)
              Rz = Rz + CellShft(3,NNN)
              R2 = Rx*Rx + Ry*Ry + Rz*Rz

              if(R2 <= Rcutoff2) then
                Nneighbor(k)  = Nneighbor(k)  + 1
                Nneighbor(kk) = Nneighbor(kk) + 1
                Nbooked(k,Nneighbor(k))   = kk
                Nbooked(kk,Nneighbor(kk)) = k
                exit inn
              end if

            end do
          end do inn

        end do

      end do

! ## attibution is reset as isolated

      Isolate = NMlip + 1
      do k = 1, NMlip
        Nattrib(k) = Isolate
      end do

      do k = 1, NMlip

        iadj = Nneighbor(k)
        minatr = Nattrib(k)

! ### determine minatr 
        do l = 1, iadj
          kk   = Nbooked(k,l)  ! counterpart molecular number
          jatr = Nattrib(kk)   ! 
          if(jatr < minatr) minatr = jatr
        end do

        if(minatr == Isolate) then
          minatr = k
        end if
        Nattrib(k) = minatr ! Nattrib: root number 

! -- anyway, here I got the minimum number of the neighboring atoms
! -- we know minatr

        do l = 1, iadj
          kk   = Nbooked(k,l)
          jatr = Nattrib(kk)
          if(jatr == Isolate) then
            Nattrib(kk) = minatr
          end if
        end do

! -- found minatr among the minatr all neighboring atoms have

        minrot = klass(minatr)

        do l = 1, iadj
          kk = Nbooked(k,l)
          jatr = Nattrib(kk)
          jrot = klass(jatr)
          if(jrot < minrot) minrot = jrot
        end do

        do l = 1, iadj
          kk = Nbooked(k,l)
          jatr = Nattrib(kk)
          jrot = klass(jatr)
          Nattrib(jrot) = minrot
        end do

      end do

! ## giving the final root number

      do k = 1, NMlip
        kk = Nattrib(k)
        Nattrib(k) = klass(kk)
      end do

! ## cluster size
      do k = 1, NMlip
        nsize(k) = 0
        ndist(k) = 0
      end do
      do k = 1, NMlip
        iatr = Nattrib(k)
        nsize(iatr) = nsize(iatr) + 1
      end do

      ncls = 0
      ncls2 = 0
      maxsz = 0
      minsz = NMlip
      sum1 = 0.d0
      sum2 = 0.d0
      do k = 1, NMlip
        isize = nsize(k)
        if(isize/=0) then
          ncls = ncls + 1
          sum1 = sum1 + dble(isize)
          sum2 = sum2 + dble(isize**2)
          ndist(isize) = ndist(isize) + 1
          if(isize>=2) then
            ncls2 = ncls2 + 1
          end if
          if(isize>maxsz) then
            maxsz = isize
          end if
          if(isize<minsz) then
            minsz = isize
          end if
        end if
      end do

      avsz1 = sum1 / dble(ncls)
      avsz2 = sum2 / sum1

      write(13,'(3i10,2f10.1,2i10)') totalstep, ncls, ncls2, avsz1, avsz2, maxsz, minsz

      if(totalstep<10) then
        write(Filename,'(a,i1,a)') './Analy/distcls0000',totalstep,'.dat'
      else if(totalstep<100) then
        write(Filename,'(a,i2,a)') './Analy/distcls000',totalstep,'.dat'
      else if(totalstep<1000) then
        write(Filename,'(a,i3,a)') './Analy/distcls00',totalstep,'.dat'
      else if(totalstep<10000) then
        write(Filename,'(a,i4,a)') './Analy/distcls0',totalstep,'.dat'
      else if(totalstep<100000) then
        write(Filename,'(a,i5,a)') './Analy/distcls',totalstep,'.dat'
      else
        write(*,*) "error : totalstep is too big"
        call Finalize
      end if
      open(14,file=trim(Filename))
      do isz = 1, NMlip
        if(ndist(isz)/=0) then
          write(14,'(2i7)') isz, ndist(isz)
        end if
      end do
      close(14)

    end do

   end do

   close(13)

Contains

  function klass(m)

  integer :: m, ms, klass

1    ms=m
     m=Nattrib(m)
     if(ms.ne.m) go to 1
     klass = m

  end function klass

end subroutine IDmicelle
