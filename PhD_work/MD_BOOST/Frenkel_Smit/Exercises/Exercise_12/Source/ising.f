      Program Ising
      Implicit None

Ccccccccccccccccccccccccccc
C     2d Ising Model      C
C     No Pbc              C
C     4 Neighbors         C
Ccccccccccccccccccccccccccc

      Integer Maxlat,Nsize,Iup(4),Jup(4),Magnet,Energy,
     &     I,J,Inew,Jnew,Lnew,Lold,Diff,Sstmm,Ilat,Jlat,
     &     Kk,Kkk,Ncycle,Ninit,Maxmag,Mnew,Mold

      Parameter (Maxlat = 32   )
      Parameter (Maxmag = 1024 )

      Integer Lattice(0:Maxlat+1,0:Maxlat+1)

      Double Precision M1,Ran_Uniform,Beta,Av1,
     &     Dist2,Weight,Move1,Move2,Dist3,Av2,
     &     W(0:Maxmag),Dist1(-Maxmag:Maxmag),
     &     Ttime,Tstart

Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1) 
 
Ccccccccccccccccccccccccc
C     Read Info         C
Ccccccccccccccccccccccccc

      Write(6,*) '2d Ising Model, 4 Neightbors, No Pbc'
      Write(6,*)
      Write(6,*) 'Written By Thijs J.H. Vlugt Nov. 1999'
      Write(6,*)
      Write(6,*)

      Write(6,*) 'Lattice Size                     ? '
      Read(*,*) Nsize

      If(Mod(Nsize,2).Ne.0.Or.Mod(Maxmag,2).Ne.0.Or.
     &     Mod(Maxlat,2).Ne.0) Stop
      
      Write(6,*) 'Beta                             ? '
      Read(*,*) Beta

      Write(6,*) 'Number Of Cycles                 ? '
      Read(*,*) Ncycle

      If(Nsize.Lt.3.Or.Nsize.Gt.Maxlat) Stop
      If(Ncycle.Lt.10) Stop

      Do I=0,Maxmag
         W(I) = 1.0d0
      Enddo
      
Cccccccccccccccccccccccccccccccccc
C     Multicanonical Algorithm   C
Cccccccccccccccccccccccccccccccccc

      Open(21,File="w.dat")

      Do I=0,Maxmag,2
         Read(21,*) Kk,Av1
         
         If(Abs(Kk).Gt.Maxmag.Or.Mod(Kk,2).Ne.0) Then
            Write(6,*) 'I Out Of Range !!'
            Stop
         Endif

         W(Kk) = Av1
      Enddo

      Close(21)
      
      Ninit  = Ncycle/3
      Av1    = 0.0d0
      Av2    = 0.0d0
      Move1  = 0.0d0
      Move2  = 0.0d0
      Tstart = Ttime()
      
      If(Ninit.Gt.100) Ninit = 100

      Do I=-Maxmag,Maxmag
         Dist1(I) = 0.0d0
      Enddo

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Lattice                     C
C     Random Lattice                         C
C                                            C
C     Add One Boundary Layer With Spin=0     C
C     This Is A Way To Avoid If-Statements   C
C     In The Computation Of The Energy       C
C                                            C
C     Nsize        = Size Of The Lattice     C
C     Lattice(I,J) = Spin Of Site (I,J)      C
Cccccccccccccccccccccccccccccccccccccccccccccc

      Do I=0,Maxlat+1
         Do J=0,Maxlat+1
            Lattice(I,J) = 0
         Enddo
      Enddo

      Do I=1,Nsize
         Do J=1,Nsize
            If(Ran_Uniform().Lt.0.5d0) Then
               Lattice(I,J) =  1
            Else
               Lattice(I,J) = -1
            Endif
         Enddo
      Enddo

Ccccccccccccccccccccccccccccccc
C     Initialize Neighbours   C
Ccccccccccccccccccccccccccccccc

      Iup(1) =  1
      Jup(1) =  0
      
      Iup(2) = -1
      Jup(2) =  0

      Iup(3) =  0
      Jup(3) =  1
      
      Iup(4) =  0
      Jup(4) = -1

Ccccccccccccccccccccccccccccccccccccccc
C     Calculate Initial Energy        C
C                                     C
C     Magnet = Total Magnetisation    C
Ccccccccccccccccccccccccccccccccccccccc

      Energy = 0
      Magnet = 0

      Do I=1,Nsize
         Do J=1,Nsize

            Magnet = Magnet + Lattice(I,J)

            Do Kk=1,4
               Inew   = I      + Iup(Kk)
               Jnew   = J      + Jup(Kk)
               Energy = Energy - 
     &              Lattice(I,J)*Lattice(Inew,Jnew)
            Enddo
         Enddo
      Enddo

      Energy = Energy/2

      Write(6,*)
      Write(6,*) 'Lattice Size               : ',Nsize
      Write(6,*) 'Beta                       : ',Beta
      Write(6,*) 'Number Of Cycles           : ',Ncycle
      Write(6,*) 'Number Of Init             : ',Ninit
      Write(6,*)
      Write(6,*)
      Write(6,*) 'Initial Energy             : ',Energy
      Write(6,*) 'Initial Magnetization      : ',Magnet
      Write(6,*)
      
      If(Abs(Magnet).Gt.Maxmag) Stop

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do Kk=1,Ncycle
         Do Kkk=1,100*Nsize*Nsize
            
Ccccccccccccccccccccccccccccccccc
C     Flip A Single Spin        C
C     Metropolis Algorithm      C
Ccccccccccccccccccccccccccccccccc

            Ilat = 1 + Idint(Ran_Uniform()*Dble(Nsize))
            Jlat = 1 + Idint(Ran_Uniform()*Dble(Nsize))
            
            Mold  = Magnet
            Lold  = Lattice(Ilat,Jlat)
            Lnew  = -Lold
            Mnew  = Mold + Lnew - Lold
            Diff  = 0
            Move2 = Move2 + 1.0d0
            
Cccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Energy Difference    C
C     Between New And Old                C
Cccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

C     End Modification

CCcccccccccccccccccccccccccccccccccc
C     Acceptance/Rejection Rule    C
C     Use Weight Function W        C
Cccccccccccccccccccccccccccccccccccc

            If(Ran_Uniform().Lt.
     &           Dexp((-Beta*Dble(Diff)) + 
     &           W(Abs(Mnew)) - W(Abs(Mold)))) Then

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Update The Lattice/Energy/Magnetisation    C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

               Move1 = Move1  + 1.0d0

C     Start Modification

C     End Modification
               
               If(Abs(Magnet).Gt.Maxmag) Then
                  Write(6,*) 'Increase The Value Of Maxmag !!'
                  Stop
               Endif
            Endif
                  
            If(Kk.Gt.Ninit) Then
               Weight = Dexp(-W(Abs(Magnet)))
               Av1    = Av1   + Weight*Dble(Energy)
               Av2    = Av2   + Weight
               
               Dist1(Magnet) = 
     &              Dist1(Magnet) + 1.0d0
            Endif
         Enddo
      Enddo

      Write(6,*)
      Write(6,*) 'Average Energy             : ',
     &     Av1/Av2

      Write(6,*) 'Fraction Accepted Swaps    : ',
     &     Move1/Move2

      Write(6,*) 'Elapsed Time [S]           : ',
     &     Ttime()-Tstart

      Write(6,*)
      Write(6,*) 'Final Energy        (Simu) : ',Energy
      Write(6,*) 'Final Magnetization (Simu) : ',Magnet

Cccccccccccccccccccccccccccccccccc
C     Calculate Distributions    C
C     Magnetisation              C
Cccccccccccccccccccccccccccccccccc

      Dist2 = 0.0d0
      Dist3 = 0.0d0
            
      Do I=-Maxmag,Maxmag,2
         Weight = Dexp(-W(Abs(I)))

         Dist2  = Dist2 + Dist1(I)
         Dist3  = Dist3 + Dist1(I)*Weight
      Enddo

      Dist2 = 1.0d0/Dist2
      Dist3 = 1.0d0/Dist3

      Open(21,File="magnetic.dat")

      Do I=-Maxmag,Maxmag,2
         If(Dist1(I).Gt.0.5d0) Then
            Weight = Dexp(-W(Abs(I)))
            
            Write(21,*) I,Dist1(I)*Dist2,
     &           Dist1(I)*Dist3*Weight
         Endif
      Enddo   

      Close(21)

Ccccccccccccccccccccccccccccc
C     Weight Distribution   C
Ccccccccccccccccccccccccccccc

      Dist2 = 0.0d0

      Do I=0,Maxmag,2
         Dist1(I) = Max(1.0d0,Dist1(I) + Dist1(-I))*
     &        Dexp(-W(Abs(I)))

         Dist1(I) = -Dlog(Dist1(I))
         Dist2    = Min(Dist2,Dist1(I))
      Enddo

      Do I=0,Maxmag,2
         Dist1(I) = Dist1(I) - Dist2
         Dist1(I) = Min(10.0d0,1.05d0*Dist1(I))
      Enddo

      Open(21,File="w.dat")

      Do I=0,Maxmag,2
         Write(21,*) I,Dist1(I)
      Enddo

      Close(21)
      
Cccccccccccccccccccccccccccccccccc
C     Calculate Final Energy     C
C     Used To Check The Code     C
Cccccccccccccccccccccccccccccccccc

      Energy = 0
      Magnet = 0

      Do I=1,Nsize
         Do J=1,Nsize

            Magnet = Magnet + Lattice(I,J)

            Do Kk=1,4
               Inew   = I      + Iup(Kk)
               Jnew   = J      + Jup(Kk)
               Energy = Energy - 
     &              Lattice(I,J)*Lattice(Inew,Jnew)
            Enddo
         Enddo
      Enddo

      Energy = Energy/2

      Write(6,*) 'Final Energy        (Calc) : ',Energy
      Write(6,*) 'Final Magnetization (Calc) : ',Magnet
      Write(6,*)     

      Stop
      End
