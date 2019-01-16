      Program Boltzmann
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Boltzmann Distribution       C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Maxn,I,N

      Parameter (Maxn = 10000)

      Double Precision Distri(0:Maxn),Norm,Beta,Temp

Cccccccccccccccccccccccccccccccccccccccccc
C     Read Info From Disk                C
C     Maxn = Maximum Number Over Levels  C
Cccccccccccccccccccccccccccccccccccccccccc

      Write(*,*) 'Number Of Energy Levels (2-10000) ? '
      Read(*,*) N

      Write(*,*) 'Temperature ? '
      Read(*,*) Temp

      If(N.Lt.2.Or.N.Ge.Maxn.Or.
     &     Temp.Lt.1.0d-7.Or.Temp.Gt.1.0d7) Stop

      Beta = 1.0d0/Temp
      Norm = 0.0d0

Cccccccccccccccccccccccccccccc
C     Loop Over All Levels   C
Cccccccccccccccccccccccccccccc

      Do I=0,(N-1)
         Temp = Dexp(-Beta*Dble(I))

         If(Temp.Lt.1.0d-70) Temp = 0.0d0

         Norm      = Norm + Temp
         Distri(I) = Temp
      Enddo

      Norm = 1.0d0/Norm

Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Open(21,File="result.dat")

      Do I=0,(N-1)
         Write(21,'(I8,F20.10)') I,Distri(I)*Norm
      Enddo

      Close(21)

      Stop
      End
