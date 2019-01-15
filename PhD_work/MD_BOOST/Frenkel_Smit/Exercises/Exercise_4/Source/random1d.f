      Program Random1d
      Implicit None

Cccccccccccccccccccccccccc
C     Random Walk 1d     C
Cccccccccccccccccccccccccc

      Integer Sstmm,Njump,Ncycle,Maxlat,I,Ip,Kkk,Kk

      Parameter (Maxlat = 100000)

      Double Precision Dist1(-Maxlat:Maxlat),Dist2,P,M1,
     &     Ran_Uniform

Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)
      Call Sample(1,1)
 
Ccccccccccccccccccccccccccccccc
C     Read Input From Disk    C
Ccccccccccccccccccccccccccccccc

      Write(*,*) 'Number Of Cycles       ? '
      Read(*,*) Ncycle

      Write(*,*) 'Number Of Jump/Cycle   ? '
      Read(*,*) Njump

      Do I=-Maxlat,Maxlat
         Dist1(I) = 0.0d0
      Enddo

      Dist2 = 0.0d0
      
Ccccccccccccccccccccccccccc
C     Loop Over Cycles    C
Ccccccccccccccccccccccccccc

      Do Kkk=1,Ncycle
         Do Kk=1,1000
            
            Call Sample(2,1)

            Ip = 0
         
CCccccccccccccccccccccccccccccccccccc
C     Perform The Random Walk       C
C     Go Up Or Down With Prob.0.5   C
Ccccccccccccccccccccccccccccccccccccc

            Do I=1,Njump
               If(Ran_Uniform().Lt.0.5d0) Then
                  Ip = Ip + 1
               Else
                  Ip = Ip - 1
               Endif

               Call Sample(3,Ip)
            Enddo

            Dist2 = Dist2 + 1.0d0

            If(Ip.Ge.-Maxlat.Or.Ip.Le.Maxlat) 
     &           Dist1(Ip) = Dist1(Ip) + 1.0d0
         Enddo
      Enddo

Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Call Sample(4,1)

      Open(21,File="output.dat")

      Dist2 = 1.0d0/Dist2

      Do I=-Maxlat,Maxlat
         If(Dist1(I).Gt.0.5d0) Then
            Dist1(I) = Dist1(I)*Dist2
            P        = 0.5d0*Dlog(2.0d0/
     &           (Dble(Njump)*4.0d0*Datan(1.0d0))) - 
     &           ((Dble(I*I))/(2.0d0*Dble(Njump)))

            Write(21,'(I10,2e20.10)') I,Dist1(I),
     &           Dexp(P)
         Endif
      Enddo

      Close(21)

      Stop
      End
