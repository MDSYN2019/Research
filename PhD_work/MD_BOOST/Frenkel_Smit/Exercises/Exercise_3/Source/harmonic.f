      Program Harmonic
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sets Of Harmonic Oscillators; E Or T Is Constant   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Sstmm,N,Maxn,Ncycle,Maxener,Ninit,Etot,
     &     Utot,I,Level,Kk,Kkk,Iup,Jup,Ilevel,Jlevel,
     &     Choice

      Parameter (Maxn    = 100000)
      Parameter (Maxener = 100000)

      Integer Elevel(0:Maxn)

      Double Precision M1,Ran_Uniform,Beta,
     &     Dist1(0:Maxener),Dist2,Av1,Av2

Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1) 
 
Ccccccccccccccccccccccccccccc
C     Read Info From Disk   C
Ccccccccccccccccccccccccccccc

      Etot = 0
      Beta = 0.0d0

      Write(*,*) 'Number Of Oscillators  ? '
      Read(*,*) N

      Write(*,*) 'Number Of Cycles       ? '
      Read(*,*) Ncycle

      If(Ncycle.Le.3) Stop

      Write(*,*) '1. Nve Ensemble'
      Write(*,*) '2. Nvt Ensemble'
      Read(*,*) Choice

      If(Choice.Eq.1) Then
         Write(*,*) 'Total Energy           ? '
         Read(*,*) Etot

         If(N.Le.3.Or.N.Ge.Maxn.Or.
     &     Etot.Le.3.Or.Etot.Ge.Maxener) Stop
      Else
         Write(*,*) 'Beta            ? '
         Read(*,*) Beta

         If(Beta.Lt.0.0d0.Or.Beta.Gt.1.0d7) Stop

         N = 1
      Endif

Cccccccccccccccccccc
C     Initialize   C
Cccccccccccccccccccc

      Do I=1,Maxn
         Elevel(I) = 0
      Enddo

      Do I=1,Maxener
         Dist1(I) = 0.0d0
      Enddo

      Ninit = Ncycle/2
      Av1   = 0.0d0
      Av2   = 0.0d0
      Dist2 = 0.0d0
      Utot  = 0
      Level = 0
      
Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Make An Initial Distribution Of The Total   C
C     Energy Over The Levels                      C
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      If(Choice.Eq.1) Then

 10      Level = Level + 1

         If(Level.Gt.N) Level = 1

         Elevel(Level) = Elevel(Level) + 1
         Utot          = Utot          + 1

         If(Utot.Ne.Etot) Goto 10
      Endif

Cccccccccccccccccccccccccccccccccccccccc
C     Calculate Total Initial Energy   C
Cccccccccccccccccccccccccccccccccccccccc
      
      Utot = 0

      Do I=1,N
         Utot = Utot + Elevel(I)
      Enddo

      Write(6,*)
      Write(6,*) 'Initial Energy : ',Utot     

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do I=1,Ncycle
         Do Kk=1,100
            Do Kkk=1,N

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Choose 2 Different Levels At Random    C
C     Exchange Energy At Random              C
Cccccccccccccccccccccccccccccccccccccccccccccc

               If(Choice.Eq.1) Then

Cccccccccccccccccccccc
C     Nve Ensemble   C
Cccccccccccccccccccccc

 20               Ilevel = 1 + 
     &                 Idint(Dble(N)*Ran_Uniform())

                  Jlevel = 1 + 
     &                 Idint(Dble(N)*Ran_Uniform())

                  If(Ilevel.Eq.Jlevel) Goto 20
                  
                  If(Ran_Uniform().Lt.0.5d0) Then
                     Iup =  1
                     Jup = -1
                  Else
                     Iup = -1
                     Jup =  1
                  Endif

Cccccccccccccccccccccccccccccccccccccccccccc
C     Calculate New Energy Of Both Levels  C
Cccccccccccccccccccccccccccccccccccccccccccc

                  Elevel(Ilevel) = Elevel(Ilevel) + Iup
                  Elevel(Jlevel) = Elevel(Jlevel) + Jup

Cccccccccccccccccccccccc
C     Reject When E<0  C
Cccccccccccccccccccccccc

                  If(Min(Elevel(Ilevel),Elevel(Jlevel)).Lt.0) Then
                     Elevel(Ilevel) = Elevel(Ilevel) - Iup
                     Elevel(Jlevel) = Elevel(Jlevel) - Jup 
                  Endif
                  
               Else
                  
Ccccccccccccccccccccccc
C     Nvt Ensemble    C
Ccccccccccccccccccccccc

                  If(Ran_Uniform().Lt.0.5d0) Then
                     Iup =  1
                  Else
                     Iup = -1
                  Endif

                  Elevel(1) = Elevel(1) + Iup

Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Accept/Reject Using The Metropolis Algorithm   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc

                  If((Elevel(1).Lt.0).Or.
     &                 (Ran_Uniform().Gt.Dexp(-Beta*Dble(Iup)))) Then

                     Elevel(1) = Elevel(1) - Iup
                  Endif
               Endif
            Enddo

Cccccccccccccccc
C     Sample   C
Cccccccccccccccc

            If(I.Ge.Ninit) Then
               Dist1(Elevel(1)) = Dist1(Elevel(1)) + 1.0d0
               Dist2            = Dist2            + 1.0d0
               Av1              = Av1              + Dble(Elevel(1))
               Av2              = Av2              + 1.0d0
            Endif
         Enddo
      Enddo

Cccccccccccccccccccccccccccccc
C     Write Results          C
Cccccccccccccccccccccccccccccc

      Utot = 0

      Do I=1,N
         Utot = Utot + Elevel(I)
      Enddo

      Write(6,*)
      Write(6,*) 'Final Energy             : ',Utot
      Write(6,*) 'Average Energy Level 1   : ',Av1/Av2

      Dist2 = 1.0d0/Dist2

      Open(21,File='output.dat')

      Do I=0,Maxener
         If(Dist1(I).Gt.0.5d0) Then
            Write(21,'(I10,E20.10)') I,Dist1(I)*Dist2
         Endif
      Enddo

      Close(21)

      Stop
      End
