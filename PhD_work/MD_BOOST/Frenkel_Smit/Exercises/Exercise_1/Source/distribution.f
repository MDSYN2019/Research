      Program Distribution
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Divide N Particles Among P Compartments    C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Sstmm,P,N,Maxp,Maxn,
     &     Ncycle,J,Jj,I,Kkk

      Parameter (Maxp = 10)
      Parameter (Maxn = 100000)

      Double Precision M1,Dist(Maxp,0:Maxn),
     &     Ran_Uniform,Faculty

      Integer          Nop(Maxp)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1) 

Ccccccccccccccccccccccccccccccccccccccc
C     Read Info From Standard Input   C
Ccccccccccccccccccccccccccccccccccccccc
 
      Write(*,*) 'Number Of Particles    ? '
      Read(*,*) N

      Write(*,*) 'Number Of Compartments ? '
      Read(*,*) P

      Write(*,*) 'Number Of Cycles       ? '
      Read(*,*) Ncycle

      If(P.Lt.2.Or.P.Gt.Maxp.Or.
     &     N.Lt.2.Or.N.Gt.Maxn) Stop

Cccccccccccccccccccc
C     Initialize   C
Cccccccccccccccccccc

      Do Jj=0,N
         Do J=1,P
            Dist(J,Jj) = 0.0d0
            Nop(J)     = 0
         Enddo
      Enddo

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do I=1,Ncycle
         Do Kkk=1,1000

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Distribute Particles                            C
C                                                     C
C     1. Loop Over All Particles                      C
C     2. Generate A Random Compartment (1...P)        C
C     3. Put This Particle In The Compartment         C
C                                                     C
C     Nop(J) = Number Of Particles In Compartment J   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification
C     End   Modification

Cccccccccccccccccccccccccccccc
C     Make Histrogram        C
Cccccccccccccccccccccccccccccc

            Do J=1,P
               Dist(J,Nop(J)) = Dist(J,Nop(J)) + 1.0d0
               Nop(J)         = 0
            Enddo
         Enddo
      Enddo

Cccccccccccccccccccccccccccccc
C     Write Results          C
Cccccccccccccccccccccccccccccc

      Open(21,File='output.dat')

      Do Jj=0,N
         Do J=1,P
            If(Dist(J,Jj).Gt.0.5d0)
     &           Dist(J,Jj) = Dist(J,Jj)/Dble(Ncycle*1000)
         Enddo
         Write(21,*) Jj,(Dist(J,Jj),J=1,P)
      Enddo

      Close(21)

Cccccccccccccccccccccccccccccccc
C     Write Analytical Dist.   C
C     For P=2                  C
Cccccccccccccccccccccccccccccccc

      If(P.Eq.2) Then
         Open(21,File='analytical.dat')

         Do Jj=0,N
            Write(21,*) Jj, Dexp(faculty(N)-faculty(JJ)-
     &           faculty(N-Jj)-dble(N)*Dlog(dble(p)))
         Enddo

         Close(21)
      Else
         Open(21,File='analytical.dat')
         Write(21,*)
         Close(21)
      Endif

Cccccccccccccccccccccccccccc
C     End Of The Program   C
Cccccccccccccccccccccccccccc

      Stop
      End
