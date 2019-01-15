      Program Pi
      Implicit None

      Integer          I,J,Nstep,Sstmm
      Double Precision Ran_Uniform,Ll,X,Y,At1,At2,Pii,M1

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculates Pi Using The Circle/Square Problem   C
C     After Execution: xmgr results.data              C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

      Write(*,*) 'Number Of Cycles ? (Example: 1000)'
      Read(*,*)  Nstep

      Write(*,*) 'Ratio L/D        ? (Always >= 1 !)'
      Read(*,*)  Ll

      Ll = Ll*2.0d0

      If(Ll.Lt.1.0d0) Then
         Write(6,*) 'Ratio Must Be At Least 1 !!!'
         Stop
      Endif

      At1 = 0.0d0
      At2 = 0.0d0

      Open(21,File='results.data',Form='Formatted')

Ccccccccccccccccccccccccccccc
C     Loop Over All Cycles  C
Ccccccccccccccccccccccccccccc

      Do I=1,Nstep
         Do J=1,1000

Ccccccccccccccccccccccccccccccccccccc
C     Generate A Uniform Point      C
C     Check If It Is In The Circle  C
C     At2 = Number Of Points        C
C     At1 = Number In The Circle    C
Ccccccccccccccccccccccccccccccccccccc

C Start Modifications


C End Modifications

         Enddo
         If(Mod(I,10).Eq.0) Write(21,*) I,Ll*Ll*At1/At2
      Enddo

      Close(21)

      Pii = Ll*Ll*At1/At2

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     The Real Value Of Pi Can Be Calculated Using       C
C     Pi = 4.0 * Arctan (1.0)                            C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Write(6,*) 'Estimate Of Pi : ',Pii
      Write(6,*) 'Real Pi        : ',4.0d0*Datan(1.0d0)
      Write(6,*) 'Relative Error : ',
     &           Dabs(Pii-4.0d0*Datan(1.0d0))/(4.0d0*Datan(1.0d0))
      
      Stop
      End
