      Program Photon
      Implicit None

      Integer          Nstep,Ninit,New,Old,I,J,Sstmm
      Double Precision Ran_Uniform,Beta,Av1,Av2,M1
      
Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Read Info From Disk                  C
Cccccccccccccccccccccccccccccccccccccccccccc

      Write(*,*) 'How Many Cycles ?                (Example: 1000 )'
      Read(*,*)  Nstep

      Write(*,*) 'How Many Initialization Cycles ? (Example: 100  )'
      Read(*,*)  Ninit

      Write(*,*) 'Beta*Epsilon ?                   (Example: 1.0d0)'
      Read(*,*)  Beta

      New = 1
      Old = 1
      Av1 = 0.0d0
      Av2 = 0.0d0

CCcccccccccccccccccccccccccccccccccccccccccc
C     Loop Over All Cycles                 C
C                                          C
C     Old = Old Position (Integer !!)      C
C     New = New Position (Integer !!)      C
Cccccccccccccccccccccccccccccccccccccccccccc

      Do I=1,Nstep
         Do J=1,1000

C Start Modification

C End   Modification

Cccccccccccccccccccccccccccccc
C     Check For Acceptance   C
Cccccccccccccccccccccccccccccc
            
            If(Ran_Uniform().Lt.Dexp(-Beta*(Dble(New-Old)))) 
     &           Old = New

cCccccccccccccccccccccccccccccccccccccccccccc
C     Calculate Average Occupancy Result    C
Ccccccccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.Ninit) Then
               Av1 = Av1 + Dble(Old)
               Av2 = Av2 + 1.0d0
            Endif
         Enddo
      Enddo

Cccccccccccccccccccccccccccccccccccccccccccc
C     Write The Final Result               C
Cccccccccccccccccccccccccccccccccccccccccccc

      Write(6,*) 'Average Value     : ',Av1/Av2
      Write(6,*) 'Theoretical Value : ',
     &     1.0d0/(Dexp(Beta)-1.0d0)
      Write(6,*) 'Ratio             : ',
     &     Dabs((Dexp(Beta)-1.0d0)*((Av1/Av2) - 
     &     (1.0d0/(Dexp(Beta)-1.0d0))))

      Stop
      End
