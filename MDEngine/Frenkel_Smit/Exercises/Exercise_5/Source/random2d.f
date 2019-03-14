      Program Random2d
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccc
C     Random Walk 2d     C
Cccccccccccccccccccccccccc

      Integer I,J,Ii,Sstmm,Iup,Iiup,Inew,Iinew,Njump,
     &     Kkk,Nlattice,Ninit,Kk

      Double Precision M1,Ran_Uniform,Av1,Av2,Rm

Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)
      Call Sample(1)
 
Cccccccccccccccccccccccccccccc
C     Read Input From Disk   C
Cccccccccccccccccccccccccccccc

      Write(*,*) 'Number Of Jumps (/1000) ? '
      Read(*,*) Njump

      Write(*,*) 'Lattice Size           ? '
      Read(*,*) Nlattice

      If(Nlattice.Lt.10.Or.Nlattice.Ge.Maxlat) Stop

      Write(*,*) 'Number Of Particles    ? '
      Read(*,*) Npart

      If(Npart.Lt.1.Or.Npart.Ge.Maxpart) Stop

      Do I=1,Maxlat
         Do Ii=1,Maxlat
            Lattice(I,Ii) = 0
         Enddo
      Enddo

      Do I=1,Maxpart
         Ipart(I)  = 0
         Iipart(I) = 0
         Mxx(I)    = 0
         Myy(I)    = 0
      Enddo

      Av1   = 0.0d0
      Av2   = 0.0d0
      Ninit = Njump/4

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Put Particles Random On The Lattice    C
C     Look For An Empty Lattice Site         C
Cccccccccccccccccccccccccccccccccccccccccccccc

      Do J=1,Npart
 10      I  = 1 + Idint(Ran_Uniform()*Dble(Nlattice))
         Ii = 1 + Idint(Ran_Uniform()*Dble(Nlattice))

         If(Lattice(I,Ii).Ne.0) Goto 10

         Lattice(I,Ii) = J
         Ipart(J)      = I
         Iipart(J)     = Ii
         Mxx(J)        = I
         Myy(J)        = Ii 
      Enddo
         
Cccccccccccccccccccccccccccccccccc
C     Perform The Random Walk    C
Cccccccccccccccccccccccccccccccccc

      Do Kkk=1,Njump
         Do Kk=1,1000
            
Ccccccccccccccccccccccccccccccccccccccccccccc
C     Choose A Random Site (J)              C
C     Choose A Random Displacement (Rm):    C
C     - Up, Down, Left, Right               C
Ccccccccccccccccccccccccccccccccccccccccccccc

            J  = 1 + Idint(Ran_Uniform()*Dble(Npart))
            Rm = 4.0d0*Ran_Uniform()

            If(Rm.Lt.1.0d0) Then
               Iup  =  1
               Iiup =  1
            Elseif(Rm.Lt.2.0d0) Then
               Iup  = -1
               Iiup = -1
            Elseif(Rm.Lt.3.0d0) Then
               Iup  = -1
               Iiup =  1
            Else
               Iup  =  1
               Iiup = -1
            Endif 

Cccccccccccccccccccccccccccccccccccccccccccc
C     New Position                         C
C     Put Particle Back On The Lattice     C
Cccccccccccccccccccccccccccccccccccccccccccc

            Inew  = Ipart(J)  + Iup
            Iinew = Iipart(J) + Iiup

            If(Inew.Lt.1) Then
               Inew = Inew + Nlattice
            Elseif(Inew.Gt.Nlattice) Then
               Inew = Inew - Nlattice
            Endif

            If(Iinew.Lt.1) Then
               Iinew = Iinew + Nlattice
            Elseif(Iinew.Gt.Nlattice) Then
               Iinew = Iinew - Nlattice
            Endif

            Av2 = Av2 + 1.0d0

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Check If New Lattice Site Is Accepted     C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

            If(Lattice(Inew,Iinew).Eq.0) Then
               Lattice(Ipart(J),Iipart(J)) = 0
               Ipart(J)                    = Inew
               Iipart(J)                   = Iinew
               Lattice(Inew,Iinew)         = J
               Mxx(J)                      = Mxx(J) + Iup
               Myy(J)                      = Myy(J) + Iiup
               Av1                         = Av1    + 1.0d0
            Endif

            If(Kkk.Ge.Ninit) Call Sample(2)
         Enddo
      Enddo

Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Call Sample(3)

      Write(6,*)
      Write(6,*) 'Fraction Accepted Jumps : ',Av1/Av2
      Write(6,*) 'Lattice Occupation      : ',
     &     Dble(Npart)/Dble(Nlattice*Nlattice)

      Stop
      End
