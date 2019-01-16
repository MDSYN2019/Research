      Subroutine Eneri(Xi, Yi, Zi, I, Jb, En, Vir)
C
C    Calculates The Energy Of Particle I With Particles J=Jb,Npart
C
C     Xi (Input)    X Coordinate Particle I
C     Yi (Input)    Y Coordinate Particle I
C     Zi (Input)    Z Coordinate Particle I
C     I  (Input)    Particle Number (Excluded !!!)
C     En  (Output)  Energy Particle I
C     Vir (Output)  Virial Particle I
C
C
      Implicit None
      
      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'
 
      Double Precision Xi, Yi, Zi, En, Dx, Dy, Dz, 
     &     R2, Vir, Virij, Enij
      Integer I, J, Jb

      En  = 0.0d0
      Vir = 0.0d0

      Do J = Jb, Npart

Ccccccccccccccccccccccccccc
C     Excluse Particle I  C
Ccccccccccccccccccccccccccc

         If (J.Ne.I) Then

            Dx = Xi - X(J)
            Dy = Yi - Y(J)
            Dz = Zi - Z(J)

Ccccccccccccccccccccccccccccccccccc
C     Nearest Image Convention    C
Ccccccccccccccccccccccccccccccccccc

            If (Dx.Gt.Hbox) Then
               Dx = Dx - Box
            Elseif (Dx.Lt.-Hbox) Then
               Dx = Dx + Box
            End If

            If (Dy.Gt.Hbox) Then
               Dy = Dy - Box
            Elseif (Dy.Lt.-Hbox) Then
               Dy = Dy + Box
            End If

            If (Dz.Gt.Hbox) Then
               Dz = Dz - Box
            Elseif (Dz.Lt.-Hbox) Then
               Dz = Dz + Box
            End If

            R2 = Dx*Dx + Dy*Dy + Dz*Dz

Cccccccccccccccccccccccccccccc
C     Calculate The Energy   C
Cccccccccccccccccccccccccccccc

            Call Ener(Enij, Virij, R2)

            En  = En + Enij
            Vir = Vir + Virij
         Endif
      Enddo

      Return
      End
