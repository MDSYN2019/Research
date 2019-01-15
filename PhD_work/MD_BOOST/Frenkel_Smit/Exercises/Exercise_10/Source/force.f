      Subroutine Force
      Implicit None
 
C     Calculate The Forces And Potential Energy
 
      Include 'system.inc'
 
      Integer I,J
      Double Precision Dx,Dy,Dz,Ff,R2i,R6i
 
C     Set Forces, Potential Energy And Pressure To Zero
 
      Do I = 1,Npart
         Fxx(I) = 0.0d0
         Fyy(I) = 0.0d0
         Fzz(I) = 0.0d0
      Enddo
 
      Upot = 0.0d0
      Press = 0.0d0
 
C     Loop Over All Particle Pairs
 
      Do I = 1,Npart - 1
         Do J = I + 1,Npart
 
C     Calculate Distance And Perform Periodic
C     Boundary Conditions
 
            Dx = Rxx(I) - Rxx(J)
            Dy = Ryy(I) - Ryy(J)
            Dz = Rzz(I) - Rzz(J)
 
            If (Dx.Gt.Hbox) Then
               Dx = Dx - Box
            Elseif (Dx.Lt. - Hbox) Then
               Dx = Dx + Box
            Endif
 
            If (Dy.Gt.Hbox) Then
               Dy = Dy - Box
            Elseif (Dy.Lt. - Hbox) Then
               Dy = Dy + Box
            Endif
 
            If (Dz.Gt.Hbox) Then
               Dz = Dz - Box
            Elseif (Dz.Lt. - Hbox) Then
               Dz = Dz + Box
            Endif
 
            R2i = Dx*Dx + Dy*Dy + Dz*Dz

C     Check If The Distance Is Within The Cutoff Radius
 
            If (R2i.Lt.Rcutsq) Then
               R2i = 1.0d0/R2i
               R6i = R2i*R2i*R2i
 
               Upot  = Upot + 4.0d0*R6i*(R6i - 1.0d0) - Ecut
               Ff    = 48.0d0*R6i*(R6i - 0.5d0)
               Press = Press + Ff
               Ff    = Ff*R2i
 
               Fxx(I) = Fxx(I) + Ff*Dx
               Fyy(I) = Fyy(I) + Ff*Dy
               Fzz(I) = Fzz(I) - Ff*Dz
 
               Fxx(J) = Fxx(J) - Ff*Dx
               Fyy(J) = Fyy(J) - Ff*Dy
               Fzz(J) = Fzz(J) - Ff*Dz
 
            Endif
         Enddo
      Enddo
 
C     Scale The Pressure
 
      Press = Press/(3.0d0*Box*Box*Box)
 
      Return
      End
