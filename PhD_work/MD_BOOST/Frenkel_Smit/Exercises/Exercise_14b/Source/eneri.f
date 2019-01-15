      Subroutine Eneri(Xi, Yi, Zi, I, Jb, En, Vir, Ib)
      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'
 
      Double Precision Xi, Yi, Zi, En, Dx, Dy, Dz, R2, Vir, Virij, Enij
      Integer I, J, Jb, Ib

      En = 0.D0
      Vir = 0.D0
      Do J = Jb, Npart
         If (Id(J).Eq.Ib) Then
            If (J.Ne.I) Then
               Dx = Xi - X(J)
               Dy = Yi - Y(J)
               Dz = Zi - Z(J)
               If (Dx.Gt.Hbox(Ib)) Then
                  Dx = Dx - Box(Ib)
               Else
                  If (Dx.Lt.-Hbox(Ib)) Dx = Dx + Box(Ib)
               End If
               If (Dy.Gt.Hbox(Ib)) Then
                  Dy = Dy - Box(Ib)
               Else
                  If (Dy.Lt.-Hbox(Ib)) Dy = Dy + Box(Ib)
               End If
               If (Dz.Gt.Hbox(Ib)) Then
                  Dz = Dz - Box(Ib)
               Else
                  If (Dz.Lt.-Hbox(Ib)) Dz = Dz + Box(Ib)
               End If
               R2 = Dx*Dx + Dy*Dy + Dz*Dz
               Call Ener(Enij, Virij, R2, Ib)
               En = En + Enij
               Vir = Vir + Virij
            End If
         End If
      End Do

      Return
      End
