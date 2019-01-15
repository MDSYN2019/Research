      Subroutine Mcloop
      Implicit None

      Include 'system.inc'

      Integer I,J,K

      Double Precision Xst(Maxchain),Yst(Maxchain),
     &     Zst(Maxchain),Weight,Weiold,Ubnew,
     &     Ubold,Unnew,Unold,R2,Ran_Uniform,
     &     Bn,Bs,Av1,Av2,Av3,Eb,En

Ccccccccccccccccccccccccccccccc
C     Generate A First Chain  C
Ccccccccccccccccccccccccccccccc

      Weight = 1.0d0
      Ubnew  = 0.0d0
      Unnew  = 0.0d0
      Weiold = 1.0d0
      Ubold  = 0.0d0
      Unold  = 0.0d0
      Bn     = 0.0d0
      Bs     = 0.0d0
      Av1    = 0.0d0
      Av2    = 0.0d0
      Av3    = 0.0d0
      
      Call Grow(.False.,Weiold,Ubold,Unold,Xst,Yst,Zst)

      Do J=1,Nuall
         Xpos(J) = Xst(J)
         Ypos(J) = Yst(J)
         Zpos(J) = Zst(J)
      Enddo

      Eb = Ubold
      En = Unold

      Call Sample(1,1.0d0,1.0d0)

      Do I=1,Nstep
         Do K = 1,1000

Ccccccccccccccccccccccc
C     Static Scheme   C
Ccccccccccccccccccccccc

            If(Lstatic) Then
               Call Grow(.False.,Weight,Ubnew,Unnew,Xst,Yst,Zst)

               Do J=1,Nuall
                  Xpos(J) = Xst(J)
                  Ypos(J) = Yst(J)
                  Zpos(J) = Zst(J)
               Enddo

               Av1 = Av1 + Ubnew*Weight
               Av2 = Av2 + Unnew*Weight
               Av3 = Av3 + Weight

            Else

Cccccccccccccccccccccccc
C     Dynamic Scheme   C
Cccccccccccccccccccccccc

               Bn = Bn + 1.0d0

               Call Grow(.True. ,Weiold,Ubold,Unold,Xst,Yst,Zst)
               Call Grow(.False.,Weight,Ubnew,Unnew,Xst,Yst,Zst)

Cccccccccccccccccccccccccccccccccccccccccc
C     Accepted ? Then Copy Coordinates   C
Cccccccccccccccccccccccccccccccccccccccccc

               If(Ran_Uniform().Lt.(Weight/Weiold)) Then

                  Bs = Bs + 1.0d0
                  En = Unnew
                  Eb = Ubnew

                  Do J=1,Nuall
                     Xpos(J) = Xst(J)
                     Ypos(J) = Yst(J)
                     Zpos(J) = Zst(J)
                  Enddo
               Endif

               Weight = 1.0d0
            Endif

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sample End-To-End Distance And Average Energies     C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            R2 = Dsqrt((Xpos(Nuall)-Xpos(1))**2 +
     &                 (Ypos(Nuall)-Ypos(1))**2 +
     &                 (Zpos(Nuall)-Zpos(1))**2)

            If(I.Gt.Ninit) Then
               Call Sample(2,R2,Weight)

               If(.Not.Lstatic) Then
                  Av1 = Av1 + Eb
                  Av2 = Av2 + En
                  Av3 = Av3 + 1.0d0
               Endif

               If(K.Eq.1.And.Mod(I,5).Eq.0) Call Writepdb
            Endif
         Enddo
      Enddo

      Call Sample(3,1.0d0,1.0d0)

      If(.Not.Lstatic) Write(6,*) 'Fraction Accepted : ',Bs/Bn

      Write(6,*) 'Average U-Bonded  : ',Av1/Av3
      Write(6,*) 'Average U-Nonb    : ',Av2/Av3

      Return
      End
