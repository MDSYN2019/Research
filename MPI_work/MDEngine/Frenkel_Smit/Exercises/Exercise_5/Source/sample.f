      Subroutine Sample(Switch)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculates Mean Square Displacement                    C
C                                                            C
C     Tmax  = Maximum Timesteps For The Correlation Time     C
C     T0max = Maximum Number Of Time Origins                 C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer   Tmax,T0max

      Parameter (Tmax  = 5000)
      Parameter (T0max = 500)

      Integer Ttel,Tt0(T0max),T0,Switch,I,T,Tvacf,Dt,
     &     Rx0(Maxpart,T0max),Ry0(Maxpart,T0max)

      Double Precision R2tx(Tmax),R2ty(Tmax),Nvacf(Tmax),Dtime
       
      Save Nvacf,R2tx,R2ty,Rx0,Ry0,Tt0,T0,Tvacf
       
      If (Switch.Eq.1) Then

Ccccccccccccccccccccccccccccccccc
C     Initialize Everything     C
Ccccccccccccccccccccccccccccccccc

         Tvacf = 0
         T0    = 0

         Do I = 1,Tmax
            R2tx(I)  = 0.0d0
            R2ty(I)  = 0.0d0
            Nvacf(I) = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

         Tvacf = Tvacf + 1

         If (Mod(Tvacf,50).Eq.0) Then

Cccccccccccccccccccccccccccc
C     New Time Origin      C
Cccccccccccccccccccccccccccc

            T0         = T0 + 1
            Ttel       = Mod(T0-1, T0max) + 1
            Tt0(Ttel)  = Tvacf

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Store Particle Positions/Velocities       C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

            Do I = 1, Npart
               Rx0(I,Ttel)  = Mxx(I)
               Ry0(I,Ttel)  = Myy(I)
            Enddo
         Endif

Cccccccccccccccccccccccccccccccccccccc
C     Loop Over All Time Origins     C
Cccccccccccccccccccccccccccccccccccccc

         Do T = 1, Min(T0, T0max)
            Dt = Tvacf - Tt0(T) + 1

            If (Dt.Lt.Tmax) Then
               Nvacf(Dt) = Nvacf(Dt) + 1.0d0
                  
               Do I = 1, Npart
                  R2tx(Dt) = R2tx(Dt) + (Dble(Mxx(I)-Rx0(I,T)))**2
                  R2ty(Dt) = R2ty(Dt) + (Dble(Myy(I)-Ry0(I,T)))**2
               Enddo
            Endif
         Enddo
                
      Else

Cccccccccccccccccccccccccccccccccccc
C     Write Everything To Disk     C
Cccccccccccccccccccccccccccccccccccc

         Open(25,File="rms.dat")

         Do I = 1, Tmax-1
            Dtime = Dble(I-1)

            If(Nvacf(I).Ge.0.5d0) Then
               R2tx(I) = R2tx(I)/Nvacf(I)
               R2ty(I) = R2ty(I)/Nvacf(I)
            Else
               R2tx(I) = 0.0d0
               R2ty(I) = 0.0d0
            Endif

            Write(25,'(3E20.10)')
     &           Dtime,R2tx(I),R2ty(I)
         Enddo

         Close(25)
      Endif

      Return
      End
