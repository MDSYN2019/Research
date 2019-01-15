      Program Cbmc
      Implicit None

      Integer Sstmm,Nchoi,Maxtrial,I,J,Kkk,
     &        Nstep,Ichoi
      
      Parameter (Maxtrial = 100)

      Double Precision M1,Uold,Cumw,Ws,Sumnew,
     &                 Sumold,Ran_Uniform,Av3,X1,
     &                 X2,X3,X1n,X2n,X3n,X1t(Maxtrial),
     &                 X2t(Maxtrial),X3t(Maxtrial),
     &                 Ut(Maxtrial),Deltax,Av1,Av2,
     &                 Bv1,Bv2
 
Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Written By Thijs Vlugt On 6-10-1998      C
C                                              C
C     Rng Initialization                       C
Cccccccccccccccccccccccccccccccccccccccccccccccc
 
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

Ccccccccccccccccccccccccccccccccccc
C     Read Input From Std Input   C
Ccccccccccccccccccccccccccccccccccc

      Write(*,*) 'How Many Cycles (1000)    ? '
      Read(*,*) Nstep

      Write(*,*) 'How Many Trial Directions ? '
      Read(*,*) Nchoi

      If(Nchoi.Gt.Maxtrial) Stop

      Write(*,*) 'Maximum Displacement      ? '
      Read(*,*) Deltax
 
      X1   = 0.0d0
      X2   = 0.0d0
      X3   = 0.0d0
      Av1  = 0.0d0
      Av2  = 0.0d0
      Av3  = 0.0d0
      Bv1  = 0.0d0
      Bv2  = 0.0d0

Cccccccccccccccccccccccccccccccccc
C     Start The Simulation       C
C     Calculate Initial Energy   C
Cccccccccccccccccccccccccccccccccc

      Uold = X1**2 + X2**2 + X3**2

      Do I=1,Nstep
         Do J=1,1000
                        
Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Cbmc; Generate Nchoi Displacewments     C
C     Selectone According To Its Boltzmann    C
C     Factor...                               C
C                                             C
C     Old Coordinates : X1, X2, X3            C
C     New Coordinates : X1t, Y1t, Z1t         C
C     Nchoi           : Number Of Trials      C
C     Sumnew          : Weight New Config.    C
Ccccccccccccccccccccccccccccccccccccccccccccccc

            Bv2    = Bv2 + 1.0d0
            Sumnew = 0.0d0

            Do Ichoi=1,Nchoi
               X1t(Ichoi) = X1 + 
     &              (Ran_Uniform()-0.5d0)*Deltax
               X2t(Ichoi) = X2 + 
     &              (Ran_Uniform()-0.5d0)*Deltax
               X3t(Ichoi) = X3 + 
     &              (Ran_Uniform()-0.5d0)*Deltax
 
               Ut(Ichoi)  = Dexp(-(X1t(Ichoi)**2 + 
     &                             X2t(Ichoi)**2 +
     &                             X3t(Ichoi)**2))

               Sumnew = Sumnew + Ut(Ichoi)
            Enddo

Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Select One Of The Trial Directions...   C
Ccccccccccccccccccccccccccccccccccccccccccccccc

            Ws    = Ran_Uniform()*Sumnew
            Cumw  = Ut(1)
            Ichoi = 1

            If(Nchoi.Ne.1) Then
               Do While(Cumw.Lt.Ws)
                  Ichoi = Ichoi + 1
                  Cumw  = Cumw  + Ut(Ichoi)
               Enddo
            Endif

Ccccccccccccccccccccccccccccccccccccccccccccc
C     Old Configuration                     C
C     Generate Nchoi-1 Positions Around     C
C     The New Config; The First One Is The  C
C     Old Configuration...                  C
C                                           C
C     Sumold = Weight Of Old Configuration  C
Ccccccccccccccccccccccccccccccccccccccccccccc

            Sumold = 0.0d0

C     Start Modification

C     End   Modification

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Accept Or Reject This Configuration    C
Cccccccccccccccccccccccccccccccccccccccccccccc

            If(Ran_Uniform().Lt.(Sumnew/Sumold)) Then
               X1 = X1t(Ichoi)
               X2 = X2t(Ichoi)
               X3 = X3t(Ichoi)
               
               Uold = X1**2 + X2**2 + X3**2
               Bv1  = Bv1 + 1.0d0
            Endif

            Av1 = Av1 + Uold
            Av3 = Av3 + Uold**2
            Av2 = Av2 + 1.0d0

         Enddo
      Enddo
 
Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Write(6,*) 'Average Energy          : ',Av1/Av2
      Write(6,*) 'Sigma <E>               : ',
     &     Dsqrt((Av3/Av2)-(Av1/Av2)**2)
      Write(6,*) 'Fraction Accepted Moves : ',Bv1/Bv2

      Stop
      End
