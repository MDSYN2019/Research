      Program Scaling
      Implicit None

      Double Precision M1,Ran_Uniform,
     &     Dist(100),Dis,C,Xold,Xnew,
     &     Av1,Av2,Maxdpl

      Logical Lscale
      Integer I,J,Sstmm,Kk,Ncycle

Cccccccccccccccccccccccccccccccccccccccccc
C     Initialize Rng                     C
Cccccccccccccccccccccccccccccccccccccccccc

      J  = Sstmm()
      M1 = 0.001d0*Dble(Mod((10+10*J),1000))
      
      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0 
      
      Call Genrand(M1)

      J = 10 + Idint(Ran_Uniform()*1000.0d0)

      Do I=1,J
         M1 = Ran_Uniform()
      Enddo

CCcccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Everything                              C
C                                                        C
C     Xold = Old Position                                C
C     Xnew = New Position                                C
C     Lscale = Do We Use Scaling ? (.True. Or .False.)   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Do I=1,100
         Dist(I) = 0.0d0
      Enddo

      Xold = 0.1d0 + 0.7d0*Ran_Uniform()
      Dis  = 0.0d0
      Av1  = 0.0d0
      Av2  = 0.0d0

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Do We Use Normal Displacements Or Scaling  ???   C
C                                                      C
C     Maxdpl = Maximum Displacement                    C
C     Ncycle = Number Of Cycles                        C
C                                                      C
C     Change The Input Parameters Here !!!!!           C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Lscale = .True.
      Maxdpl = 0.7d0
      Ncycle = 10000

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Start The Simulation                             C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Do I=1,Ncycle
         Do J=1,10000

            Av2 = Av2 + 1.0d0

            If(Lscale) Then

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Scale The Position By A Factor C                 C
C     Select At Random To Divide Or Multiply By C      C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               C = 1.0d0 + Ran_Uniform()

               If(Ran_Uniform().Lt.0.5d0) C = 1.0d0/C

               Xnew = C*Xold

C     Start Modification

C     End   Modification

            Else

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Random Displacement                              C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               Xnew = Xold + 
     &              2.0d0*(Ran_Uniform()-0.5d0)*
     &              Maxdpl

C     Start Modification

C     End   Modification

            Endif

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sample Distribution                              C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            Kk = 1 + Idint(Xold*100.D0)

            If(Kk.Ge.1.And.Kk.Le.100)
     &           Dist(Kk) = Dist(Kk) + 1.0d0
            
            Dis = Dis + 1.0d0

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Some Dummy Random Numbers     C
C     Does This Help ???!!!???               C
Cccccccccccccccccccccccccccccccccccccccccccccc

            Do Kk = 1,5
               C = Ran_Uniform()
            Enddo
         Enddo
      Enddo

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Write Results To Disk                            C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Write(6,*) 'Fraction Accepted Trialmoves : ',
     &     Av1/Av2

      Dis = 1.0d0/Dis

      Open(21,File="distri.dat")

      Do I=1,100
         Write(21,*) 0.01d0*Dble(I),Dist(I)*Dis
      Enddo
      
      Close(21)

      Stop
      End
