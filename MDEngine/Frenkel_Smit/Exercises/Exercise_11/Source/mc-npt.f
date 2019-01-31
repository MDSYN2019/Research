      Program Mcnpt
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccc
C     Npt Simulation Of Hard-Spheres      C
Ccccccccccccccccccccccccccccccccccccccccccc

      Integer          Nstep,I,J,K,Iii,Kkk,Ipart,Ninit,Sstmm
      Logical          Lready

      Double Precision Ran_Uniform,R2,Att1,Att2,Xn,Yn,Zn,
     &     Dispmax,Dvol,Avv1,Avv2,Xold(1000),Yold(1000),
     &     Zold(1000),Press,Vnew,Vold,Ibox,Iboxold,Boxold,
     &     Box,Dx,Dy,Dz,Mvo1,Mvo2,M1
      
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

      Write(*,*) 
     &     'How Many Cycles       ?        (Example: 200           )'
      Read(*,*)  Nstep

      Write(*,*) 
     &     'How Many Init. Cycles ?        (Example: 10            )'
      Read(*,*)  Ninit
      
      Write(*,*) 
     &     'How Many Particles    ?        (Always 100 < I < 700   )'
      Read(*,*)  Npart

      Write(*,*) 
     &     'Maximum Displacement  ?        (Disp > 0 (Ex. 1.0d0    )'
      Read(*,*)  Dispmax

      Write(*,*) 
     &     'Maximum Vol. Change   ?        (Dvol > 0 (Ex. 0.01d0   )'
      Read(*,*)  Dvol

      Write(*,*) 
     &     'Product Beta*P                 (Beta*P > 0 (Ex 1.0d0)  )'
      Read(*,*)  Press


      Write(6,*)
      Write(6,*)

      If(Npart.Gt.700.Or.Npart.Lt.100) Stop

Cccccccccccccccccccccccccc
C     Initialize Stuff   C
Cccccccccccccccccccccccccc

      Att1 = 0.0d0
      Att2 = 0.0d0
      Avv1 = 0.0d0
      Avv2 = 0.0d0   
      Mvo1 = 0.0d0
      Mvo2 = 0.0d0
      Box  = 10.0d0
      Ibox = 0.1d0
      Kkk  = 0

Cccccccccccccccccccccccccccccccccccccccccc
C     Put Particles On A Lattice         C
C     Number Of Particles Always Have    C
C     To Be Smaller Than 725             C
Cccccccccccccccccccccccccccccccccccccccccc

      Do I=1,9
         Do J=1,9
            Do K=1,9

               Kkk    = Kkk + 1
               X(Kkk) = Dble(I)*1.1d0
               Y(Kkk) = Dble(J)*1.1d0
               Z(Kkk) = Dble(K)*1.1d0
            Enddo
         Enddo
      Enddo

Ccccccccccccccccccccccccccccccccc
C     Start Of The Simulation   C
Ccccccccccccccccccccccccccccccccc

      Call Writepdb

      Do Iii=1,Nstep
         Do Kkk=1,1000
            
Cccccccccccccccccccccccccccccccccccccccccccc
C     Select Particle And Move At Random   C
Cccccccccccccccccccccccccccccccccccccccccccc

            Ipart = 1 + Idint(Ran_Uniform()*Dble(Npart+1))

            If(Ipart.Gt.Npart) Then

Ccccccccccccccccccccccccccccccc
C     Volume Change           C
C     Random Walk In Ln(V)    C
Ccccccccccccccccccccccccccccccc
               
               Avv1    = Avv1 + 1.0d0
               Vold    = Box**3
               Lready  = .False.
               Vnew    = Dexp(Dlog(Vold) + 
     &              (Ran_Uniform()-0.5d0)*Dvol)
               Boxold  = Box
               Iboxold = Ibox
               Box     = Vnew**(1.0d0/3.0d0)
               Ibox    = 1.0d0/Box

Ccccccccccccccccccccccccccccccc
C     Transform Coordinates   C
Ccccccccccccccccccccccccccccccc

               Do I=1,Npart
                  Xold(I) = X(I)
                  Yold(I) = Y(I)
                  Zold(I) = Z(I)

                  X(I) = X(I)*Box*Iboxold
                  Y(I) = Y(I)*Box*Iboxold
                  Z(I) = Z(I)*Box*Iboxold
               Enddo

Ccccccccccccccccccccccccccccc
C     Check For Overlaps    C
Ccccccccccccccccccccccccccccc

               Do I=1,(Npart-1)
                  Do J=(I+1),Npart

                     Dx = X(I) - X(J)
                     Dy = Y(I) - Y(J)
                     Dz = Z(I) - Z(J)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
 
                     Dx = Dx - Box*
     &                    Dble(Idint(Dx*Ibox + 9999.5d0) - 9999)
                     Dy = Dy - Box*
     &                    Dble(Idint(Dy*Ibox + 9999.5d0) - 9999)
                     Dz = Dz - Box*
     &                    Dble(Idint(Dz*Ibox + 9999.5d0) - 9999)

                     R2 = Dx**2 + Dy**2 + Dz**2 

                     If(R2.Lt.1.0d0) Goto 110
                  Enddo
               Enddo

Ccccccccccccccccccccccccccccccccccccccccccc
C     No Overlap... Use Acceptance Rule   C
Ccccccccccccccccccccccccccccccccccccccccccc
 
               If(Ran_Uniform().Lt.Dexp(-Press*(Vnew - Vold) + 
     &              Dlog(Vnew/Vold)*Dble(Npart+1)))
     &              Lready = .True.
               
 110           If(Lready) Then

Cccccccccccccccccccccc
C     Accepted !!!!  C
Cccccccccccccccccccccc

                  Avv2 = Avv2 + 1.0d0
               Else

Ccccccccccccccccccccccccccccc
C     Rejected !!!          C
C     Restore Coordinates   C
Ccccccccccccccccccccccccccccc

                  Do I=1,Npart
                     X(I) = Xold(I)
                     Y(I) = Yold(I)
                     Z(I) = Zold(I)
                  Enddo

                  Box  = Boxold
                  Ibox = Iboxold

               Endif

            Else

Cccccccccccccccccccccc
C     Displacement   C
Cccccccccccccccccccccc

               Att1 = Att1 + 1.0d0

               Xn = X(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax
               Yn = Y(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax
               Zn = Z(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax

Cccccccccccccccccccccccccccccccccccccc
C     Put Particle Back In The Box   C
Cccccccccccccccccccccccccccccccccccccc

               Xn = Xn - Box*
     &              (Dble(Idint(Xn*Ibox + 9999.0d0) - 9999))
               Yn = Yn - Box*
     &              (Dble(Idint(Yn*Ibox + 9999.0d0) - 9999))
               Zn = Zn - Box*
     &              (Dble(Idint(Zn*Ibox + 9999.0d0) - 9999))
               
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     See If There Is An Overlap With Any Other Particle  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               Do I=1,Npart
                  If(I.Ne.Ipart) Then

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                     Dx = Xn - X(I)
                     Dy = Yn - Y(I)
                     Dz = Zn - Z(I)

                     Dx = Dx - Box*
     &                    Dble(Idint(Dx*Ibox + 9999.5d0) - 9999)
                     Dy = Dy - Box*
     &                    Dble(Idint(Dy*Ibox + 9999.5d0) - 9999)
                     Dz = Dz - Box*
     &                    Dble(Idint(Dz*Ibox + 9999.5d0) - 9999)

                     R2 = Dx**2 + Dy**2 + Dz**2 

                     If(R2.Lt.1.0d0) Goto 111
                  Endif
               Enddo

Cccccccccccccccccccccccccccccccccc
C     No Overlaps, So Accepted   C
Cccccccccccccccccccccccccccccccccc

               Att2     = Att2 + 1.0d0
               X(Ipart) = Xn
               Y(Ipart) = Yn
               Z(Ipart) = Zn

 111        Endif
         Enddo

Ccccccccccccccccccccccccccccccc
C     Make A Movie File       C
C     Print Current Volume    C
Ccccccccccccccccccccccccccccccc

         If(Mod(Iii,10).Eq.0) Then
            Call Writepdb

            Write(6,*) 'Volume / Box             : ',Box**3,Box
         Endif            

Cccccccccccccccccccccccc
C     Sample Volume    C
Cccccccccccccccccccccccc

         If(Iii.Gt.Ninit) Then
            Mvo1 = Mvo1 + 1.0d0
            Mvo2 = Mvo2 + (Box**3)
         Endif
      Enddo

Ccccccccccccccccccccccc
C     Print Results   C
Ccccccccccccccccccccccc

      Write(6,*)
      Write(6,*)
      Write(6,*) 'Fraction Succes (Displ.) : ',Att2/Att1
      Write(6,*) 'Fraction Succes (Volume) : ',Avv2/Avv1
      Write(6,*) 'Average Volume           : ',Mvo2/Mvo1
      Write(6,*) 'Average Density          : ',Dble(Npart)/(Mvo2/Mvo1)
      Write(6,*)

      Stop
      End
