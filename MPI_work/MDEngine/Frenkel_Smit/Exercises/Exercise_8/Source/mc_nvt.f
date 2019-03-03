      Program Mc_Nvt
      Implicit None
C____________________________________________________C
C                                                    C
C                                                    C
C   Equation Of State Of The Lennard-Jones Fluid     C
C                                                    C
C____________________________________________________C
 
      Include 'potential.inc'
      Include 'parameter.inc'
      Include 'system.inc'
      Include 'conf.inc'

      Integer Equil, Prod, Nsamp, Ii, Icycl, Ndispl, Attempt, 
     &        Nacc, Ncycl, Nmoves, Imove, Kkk
      Double Precision En, Ent, Vir, Virt, Dr, Av1, Av2, Press,
     &                 Bv1,Bv2
 
      Write (6, *) '**************** Mc_Nvt ***************'

C     ---Initialize Sysem

      Call Readdat(Equil, Prod, Nsamp, Ndispl, Dr)

      Nmoves = Ndispl
      Av1    = 0.0d0
      Av2    = 0.0d0
      Bv1    = 0.0d0
      Bv2    = 0.0d0

C     ---Total Energy Of The System

      Call Toterg(En, Vir)
      Write (6, 99001) En, Vir

C     ---Start Mc-Cycle

      Do Ii = 1, 2

C        --- Ii=1 Equilibration
C        --- Ii=2 Production

         If (Ii.Eq.1) Then
            Ncycl = Equil
            If (Ncycl.Ne.0) Write (6, *) ' Start Equilibration '
         Else
            If (Ncycl.Ne.0) Write (6, *) ' Start Production '
            Ncycl = Prod
         End If
        
         Attempt = 0
         Nacc    = 0

C        ---Intialize The Subroutine That Adjust The Maximum Displacement

         Call Adjust(Attempt, Nacc, Dr)

         Do Icycl = 1, Ncycl

            Do Imove = 1, Nmoves

C              ---Attempt To Displace A Particle

               Call Mcmove(En, Vir, Attempt, Nacc, Dr)
            End Do

            If (Ii.Eq.2) Then

C              ---Sample Averages

               If (Mod(Icycl,Nsamp).Eq.0) Then
                  Call Sample(Icycl, En, Vir, Press)
                  Av1 = Av1 + Press
                  Av2 = Av2 + 1.0d0

Ccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Chemical Potential  C
C     Do 10 Trial Chains                C
C     Calculate The Average Of          C
C     [Exp(-Beta*Energy)]               C
C     You Can Use The Subroutine Eneri  C
C     For This. Good Luck !!!           C
Ccccccccccccccccccccccccccccccccccccccccc

                  Do Kkk=1,10

C     Start Modification

                     Bv2 = Bv2 + 1.0d0

C     End   Modification

                  Enddo
               Endif
            Endif

            If(Mod(Icycl,20).Eq.0) Call Writepdb

            If (Mod(Icycl,Ncycl/5).Eq.0) Then
               Write (6, *) '======>> Done ', Icycl, ' Out Of ', Ncycl

C              ---Write Intermediate Configuration To File

               Call Store(8, Dr)

C              ---Adjust Maximum Displacements

               Call Adjust(Attempt, Nacc, Dr)
            End If
         End Do
         If (Ncycl.Ne.0) Then
            If (Attempt.Ne.0) Write (6, 99003) Attempt, Nacc, 
     &                               100.0d0*Dble(Nacc)/Dble(Attempt)

C           ---Test Total Energy

            Call Toterg(Ent, Virt)
            If (Abs(Ent-En).Gt.1.D-6) Then
               Write (6, *) 
     &                    ' ######### Problems Energy ################ '
            End If
            If (Abs(Virt-Vir).Gt.1.D-6) Then
               Write (6, *) 
     &              ' ######### Problems Virial ################ '
            End If
            Write (6, 99002) Ent, En, Ent - En, Virt, Vir, Virt - Vir
            Write (6,*)

Ccccccccccccccccccccccccccccccccccccccccccccccc
Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Print Chemical Potential And Pressure   C
Ccccccccccccccccccccccccccccccccccccccccccccccc
Ccccccccccccccccccccccccccccccccccccccccccccccc

            If(Ii.Eq.2) Then
               Write (6,*) 'Average Pressure                  : ',
     &              Av1/Av2
               Write (6,*) 'Chemical Potential                : ',
     &              -Log((Bv1/Bv2)*(Box*Box*Box/
     &              Dble(Npart)))/Beta
            End If
         End If
      End Do
      Call Store(21, Dr)
      Stop
 
99001 Format (' Total Energy Initial Configuration: ', F12.5, /, 
     &        ' Total Virial Initial Configuration: ', F12.5)
99002 Format (' Total Energy End Of Simulation    : ', F12.5, /, 
     &        '       Running Energy              : ', F12.5, /, 
     &        '       Difference                  :  ', E12.5, /, 
     &        ' Total Virial End Of Simulation    : ', F12.5, /, 
     &        '       Running Virial              : ', F12.5, /, 
     &        '       Difference                  :  ', E12.5)
99003 Format (' Number Of Att. To Displ. A Part.  : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
      End
