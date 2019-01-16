      Program Gibbs
      Implicit None

C     ---Gibbs-Ensemble Simulation Of The Lennard-Joned Fluid

      Integer Equil, Prod, Nsamp, Ii, Icycl, Ndispl, Attempt, 
     &        Nacc, Ncycl, Nmoves, Imove, Nvol, Accv, Attv, Ib, Nswap, 
     &        Accsw, Attsw, Sstmm
      Double Precision En(2), Ent(2), Vir(2), Virt(2), Vmax, Dr, 
     &                 Ran, Succ,Ran_Uniform,M1

      Include 'parameter.inc'
      Include 'chem.inc'
      Include 'conf.inc'
      Include 'potential.inc'
      Include 'system.inc'
 
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

      Write (6, *) '**************** GIBBS ***************'

C     ---Initialize Sysem

      Call Readdat(Equil, Prod, Nsamp, Ndispl, Dr, Nvol, Vmax, Nswap, 
     &             Succ)
      Nmoves = Ndispl + Nvol + Nswap

C     ---Total Energy Of The System

      Do Ib = 1, 2
         Call Toterg(En(Ib), Vir(Ib), Ib)
         Write (6, 99001) Ib, En(Ib), Vir(Ib)
      End Do

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
         Nacc = 0
         Attv = 0
         Accv = 0
         Attsw = 0
         Accsw = 0

C        ---Initialize Calculation Chemical Potential

         Call Init_Chem(0)

C        ---Intialize The Subroutine That Adjust The Maximum Displacement

         Call Adjust(Attempt, Nacc, Dr, Attv, Accv, Vmax, Succ)

         Do Icycl = 1, Ncycl
            Do Imove = 1, Nmoves
               Ran = Ran_Uniform()*Dble(Ndispl+Nvol+Nswap)
               If (Ran.Lt.Dble(Ndispl)) Then

C                 ---Attempt To Displace A Particle

                  Call Mcmove(En, Vir, Attempt, Nacc, Dr)
               Else If (Ran.Lt.Dble(Ndispl+Nvol)) Then

C                 ---Attempt To Change The Volume

                  Call Mcvol(En, Vir, Attv, Accv, Vmax)
               Else

C                 ---Attemp To Exchange Particles

                  Call Mcswap(En, Vir, Attsw, Accsw)
 
               End If
            End Do
            If (Ii.Eq.2) Then

C              ---Sample Averages

               If (Mod(Icycl,Nsamp).Eq.0) Call Sample(Icycl, En, Vir)
            End If
            If (Mod(Icycl,Ncycl/5).Eq.0) Then
               Write (6, *) '======>> Done ', Icycl, ' Out Of ', Ncycl, 
     &                      Npbox(1), Npbox(2)

C              ---Write Intermediate Configuration To File

               Call Store(8, Dr, Vmax)

C              ---Adjust Maximum Displacements

               Call Adjust(Attempt, Nacc, Dr, Attv, Accv, Vmax, Succ)
            End If
         End Do
         If (Ncycl.Ne.0) Then
            If (Attempt.Ne.0) Write (6, 99003) Attempt, Nacc, 
     &                               100.*Float(Nacc)/Float(Attempt)
            If (Attv.Ne.0) Write (6, 99004) Attv, Accv, 100.*Float(Accv)
     &                            /Float(Attv)
            If (Attsw.Ne.0) Write (6, 99005) Attsw, Accsw, 
     &                             100.*Float(Accsw)/Float(Attsw)
            Do Ib = 1, 2

C              ---Test Total Energy

               Call Toterg(Ent(Ib), Virt(Ib), Ib)
               If (Abs(Ent(Ib)-En(Ib)).Gt.1.D-6) Then
                  Write (6, *) 
     &                    ' ######### Problems Energy ################ '
               End If
               If (Abs(Virt(Ib)-Vir(Ib)).Gt.1.D-6) Then
                  Write (6, *) 
     &                    ' ######### Problems Virial ################ '
               End If
               Write (6, 99002) Ib, Ent(Ib), En(Ib), Ent(Ib) - En(Ib), 
     &                          Virt(Ib), Vir(Ib), Virt(Ib) - Vir(Ib)
            End Do
C           ---Calculation Chemical Potential
            Call Init_Chem(2)
         End If
      End Do
      Call Store(21, Dr, Vmax)

 
99001 Format (' Box : ', I3, /, 
     &        ' Total Energy Initial Configuration : ', F12.5, /, 
     &        ' Total Virial Initial Configuration : ', F12.5)
99002 Format (' Box : ', I3, /, ' Total Energy End Of Simulation   : ', 
     &        F12.5, /, '       Running Energy             : ', F12.5, 
     &        /, '       Difference                 :  ', E12.5, /, 
     &        ' Total Virial End Of Simulation   : ', F12.5, /, 
     &        '       Running Virial             : ', F12.5, /, 
     &        '       Difference                 :  ', E12.5)
99003 Format (' Number Of Att. To Displ. A Part.  : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
99004 Format (' Number Of Att. To Chan. Volume    : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
99005 Format (' Number Of Att. To Exchange Part.  : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
      
      Stop
      End
