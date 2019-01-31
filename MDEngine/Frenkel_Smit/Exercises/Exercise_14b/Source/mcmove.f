      Subroutine Mcmove(En, Vir, Attempt, Nacc, Dr)

C     Attempts To Displace A Randomly Selected Particle

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Double Precision Enn, Eno, En, Xn, Yn, Zn, Viro, Virn, Dr, 
     &                 Vir,Ran_Uniform

      Dimension En(*), Vir(*)
      Integer O, Attempt, Nacc, Jb, Ido
 
      Attempt = Attempt + 1
      Jb = 1

C     ---Select A Particle At Random

      O = Int(Dble(Npart)*Ran_Uniform()) + 1
      Ido = Id(O)

C     ---Calculate Energy Old Configuration

      Call Eneri(X(O), Y(O), Z(O), O, Jb, Eno, Viro, Ido)

C     ---Give Particle A Random Displacement

      Xn = X(O) + (Ran_Uniform()-0.5d0)*Dr
      Yn = Y(O) + (Ran_Uniform()-0.5d0)*Dr
      Zn = Z(O) + (Ran_Uniform()-0.5d0)*Dr

C     ---Calculate Energy New Configuration:

      Call Eneri(Xn, Yn, Zn, O, Jb, Enn, Virn, Ido)

C     ---Acceptance Test

      If (Ran_Uniform().Lt.Exp(-Beta*(Enn-Eno))) Then

C        --Accepted

         Nacc     = Nacc     + 1
         En(Ido)  = En(Ido)  + (Enn-Eno)
         Vir(Ido) = Vir(Ido) + (Virn-Viro)

C        ---Put Particle In Simulation Box

         If (Xn.Lt.0)        Xn = Xn + Box(Ido)
         If (Xn.Gt.Box(Ido)) Xn = Xn - Box(Ido)
         If (Yn.Lt.0)        Yn = Yn + Box(Ido)
         If (Yn.Gt.Box(Ido)) Yn = Yn - Box(Ido)
         If (Zn.Lt.0)        Zn = Zn + Box(Ido)
         If (Zn.Gt.Box(Ido)) Zn = Zn - Box(Ido)
         X(O) = Xn
         Y(O) = Yn
         Z(O) = Zn
      Endif
      Return
      End
