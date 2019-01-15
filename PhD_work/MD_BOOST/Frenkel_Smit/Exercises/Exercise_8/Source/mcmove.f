      Subroutine Mcmove(En, Vir, Attempt, Nacc, Dr)

C
C     Attempts To Displace A Randomly Selected Particle
C
C
C     Ener   (Input/Output) : Total Energy
C     Vir    (Input/Output) : Total Virial
C     Attemp (Input/Output) : Number Of Attemps That Have Been
C                             Performed To Displace A Particle
C     Nacc   (Input/Output) : Number Of Successful Attemps
C                             To Displace A Particle
C     Dr     (Input)        : Maximum Displacement
C
      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Double Precision Enn, Eno, En, Ran_Uniform, Xn, Yn, 
     &     Zn, Viro, Virn, Vir, Dr
      Integer O, Attempt, Nacc, Jb
 
      Attempt = Attempt + 1
      Jb = 1

C     ---Select A Particle At Random

      O = Int(Dble(Npart)*Ran_Uniform()) + 1

C     ---Calculate Energy Old Configuration

      Call Eneri(X(O), Y(O), Z(O), O, Jb, Eno, Viro)

C     ---Give Particle A Random Displacement

      Xn = X(O) + (Ran_Uniform()-0.5d0)*Dr
      Yn = Y(O) + (Ran_Uniform()-0.5d0)*Dr
      Zn = Z(O) + (Ran_Uniform()-0.5d0)*Dr

C     ---Calculate Energy New Configuration:

      Call Eneri(Xn, Yn, Zn, O, Jb, Enn, Virn)

C     ---Acceptance Test

      If (Ran_Uniform().Lt.Exp(-Beta*(Enn-Eno))) Then

C        --Accepted

         Nacc = Nacc + 1
         En = En + (Enn-Eno)
         Vir = Vir + (Virn-Viro)

C        ---Put Particle In Simulation Box

         If (Xn.Lt.0.0d0) Then
            Xn = Xn + Box
         Elseif (Xn.Gt.Box) Then
            Xn = Xn - Box
         Endif

         If (Yn.Lt.0.0d0) Then
            Yn = Yn + Box
         Elseif (Yn.Gt.Box) Then
            Yn = Yn - Box
         Endif

         If (Zn.Lt.0.0d0) Then
            Zn = Zn + Box
         Elseif (Zn.Gt.Box) Then
            Zn = Zn - Box
         Endif

         X(O) = Xn
         Y(O) = Yn
         Z(O) = Zn
      End If
      Return
      End
