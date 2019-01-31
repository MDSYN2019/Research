      SUBROUTINE MCSWAP(En, Vir, Attempt, Acc)
      Implicit None

C     ---Exchange A Particle Bewteen The Two Boxes
 
      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'
      Include 'chem.inc'
 
      Double Precision En, Vir, Xn, Yn, Zn, Enn, Virn, Eno, Viro, 
     &                 Arg, Vola, Vold, 
     &                 Xo, Yo, Zo, Dele,Ran_Uniform

      Integer Attempt, O, Iadd, Idel, Jb, Idi, Acc
      Dimension En(*), Vir(*)
 
 
      Attempt = Attempt + 1

C     ===Select A Box At Random

      If (Ran_Uniform().Lt.0.5d0) Then
         Iadd = 1
         Idel = 2
      Else
         Iadd = 2
         Idel = 1
      End If

      Vola = Box(Iadd)**3
      Vold = Box(Idel)**3

C     ---Add A Particle To Box Iadd

      Xn = Box(Iadd)*Ran_Uniform()
      Yn = Box(Iadd)*Ran_Uniform()
      Zn = Box(Iadd)*Ran_Uniform()

C     ---Calculate Energy Of This Particle

      Jb = 1
      O  = Npart + 1
      
      Call Eneri(Xn, Yn, Zn, O, Jb, Enn, Virn, Iadd)

C     ---Calculate Contibution To The Chemical Potential:

      Arg = -Beta*Enn
      Chp(Iadd) = Chp(Iadd) + Vola*Exp(Arg)/Dble(Npbox(Iadd)+1)
      If (Npbox(Iadd).Eq.Npart) Chp(Iadd) = Chp(Iadd) + Vola*Exp(Arg)
     &    /Dble(Npbox(Iadd)+1)
      Ichp(Iadd) = Ichp(Iadd) + 1
 

C     ---Delete Particle From Box B:

      If (Npbox(Idel).Eq.0) Then
         Return
      End If
      Idi = 0
      Do While (Idi.Ne.Idel)
         O = Int(Dble(Npart)*Ran_Uniform()) + 1
         Idi = Id(O)
      End Do
      Xo = X(O)
      Yo = Y(O)
      Zo = Z(O)
      Call Eneri(Xo, Yo, Zo, O, Jb, Eno, Viro, Idel)
 
C     ---Acceptence Test:

      Dele = Enn - Eno + Dlog(Vold*Dble((Npbox(Iadd)+1))/
     &     (Vola*Dble(Npbox(Idel))))
     &       /Beta

      If (Ran_Uniform().Lt.Exp(-Beta*Dele)) Then

C        ---Accepted:

         Acc = Acc + 1
         Npbox(Iadd) = Npbox(Iadd) + 1
         X(O) = Xn
         Y(O) = Yn
         Z(O) = Zn
         Id(O) = Iadd
         En(Iadd) = En(Iadd) + Enn
         Vir(Iadd) = Vir(Iadd) + Virn
         Npbox(Idel) = Npbox(Idel) - 1
         En(Idel) = En(Idel) - Eno
         Vir(Idel) = Vir(Idel) - Viro
      End If
      Return
      End
