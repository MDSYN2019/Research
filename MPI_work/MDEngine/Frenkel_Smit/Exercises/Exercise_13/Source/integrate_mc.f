      Subroutine Integrate_Mc
      Implicit None
 
      Include 'system.inc'
 
Cccccccccccccccccccccc
C     Mc Simulation  C
Cccccccccccccccccccccc

      Double Precision Uold,Unew,F,Xnew,Ran_Uniform
 
      Call Force(Xpos,Uold,F)

      Xnew = Xpos + 2.5d0*(Ran_Uniform()-0.5d0)

      Call Force(Xnew,Unew,F)

      If(Ran_Uniform().Lt.Dexp(-(Unew-Uold)/Temp)) Xpos = Xnew

      Oldf = 0.0d0
      Vpos = 0.0d0
      Cons = 0.0d0
 
      Return
      End
