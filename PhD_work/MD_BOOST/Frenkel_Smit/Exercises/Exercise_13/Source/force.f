      Subroutine Force(X,U,F)
      Implicit None
 
C     Calculate The Forces And Potential Energy
 
      Include 'system.inc'
 
      Double Precision X,U,F

      If(X.Le.0.0d0) Then
         U =  2.0d0*Onepi*Onepi*X*X
         F = -4.0d0*Onepi*Onepi*X
      Elseif(X.Ge.1.0d0) Then
         U =  2.0d0*Onepi*Onepi*(X - 1.0d0)*(X - 1.0d0)
         F = -4.0d0*Onepi*Onepi*(X - 1.0d0)
      Else
         U =  1.0d0 - Dcos(2.0d0*Onepi*X)
         F = -2.0d0*Onepi*Dsin(2.0d0*Onepi*X)
      Endif
 
      Return
      End
