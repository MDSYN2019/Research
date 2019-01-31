      Subroutine Toterg(Ener, Vir)
C
C     Calculates Total Energy Of The System
C     Only Used In The Beginning Or At The End Of The Program
C
C     Ener (Output) : Total Energy
C     Vir  (Output) : Total Virial
C
 
      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
 
      Double Precision Xi, Yi, Zi, Ener, Eni, Viri, Vir
      Integer I, Jb
 
      Ener = 0.0d0
      Vir  = 0.0d0

      Do I = 1, Npart - 1
         Xi = X(I)
         Yi = Y(I)
         Zi = Z(I)
         Jb = I + 1

         Call Eneri(Xi, Yi, Zi, I, Jb, Eni, Viri)

         Ener = Ener + Eni
         Vir  = Vir + Viri
      End Do

      Return
      End
