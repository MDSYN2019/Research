      Subroutine Ener(En, Vir, R2)
C
C     Calculate Energy Between A (Single) Pair Of Atoms
C
C     En : (Output) Energy
C     Vir: (Output) Virial
C     R2 : (Input) Distance Squared Between Two Particles
C
      Implicit None
      Double Precision R2, R2i, R6i, En, Vir

      Include 'potential.inc'
 
      En = 0.0d0
      Vir = 0.0d0

      If (R2.Lt.Rc2) Then
         R2i = Sig2/R2
         R6i = R2i*R2i*R2i
         En  = Eps4*(R6i*R6i-R6i)

C     Start Modification Virial

C     End   Modification Virial

      End If
      
      Return
      End
