Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Random Number Generator In Fortran77,      C 
C     Can Be Used To Replace ran_uniform.c       C
C                                                C
C     Thijs J.H. Vlugt, March 6 2000             C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine Genrand(X)
      Implicit None

      Integer I,J
      Double Precision X,Y,Ran_Uniform

      I = 1 + Idint(Dabs(X)*100.0d0)

      Do J=1,I
         Y = Ran_Uniform()
      Enddo

      Return
      End

      Function Ran_Uniform()
      Implicit None

      Integer Idum
      Double Precision Ran_Uniform,Ran2

      Save Idum
      Data Idum/-10/

      Ran_Uniform = Ran2(Idum)

      Return
      End

      Function Ran2(Idum)
      Implicit None

      Integer Idum,Im1,Im2,Imm1,Ia1,Ia2,
     &     Iq1,Iq2,Ir1,Ir2,Ntab,Ndiv

      Double Precision Ran2,Am,Eps,Rnmx

      Parameter (Im1=2147483563,Im2=2147483399,
     &     Am=1.0d0/Im1,Imm1=Im1-1,Ia1=40014,
     &     Ia2=40692,Iq1=53668,Iq2=52774,Ir1=12211,
     &     Ir2=3791,Ntab=32,Ndiv=1+Imm1/Ntab,
     &     Eps=1.2d-7,Rnmx=1.0d0-Eps)

      Integer Idum2,J,K,Iv(Ntab),Iy

      Save Iv,Iy,Idum2
      Data Idum2/123456789/, Iv/Ntab*0/, Iy/0/

      If (Idum.Le.0) Then
         Idum=Max(-Idum,1)
         Idum2=Idum
         Do J=Ntab+8,1,-1
            K=Idum/Iq1
            Idum=Ia1*(Idum-K*Iq1)-K*Ir1
            If (Idum.Lt.0) Idum=Idum+Im1
            If (J.Le.Ntab) Iv(J)=Idum
         Enddo
         Iy=Iv(1)
      Endif

      K=Idum/Iq1
      Idum=Ia1*(Idum-K*Iq1)-K*Ir1
      If (Idum.Lt.0) Idum=Idum+Im1
      K=Idum2/Iq2
      Idum2=Ia2*(Idum2-K*Iq2)-K*Ir2
      If (Idum2.Lt.0) Idum2=Idum2+Im2
      J=1+Iy/Ndiv
      Iy=Iv(J)-Idum2
      Iv(J)=Idum
      If(Iy.Lt.1)Iy=Iy+Imm1
      Ran2=Min(Am*Iy,Rnmx)

      Return
      End
