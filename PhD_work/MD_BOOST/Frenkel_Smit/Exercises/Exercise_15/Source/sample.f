      Subroutine Sample(Choise,R2,W)
      Implicit None

      Include 'system.inc'

      Integer I,Choise,Maxx
      Parameter(Maxx=500)

      Double Precision Ggt,Ggg(Maxx),Deltax,R2,W

      Save Ggt,Ggg,Deltax

      If(Choise.Eq.1) Then
         Deltax = Dble(Maxx-1)/Dble(Maxchain)
         Ggt    = 0.0d0

         Do I=1,Maxx
            Ggg(I) = 0.0d0
         Enddo
      Elseif(Choise.Eq.2) Then
         I      = 1 + Idint(R2*Deltax)
         Ggg(I) = Ggg(I) + W*R2
         Ggt    = Ggt    + W
      Else
         Ggt    = 1.0d0/Ggt
         Deltax = 1.0d0/Deltax

         Do I=1,Maxx
            Write(23,*) Dble(I-1)*Deltax,Ggg(I)*Ggt
         Enddo
      Endif

      Return
      End
