      Subroutine Init_Chem(Switch)

C     ---Initialize And Calculate Chemical Potentials
C

      Implicit None

      Integer Switch, Ib

      Include 'chem.inc'
      Include 'system.inc'

      If (Switch.Eq.0) Then
C        ---Initialize
         Do Ib = 1, 2
            Chp(Ib) = 0.0d0
            Ichp(Ib) = 0
         End Do
      Else If (Switch.Eq.2) Then
C        ---Print Final Results
         Do Ib = 1, 2
            If (Ichp(Ib).Ne.0) Then
               Chp(Ib) = -Dlog(Chp(Ib)/Dble(Ichp(Ib)))/Beta
            End If
         End Do
         Write (6, 99001) (Ichp(1)+Ichp(2))/2, Chp(1), Chp(2)
      Else
         Stop 'Error: Init_Chem'
      End If
99001 Format (' Chemical Potentials : ', /, ' Number Of Samples : ', 
     &        I12, /, ' Box 1 ', F7.3, /, ' Box 2 ', F7.3, /)
      Return
      End
