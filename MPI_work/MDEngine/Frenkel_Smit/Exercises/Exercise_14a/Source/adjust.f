      Subroutine Adjust(Attemp, Nacc, Dr)
      Implicit None
C
C     Adjusts Maximum Displacement Such That 50% Of The
C     Movels Will Be Accepted
C
C  Attemp (Input)  Number Of Attemps That Have Been Performed 
C                  To Displace A Particle
C  Nacc   (Input)  Number Of Successful Attemps To 
C                  Displace A Particle
C  Dr     (Output) New Maximum Displacement
C
 
      Include 'system.inc'
      
      Integer Attemp, Nacc, Attempp, Naccp
      Double Precision Dro, Frac, Dr
      Save Naccp, Attempp
      Data Naccp/0/
      Data Attempp/0/
 
      If (Attemp.Eq.0.Or.Attempp.Ge.Attemp) Then
         Naccp = Nacc
         Attempp = Attemp
      Else
         Frac = Dble(Nacc-Naccp)/Dble(Attemp-Attempp)
         Dro  = Dr
         Dr   = Dr*Abs(Frac/0.5d0)

C        ---Limit The Change:

         If (Dr/Dro.Gt.1.5d0) Dr = Dro*1.5d0
         If (Dr/Dro.Lt.0.5d0) Dr = Dro*0.5d0
         If (Dr.Gt.Hbox/2.D0) Dr = Hbox/2.D0
         Write (6, 99001) Dr, Dro, Frac, Attemp - Attempp, Nacc - Naccp

C        ---Store Nacc And Attemp For Next Use

         Naccp = Nacc
         Attempp = Attemp
      End If
      Return
99001 Format (' Max. Displ. Set To : ', F6.3, ' (Old : ', F6.3, ')', /, 
     &        ' Frac. Acc.: ', F4.2, ' Attempts: ', I7, ' Succes: ', I7)
      End
