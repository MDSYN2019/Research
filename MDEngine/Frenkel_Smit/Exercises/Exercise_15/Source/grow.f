      Subroutine Grow(Lold,Weight,Ubonded,Unonb,Xst,Yst,Zst)
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Grow A Chain Using Cbmc Or Random Insertion  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccc

      Include 'system.inc'

      Logical Lold,Lready
      Integer I,J,K,Ntrial,Iwalk

      Double Precision Ubonded,Unonb,Weight,Xst(Maxchain),
     &     Yst(Maxchain),Zst(Maxchain),Xtrial(Maxtrial),
     &     Ytrial(Maxtrial),Ztrial(Maxtrial),Dx1,Dx2,Dy1,
     &     Dy2,Dz1,Dz2,Ubend(Maxtrial),Uext(Maxtrial),X,Y,
     &     Z,Ws,Ran_Uniform,Cumw,Cmmm,U,Bfac(Maxtrial)

      Ubonded = 0.0d0
      Unonb   = 0.0d0
      Weight  = 1.0d0
      Ws      = 0.0d0

      If(Lcbmc) Then
         Ntrial = Nchoi
      Else
         Ntrial = 1
      Endif

      Xst(1) = 0.0d0
      Yst(1) = 0.0d0
      Zst(1) = 0.0d0

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Grow The Chain; Loop Over All Segments Exept The First One  C
C                                                                 C
C     Xpos/Ypos/Zpos       : Position Of The Old Chain            C
C     Xst/Yst/Zst          : Position Of The Trial Chain          C
C     Xtrial/Ytrial/Ztrial : Position Of A Trial Segment          C
C     Lcbmc                : Do We Use Cbmc (.True. Or .False.)   C
C     Lold                 : Do We Have The Old Configuration     C
C                            (Only For Dynamic Schemes)           C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Do I=2,Nuall

Ccccccccccccccccccccccccccccccccccccccc
C     Loop Over All Trial Positions   C
Ccccccccccccccccccccccccccccccccccccccc

         Do J=1,Ntrial

 10         Ubend(J) = 0.0d0
            Uext(J)  = 0.0d0
            Lready   = .False.
            Ws       = 0.0d0
            
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Trial Position                              C
C     J=1 And Lold : Old Position                          C
C     Otherwise    : New Position                          C
C                                                          C
C     Old Coordinates Are Stored In Xpos/Ypos/Zpos         C
C     The New Chain Is Stored In Xst/Yst/Zst               C
C     Trial Positions Are Stored In Xtrial/Ytrial/Ztrial   C
C                                                          C
C     Subroutine Ran_Sphere(X,Y,Z): Creates A Random Unit  C
C     Vector On A Sphere                                   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

C     End Modification

Cccccccccccccccccccccccccccccccccccccccccc
C     Calculate Bond-Bending Potential   C
Cccccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.2) Then
               Dx1 = Xtrial(J) - Xst(I-1)
               Dy1 = Ytrial(J) - Yst(I-1)
               Dz1 = Ztrial(J) - Zst(I-1)
                                         
               Dx2 = Xst(I-2)  - Xst(I-1)
               Dy2 = Yst(I-2)  - Yst(I-1)
               Dz2 = Zst(I-2)  - Zst(I-1)

               Call Angle(Dx1,Dy1,Dz1,Dx2,Dy2,Dz2,U)

               Ubend(J) = U
            Endif

Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Accept Or Reject The Chosen Configuration      C
C     This Is For Cmbc Only, Because Only In Cbmc    C
C     A Trial Position According To A Bond-Bending   C
C     Potential Has To Be Generated...               C
C                                                    C
C     When J=1.And.Lold : Always Accepted !! (Why?)  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc

            If(Lcbmc) Then

C     Start Modification

C     End Modification

            Else
               Lready = .True.
            Endif

            If(.Not.Lready) Goto 10

Ccccccccccccccccccccccccccccccccccccccc
C     Calculate The External Energy   C
C     Only When There Are More Than 4 C
C     Beads...                        C
Ccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.3) Then
               Do K=1,(I-3)

                  Dx1 = Xtrial(J) - Xst(K)
                  Dy1 = Ytrial(J) - Yst(K)
                  Dz1 = Ztrial(J) - Zst(K)
                  
                  Call Repulsion(Dx1,Dy1,Dz1,U)

                  Uext(J) = Uext(J) + U
               Enddo
            Endif
         Enddo

Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     End Loop Over All Trial Positions           C
C     Choose A Trial Position For Cbmc            C
C     For Cbmc And An Old Config.: 1st Config.    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

         Iwalk = 1

         If(Lcbmc) Then
            Do K=1,Ntrial
               Bfac(K) = Dexp(-Beta*Uext(K))
               Ws      = Ws + Bfac(K)
            Enddo

            Cumw  = Bfac(1)
            Cmmm  = Ws*Ran_Uniform()

            If(.Not.Lold) Then
               Do While(Cumw.Lt.Cmmm)
                  Iwalk = Iwalk + 1
                  Cumw  = Cumw  + Bfac(Iwalk)
               Enddo
            Endif
         Else
            Ws = Dexp(-Beta*(Uext(1)+Ubend(1)))
         Endif

Ccccccccccccccccccccccccccccccccccccccccccc
C     Store The Chosen Configuration      C
C     Update Energies/Rosenbluth Weight   C
Ccccccccccccccccccccccccccccccccccccccccccc

         Xst(I) = Xtrial(Iwalk)
         Yst(I) = Ytrial(Iwalk)
         Zst(I) = Ztrial(Iwalk)

         Ubonded     = Ubonded + Ubend(Iwalk)
         Unonb       = Unonb   + Uext(Iwalk)
         Weight      = Weight*Ws/Dble(Ntrial)
      Enddo

      Return
      End
