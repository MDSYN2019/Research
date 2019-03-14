      Double Precision X,Y,Z
      Integer Npart,Id,Npbox

      Common /Conf1/ X(Npmax),Y(Npmax),Z(Npmax),Id(Npmax),Npart
      Common /Conf2/ Npbox(2)
C
C     X(I),Y(I),Z(I)    : Position Particle I
C     Id(I)             : Box 1 Or 2
C     Npart             : Actual Number Of Particles
C     Npbox(I)          : Number Of Particles In Box I
