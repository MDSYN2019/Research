subroutine cfftb ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTB computes the backward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    CFFTB computes a complex periodic sequence from its Fourier coefficients.
!
!    A call of CFFTF followed by a call of CFFTB will multiply the
!    sequence by N.  In other words, the transforms are not normalized.
!
!    The array WSAVE must be initialized by CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N ) 
!        C_in(K) * exp ( sqrt ( - 1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations,
!    G. Rodrigue, editor, 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling 
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N.  
!
  implicit none
!
  integer n
!
  complex(8) c(n)
  real(8) wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftb1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTB1 is a lower-level routine used by CFFTB.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  complex(8) c(n)
  complex(8) ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nac
  integer nf
  real(8) wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passb4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passb4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passb2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passb2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido
 
      if ( na == 0 ) then
        call passb3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passb3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passb5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passb5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passb ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passb ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftf ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTF computes the forward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    CFFTF computes the Fourier coefficients of a complex periodic sequence.
!
!    The transform is not normalized.  To obtain a normalized transform,
!    the output must be divided by N.  Otherwise a call of CFFTF
!    followed by a call of CFFTB will multiply the sequence by N.
!
!    The array WSAVE must be initialized by calling CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N ) 
!        C_in(K) * exp ( - sqrt ( -1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations,
!    G. Rodrigue, editor, 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling 
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N. 
!
  implicit none
!
  integer n
!
  complex(8) c(n)
  real(8) wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTF1 is a lower level routine used by CFFTF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  complex(8) c(n)
  complex(8) ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nac
  integer nf
  real(8) wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
 
      if ( na == 0 ) then
        call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passf4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passf2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passf2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cffti ( n, wsave )
!
!*******************************************************************************
!
!! CFFTI initializes WSAVE, used in CFFTF and CFFTB. 
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the 
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber, 
!    Vectorizing the FFT's, 
!    in Parallel Computations,
!    G. Rodrigue, editor, 
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee, 
!    The SLATEC Common Math Library, 
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Output, real WSAVE(4*N+15), contains data, dependent on the value
!    of N, which is necessary for the CFFTF or CFFTB routines.  
!
  implicit none
!
  integer n
!
  real(8) wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cffti1 ( n, wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! CFFTI1 is a lower level routine used by CFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  real(8) arg
  real(8) argh
  real(8) argld
  real(8) fi
  integer i
  integer i1
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer j
  integer k1
  integer l1
  integer l2
  integer ld
  integer nf
  real(8) r_pi
  real(8) wa(2*n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * r_pi ( ) / dble ( n )
  i = 2
  l1 = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      i1 = i
      wa(i-1) = 1.0D+00
      wa(i) = 0.0D+00
      ld = ld + l1
      fi = 0.0D+00
      argld = dble ( ld ) * argh

      do ii = 4, 2*ido+2, 2
        i = i + 2
        fi = fi + 1.0D+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      if ( 5 < ip ) then
        wa(i1-1) = wa(i-1)
        wa(i1) = wa(i)
      end if

    end do

    l1 = l2

  end do

  return
end
subroutine i_factor ( n, ifac )
!
!*******************************************************************************
!
!! I_FACTOR factors an integer.
!
!
!  Modified:
!
!    14 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the number to be factored.
!
!    Output, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer i
  integer ifac(15)
  integer j
  integer n
  integer nf
  integer nl
  integer nq
  integer nr
  integer ntry
!
  ifac(1) = n

  nf = 0
  nl = n

  if ( n == 0 ) then
    nf = 1
    ifac(2) = nf
    ifac(2+nf) = 0
    return
  end if

  if ( n < 1 ) then
    nf = nf + 1
    ifac(2+nf) = -1
    nl = - n
  end if

  if ( nl == 1 ) then
    nf = nf + 1
    ifac(2) = nf
    ifac(2+nf) = 1
    return
  end if

  j = 0

  do while ( 1 < nl )

    j = j + 1
!
!  Choose a trial divisor, NTRY.
!
    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if
!
!  Divide by the divisor as many times as possible.
!
    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nl = nq
      nf = nf + 1
!
!  Make sure factors of 2 appear in the front of the list.
!
      if ( ntry /= 2 ) then

        ifac(2+nf) = ntry

      else

        do i = nf, 2, -1
          ifac(i+2) = ifac(i+1)
        end do
        ifac(3) = 2

      end if

    end do

  end do

  ifac(2) = nf

  return
end
subroutine passb ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSB is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real(8) c1(ido,l1,ip)
  real(8) c2(idl1,ip)
  real(8) cc(ido,ip,l1)
  real(8) ch(ido,l1,ip)
  real(8) ch2(idl1,ip)
  integer i
  integer idij
  integer idj
  integer idl
  integer idlj
  integer idp
  integer ik
  integer inc
  integer ipph
  integer j
  integer jc
  integer k
  integer l
  integer lc
  integer nac
  integer nt
  real(8) wa(*)
  real(8) wai
  real(8) war
!
  nt = ip * idl1
  ipph = ( ip + 1 ) / 2
  idp = ip * ido

  if ( l1 <= ido ) then

    do j = 2, ipph
      jc = ip + 2 - j
      do k = 1, l1
        ch(1:ido,k,j)  = cc(1:ido,j,k) + cc(1:ido,jc,k)
        ch(1:ido,k,jc) = cc(1:ido,j,k) - cc(1:ido,jc,k)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      do i = 1, ido
        ch(i,1:l1,j)  = cc(i,j,1:l1) + cc(i,jc,1:l1)
        ch(i,1:l1,jc) = cc(i,j,1:l1) - cc(i,jc,1:l1)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l) = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =            wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j
      idlj = idlj + inc
      if ( idp < idlj ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  nac = 1

  if ( ido == 2 ) then
    return
  end if

  nac = 0
  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1:2,1:l1,2:ip) = ch(1:2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSB2 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,2,l1)
  real(8) ch(ido,l1,2)
  integer i
  integer k
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   + cc(i,2,k)
        ti2         = cc(i,1,k)   - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 + wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 - wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSB3 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,3,l1)
  real(8) ch(ido,l1,3)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer i
  integer k
  real(8) taui
  real(8), parameter :: taur = -0.5D+00
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
  real(8) wa2(ido)
!
  taui = sqrt ( 3.0D+00 ) / 2.0D+00

  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k) - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSB4 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,4,l1)
  real(8) ch(ido,l1,4)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer i
  integer k
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,4,k) - cc(2,2,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,2,k) - cc(1,4,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k) - cc(i,3,k)
        ti2 = cc(i,1,k) + cc(i,3,k)
        ti3 = cc(i,2,k) + cc(i,4,k)
        tr4 = cc(i,4,k) - cc(i,2,k)

        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,2,k) - cc(i-1,4,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i,k,1) = ti2 + ti3
        ci3 = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 - wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 + wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 - wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 + wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 - wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 + wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSB5 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,5,l1)
  real(8) ch(ido,l1,5)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer i
  integer k
  real(8), parameter :: ti11 = 0.951056516295154D+00
  real(8), parameter :: ti12 = 0.587785252292473D+00
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8), parameter :: tr11 = 0.309016994374947D+00
  real(8), parameter :: tr12 = -0.809016994374947D+00
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
  real(8) wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 - wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 + wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 - wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 + wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine passf ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSF is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real(8) c1(ido,l1,ip)
  real(8) c2(idl1,ip)
  real(8) cc(ido,ip,l1)
  real(8) ch(ido,l1,ip)
  real(8) ch2(idl1,ip)
  integer i
  integer idij
  integer idj
  integer idl
  integer idlj
  integer idp
  integer ik
  integer inc
  integer ipph
  integer j
  integer jc
  integer k
  integer l
  integer lc
  integer nac
  integer nt
  real(8) wa(*)
  real(8) wai
  real(8) war
!
  nt = ip * idl1
  ipph = (ip+1) / 2
  idp = ip * ido

  if ( l1 <= ido ) then

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j

      idlj = idlj + inc
      if ( idp < idlj ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  if ( ido == 2 ) then
    nac = 1
    return
  end if

  nac = 0

  c2(1:idl1,1)    = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
  c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSF2 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,2,l1)
  real(8) ch(ido,l1,2)
  integer i
  integer k
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)

        ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
        ti2       = cc(i,1,k) - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSF3 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,3,l1)
  real(8) ch(ido,l1,3)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer i
  integer k
  real(8) taui
  real(8), parameter :: taur = -0.5D+00
  real(8) ti2
  real(8) tr2
  real(8) wa1(ido)
  real(8) wa2(ido)
!
  taui = - sqrt ( 3.0D+00 ) / 2.0D+00

  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSF4 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,4,l1)
  real(8) ch(ido,l1,4)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer i
  integer k
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,2,k) - cc(2,4,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,4,k) - cc(1,2,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k)   - cc(i,3,k)
        ti2 = cc(i,1,k)   + cc(i,3,k)
        ti3 = cc(i,2,k)   + cc(i,4,k)
        tr4 = cc(i,2,k)   - cc(i,4,k)
        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,4,k) - cc(i-1,2,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSF5 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real(8) cc(ido,5,l1)
  real(8) ch(ido,l1,5)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer i
  integer k
  real(8), parameter :: ti11 = -0.951056516295154D+00
  real(8), parameter :: ti12 = -0.587785252292473D+00
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
!
!  cos ( 72 ) = +0.3090
!
  real(8), parameter :: tr11 =  0.309016994374947D+00
!
!  cos ( 36 ) = +0.809016
!
  real(8), parameter :: tr12 = -0.809016994374947D+00
  real(8) wa1(ido)
  real(8) wa2(ido)
  real(8) wa3(ido)
  real(8) wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

      end do
    end do

  end if

  return
end
function r_pi ()
!
!*******************************************************************************
!
!! R_PI returns the value of pi.
!
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_PI, the value of PI.
!
  implicit none
!
  real(8) r_pi
!
  r_pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
