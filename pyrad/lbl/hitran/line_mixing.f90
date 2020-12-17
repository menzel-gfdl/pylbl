!Calculate the relaxation matrix.
subroutine create_relaxation_matrix(nlines, temperature, iso, li, lf, gamma_hwhm, &
                                    population, dk0, ji, jf, &
                                    b0pp, b0pq, b0pr, &
                                    b0qp, b0qq, b0qr, &
                                    b0rp, b0rq, b0rr, &
                                    w0pp, w0pq, w0pr, &
                                    w0qp, w0qq, w0qr, &
                                    w0rp, w0rq, w0rr, w)
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  integer(kind=int32), intent(in) :: nlines
  real(kind=real64), intent(in) :: temperature
  integer(kind=int32), intent(in) :: iso, lf, li
  integer(kind=int32), dimension(nlines), intent(in), target :: jf, ji
  real(kind=real64), dimension(nlines), intent(in) :: dk0, gamma_hwhm, population
  real(kind=real64), dimension(0:9,0:9,0:130,0:130), intent(in) :: b0pp, b0pq, b0pr, &
                                                                   b0qp, b0qq, b0qr, &
                                                                   b0rp, b0rq, b0rr, &
                                                                   w0pp, w0pq, w0pr, &
                                                                   w0qp, w0qq, w0qr, &
                                                                   w0rp, w0rq, w0rr
  real(kind=real64), dimension(nlines,nlines), intent(inout) :: w

  integer(kind=int32) :: i, j, lf_, li_
  integer(kind=int32), dimension(:), pointer :: jf_, ji_
  logical :: iso_flag
  real(kind=real64) :: b0, logt, sumlw, sumup, w0, ycal
  real(kind=real64), parameter :: t0 = 296._real64

  w(:,:) = 0._real64
  logt = log(t0/temperature)
  iso_flag = iso .gt. 2 .and. iso .ne. 7 .and. iso .ne. 10

  !Define the transition so it is "downward".
  li_ = min(li, lf)
  lf_ = max(li, lf)
  if (li .le. lf) then
    ji_ => ji
    jf_ => jf
  else
    ji_ => jf
    jf_ => ji
  endif

  do i = 1, nlines
    do j = 1, nlines
      if (ji_(j) .gt. ji_(i)) cycle
      if (iso_flag .and. mod(abs(ji(i) - ji(j)), 2) .ne. 0) cycle
      if (ji_(i) .gt. jf_(i) .and. ji_(j) .gt. jf_(j)) then
        w0 = w0pp(li_,lf_,ji_(i),ji_(j))
        b0 = b0pp(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .gt. jf_(i) .and. ji_(j) .eq. jf_(j))then
        w0 = w0pq(li_,lf_,ji_(i),ji_(j))
        b0 = b0pq(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .gt. jf_(i) .and. ji_(j) .lt. jf_(j))then
        w0 = w0pr(li_,lf_,ji_(i),ji_(j))
        b0 = b0pr(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .lt. jf_(i) .and. ji_(j) .gt. jf_(j))then
        w0 = w0rp(li_,lf_,ji_(i),ji_(j))
        b0 = b0rp(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .lt. jf_(i) .and. ji_(j) .eq. jf_(j))then
        w0 = w0rq(li_,lf_,ji_(i),ji_(j))
        b0 = b0rq(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .lt. jf_(i) .and. ji_(j) .lt. jf_(j))then
        w0 = w0rr(li_,lf_,ji_(i),ji_(j))
        b0 = b0rr(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .eq. jf_(i) .and. ji_(j) .gt. jf_(j))then
        w0 = w0qp(li_,lf_,ji_(i),ji_(j))
        b0 = b0qp(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .eq. jf_(i) .and. ji_(j) .eq. jf_(j))then
        w0 = w0qq(li_,lf_,ji_(i),ji_(j))
        b0 = b0qq(li_,lf_,ji_(i),ji_(j))
      elseif(ji_(i) .eq. jf_(i) .and. ji_(j) .lt. jf_(j))then
        w0 = w0qr(li_,lf_,ji_(i),ji_(j))
        b0 = b0qr(li_,lf_,ji_(i),ji_(j))
      endif
      ycal = exp(w0 - b0*logt)
      w(j,i) = ycal
      w(i,j) = ycal*population(i)/population(j)
    enddo
  enddo

  do i = 1, nlines
    do j = 1, nlines
      if (i .ne. j) then
        w(i,j) = -abs(w(i,j))
      endif
    enddo
  enddo

  do i = 1, nlines
    w(i,i) = gamma_hwhm(i)
  enddo

  do i = 1, nlines
    sumlw = 0._real64
    sumup = 0._real64
    do j = 1, nlines
      if (iso_flag .and. mod(abs(ji(i) - ji(j)), 2) .ne. 0) cycle
      if (j .gt. i) then
        sumlw = sumlw + abs(dk0(j))*w(j,i)
      else
        sumup = sumup + abs(dk0(j))*w(j,i)
      endif
    enddo
    do j = i + 1, nlines
      if (sumlw .eq. 0._real64) then
        w(j,i) = 0._real64
        w(i,j) = 0._real64
      else
        w(j,i) = w(j,i)*(-sumup/sumlw)
        w(i,j) = w(j,i)*population(i)/population(j)
      endif
    enddo
  enddo
end subroutine create_relaxation_matrix


!Calculate the first-order line-mixing coefficients.
subroutine calculate_coefficients(nlines, iso, dipole, ji, v, w, y)
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  integer(kind=int32), intent(in) :: iso, nlines
  integer(kind=int32), dimension(nlines), intent(in) :: ji
  real(kind=real64), dimension(nlines), intent(in) :: dipole, v
  real(kind=real64), dimension(nlines,nlines), intent(in) :: w
  real(kind=real64), dimension(nlines), intent(inout) :: y

  logical :: iso_flag
  integer(kind=int32) :: i, j
  real(kind=real64) :: dv

  iso_flag = iso .gt. 2 .and. iso .ne. 7 .and. iso .ne. 10
  do i = 1, nlines
    y(i) = 0._real64
    do j = 1, nlines
      if (j .eq. i) cycle
      if (iso_flag .and. mod(abs(ji(i) - ji(j)), 2) .ne. 0) cycle
      dv = v(i) - v(j)
      if (abs(dv) .lt. 1.d-4) dv = 1.d-4
      y(i) = y(i) + 2._real64*abs(dipole(j))/abs(dipole(i))*(1._real64/dv)*w(j,i)
    enddo
  enddo
end subroutine calculate_coefficients
