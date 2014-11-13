PROGRAM gen
  USE omp_lib
  IMPLICIT NONE

  CHARACTER(LEN=256)::program_name, str_number_of_fibers

  INTEGER::N
  REAL*4::min_distance, min_distance_2, avg_distance, step
  REAL*4::domain

  REAL*4,SAVE,DIMENSION(:,:),ALLOCATABLE::positions
  REAL*4,SAVE,DIMENSION(:,:),ALLOCATABLE::orientations

  LOGICAL::optimal
  INTEGER::i,j
  REAL*4,DIMENSION(3)::p_i,p_j
  REAL*4::r
  REAL*4::distance, nearest_distance

  CALL init_random_seed()

  CALL GETARG(0, program_name)
  CALL GETARG(1, str_number_of_fibers)
  READ(str_number_of_fibers,'(I10)') N

  min_distance = 0.2
  min_distance_2 = min_distance**2
  avg_distance = min_distance * 1.7 !0.4

  domain = (N-1)**(1.0/3.0) * avg_distance

  step = 0.1 * min_distance

  ALLOCATE(positions(N,3))
  ALLOCATE(orientations(N,3))

  PRINT *, "N:", N
  PRINT *, "DOMAIN:", -domain, domain

  CALL RANDOM_NUMBER(positions)
  CALL RANDOM_NUMBER(orientations)

  positions = positions * (domain - (-domain)) + (-domain)
  orientations = orientations * (1 - (-1)) + (-1)

  PRINT *, "BEFORE:"
  CALL stats(N, positions)

  optimal = .FALSE.

  DO WHILE (optimal .NEQV. .TRUE.)
    optimal = .TRUE.

    DO i = 1, N

      p_i = positions(i,:)

      ! BOUNDARY
      IF (p_i(1) > domain) THEN
        CALL RANDOM_NUMBER(r)
        p_i(1) = p_i(1) - r * step
        optimal = optimal .AND. .FALSE.
      ELSE IF (p_i(1) < -domain) THEN
        CALL RANDOM_NUMBER(r)
        p_i(1) = p_i(1) + r * step
        optimal = optimal .AND. .FALSE.
      END IF

      IF (p_i(2) > domain) THEN
        CALL RANDOM_NUMBER(r)
        p_i(2) = p_i(2) - r * step
        optimal = optimal .AND. .FALSE.
      ELSE IF (p_i(2) < -domain) THEN
        CALL RANDOM_NUMBER(r)
        p_i(2) = p_i(2) + r * step
        optimal = optimal .AND. .FALSE.
      END IF

      IF (p_i(3) > domain) THEN
        CALL RANDOM_NUMBER(r)
        p_i(3) = p_i(3) - r * step
        optimal = optimal .AND. .FALSE.
      ELSE IF (p_i(3) < -domain) THEN
        CALL RANDOM_NUMBER(r)
        p_i(3) = p_i(3) + r * step
        optimal = optimal .AND. .FALSE.
      END IF

      nearest_distance = 99999.0
      DO j = 1, N

        IF (i /= j) THEN

          p_j = positions(j,:)

          distance = SUM((p_i - p_j)**2)

          IF (distance < nearest_distance) THEN
            nearest_distance = distance
          END IF

        END IF

      END DO

      IF (nearest_distance < min_distance_2) THEN
        CALL RANDOM_NUMBER(r)
        p_i(1) = p_i(1) + r * (step - (-step)) + (-step)
        CALL RANDOM_NUMBER(r)
        p_i(2) = p_i(2) + r * (step - (-step)) + (-step)
        CALL RANDOM_NUMBER(r)
        p_i(3) = p_i(3) + r * (step - (-step)) + (-step)

        optimal = optimal .AND. .FALSE.
      END IF

      positions(i,:) = p_i

    END DO
  END DO

  PRINT *, "AFTER:"
  CALL stats(N, positions)

  OPEN(10,file="XcT_gen"//TRIM(str_number_of_fibers)//".in")
  WRITE(10,'(*(I10))') (N)
  DO i=1,N
    WRITE(10,'(*(F32.16))') (positions(i,:))
    WRITE(10,'(*(F32.16))') (orientations(i,:) / SQRT(DOT_PRODUCT(orientations(i,:),orientations(i,:))))
  END DO
  CLOSE(10)

END PROGRAM gen

SUBROUTINE stats(N, positions)
  IMPLICIT NONE

  INTEGER,INTENT(IN)::N
  REAL*4,INTENT(IN),DIMENSION(N,3)::positions

  INTEGER::i,j,count

  REAL*4::distance, nearest_distance, min_dist, total
  REAL*4,DIMENSION(3)::p_i,p_j

  count = 0
  total = 0.0
  min_dist = 99999.0

  DO i = 1, N

    p_i = positions(i,:)

    nearest_distance = 99999.0
    DO j= 1, N

      IF (i /= j) THEN

        p_j = positions(j,:)

        distance = SQRT(SUM((p_i - p_j)**2))

        IF (distance < nearest_distance) THEN
          nearest_distance = distance
        END IF

      END IF

    END DO

    total = total + nearest_distance
    count = count + 1

    IF (nearest_distance < min_dist) THEN
      min_dist = nearest_distance
    END IF

  END DO

  PRINT *, "MIN:", min_dist
  PRINT *, "AVG:", total/count

END SUBROUTINE stats

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
