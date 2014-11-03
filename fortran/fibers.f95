#include "constants.incl"

#define PLUS_EQUALS(A, B) A = A + B
#define MINUS_EQUALS(A, B) A = A - B
#define MULTIPLY_EQUALS(A, B) A = A * B
#define DIVIDE_EQUALS(A, B) A = A / B

PROGRAM fibers
  USE omp_lib
  IMPLICIT NONE

  CHARACTER(LEN=256)::program_name, data_name

  INTEGER total_count_rate,total_count_max,total_count1,total_count2
  REAL*8 total_CPU_p
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p

  !--------------------------------------------------
  ! Initalize Memory
  !--------------------------------------------------
  REAL*4,SAVE,DIMENSION(DIMENSIONS*NUMBER_OF_FIBERS)::previous_positions, current_positions, next_positions
  REAL*4,SAVE,DIMENSION(DIMENSIONS*NUMBER_OF_FIBERS)::previous_orientations, current_orientations, next_orientations
  REAL*4,SAVE,DIMENSION(DIMENSIONS*NUMBER_OF_FIBERS)::previous_translational_velocities, current_translational_velocities
  REAL*4,SAVE,DIMENSION(DIMENSIONS*NUMBER_OF_FIBERS)::previous_rotational_velocities, current_rotational_velocities

  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_ROWS,TOTAL_NUMBER_OF_ROWS)::a_matrix
  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_ROWS)::b_vector

  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_points
  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_weights
  REAL*4,SAVE,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS,NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::legendre_polynomials

  REAL*4,SAVE,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::lambda
  REAL*4,SAVE,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::eigen

  REAL*4,SAVE,DIMENSION(DIMENSIONS)::external_force

  INTEGER::i,j,force_index,force_index_i, force_index_j,quadrature_index_i,quadrature_index_j
  REAL*4,DIMENSION(DIMENSIONS)::position_i, orientation_i
  REAL*4,DIMENSION(DIMENSIONS)::position_j, orientation_j
  REAL*4::quadrature_weight, legendre_polynomial

  REAL*4,DIMENSION(6)::T
  REAL*4,DIMENSION(3)::Q
  REAL*4,DIMENSION(3)::TF, TFA0, TFA1, TFA1_TMP
  REAL*4::QF

  REAL*4,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS, 6)::G
  REAL*4,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS, 3)::GF
  REAL*4,DIMENSION(DIMENSIONS)::position_on_fiber_i,position_on_fiber_j,difference,difference2, oriented_force
  REAL*4::invDistance,invDistance3,invDistance5
  REAL*4,DIMENSION(6)::K
  REAL*4,DIMENSION(DIMENSIONS)::force_on_fiber_j

  INTEGER::x_row_index,y_row_index,z_row_index
  INTEGER::x_col_index,y_col_index,z_col_index
  REAL*4::c,d,e,cc,D1,gamma

  INTEGER,SAVE,DIMENSION(TOTAL_NUMBER_OF_ROWS)::IPIV
  INTEGER::INFO

  INTEGER::IDUMMY,ind

  c = LOG(SLENDERNESS * SLENDERNESS * EXP(1.0))
  d = -c
  e = 2.0
  cc = 1.0
  D1 = 0.75 / (d - 2.0 * cc)

  a_matrix = 0.0
  b_vector = 0.0
  external_force = (/ 0.0, 0.0, -0.5 /)

  !--------------------------------------------------
  ! Load positions and orientations
  !--------------------------------------------------
  !me:  Intializes the positons and @todo tVecs with data from the specified
  !     input file.
  CALL GETARG(0, program_name)
	CALL GETARG(1, data_name)

  PRINT *, TRIM(data_name)
  PRINT *, NUMBER_OF_FIBERS

  OPEN(10,file=TRIM(data_name));
  READ(10,*) IDUMMY
  DO i=1,NUMBER_OF_FIBERS
     ind=(i-1)*3
     READ(10,*) current_positions(ind+1),current_positions(ind+2),current_positions(ind+3),current_orientations(ind+1),current_orientations(ind+2),current_orientations(ind+3)
  END DO
  CLOSE(10)
  PRINT *,"Read initial data from file. "
  !PRINT '(*(F16.8))', current_positions

  !--------------------------------------------------
  ! Precompute constants
  !--------------------------------------------------
  CALL precomputeLegendrePolynomials(quadrature_points, quadrature_weights, legendre_polynomials)
  CALL precomputeLambda(lambda, eigen)

  !--------------------------------------------------
  ! Simulation Step
  !--------------------------------------------------
  CALL SYSTEM_CLOCK(total_count1, total_count_rate, total_count_max)

    !--------------------------------------------------
    ! 1. Assemble System
    !--------------------------------------------------
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)

    b_vector = 0.0

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(quadrature_points, quadrature_weights, legendre_polynomials, lambda, eigen, a_matrix, b_vector, external_force, c, d, e, cc, D1)
    DO i = 0, NUMBER_OF_FIBERS-1

      position_i = current_positions(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS)
      orientation_i = current_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS)

      DO j = 0, NUMBER_OF_FIBERS-1

        IF (i==j) THEN

          DO force_index = 0, NUMBER_OF_TERMS_IN_FORCE_EXPANSION-1
            a_matrix( &
              i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1, &
              j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 1) = &
                1.0
            a_matrix( &
              i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2, &
              j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 2) = &
                1.0
            a_matrix( &
              i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 3, &
              j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index * DIMENSIONS + 3) = &
                1.0
          END DO

        ELSE

          position_j = current_positions(j*DIMENSIONS + 1:j*DIMENSIONS + DIMENSIONS)
          orientation_j = current_orientations(j*DIMENSIONS + 1:j*DIMENSIONS + DIMENSIONS)

          DO force_index_j = 0, NUMBER_OF_TERMS_IN_FORCE_EXPANSION-1

            DO quadrature_index_i = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1

              G(quadrature_index_i+1,:) = 0.0
              GF(quadrature_index_i+1,:) = 0.0

              position_on_fiber_i = position_i + quadrature_points(quadrature_index_i+1) * orientation_i

              DO quadrature_index_j = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
                position_on_fiber_j = position_j + quadrature_points(quadrature_index_j+1) * orientation_j

                difference = position_on_fiber_i - position_on_fiber_j
                difference2 = difference**2

                invDistance = 1.0/SQRT(difference2(1)+difference2(2)+difference2(3))
                invDistance3 = invDistance * invDistance * invDistance
                invDistance5 = invDistance3 * invDistance * invDistance

                K(1) = invDistance + invDistance3 * difference2(1) + 2.0 * SLENDERNESS * SLENDERNESS * (invDistance3 - 3.0 * invDistance5 * difference2(1))
                K(2) = invDistance + invDistance3 * difference2(2) + 2.0 * SLENDERNESS * SLENDERNESS * (invDistance3 - 3.0 * invDistance5 * difference2(2))
                K(3) = invDistance + invDistance3 * difference2(3) + 2.0 * SLENDERNESS * SLENDERNESS * (invDistance3 - 3.0 * invDistance5 * difference2(3))
                K(4) = invDistance3 * difference(1) * difference(2) + 2.0 * SLENDERNESS * SLENDERNESS * (-3.0) * invDistance5 * difference(1) * difference(2)
                K(5) = invDistance3 * difference(1) * difference(3) + 2.0 * SLENDERNESS * SLENDERNESS * (-3.0) * invDistance5 * difference(1) * difference(3)
                K(6) = invDistance3 * difference(2) * difference(3) + 2.0 * SLENDERNESS * SLENDERNESS * (-3.0) * invDistance5 * difference(2) * difference(3)

                quadrature_weight = quadrature_weights(quadrature_index_j+1)
                legendre_polynomial = legendre_polynomials(quadrature_index_j+1, force_index_j+1)

                PLUS_EQUALS(G(quadrature_index_i+1,1),quadrature_weight * K(1) * legendre_polynomial)
                PLUS_EQUALS(G(quadrature_index_i+1,2),quadrature_weight * K(2) * legendre_polynomial)
                PLUS_EQUALS(G(quadrature_index_i+1,3),quadrature_weight * K(3) * legendre_polynomial)
                PLUS_EQUALS(G(quadrature_index_i+1,4),quadrature_weight * K(4) * legendre_polynomial)
                PLUS_EQUALS(G(quadrature_index_i+1,5),quadrature_weight * K(5) * legendre_polynomial)
                PLUS_EQUALS(G(quadrature_index_i+1,6),quadrature_weight * K(6) * legendre_polynomial)

                PLUS_EQUALS(GF(quadrature_index_i+1,1), quadrature_weight * (K(1) * external_force(1) + K(4) * external_force(2) + K(5) * external_force(3)))
                PLUS_EQUALS(GF(quadrature_index_i+1,2), quadrature_weight * (K(4) * external_force(1) + K(2) * external_force(2) + K(6) * external_force(3)))
                PLUS_EQUALS(GF(quadrature_index_i+1,3), quadrature_weight * (K(5) * external_force(1) + K(6) * external_force(2) + K(3) * external_force(3)))

              END DO
            END DO

            force_index_i = 0
            T = 0.0
            TF = 0.0

            DO quadrature_index_i = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
              quadrature_weight = quadrature_weights(quadrature_index_i+1)
              legendre_polynomial = legendre_polynomials(quadrature_index_i+1, 0 + 1)

              PLUS_EQUALS(T,quadrature_weight * G(quadrature_index_i+1,:) * legendre_polynomial)

              IF (force_index_j == 0) THEN
                PLUS_EQUALS(TF, quadrature_weight * GF(quadrature_index_i+1,:) * legendre_polynomial)
              END IF
            END DO

            Q(1) = T(1) * orientation_i(1) + T(4) * orientation_i(2) + T(5) * orientation_i(3)
            Q(2) = T(4) * orientation_i(1) + T(2) * orientation_i(2) + T(6) * orientation_i(3)
            Q(3) = T(5) * orientation_i(1) + T(6) * orientation_i(2) + T(3) * orientation_i(3)

            x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 1
            y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 2
            z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 3

            x_col_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 1
            y_col_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 2
            z_col_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 3

            a_matrix(x_row_index, x_col_index) = D1 * orientation_i(1) * Q(1);
            a_matrix(x_row_index, y_col_index) = D1 * orientation_i(1) * Q(2);
            a_matrix(x_row_index, z_col_index) = D1 * orientation_i(1) * Q(3);
            a_matrix(y_row_index, x_col_index) = D1 * orientation_i(2) * Q(1);
            a_matrix(y_row_index, y_col_index) = D1 * orientation_i(2) * Q(2);
            a_matrix(y_row_index, z_col_index) = D1 * orientation_i(2) * Q(3);
            a_matrix(z_row_index, x_col_index) = D1 * orientation_i(3) * Q(1);
            a_matrix(z_row_index, y_col_index) = D1 * orientation_i(3) * Q(2);
            a_matrix(z_row_index, z_col_index) = D1 * orientation_i(3) * Q(3);

            IF (force_index_j == 0) THEN
              QF = TF(1) * orientation_i(1) + TF(2) * orientation_i(2) + TF(3) * orientation_i(3)

              MINUS_EQUALS(b_vector(x_row_index), D1 * orientation_i(1) * QF)
              MINUS_EQUALS(b_vector(y_row_index), D1 * orientation_i(2) * QF)
              MINUS_EQUALS(b_vector(z_row_index), D1 * orientation_i(3) * QF)

            END IF

            DO force_index_i = 1, NUMBER_OF_TERMS_IN_FORCE_EXPANSION-1
              gamma = 0.5 * (2.0 * (force_index_i + 1) + 1.0) / (d + e - cc * lambda(force_index_i+1));

              T = 0.0
              TF = 0.0

              DO quadrature_index_i = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
                quadrature_weight = quadrature_weights(quadrature_index_i+1)
                legendre_polynomial = legendre_polynomials(quadrature_index_i+1, force_index_i + 1)

                PLUS_EQUALS(T,quadrature_weight * G(quadrature_index_i+1,:) * legendre_polynomial)

                IF (force_index_j == 0) THEN
                  PLUS_EQUALS(TF, quadrature_weight * GF(quadrature_index_i+1,:) * legendre_polynomial)
                END IF
              END DO

              Q(1) = T(1) * orientation_i(1) + T(4) * orientation_i(2) + T(5) * orientation_i(3)
              Q(2) = T(4) * orientation_i(1) + T(2) * orientation_i(2) + T(6) * orientation_i(3)
              Q(3) = T(5) * orientation_i(1) + T(6) * orientation_i(2) + T(3) * orientation_i(3)

              x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 1
              y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 2
              z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 3

              x_col_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 1
              y_col_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 2
              z_col_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 3

              a_matrix(x_row_index,x_col_index) = gamma * (T(1) - eigen(force_index_i+1) * orientation_i(1) * Q(1));
              a_matrix(x_row_index,y_col_index) = gamma * (T(4) - eigen(force_index_i+1) * orientation_i(1) * Q(2));
              a_matrix(x_row_index,z_col_index) = gamma * (T(5) - eigen(force_index_i+1) * orientation_i(1) * Q(3));
              a_matrix(y_row_index,x_col_index) = gamma * (T(4) - eigen(force_index_i+1) * orientation_i(2) * Q(1));
              a_matrix(y_row_index,y_col_index) = gamma * (T(2) - eigen(force_index_i+1) * orientation_i(2) * Q(2));
              a_matrix(y_row_index,z_col_index) = gamma * (T(6) - eigen(force_index_i+1) * orientation_i(2) * Q(3));
              a_matrix(z_row_index,x_col_index) = gamma * (T(5) - eigen(force_index_i+1) * orientation_i(3) * Q(1));
              a_matrix(z_row_index,y_col_index) = gamma * (T(6) - eigen(force_index_i+1) * orientation_i(3) * Q(2));
              a_matrix(z_row_index,z_col_index) = gamma * (T(3) - eigen(force_index_i+1) * orientation_i(3) * Q(3));

              IF (force_index_j == 0) THEN
                QF = TF(1) * orientation_i(1) + TF(2) * orientation_i(2) + TF(3) * orientation_i(3)

                MINUS_EQUALS(b_vector(x_row_index), gamma * (TF(1) - eigen(force_index_i+1) * orientation_i(1) * QF))
                MINUS_EQUALS(b_vector(y_row_index), gamma * (TF(2) - eigen(force_index_i+1) * orientation_i(2) * QF))
                MINUS_EQUALS(b_vector(z_row_index), gamma * (TF(3) - eigen(force_index_i+1) * orientation_i(3) * QF))
              END IF
            END DO
          END DO
        END IF
      END DO
    END DO
    !$OMP END PARALLEL DO

    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:assemble_matrix:", CPU_p

    ! OPEN(10,file="AMat.out");
    ! DO i=1,TOTAL_NUMBER_OF_ROWS
    !   WRITE(10,'(*(F16.8))') (a_matrix(i,j),j=1,TOTAL_NUMBER_OF_ROWS)
    ! END DO
    ! CLOSE(10)
    !
    ! OPEN(10,file="BVec.out");
    ! DO i=1,TOTAL_NUMBER_OF_ROWS
    !   WRITE(10,'(*(F16.8))') (b_vector(i))
    ! END DO
    ! CLOSE(10)

    !--------------------------------------------------
    ! 2. Solve System
    !--------------------------------------------------
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)
    CALL sgesv(TOTAL_NUMBER_OF_ROWS, 1, a_matrix, TOTAL_NUMBER_OF_ROWS, IPIV, b_vector, TOTAL_NUMBER_OF_ROWS, INFO)
    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:solve_system:", CPU_p

    ! OPEN(10,file="XVec.out");
    ! DO i=1,TOTAL_NUMBER_OF_ROWS
    !   WRITE(10,'(*(F16.8))') (b_vector(i))
    ! END DO
    ! CLOSE(10)

    !--------------------------------------------------
    ! 3. Update System
    !--------------------------------------------------
    !--------------------------------------------------
    ! 3.1. Update Velocity
    !--------------------------------------------------
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(quadrature_points, quadrature_weights, legendre_polynomials, lambda, eigen, current_translational_velocities, current_rotational_velocities, external_force, c, d, e, cc, D1)
    DO i = 0, NUMBER_OF_FIBERS-1

      position_i = current_positions(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS)
      orientation_i = current_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS)

      oriented_force = sum(orientation_i * (2.0 * external_force)) * orientation_i;
      current_translational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS) = (d + 2.0) * 2.0 * external_force + (d - 2.0) * oriented_force
      current_rotational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS) = 0.0

      DO j = 0, NUMBER_OF_FIBERS-1

        IF (i /= j) THEN

          position_j = current_positions(j*DIMENSIONS + 1:j*DIMENSIONS + DIMENSIONS)
          orientation_j = current_orientations(j*DIMENSIONS + 1:j*DIMENSIONS + DIMENSIONS)

          DO quadrature_index_i = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1

            GF(quadrature_index_i+1,:) = 0.0

            position_on_fiber_i = position_i + quadrature_points(quadrature_index_i+1) * orientation_i

            DO quadrature_index_j = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
              position_on_fiber_j = position_j + quadrature_points(quadrature_index_j+1) * orientation_j

              difference = position_on_fiber_i - position_on_fiber_j
              difference2 = difference**2

              invDistance = 1.0/SQRT(difference2(1)+difference2(2)+difference2(3))
              invDistance3 = invDistance * invDistance * invDistance
              invDistance5 = invDistance3 * invDistance * invDistance

              K(1) = invDistance + invDistance3 * difference2(1) + 2.0 * SLENDERNESS * SLENDERNESS * (invDistance3 - 3.0 * invDistance5 * difference2(1))
              K(2) = invDistance + invDistance3 * difference2(2) + 2.0 * SLENDERNESS * SLENDERNESS * (invDistance3 - 3.0 * invDistance5 * difference2(2))
              K(3) = invDistance + invDistance3 * difference2(3) + 2.0 * SLENDERNESS * SLENDERNESS * (invDistance3 - 3.0 * invDistance5 * difference2(3))
              K(4) = invDistance3 * difference(1) * difference(2) + 2.0 * SLENDERNESS * SLENDERNESS * (-3.0) * invDistance5 * difference(1) * difference(2)
              K(5) = invDistance3 * difference(1) * difference(3) + 2.0 * SLENDERNESS * SLENDERNESS * (-3.0) * invDistance5 * difference(1) * difference(3)
              K(6) = invDistance3 * difference(2) * difference(3) + 2.0 * SLENDERNESS * SLENDERNESS * (-3.0) * invDistance5 * difference(2) * difference(3)

              force_on_fiber_j = external_force

              DO force_index_j = 0, NUMBER_OF_TERMS_IN_FORCE_EXPANSION-1
                legendre_polynomial = legendre_polynomials(quadrature_index_j+1, force_index_j+1)

                x_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1
                y_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2
                z_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 3

                PLUS_EQUALS(force_on_fiber_j(1), b_vector(x_row_index) * legendre_polynomial)
                PLUS_EQUALS(force_on_fiber_j(2), b_vector(y_row_index) * legendre_polynomial)
                PLUS_EQUALS(force_on_fiber_j(3), b_vector(z_row_index) * legendre_polynomial)
              END DO

              quadrature_weight = quadrature_weights(quadrature_index_j+1)

              PLUS_EQUALS(GF(quadrature_index_i+1,1), quadrature_weight * (K(1) * force_on_fiber_j(1) + K(4) * force_on_fiber_j(2) + K(5) * force_on_fiber_j(3)))
              PLUS_EQUALS(GF(quadrature_index_i+1,2), quadrature_weight * (K(4) * force_on_fiber_j(1) + K(2) * force_on_fiber_j(2) + K(6) * force_on_fiber_j(3)))
              PLUS_EQUALS(GF(quadrature_index_i+1,3), quadrature_weight * (K(5) * force_on_fiber_j(1) + K(6) * force_on_fiber_j(2) + K(3) * force_on_fiber_j(3)))
            END DO

          END DO

          TFA0 = 0.0
          TFA1 = 0.0

          DO quadrature_index_i = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
            quadrature_weight = quadrature_weights(quadrature_index_i+1)
            legendre_polynomial = legendre_polynomials(quadrature_index_i+1, 0 + 1)

            PLUS_EQUALS(TFA0, quadrature_weight * GF(quadrature_index_i+1,:))
            PLUS_EQUALS(TFA1, quadrature_weight * GF(quadrature_index_i+1,:) * legendre_polynomial)
          END DO

          PLUS_EQUALS(current_translational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS), TFA0)

          TFA1_TMP = TFA1 * orientation_i

          PLUS_EQUALS(current_rotational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS), TFA1 - (orientation_i * TFA1_TMP(1) + orientation_i * TFA1_TMP(2) + orientation_i * TFA1_TMP(3)))

        END IF

      END DO

      MULTIPLY_EQUALS(current_translational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS), (0.5 / d))
      MULTIPLY_EQUALS(current_rotational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS), (1.5 / d))

    END DO
    !$OMP END PARALLEL DO

    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:update_velocities:", CPU_p

    ! OPEN(10,file="TRANSVel.out");
    ! DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
    !  WRITE(10,'(*(F16.8))') (current_translational_velocities(i))
    ! END DO
    ! CLOSE(10)
    ! OPEN(10,file="ROTVel.out");
    ! DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
    !  WRITE(10,'(*(F16.8))') (current_rotational_velocities(i))
    ! END DO
    ! CLOSE(10)

    !--------------------------------------------------
    ! 3.2. Update Fibers
    !--------------------------------------------------
    CALL SYSTEM_CLOCK(count1, count_rate, count_max)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    DO i = 0, NUMBER_OF_FIBERS-1
      next_positions(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS) = &
        current_positions(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS) + TIMESTEP * current_translational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS)
      next_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS) = &
        current_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS) + TIMESTEP * current_rotational_velocities(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS)

      DIVIDE_EQUALS(next_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS), SQRT(DOT_PRODUCT(next_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS), next_orientations(i*DIMENSIONS + 1:i*DIMENSIONS + DIMENSIONS))))
    END DO
    !$OMP END PARALLEL DO

    ! OPEN(10,file="POS.out");
    ! DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
    !   WRITE(10,'(*(F16.8))') (next_positions(i))
    ! END DO
    ! CLOSE(10)
    ! OPEN(10,file="ORIENT.out");
    ! DO i=1,NUMBER_OF_FIBERS * DIMENSIONS
    !   WRITE(10,'(*(F16.8))') (next_orientations(i))
    ! END DO
    ! CLOSE(10)

    CALL SYSTEM_CLOCK(count2, count_rate, count_max)
    CPU_p = real(count2-count1)/count_rate
    PRINT *,"BENCHMARK:update_fibers:", CPU_p

  CALL SYSTEM_CLOCK(total_count2, total_count_rate, total_count_max)
  total_CPU_p = real(total_count2-total_count1)/total_count_rate
  PRINT *,"BENCHMARK:$total:", total_CPU_p

END PROGRAM fibers

SUBROUTINE precomputeLegendrePolynomials(quadrature_points, quadrature_weights, legendre_polynomials)

  REAL*8::p0,p1,p2
  REAL*8::w0,w1,w2

  REAL*8::lower_bound
  REAL*8::interval_size

  REAL*4,INTENT(OUT),DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_points
  REAL*4,INTENT(OUT),DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::quadrature_weights
  REAL*4,INTENT(OUT),DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS,NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::legendre_polynomials
  REAL*8,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::internal_quadrature_points
  REAL*8,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS)::internal_quadrature_weights
  REAL*8,DIMENSION(TOTAL_NUMBER_OF_QUADRATURE_POINTS,NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::internal_legendre_polynomials

  INTEGER::interval_index, force_index, point_index

  p0 = -SQRT(15.0d0) / 5.0d0
  p1 = 0.0d0
  p2 = SQRT(15.0d0) / 5.0d0

  w0 = 5.0d0 / 9.0d0
  w1 = 8.0d0 / 9.0d0
  w2 = 5.0d0 / 9.0d0

  lower_bound = -1.0d0

  interval_size = 2.0d0 / NUMBER_OF_QUADRATURE_INTERVALS

  DO interval_index = 0, NUMBER_OF_QUADRATURE_INTERVALS-1
    internal_quadrature_points(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 1) = &
      (2.0d0 * lower_bound + interval_size + p0 * interval_size) / 2.0d0
    internal_quadrature_points(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 2) = &
      (2.0d0 * lower_bound + interval_size + p1 * interval_size) / 2.0d0
    internal_quadrature_points(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 3) = &
      (2.0d0 * lower_bound + interval_size + p2 * interval_size) / 2.0d0

    internal_quadrature_weights(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 1) = &
      w0 / NUMBER_OF_QUADRATURE_INTERVALS
    internal_quadrature_weights(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 2) = &
      w1 / NUMBER_OF_QUADRATURE_INTERVALS
    internal_quadrature_weights(interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL + 3) = &
      w2 / NUMBER_OF_QUADRATURE_INTERVALS

    lower_bound = lower_bound + interval_size
  END DO

  quadrature_points = internal_quadrature_points
  quadrature_weights = internal_quadrature_weights

  DO force_index = 0, NUMBER_OF_TERMS_IN_FORCE_EXPANSION-1
    DO point_index = 0, TOTAL_NUMBER_OF_QUADRATURE_POINTS-1
      internal_legendre_polynomials(point_index + 1,force_index + 1) = &
        calculateLegendrePolynomial(internal_quadrature_points(point_index+1), force_index+1)
    END DO
  END DO

  legendre_polynomials = internal_legendre_polynomials

END SUBROUTINE precomputeLegendrePolynomials

SUBROUTINE precomputeLambda(lambda, eigen)

  REAL*4,INTENT(OUT),DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::lambda
  REAL*4,INTENT(OUT),DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::eigen
  REAL*8,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::internal_lambda
  REAL*8,DIMENSION(NUMBER_OF_TERMS_IN_FORCE_EXPANSION)::internal_eigen

  INTEGER::force_index

  REAL*8::c, d, e, cc
  c = LOG(SLENDERNESS * SLENDERNESS * EXP(1.0d0))
  d = -c
  e = 2.0d0
  cc = 1.0d0

  internal_lambda(1) = 2.0d0
  internal_eigen(1) = ((d - e - cc * internal_lambda(1)) / 2.0d0) / (d - cc * internal_lambda(1))

  DO force_index = 2, NUMBER_OF_TERMS_IN_FORCE_EXPANSION
    internal_lambda(force_index) = internal_lambda(force_index - 1) + 2.0d0 / force_index;
    internal_eigen(force_index) = ((d - e - cc * internal_lambda(force_index)) / 2.0d0) / (d - cc * internal_lambda(force_index))
  END DO

  lambda = internal_lambda
  eigen = internal_eigen

END SUBROUTINE precomputeLambda

FUNCTION calculateLegendrePolynomial(x, n)

  REAL*8,INTENT(IN)::x
  INTEGER,INTENT(IN)::n

  IF (n == 0) THEN
    calculateLegendrePolynomial = 1.0d0
  ELSEIF (n == 1) THEN
    calculateLegendrePolynomial = x
  ELSEIF (n == 2) THEN
    calculateLegendrePolynomial = (1.0d0 / 2.0d0) * &
      (3.0d0 * x**2 - 1.0d0)
  ELSEIF (n == 3) THEN
    calculateLegendrePolynomial = (1.0d0 / 2.0d0) * &
      (5.0d0 * x**3 - 3.0d0 * x)
  ELSEIF (n == 4) THEN
    calculateLegendrePolynomial = (1.0d0 / 8.0d0) * &
      (35.0d0 * x**4 - 30.0d0 * x**2 + 3.0d0)
  ELSEIF (n == 5) THEN
    calculateLegendrePolynomial = (1.0d0 / 8.0d0) * &
      (63.0d0 * x**5 - 70.0d0 * x**3 + 15.0d0 * x)
  ELSEIF (n == 6) THEN
    calculateLegendrePolynomial = (1.0d0 / 16.0d0) * &
      (231.0d0 * x**6 - 315.0d0 * x**4 + 105.0d0 * x**2 - 5.0d0)
  ELSEIF (n == 7) THEN
    calculateLegendrePolynomial = (1.0d0 / 16.0d0) * &
      (429.0d0 * x**7 - 693.0d0 * x**5 + 315.0d0 * x**3 - 35.0d0 * x)
  ELSEIF (n == 8) THEN
    calculateLegendrePolynomial = (1.0d0 / 128.0d0) * &
      (6435.0d0 * x**8 - 12012.0d0 * x**6 + 6930.0d0 * x**4 - 1260.0d0 * x**2 + 35.0d0)
  ELSE
    PRINT *,"Could not precompute legendre polynomials - n not in range [1..8]: "
  ENDIF

END FUNCTION calculateLegendrePolynomial
