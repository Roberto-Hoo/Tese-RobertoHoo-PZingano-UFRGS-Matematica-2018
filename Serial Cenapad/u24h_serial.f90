!
! **************************************************************************
!
!   NUMERICAL SOLUTION OF  u  =  ( F(x,y) u_x )  +  ( G(x,y) u_y )
!                           t                  x                 y
!             where  F(x,y) = G(x,y) = | grad u |^(p-2)
!
! **************************************************************************
!
!  Code INITIAL_APPROXIMATION must have already been run !!! <---IMPORTANT
!
!          The code runs from t = t0  to t = tF
!
!  INPUT DATA: this code reads data from file "u24h_serial_input.dat" (*)
!
!  OUTPUT FILES: this code generates the following data files on disk:
!
!  (1) u24h_serial_input_previous.dat : a safe copy of the latest input file
!                                    used by this code (i.e., file (*) above)
!  (2) u24h_serial_INPUT.dat : output file with all major quantities computed
!                  by this code PLUS all previous input data contained in (*)
!  (3) u24h_serial_LOG.txt : a log file with running statistics and a summary
!              of major output results computed by this and previous runnings
!  (4) u24h_serial_ERROR.txt : an ERROR log file explaining some instances of
!           program error (inconsistencies or lack of proper input data file)
!
! ***************************************************************************
!

PROGRAM u24h_serial   ! <--- to be executed ONLY after program
!                        "initial_approximation.f90" has been run
IMPLICIT NONE

INTEGER, PARAMETER ::  M = 150, N = 150;   ! <--- no. of grid points

INTEGER  M_read, N_read   ! <--- M_read MUST be equal to M
! <--- N_read MUST be equal to N
INTEGER  nn   ! <--- nn is the space dimension
REAL*8   p    ! <--- p must be greater than nn

REAL*8   h    ! <--- spatial grid spacing

REAL*8,  PARAMETER ::  cfl = 0.01D0   ! <--- Courant-Friedrichs-Lewy number

REAL*8   dt   ! <--- timestep length

REAL*8   x(-M:M), y(-N:N)      ! <--- grid points
REAL*8   u(-M-1:M+1,-N-1:N+1)  ! <--- extended solution
REAL*8   v(-M-1:M+1,-N-1:N+1)  ! <--- (temporary) solution at new time level
REAL*8   u_previous(-M-1:M+1,-N-1:N+1) ! <--- solution at the previous statistics dumping

REAL*8   F(-M:M+1,-N:N)
REAL*8   G(-M:M,-N:N+1)

REAL*8   t0, tF
REAL*8   t

REAL*8   dt_dump           ! <--- specified by the program: initial_approximation.f90
REAL*8   time_next_dump    !           (usually defined to be: dt_dump = 0.01)

REAL*8   x1, x2, x3
REAL*8   y1, y2, y3
INTEGER  i1, i2, i3
INTEGER  j1, j2, j3

REAL*8   b1, b2, b3
REAL*8   b
REAL*8   m1, m2, m3

INTEGER  i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
INTEGER  j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
INTEGER  length_i1, length_i2, length_i3
INTEGER  length_j1, length_j2, length_j3

REAL*8   t1_min, t2_min, t3_min
REAL*8   t1_max, t2_max, t3_max

REAL*8   ff1_value_u0, ff2_value_u0, ff3_value_u0   ! <--- far field values of initial data

REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff1_value_u,  ff2_value_u,  ff3_value_u
REAL*8,  DIMENSION(:), ALLOCATABLE ::  min_u, max_u
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_mass, variation_sup
REAL*8,  DIMENSION(:), ALLOCATABLE ::  time_per_cycle

REAL*8   total_elapsed_time

INTEGER  count

REAL*8   min_u0, max_u0
REAL*8   mass_u0         ! <--- not used

INTEGER  no_runs

! PREVIOUS data (already computed) to be retained:

REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff1_value_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff2_value_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff3_value_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  min_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  max_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_mass_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_sup_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  time_per_cycle_previous

INTEGER, DIMENSION(:,:), ALLOCATABLE :: DATE_start_LOG_previous
INTEGER, DIMENSION(:,:), ALLOCATABLE :: DATE_finish_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff1_value_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff2_value_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff3_value_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_mass_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_sup_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  min_u_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  max_u_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  tF_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  elapsed_time_LOG_previous

REAL*8   variation_mass_0
REAL*8   variation_sup_0

! Local variables:

INTEGER  i, j, numLoops

REAL*8   h_to_pm2, q

LOGICAL  file_present

INTEGER  DATE_start(8)
INTEGER  DATE_finish(8)
INTEGER  values
REAL*8   time_start, time_finish   ! <--- used by CPU_TIME timing routine

REAL*8   tF_new

INTEGER  previous_length
INTEGER  new_length

INTEGER  no_dumps

INTEGER  new_count

REAL*8   FP_count, FLOPs

! Intrinsic functions used:

REAL*8,  INTRINSIC :: sum
REAL*8,  INTRINSIC :: abs
!  REAL*8,  INTRINSIC :: min, max
REAL*8,  INTRINSIC :: minval, maxval
INTEGER, INTRINSIC :: size, nint    ! <--- NINT rounds REALs to nearest INTEGERs

INTRINSIC DATE_AND_TIME
INTRINSIC CPU_TIME

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------
!
!  READING INPUT DATA FILE  u24h_serial_INPUT.dat:
!   (created by program previous run of "u24h_serial.f90"
!    or by program "initial_approximation.f90" on first run)
!
!--------------------------------------------------------


CALL DATE_AND_TIME(VALUES = DATE_start)

!  READING INPUT DATA:

INQUIRE( FILE = 'u24h_serial_INPUT.dat', EXIST = file_present )

IF ( file_present .EQV. .FALSE. ) THEN

    OPEN(11, FILE = 'u24h_serial_ERROR.txt', ACTION = 'WRITE', &
        STATUS = 'REPLACE')
    write(11,*)
    write(11,*) '*******************************************************'
    write(11,*) ' ERROR ***** ERROR ***** ERROR ***** ERROR ***** ERROR '
    write(11,*) '*******************************************************'
    write(11,*)
    write(11,51) DATE_start(3), DATE_start(2), DATE_start(1)
    51     FORMAT(1X,'DAY :  ',I2.2,' / ',I2.2,' / ',I4)
    write(11,52) DATE_start(5:7)
    52     FORMAT(1X,'TIME:  ',I2.2,' : ',I2.2,' : ',I2.2)
    write(11,*)
    write(11,*) '*******************************************************'
    write(11,*) ' ERROR CONDITION:'
    write(11,*) ' FILE  u24h_serial_INPUT.dat  NOT FOUND'
    write(11,*)
    write(11,*) ' EXECUTION WILL TERMINATE'
    write(11,*)
    write(11,*) '*******************************************************'

    CLOSE(11)

    write(6,*) '*****************************************'
    write(6,*) ' ERROR: see file "u24h_serial_ERROR.txt" '
    write(6,*) '        for explanation                  '
    write(6,*) ' ***** EXECUTION HAS BEEN ABORTED *****  '
    write(6,*) '*****************************************'

    STOP

ENDIF


!*****************************************************************************
!*****************************************************************************
!**              EXECUTION of  u24h_serial.f90  STARTS HERE:                **
!*****************************************************************************
!*****************************************************************************
do numLoops = 1, 2000

    previous_length = 0;

    OPEN(10, FILE = 'u24h_serial_INPUT.dat', ACTION = 'READ')

    510 FORMAT( 10(1X,I10) )
    520 FORMAT( 10(1X,E25.16) )

    read(10,510) M_read, N_read
    read(10,510) nn

    IF ( (M_read .NE. M) .OR. (N_read .NE. N) ) THEN

        OPEN(11, FILE = 'u24h_serial_ERROR.txt', ACTION = 'WRITE', &
            STATUS = 'REPLACE')
        write(11,*)
        write(11,*) '*******************************************************'
        write(11,*) ' ERROR ***** ERROR ***** ERROR ***** ERROR ***** ERROR '
        write(11,*) '*******************************************************'
        write(11,*)
        write(11,51) DATE_start(3), DATE_start(2), DATE_start(1)
        !51    FORMAT(1X,'DAY :  ',I2.2,' / ',I2.2,' / ',I4)
        write(11,52) DATE_start(5:7)
        !52    FORMAT(1X,'TIME:  ',I2.2,' : ',I2.2,' : ',I2.2)
        write(11,*)
        write(11,*) '*******************************************************'
        write(11,*) ' ERROR CONDITION:'
        write(11,*) ' INPUT VALUES of  M, N (in file: u24h_serial_input.dat)'
        write(11,*) ' ARE NOT THE EXPECTED VALUES: M = ', M, ' N = ', N
        write(11,*)
        write(11,*) ' EXECUTION WILL TERMINATE'
        write(11,*)
        write(11,*) '*******************************************************'

        CLOSE(11)

        write(6,*) '*****************************************'
        write(6,*) ' ERROR: see file "u24h_serial_ERROR.txt" '
        write(6,*) '        for explanation                  '
        write(6,*) ' ***** EXECUTION HAS BEEN ABORTED *****  '
        write(6,*) '*****************************************'

        STOP

    ENDIF

    read(10,520) p
    read(10,520) h

    read(10,520) t0, tF, dt_dump

    read(10,520) x1, x2, x3
    read(10,520) y1, y2, y3
    read(10,510) i1, i2, i3
    read(10,510) j1, j2, j3

    read(10,520) b1, b2, b3
    read(10,520) b

    read(10,510) count
    read(10,510) no_runs

    read(10,510) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
    read(10,510) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
    read(10,510) length_i1, length_i2, length_i3
    read(10,510) length_j1, length_j2, length_j3

    read(10,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

    read(10,520) min_u0, max_u0
    read(10,520) mass_u0

    read(10,520) t1_min, t2_min, t3_min
    read(10,520) t1_max, t2_max, t3_max

    read(10,520)

    ALLOCATE( DATE_start_LOG_previous (0:no_runs,1:8) )
    ALLOCATE( DATE_finish_LOG_previous(0:no_runs,1:8) )
    ALLOCATE( tF_LOG_previous(0:no_runs) )
    ALLOCATE( ff1_value_LOG_previous(0:no_runs) )
    ALLOCATE( ff2_value_LOG_previous(0:no_runs) )
    ALLOCATE( ff3_value_LOG_previous(0:no_runs) )
    ALLOCATE( variation_mass_LOG_previous(0:no_runs) )
    ALLOCATE( variation_sup_LOG_previous (0:no_runs) )
    ALLOCATE( min_u_LOG_previous(0:no_runs) )
    ALLOCATE( max_u_LOG_previous(0:no_runs) )
    ALLOCATE( elapsed_time_LOG_previous(0:no_runs) )

    DO i = 0, no_runs
        read(10,510) (DATE_start_LOG_previous (i,j), j = 1,8)
    ENDDO
    DO i = 0, no_runs
        read(10,510) (DATE_finish_LOG_previous(i,j), j = 1,8)
    ENDDO

    DO i = 0, no_runs
        read(10,520) tF_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        read(10,520) ff1_value_LOG_previous(i)
    ENDDO
    DO i = 0, no_runs
        read(10,520) ff2_value_LOG_previous(i)
    ENDDO
    DO i = 0, no_runs
        read(10,520) ff3_value_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        read(10,520) variation_mass_LOG_previous(i)
    ENDDO
    DO i = 0, no_runs
        read(10,520) variation_sup_LOG_previous (i)
    ENDDO

    DO i = 0, no_runs
        read(10,520) min_u_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        read(10,520) max_u_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        read(10,520) elapsed_time_LOG_previous(i)
    ENDDO

    read(10,520)

    read(10,520) (x(i), i = -M, M)
    read(10,520)
    read(10,520) (y(j), j = -N, N)

    read(10,520)

    DO i = -M, M
        read(10,520) (u(i,j), j = -N, N)   ! <--- initial state for the present running
    ENDDO

    IF ( no_runs > 0 ) THEN

        read(10, 510)
        read(10, 510) previous_length    ! <--- previous_length should
        read(10, 510)                    !      be equal to: count

        if ( previous_length .NE. count ) then

            OPEN(11, FILE = 'u24h_serial_ERROR.txt', ACTION = 'WRITE', &
                STATUS = 'REPLACE')
            write(11,*)
            write(11,*) '*******************************************************'
            write(11,*) ' ERROR ***** ERROR ***** ERROR ***** ERROR ***** ERROR '
            write(11,*) '*******************************************************'
            write(11,*)
            write(11,51) DATE_start(3), DATE_start(2), DATE_start(1)
            !51         FORMAT(1X,'DAY :  ',I2.2,' / ',I2.2,' / ',I4)
            write(11,52) DATE_start(5:7)
            !52         FORMAT(1X,'TIME:  ',I2.2,' : ',I2.2,' : ',I2.2)
            write(11,*)
            write(11,*) '*******************************************************'
            write(11,*) ' ERROR CONDITION:'
            write(11,*) ' variable "count" should be equal to '
            write(11,*) ' variable "previous_length"'
            write(11,*)
            write(11,*) ' count = ', count
            write(11,*) ' previous_length = ', previous_length
            write(11,*)
            write(11,*) ' EXECUTION WILL TERMINATE'
            write(11,*)
            write(11,*) '*******************************************************'

            CLOSE(11)

            write(6,*) '*****************************************'
            write(6,*) ' ERROR: see file "u24h_serial_ERROR.txt" '
            write(6,*) '        for explanation '
            write(6,*)
            write(6,*) 'count = ', count
            write(6,*) 'previous_length = ', previous_length
            write(6,*)
            write(6,*) ' ***** EXECUTION HAS BEEN ABORTED *****  '
            write(6,*) '*****************************************'

            STOP

        endif

        ALLOCATE( ff1_value_u_previous(previous_length) )
        ALLOCATE( ff2_value_u_previous(previous_length) )
        ALLOCATE( ff3_value_u_previous(previous_length) )
        ALLOCATE( min_u_previous(previous_length) )
        ALLOCATE( max_u_previous(previous_length) )
        ALLOCATE( variation_mass_previous(previous_length) )
        ALLOCATE( variation_sup_previous (previous_length) )
        ALLOCATE( time_per_cycle_previous(previous_length) )

        DO i = 1, previous_length
            read(10, 520) ff1_value_u_previous(i), ff2_value_u_previous(i), &
                ff1_value_u_previous(i)
        ENDDO

        read(10,520)
        DO i = 1, previous_length
            read(10, 520) variation_mass_previous(i), variation_sup_previous(i)
        ENDDO

        read(10,520)
        DO i = 1, previous_length
            read(10, 520) min_u_previous(i), max_u_previous(i)
        ENDDO

        read(10,520)
        DO i = 1, previous_length
            read(10, 520) time_per_cycle_previous(i)
        ENDDO

    ENDIF

    CLOSE(10)

    !-------------------------------------------------------------------------------
    !  CREATING SAFE COPY  u24h_serial_input_previous.dat  OF THE INPUT FILE ABOVE:
    !-------------------------------------------------------------------------------

    OPEN(20, FILE = 'u24h_serial_input_previous.dat', ACTION = 'WRITE', &
        STATUS = 'REPLACE')

    write(20,510) M_read, N_read
    write(20,510) nn

    write(20,520) p
    write(20,520) h

    write(20,520) t0, tF, dt_dump

    write(20,520) x1, x2, x3
    write(20,520) y1, y2, y3
    write(20,510) i1, i2, i3
    write(20,510) j1, j2, j3

    write(20,520) b1, b2, b3
    write(20,520) b

    write(20,510) count
    write(20,510) no_runs

    write(20,510) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
    write(20,510) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
    write(20,510) length_i1, length_i2, length_i3
    write(20,510) length_j1, length_j2, length_j3

    write(20,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

    write(20,520) min_u0, max_u0
    write(20,520) mass_u0

    write(20,520) t1_min, t2_min, t3_min
    write(20,520) t1_max, t2_max, t3_max

    write(20,520)

    DO i = 0, no_runs
        write(20,510) (DATE_start_LOG_previous (i,j), j = 1,8)
    ENDDO
    DO i = 0, no_runs
        write(20,510) (DATE_finish_LOG_previous(i,j), j = 1,8)
    ENDDO

    DO i = 0, no_runs
        write(20,520) tF_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        write(20,520) ff1_value_LOG_previous(i)
    ENDDO
    DO i = 0, no_runs
        write(20,520) ff2_value_LOG_previous(i)
    ENDDO
    DO i = 0, no_runs
        write(20,520) ff3_value_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        write(20,520) variation_mass_LOG_previous(i)
    ENDDO
    DO i = 0, no_runs
        write(20,520) variation_sup_LOG_previous (i)
    ENDDO

    DO i = 0, no_runs
        write(20,520) min_u_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        write(20,520) max_u_LOG_previous(i)
    ENDDO

    DO i = 0, no_runs
        write(20,520) elapsed_time_LOG_previous(i)
    ENDDO

    write(20,520)

    write(20,520) (x(i), i = -M_read, M_read)
    write(20,520)
    write(20,520) (y(j), j = -N_read, N_read)

    write(20,520)

    DO i = -M_read, M_read
        write(20,520) (u(i,j), j = -N_read, N_read)   ! <--- current initial state
    ENDDO

    IF ( no_runs > 0 ) THEN

        write(20, 510)
        write(20, 510) previous_length
        write(20, 510)

        DO i = 1, previous_length
            write(20,520) ff1_value_u_previous(i), ff2_value_u_previous(i), &
                ff3_value_u_previous(i)
        ENDDO

        write(20,520)
        DO i = 1, previous_length
            write(20,520) variation_mass_previous(i), variation_sup_previous(i)
        ENDDO

        write(20,520)
        DO i = 1, previous_length
            write(20,520) min_u_previous(i), max_u_previous(i)
        ENDDO

        write(20,520)
        DO i = 1, previous_length
            write(20, 520) time_per_cycle_previous(i)
        ENDDO

    ENDIF

    CLOSE(20)

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    !----------------------------------------------------------
    !              Setting new t0, tF:
    !----------------------------------------------------------

    tF_new = tF + 100D0;
    if ( tF < 50.1D0 ) tF_new = 100D0;
    if ( tF < 20.1D0 ) tF_new =  50D0;
    if ( tF < 10.1D0 ) tF_new =  20D0;
    if ( tF <  5.1D0 ) tF_new =  10D0;
    if ( tF <  1.1D0 ) tF_new =   5D0;
    if ( tF <  0.1D0 ) tF_new =   1D0;

    !  if ( tF < 0.2D0 ) tF_new = 0.20D0;  ! <------ eliminate
    !  if ( tF < 0.1D0 ) tF_new = 0.10D0;  ! <------ eliminate

    t0 = tF;
    tF = tF_new;

    !----------------------------------------------------------
    !     ALLOCATING the new solution-statistics arrays:
    !----------------------------------------------------------

    no_dumps = NINT( (tF-t0)/dt_dump );
    new_length = previous_length + no_dumps;

    ALLOCATE( ff1_value_u(1:no_dumps) )
    ALLOCATE( ff2_value_u(1:no_dumps) )
    ALLOCATE( ff3_value_u(1:no_dumps) )

    ALLOCATE( variation_mass(1:no_dumps) )
    ALLOCATE( variation_sup (1:no_dumps) )

    ALLOCATE( min_u(1:no_dumps) )
    ALLOCATE( max_u(1:no_dumps) )

    ALLOCATE( time_per_cycle(1:no_dumps) )


    !----------------------------------------------------------
    !      Extending solution values to extended grid:
    !----------------------------------------------------------

    u(-M-1,-N:N) = u(-M,-N:N);
    u( M+1,-N:N) = u( M,-N:N);
    u(-M-1:M+1,-N-1) = u(-M-1:M+1,-N);
    u(-M-1:M+1, N+1) = u(-M-1:M+1, N);

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    !
    ! Computing steady-state limit of solution u(.,t):
    !

    dt = cfl*h*h;  ! <--- time step length

    ! --------------------------------------------------------------------

    write(6,*) ' '
    write(6,*) '***************************************************************'
    write(6,*) 'Computing ...'
    write(6,*) '***************************************************************'
    write(6,*) ' '
    write(6,*) '---------------------------------------------------------------'
    write(6,*) '  variation mass   variation (sup)  far field value 2 '
    write(6,*) '---------------------------------------------------------------'

    q = (p - 2D0)/2D0;

    h_to_pm2 = h**(p-2);

    time_next_dump = t0 + dt_dump;

    t = t0;

    u_previous = u;

    new_count = 0;

    no_runs = no_runs + 1;

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !                    MAIN COMPUTATION SECTION
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------


    DO WHILE ( t < tF )


        CALL CPU_TIME(time_start)


        DO WHILE ( t < time_next_dump )


            t = t + dt;


            ! ****************************************************************** !
            !                computation of F values on the grid:                !
            ! ****************************************************************** !
            ! ******  F(i,j) for -M <= i <= M + 1 and -N <= j <= N:

            F(-M:M+1,-N:N) = ( (u(-M:M+1,-N:N) - u(-M-1:M,-N:N))**2 + &
                ( (u(-M-1:M,-N+1:N+1)+u(-M:M+1,-N+1:N+1))-(u(-M-1:M,-N-1:N-1)+u(-M:M+1,-N-1:N-1)) )**2 /16 )**q / h_to_pm2;

            ! ****************************************************************** !
            !                computation of G values on the grid:                !
            ! ****************************************************************** !
            ! ******  G(i,j) for -M <= i <= M and -N <= j <= N + 1:

            G(-M:M,-N:N+1) = ( (u(-M:M,-N:N+1) - u(-M:M,-N-1:N))**2 + &
                ( (u(-M+1:M+1,-N-1:N)+u(-M+1:M+1,-N:N+1))-(u(-M-1:M-1,-N-1:N)+u(-M-1:M-1,-N:N+1)) )**2 /16 )**q / h_to_pm2;


            ! ****************************************************************** !
            !      computation of v = [ new u values at the new time level ]:    !
            ! ****************************************************************** !
            ! ******  v(i,j) for -M <= i <= M and -N <= j <= N:

            v(-M:M,-N:N) = u(-M:M,-N:N) + &
                cfl*( F(-M+1:M+1,-N:N)*(u(-M+1:M+1,-N:N)-u(-M:M,-N:N)) - F(-M:M,-N:N)*(u(-M:M,-N:N)-u(-M-1:M-1,-N:N)) ) + &
                cfl*( G(-M:M,-N+1:N+1)*(u(-M:M,-N+1:N+1)-u(-M:M,-N:N)) - G(-M:M,-N:N)*(u(-M:M,-N:N)-u(-M:M,-N-1:N-1)) );

            ! ******  extending  v  to the extra grid points:

            v(-M-1,-N:N) = v(-M,-N:N);
            v( M+1,-N:N) = v( M,-N:N);
            v(-M-1:M+1,-N-1) = v(-M-1:M+1,-N);
            v(-M-1:M+1, N+1) = v(-M-1:M+1, N);

            ! ****************************************************************** !
            !            computation of new u values completed!                  !
            ! ****************************************************************** !


            u = v;    ! <--- updating u (at the new time level)

            ! ******  correcting u values at the singular points:

            u(i1,j1) = b1;

            u(i2,j2) = b2;

            u(i3,j3) = b3;

        ENDDO

        !--- time to save solution statistics:

        CALL CPU_TIME(time_finish)

        new_count = new_count + 1;

        time_per_cycle(new_count) = time_finish - time_start;

        !--- solution far-field values:

        ff1_value_u(new_count) = ( sum( u(i_ff1a:i_ff1a+10,j_ff1a:j_ff1b) )/(11*length_j1) &
            + sum( u(i_ff1b-10:i_ff1b,j_ff1a:j_ff1b) )/(11*length_j1) &
            + sum( u(i_ff1a:i_ff1b,j_ff1a:j_ff1a+10) )/(11*length_i1) &
            + sum( u(i_ff1a:i_ff1b,j_ff1b-10:j_ff1b) )/(11*length_i1) )/4D0;

        ff2_value_u(new_count) = ( sum( u(i_ff2a:i_ff2a+10,j_ff2a:j_ff2b) )/(11*length_j2) &
            + sum( u(i_ff2b-10:i_ff2b,j_ff2a:j_ff2b) )/(11*length_j2) &
            + sum( u(i_ff2a:i_ff2b,j_ff2a:j_ff2a+10) )/(11*length_i2) &
            + sum( u(i_ff2a:i_ff2b,j_ff2b-10:j_ff2b) )/(11*length_i2) )/4D0;

        ff3_value_u(new_count) = ( sum( u(i_ff3a:i_ff3a+10,j_ff3a:j_ff3b) )/(11*length_j3) &
            + sum( u(i_ff3b-10:i_ff3b,j_ff3a:j_ff3b) )/(11*length_j3) &
            + sum( u(i_ff3a:i_ff3b,j_ff3a:j_ff3a+10) )/(11*length_i3) &
            + sum( u(i_ff3a:i_ff3b,j_ff3b-10:j_ff3b) )/(11*length_i3) )/4D0;

        !--------------------------------------------------

        ! saving minimum and maximum solution values at current time:

        min_u(new_count) = minval( u );
        max_u(new_count) = maxval( u );

        ! computing supnorm of solution variation from previous u:

        variation_sup(new_count) = maxval( abs( u - u_previous ) );

        ! computing mass of solution variation from previous u:

        variation_mass(new_count) = sum( u - u_previous )*h**2;

        ! current solution values become previous solution values
        ! for the next iteration cycle:

        u_previous = u;

        ! printing out some statistics before going to next cycle:

        write(6,50) (t-t0)/(tF-t0)
        50  FORMAT( 3E15.6 )
        write(6,50) variation_mass(new_count), variation_sup(new_count), ff2_value_u(new_count)

        ! updating time of next dumping for solution statistics:

        time_next_dump = time_next_dump + dt_dump;

    ENDDO

    total_elapsed_time = sum( time_per_cycle );

    count = count + new_count;   ! <--- updating variable "count"

    CALL DATE_AND_TIME(VALUES = DATE_finish)


    write(6,*) '------------------------------------------------------------'
    write(6,*) 'TOTAL elapsed time:'
    write(6,*)  total_elapsed_time
    write(6,*) '------------------------------------------------------------'
    write(6,*) '             SUCCESSFUL EXECUTION                           '
    write(6,*) '------------------------------------------------------------'
    write(6,*) '------------------------------------------------------------'

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !                       SAVING OUTPUT DATA TO DISK
    !      (this data will serve as new input data for the program
    !                "u24h_serial.f90" on its next execution)
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    OPEN(50, FILE = 'u24h_serial_INPUT.dat', ACTION = 'WRITE', &
        STATUS = 'REPLACE')
    530 FORMAT( 10(1X,I10) )
    540 FORMAT( 10(1X,E25.16) )

    write(50,530) M, N
    write(50,530) nn

    write(50,540) p
    write(50,540) h

    write(50,540) t0, tF, dt_dump

    write(50,540) x1, x2, x3
    write(50,540) y1, y2, y3
    write(50,530) i1, i2, i3
    write(50,530) j1, j2, j3

    write(50,540) b1, b2, b3
    write(50,540) b

    write(50,530) count
    write(50,530) no_runs

    write(50,530) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
    write(50,530) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
    write(50,530) length_i1, length_i2, length_i3
    write(50,530) length_j1, length_j2, length_j3

    write(50,540) ff1_value_u0, ff2_value_u0, ff3_value_u0

    write(50,540) min_u0, max_u0
    write(50,540) mass_u0

    write(50,540) t1_min, t2_min, t3_min
    write(50,540) t1_max, t2_max, t3_max

    write(50,540)

    ! Saving LOG variables: VVVVV

    DO i = 0, no_runs-1
        write(50,530) (DATE_start_LOG_previous (i,j), j = 1,8)
    ENDDO
    write(50,530) (DATE_start (j), j = 1,8)

    DO i = 0, no_runs-1
        write(50,530) (DATE_finish_LOG_previous(i,j), j = 1,8)
    ENDDO
    write(50,530) (DATE_finish(j), j = 1,8)

    DO i = 0, no_runs-1
        write(50,540) tF_LOG_previous(i)
    ENDDO
    write(50,540) tF

    DO i = 0, no_runs-1
        write(50,540) ff1_value_LOG_previous(i)
    ENDDO
    write(50,540) ff1_value_u(new_count)

    DO i = 0, no_runs-1
        write(50,540) ff2_value_LOG_previous(i)
    ENDDO
    write(50,540) ff2_value_u(new_count)

    DO i = 0, no_runs-1
        write(50,540) ff3_value_LOG_previous(i)
    ENDDO
    write(50,540) ff3_value_u(new_count)

    DO i = 0, no_runs-1
        write(50,540) variation_mass_LOG_previous(i)
    ENDDO
    write(50,540) variation_mass(new_count)

    DO i = 0, no_runs-1
        write(50,540) variation_sup_LOG_previous (i)
    ENDDO
    write(50,540) variation_sup(new_count)

    DO i = 0, no_runs-1
        write(50,540) min_u_LOG_previous(i)
    ENDDO
    write(50,540) min_u(new_count)

    DO i = 0, no_runs-1
        write(50,540) max_u_LOG_previous(i)
    ENDDO
    write(50,540) max_u(new_count)

    DO i = 0, no_runs-1
        write(50,540) elapsed_time_LOG_previous(i)
    ENDDO
    write(50,540) total_elapsed_time

    ! LOG variables saved ^^^^^

    write(50,540)

    write(50,540) (x(i), i = -M, M)
    write(50,540)
    write(50,540) (y(j), j = -N, N)

    write(50,540)

    DO i = -M, M
        write(50,540) (u(i,j), j = -N, N)   ! <--- current final solution
    ENDDO

    write(50,530)
    write(50,530) new_length   ! <--- new+length = previous_length + no_dumps
    write(50,530)

    DO i = 1, previous_length
        write(50,540) ff1_value_u_previous(i), ff2_value_u_previous(i), &
            ff3_value_u_previous(i)
    ENDDO
    DO i = 1, no_dumps
        write(50,540) ff1_value_u(i), ff2_value_u(i), ff3_value_u(i)
    ENDDO

    write(50,540)
    DO i = 1, previous_length
        write(50,540) variation_mass_previous(i), variation_sup_previous(i)
    ENDDO
    DO i = 1, no_dumps
        write(50,540) variation_mass(i), variation_sup(i)
    ENDDO

    write(50,540)
    DO i = 1, previous_length
        write(50,540) min_u_previous(i), max_u_previous(i)
    ENDDO
    DO i = 1, no_dumps
        write(50,540) min_u(i), max_u(i)
    ENDDO

    write(50,540)
    DO i = 1, previous_length
        write(50,540) time_per_cycle_previous(i)
    ENDDO
    DO i = 1, no_dumps
        write(50,540) time_per_cycle(i)
    ENDDO

    CLOSE(50)


    !-------------------------------------------------------
    !    Generating  LOG FILE  for current run:
    !-------------------------------------------------------

    OPEN(70, FILE = 'u24h_serial_LOG.txt', ACTION = 'WRITE', &
        STATUS = 'REPLACE')
    write(70,*)
    write(70,*)
    write(70,*) '*******************************************************'
    write(70,*) '*           LOG DATA for u24h_serial.f90              *'
    write(70,*) '*******************************************************'
    write(70,*)
    write(70,*)

    !-----------------------------------------------------------------------
    !  TABLE 1: Run no., DATE/TIME (start & finish), tF, FP count, FLOP/s:
    !-----------------------------------------------------------------------

    701 FORMAT(T2,'|-----------','|---------------------',3('|----------------'),'|')
    702 FORMAT(T2,'|', T14,'|', T36,'|', T53,'|', T70,'|', T87,'|')       ! <--- blank line
    703 FORMAT(T2,'|', T5,'Run no.', T14,'|', T20,'DATE / TIME', T36,'|', &
        T44,'tF', T53,'|', T58,'FP count', T70,'|', T76,'FLOP/s', T87,'|')
    704 FORMAT(T2,'|', T14,'|', T16,I2.2,'/',I2.2,'/',I4,1X,I2.2,':',I2.2,':',I2.2,&
        T36,'|', T53,'|', T70,'|', T87,'|')
    705 FORMAT(T2,'|', T6,I3, T14,'|', T16,I2.2,'/',I2.2,'/',I4,1X,I2.2,':',I2.2,':',I2.2,&
        T36,'|', T39,E12.6, T53,'|', T56,E12.6, T70,'|', T73,E12.6, T87,'|')

    write(70,701)
    write(70,702)
    write(70,703)
    write(70,702)
    write(70,701)
    i = 0;
    FP_count = 20*(2*M+1)*(2*N+1)* dt_dump/dt;
    FLOPs = FP_count / elapsed_time_LOG_previous(0);
    write(70,704) DATE_start_LOG_previous(i,3), DATE_start_LOG_previous(i,2), &
        DATE_start_LOG_previous(i,1), DATE_start_LOG_previous(i,5:7)
    write(70,705) i, DATE_finish_LOG_previous(i,3), DATE_finish_LOG_previous(i,2), &
        DATE_finish_LOG_previous(i,1), DATE_finish_LOG_previous(i,5:7), &
        tF_LOG_previous(i), FP_count, FLOPs
    write(70,701)
    IF ( no_runs > 1 ) THEN
        do i = 1, no_runs-1
            FP_count = 20*(2*M+1)*(2*N+1)*(tF_LOG_previous(i) - tF_LOG_previous(i-1))/dt;
            FLOPs = FP_count / elapsed_time_LOG_previous(i);
            write(70,704) DATE_start_LOG_previous(i,3), DATE_start_LOG_previous(i,2), &
                DATE_start_LOG_previous(i,1), DATE_start_LOG_previous(i,5:7)
            write(70,705) i, DATE_finish_LOG_previous(i,3), DATE_finish_LOG_previous(i,2), &
                DATE_finish_LOG_previous(i,1), DATE_finish_LOG_previous(i,5:7), &
                tF_LOG_previous(i), FP_count, FLOPs
            write(70,701)
        enddo
    ENDIF
    i = no_runs;
    FP_count = 20*(2*M+1)*(2*N+1)*(tF - tF_LOG_previous(i-1))/dt;
    FLOPs = FP_count / total_elapsed_time;
    write(70,704) DATE_start(3), DATE_start(2), DATE_start(1), DATE_start(5:7)
    write(70,705) i, DATE_finish(3), DATE_finish(2), DATE_finish(1), DATE_finish(5:7), &
        tF, FP_count, FLOPs
    write(70,701)

    write(70,*)
    write(70,*)
    write(70,*)

    !-----------------------------------------------------------------------
    !  TABLE 2: Run no., ff1_value, ff2_value, ff3_value, elapsed time:
    !-----------------------------------------------------------------------

    711 FORMAT(T2,'|-----------',4('|----------------'),'|')
    712 FORMAT(T2,'|', T14,'|', T31,'|', T48,'|', T65,'|', T82,'|')       ! <--- blank line
    713 FORMAT(T2,'|', T5,'Run no.', T14,'|', T21,'ffv1', T31,'|', &
        T38,'ffv2', T48,'|', T55,'ffv3', T65,'|', T68,'elapsed time', T82,'|')
    714 FORMAT(T2,'|', T14,'|', T31,'|', T48,'|', T65,'|', T69, '(seconds)', T82,'|')
    715 FORMAT(T2,'|', T6,I3, T14,'|', T17,E12.6, T31,'|', T34,E12.6, T48,'|', &
        T51,E12.6, T65,'|', T68,E12.6, T82,'|')

    write(70,711)
    write(70,712)
    write(70,713)
    write(70,714)
    write(70,711)

    do i = 0, no_runs-1
        write(70,715) i, ff1_value_LOG_previous(i), ff2_value_LOG_previous(i), &
            ff3_value_LOG_previous(i), elapsed_time_LOG_previous(i)
        write(70,711)
    enddo
    i = no_runs;
    write(70,715) i, ff1_value_u(new_count), ff2_value_u(new_count), &
        ff3_value_u(new_count), total_elapsed_time
    write(70,711)

    write(70,*)
    write(70,*)
    write(70,*)

    !-----------------------------------------------------------------------
    !  TABLE 3: Run no., variation_mass, variation_sup, min_u, max_u:
    !-----------------------------------------------------------------------

    721 FORMAT(T2,'|-----------',4('|----------------'),'|')
    722 FORMAT(T2,'|', T14,'|', T31,'|', T48,'|', T65,'|', T82,'|')       ! <--- blank line
    723 FORMAT(T2,'|', T5,'Run no.', T14,'|', T16,'variation_mass', T31,'|', &
        T34,'variation_sup', T48,'|', T54,'min_u', T65,'|', T71,'max_u', T82,'|')
    724 FORMAT(T2,'|', T6,I3, T14,'|', T18,'----------', T31,'|', T35,'----------', &
        T48,'|', T51,E12.6, T65,'|', T68,E12.6, T82,'|')
    725 FORMAT(T2,'|', T6,I3, T14,'|', T16,E13.6, T31,'|', T34,E12.6, T48,'|', &
        T51,E12.6, T65,'|', T68,E12.6, T82,'|')

    write(70,721)
    write(70,722)
    write(70,723)
    write(70,722)
    write(70,721)

    i = 0;
    ! write(70,724) i, min_u_LOG_previous(0), max_u_LOG_previous(0)
    ! write(70,721)
    do i = 0, no_runs-1
        write(70,725) i, variation_mass_LOG_previous(i), variation_sup_LOG_previous(i), &
            min_u_LOG_previous(i), max_u_LOG_previous(i)
        write(70,721)
    enddo
    i = no_runs;
    write(70,725) i, variation_mass(new_count), variation_sup(new_count), &
        min_u(new_count), max_u(new_count)
    write(70,721)


    CLOSE(70)
 !--------------------------------------------------------------------------
 !--------------------------------------------------------------------------

    !-----------------------------------------
    !  Deallocating variables from memory:
    !-----------------------------------------

    DEALLOCATE( x, y )
    DEALLOCATE( u, v )
    DEALLOCATE( u_previous )
    DEALLOCATE( F, G )
    DEALLOCATE( ff1_value_u, ff2_value_u, ff3_value_u )
    DEALLOCATE( min_u, max_u )
    DEALLOCATE( variation_mass, variation_sup )
    DEALLOCATE( time_per_cycle )

    DEALLOCATE( ff1_value_u_previous )
    DEALLOCATE( ff2_value_u_previous )
    DEALLOCATE( ff3_value_u_previous )
    DEALLOCATE( min_u_previous )
    DEALLOCATE( max_u_previous )
    DEALLOCATE( variation_mass_previous )
    DEALLOCATE( variation_sup_previous )
    DEALLOCATE( time_per_cycle_previous )

    DEALLOCATE( DATE_start_LOG_previous )
    DEALLOCATE( DATE_finish_LOG_previous )
    DEALLOCATE( ff1_value_LOG_previous )
    DEALLOCATE( ff2_value_LOG_previous )
    DEALLOCATE( ff3_value_LOG_previous )
    DEALLOCATE( variation_mass_LOG_previous )
    DEALLOCATE( variation_sup_LOG_previous )
    DEALLOCATE( min_u_LOG_previous )
    DEALLOCATE( max_u_LOG_previous )
    DEALLOCATE( tF_LOG_previous )
    DEALLOCATE( elapsed_time_LOG_previous )

    !-----------------------------------------
    !       Deallocation completed
    !-----------------------------------------
end do !do numLoops = 1, 2000 linha mais ou menos 170
END PROGRAM
