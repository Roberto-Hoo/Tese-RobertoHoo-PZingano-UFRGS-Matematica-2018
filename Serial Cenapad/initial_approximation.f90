
 !--------------------------------------------------------------------------
 !--------------------------------------------------------------------------
 
 PROGRAM initial_approximation

    IMPLICIT NONE

    INTEGER, PARAMETER ::  M = 150      ! x-grid is: x(-M),...,x(M)
    INTEGER, PARAMETER ::  N = 150      ! y-grid is: y(-N),...,x(N)
    
    INTEGER, PARAMETER ::  nn = 2       ! space dimension
    REAL*8,  PARAMETER ::  p = 2.1D0    ! p must be greater than nn

    REAL*8,  PARAMETER ::  h = 0.1D0    ! spatial meshgrid spacing
    
    REAL*8   x(-M:M), y(-N:N)           ! grid point coordinates   
    REAL*8   u0(-M:M,-N:N)              ! initial solution values

    REAL*8,  PARAMETER ::  t0 = 0D0     ! initial time
    REAL*8,  PARAMETER ::  tF = 0D0     ! <--- new solution in NOT computed here
    
    INTEGER, PARAMETER ::  no_runs = 0  ! <--- again: new solution is NOT computed here
    INTEGER, PARAMETER ::  count = 0    ! <--- again: new solution is NOT computed here
    
    REAL*8,  PARAMETER ::  dt_dump = 0.01D0   ! time between consecutive dumpings 
                                              ! of solution statistics
    REAL*8   x1, x2, x3    ! x-coordinates of singularities
    REAL*8   y1, y2, y3    ! y-coordinates of singularities
    INTEGER  i1, i2, i3    ! i-indices of singularities
    INTEGER  j1, j2, j3    ! j-indices of singularities

    REAL*8   b1, b2, b3    ! solution values at the singularities
    REAL*8   b             ! initial guess of the far-field values
    REAL*8   m1, m2, m3    ! mass of initial cones relative to height b

    INTEGER  i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b    ! i-indices of far-field zones
    INTEGER  j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b    ! j-indices of far-field zones
    INTEGER  length_i1, length_i2, length_i3                   ! i-length of far-field zones
    INTEGER  length_j1, length_j2, length_j3                   ! j-length of far-field zones

    REAL*8   ff1_value_u0, ff2_value_u0, ff3_value_u0          ! average values of u0 in far-field zones

    REAL*8   variation_mass_0, variation_sup_0
    
    REAL*8   min_u0, max_u0    ! minimum and maximum values of initial solution over the grid
    REAL*8   mass_u0           ! mass of initial solution u0 relative to reference height b

 ! Extra variables needed to compute 
 ! variation_mass_0, variation_sup_0:
 
    REAL*8,  PARAMETER ::  cfl = 0.01D0

    REAL*8   u(-M-1:M+1,-N-1:N+1)
    REAL*8   v(-M-1:M+1,-N-1:N+1)
    REAL*8   F(-M:M+1,-N:N)
    REAL*8   G(-M:M,-N:N+1)

    REAL*8   t
    REAL*8   dt
    REAL*8   time_next_dump

    REAL*8   q
    REAL*8   h_to_pm2
    
 ! Local values of some additional mesh parameters:

    REAL*8   Rx, Ry

    REAL*8   minimum_distance_to_travel_1
    REAL*8   minimum_distance_to_travel_2
    REAL*8   minimum_distance_to_travel_3

    REAL*8   maximum_distance_to_travel_1
    REAL*8   maximum_distance_to_travel_2
    REAL*8   maximum_distance_to_travel_3

    REAL*8   t1_min, t2_min, t3_min
    REAL*8   t1_max, t2_max, t3_max

    REAL*8   d1, d2, d3 

    REAL*8   r12, r13, r23, r0

    INTEGER  DATE_start (8)
    INTEGER  DATE_finish(8)
    INTEGER  values(8)
    
    REAL*8   time_start, time_finish
    REAL*8   total_elapsed_time
    
 ! Local variables:
 
    INTEGER  i, j

    REAL*8,  PARAMETER ::  pi = 3.141592653589793238D0;

    REAL*8   expoente, common_factor
    REAL*8   C1, C2, C3, cc1, cc2, cc3

 ! Intrinsic functions used:

    REAL*8, INTRINSIC :: sum
    REAL*8, INTRINSIC :: abs, sqrt
    REAL*8, INTRINSIC :: min, max
    REAL*8, INTRINSIC :: minval, maxval
    REAL*8, INTRINSIC :: gamma 

    INTRINSIC  DATE_AND_TIME
 !  INTRINSIC  CPU_TIME
    INTRINSIC  DBLE

 
 !******************************************************************
 !*                                                                *
 !*                  MAIN PROGRAM starts here:                     *
 !*                                                                *
 !******************************************************************

    
    CALL DATE_AND_TIME(VALUES = DATE_start)

    time_start = dble(DATE_start(7)) + dble(DATE_start(8))/1000D0;   
 !  CALL CPU_TIME(time_start)
    
 !--------------------------------------------------------------
 ! Location of singularities and corresponding solution values:
 !--------------------------------------------------------------

    x1 = +1.0D0;   ! <--- x-coordinate of 1st singularity
    y1 = +0.5D0;   ! <--- y-coordinate of 1st singularity

    x2 = -1.5D0;   ! <--- x-coordinate of 2nd singularity
    y2 = +1.0D0;   ! <--- y-coordinate of 2nd singularity

    x3 = +0.5D0;   ! <--- x-coordinate of 3rd singularity
    y3 = -1.5D0;   ! <--- y-coordinate of 3rd singularity


    b1 = 6D0;

    b2 = 5D0;

    b3 = 1D0;

    b = (b1 + b2 + b3)/3D0 + 0.5D0;   ! expected far-field solution value is:
                                      ! (b1 + b2 + b3)/3

 !--------------------------------------------------------------


    i1 = -M; i2 = -M; i3 = -M;
    j1 = -N; j2 = -N; j3 = -N;

    DO i = -M, M
       x(i) = i*h;
       IF ( x(i) < x1 + h/100D0 ) i1 = i;    ! <--- i-index of 1st singularity
       IF ( x(i) < x2 + h/100D0 ) i2 = i;    ! <--- i-index of 2nd singularity
       IF ( x(i) < x3 + h/100D0 ) i3 = i;    ! <--- i-index of 3rd singularity
    ENDDO

    DO j = -N, N
       y(j) = j*h;
       IF ( y(j) < y1 + h/100D0 ) j1 = j;    ! <--- j-index of 1st singularity
       IF ( y(j) < y2 + h/100D0 ) j2 = j;    ! <--- j-index of 2nd singularity
       IF ( y(j) < y3 + h/100D0 ) j3 = j;    ! <--- j-index of 3rd singularity
    ENDDO

 !--------------------------------------------------
 !
 !             INITIAL APPROXIMATION
 !                     to
 !             STEADY STATE SOLUTION
 !
 !--------------------------------------------------

    r12 = sqrt( (x1 - x2)**2 + (y1 - y2)**2 );
    r13 = sqrt( (x1 - x3)**2 + (y1 - y3)**2 );
    r23 = sqrt( (x2 - x3)**2 + (y2 - y3)**2 );

    r0 = 1D0/2D0 * min(r12, r13, r23);

 !--------------------------------------------------

 ! Computing the mass of each initial cone
 ! relatively to the base plane z = b:

    m1 = 1D0/3D0 * pi*r0**2 * abs( b1 - b );
    m2 = 1D0/3D0 * pi*r0**2 * abs( b2 - b );
    m3 = 1D0/3D0 * pi*r0**2 * abs( b3 - b );

 ! The radius of the support of each initial cone
 ! for large t is asymptotically given by
 ! Rj = ccj * t^(1/(n*(p-2)+p))
 ! where ccj is given below (using Barenblatt solution):

    expoente = p/(p + nn*(p-2)) * (p-2)/(p-1);
    common_factor = ( 1/(nn*(p-2)+p)**(nn/p) * p/(p-1) / (2*pi) * &
        ((p-2)/p)**(nn*(p-1)/p) / beta( nn*(p-1)/p, (2*p-3)/(p-2) ) )**expoente;

    C1 = common_factor * m1**expoente;
    C2 = common_factor * m2**expoente;
    C3 = common_factor * m3**expoente;
    
    cc1 = 1/(nn*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C1**((p-1)/p); 
    cc2 = 1/(nn*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C2**((p-1)/p); 
    cc3 = 1/(nn*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C3**((p-1)/p); 

 ! Estimating time need to reach the boundary of
 ! computational region [-Rx, Rx] x [-Ry, Ry]:

    Rx = M*h; 
    Ry = N*h;
    
    minimum_distance_to_travel_1 = min(Rx-x1, Rx+x1, Ry-y1, Ry+y1);
    minimum_distance_to_travel_2 = min(Rx-x2, Rx+x2, Ry-y2, Ry+y2);
    minimum_distance_to_travel_3 = min(Rx-x3, Rx+x3, Ry-y3, Ry+y3);

    maximum_distance_to_travel_1 = max(Rx-x1, Rx+x1, Ry-y1, Ry+y1);
    maximum_distance_to_travel_2 = max(Rx-x2, Rx+x2, Ry-y2, Ry+y2);
    maximum_distance_to_travel_3 = max(Rx-x3, Rx+x3, Ry-y3, Ry+y3);
    
    t1_min = ( cc1 * minimum_distance_to_travel_1 )**(nn*(p-2)+p);
    t2_min = ( cc2 * minimum_distance_to_travel_2 )**(nn*(p-2)+p);
    t3_min = ( cc3 * minimum_distance_to_travel_3 )**(nn*(p-2)+p);

    t1_max = ( cc1 * maximum_distance_to_travel_1 )**(nn*(p-2)+p);
    t2_max = ( cc2 * maximum_distance_to_travel_2 )**(nn*(p-2)+p);
    t3_max = ( cc3 * maximum_distance_to_travel_3 )**(nn*(p-2)+p);

    write(6,*) '-----------------------------------------------------------'
    write(6,*) 'Estimated minimum time to reach the computational boundary:' 
    write(6,100) t1_min, t2_min, t3_min
    
100 FORMAT( 1X, E15.6 )
   
    write(6,*) 'Estimated time to pass the computational boundary entirely:' 
    write(6,100) t1_max, t2_max, t3_max
    write(6,*) '-----------------------------------------------------------'

 !--------------------------------------------------
 !  Setting the values of initial solution u0:
 !--------------------------------------------------
 
    u0 = b;

    DO i = -M, M
       DO j = -N, N
        
        d1 = sqrt( (x(i) - x1)**2 + (y(j) - y1)**2 );
        
        d2 = sqrt( (x(i) - x2)**2 + (y(j) - y2)**2 );
        
        d3 = sqrt( (x(i) - x3)**2 + (y(j) - y3)**2 );
        
        if ( d1 < r0 ) then
           u0(i,j) = b1 - (b1 - b)*d1/r0;
        endif
        
        if ( d2 < r0 ) then
           u0(i,j) = b2 - (b2 - b)*d2/r0;
        endif  
        
        if ( d3 < r0 ) then
           u0(i,j) = b3 - (b3 - b)*d3/r0;
        endif 
        
       ENDDO
    
    ENDDO   
       
    u0(i1,j1) = b1;
    u0(i2,j2) = b2;
    u0(i3,j3) = b3;

 !--------------------------------------------------
 !  Finding the indices of far-field zones:
 !--------------------------------------------------

    i_ff1a = - M; i_ff1b = - M;
    i_ff2a = - M; i_ff2b = - M;
    i_ff3a = - M; i_ff3b = - M;

    DO i = -M, M
       if ( x(i) < -0.50*Rx + h/1000D0 ) i_ff1a = i;
       if ( x(i) < -0.70*Rx + h/1000D0 ) i_ff2a = i;
       if ( x(i) < -0.90*Rx + h/1000D0 ) i_ff3a = i;
       if ( x(i) < +0.50*Rx + h/1000D0 ) i_ff1b = i;
       if ( x(i) < +0.70*Rx + h/1000D0 ) i_ff2b = i;
       if ( x(i) < +0.90*Rx + h/1000D0 ) i_ff3b = i;
    ENDDO

    j_ff1a = - N; j_ff1b = - N;
    j_ff2a = - N; j_ff2b = - N;
    j_ff3a = - N; j_ff3b = - N;

    DO j = -N, N
       if ( y(j) < -0.50*Ry + h/1000D0 ) j_ff1a = j;
       if ( y(j) < -0.70*Ry + h/1000D0 ) j_ff2a = j;
       if ( y(j) < -0.90*Ry + h/1000D0 ) j_ff3a = j;
       if ( y(j) < +0.50*Ry + h/1000D0 ) j_ff1b = j;
       if ( y(j) < +0.70*Ry + h/1000D0 ) j_ff2b = j;
       if ( y(j) < +0.90*Ry + h/1000D0 ) j_ff3b = j;
    ENDDO

    length_i1 = i_ff1b - i_ff1a + 1;
    length_i2 = i_ff2b - i_ff2a + 1;
    length_i3 = i_ff3b - i_ff3a + 1;
    
    length_j1 = j_ff1b - j_ff1a + 1;
    length_j2 = j_ff2b - j_ff2a + 1;
    length_j3 = j_ff3b - j_ff3a + 1;


 !-----------------------------------------------------
 !  Computing initial values of u0 on far-field zones:
 !-----------------------------------------------------

    ff1_value_u0  = ( sum( u0(i_ff1a:i_ff1a+10,j_ff1a:j_ff1b) )/(11*length_j1) &
                    + sum( u0(i_ff1b-10:i_ff1b,j_ff1a:j_ff1b) )/(11*length_j1) &
                    + sum( u0(i_ff1a:i_ff1b,j_ff1a:j_ff1a+10) )/(11*length_i1) &
                    + sum( u0(i_ff1a:i_ff1b,j_ff1b-10:j_ff1b) )/(11*length_i1) )/4D0; 

    ff2_value_u0  = ( sum( u0(i_ff2a:i_ff2a+10,j_ff2a:j_ff2b) )/(11*length_j2) &
                    + sum( u0(i_ff2b-10:i_ff2b,j_ff2a:j_ff2b) )/(11*length_j2) &
                    + sum( u0(i_ff2a:i_ff2b,j_ff2a:j_ff2a+10) )/(11*length_i2) &
                    + sum( u0(i_ff2a:i_ff2b,j_ff2b-10:j_ff2b) )/(11*length_i2) )/4D0;

    ff3_value_u0  = ( sum( u0(i_ff3a:i_ff3a+10,j_ff3a:j_ff3b) )/(11*length_j3) &
                    + sum( u0(i_ff3b-10:i_ff3b,j_ff3a:j_ff3b) )/(11*length_j3) &
                    + sum( u0(i_ff3a:i_ff3b,j_ff3a:j_ff3a+10) )/(11*length_i3) &
                    + sum( u0(i_ff3a:i_ff3b,j_ff3b-10:j_ff3b) )/(11*length_i3) )/4D0;

 !--------------------------------------------------

    min_u0 = minval( u0 );
    max_u0 = maxval( u0 );

    mass_u0 = sum( u0 - b )*h**2;


 !******************************************************************
 !*                                                                *
 !*     The following computation determines a starting value      *
 !*     for the quantities  variation_mass, variation_sup:         *
 !*                                                                *
 !******************************************************************     

 !----------------------------------------------------------
 !      Extending solution values to extended grid:
 !----------------------------------------------------------

    u(-M:M,-N:N) = u0; 
    
    u(-M-1,-N:N) = u(-M,-N:N);
    u( M+1,-N:N) = u( M,-N:N);
    u(-M-1:M+1,-N-1) = u(-M-1:M+1,-N);
    u(-M-1:M+1, N+1) = u(-M-1:M+1, N);

 !----------------------------------------------------------
 
    q = (p - 2D0)/2D0;
 
    h_to_pm2 = h**(p-2); 
    
    time_next_dump = dt_dump;

    dt = cfl*h*h;
    
    t = 0;
    
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

 !******************************************************************
 !*                                                                *
 !*               Computations of u COMPLETED                      *
 !*                                                                *
 !******************************************************************
 
 !--- Computing starting values of  variation_mass, variation_sup:
   
    variation_mass_0 = sum( u(-M:M,-N:N) - u0 )*h**2;

    variation_sup_0 = maxval( abs( u(-M:M,-N:N) - u0 ) );
   
 !--- Final computations:

    CALL DATE_AND_TIME(VALUES = DATE_finish)

    time_finish = dble(DATE_finish(7)) + dble(DATE_finish(8))/1000D0;
 !  CALL CPU_TIME(time_finish)

    total_elapsed_time = time_finish - time_start;   ! <--- unit: seconds
    
    
    write(6,*) '-----------------------------------------------------------'
    write(6,*) 'este eh o NOVO initial_approximation.f90'
    write(6,*) '-----------------------------------------------------------'
 !--------------------------------------------------

 !--------------------------------------------------
 !  Saving important data to output file
 ! (to be used by code u24h_serial.f90):
 !--------------------------------------------------

    OPEN(10, FILE = 'u24h_serial_INPUT.dat', ACTION = 'WRITE')
    
510 FORMAT( 10(1X,I10) )
520 FORMAT( 10(1X,E25.16) )

    write(10,510) M, N
    write(10,510) nn
    
    write(10,520) p
    write(10,520) h
    
    write(10,520) t0, tF, dt_dump

    write(10,520) x1, x2, x3
    write(10,520) y1, y2, y3
    write(10,510) i1, i2, i3
    write(10,510) j1, j2, j3

    write(10,520) b1, b2, b3
    write(10,520) b

    write(10,510) count
    write(10,510) no_runs
    
    write(10,510) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
    write(10,510) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
    write(10,510) length_i1, length_i2, length_i3
    write(10,510) length_j1, length_j2, length_j3

    write(10,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

    write(10,520) min_u0, max_u0
    write(10,520) mass_u0

    write(10,520) t1_min, t2_min, t3_min
    write(10,520) t1_max, t2_max, t3_max

    write(10,520)

 ! Saving LOG variable values:

    write(10,510) (DATE_start (i), i = 1,8)
    write(10,510) (DATE_finish(i), i = 1,8)
    
    write(10,520) tF

    write(10,520) ff1_value_u0
    write(10,520) ff2_value_u0
    write(10,520) ff3_value_u0

    write(10,520) variation_mass_0
    write(10,520) variation_sup_0
    
    write(10,520) min_u0
    write(10,520) max_u0
    
    write(10,520) total_elapsed_time
    
    write(10,520)
    
    write(10,520) (x(i), i = -M, M)
    write(10,520)
    write(10,520) (y(j), j = -N, N)

    write(10,520)
    
    DO i = -M, M
       write(10,520) (u0(i,j), j = -N, N)
    ENDDO
    
    CLOSE(10)
    
 !--------------------------------------------------------------------------
 !--------------------------------------------------------------------------

 CONTAINS 

     FUNCTION beta(x,y)
     
       real*8  x, y, beta

       beta = gamma(x)*gamma(y)/gamma(x+y);

     END FUNCTION beta

 END PROGRAM initial_approximation

 !--------------------------------------------------------------------------
 !--------------------------------------------------------------------------
