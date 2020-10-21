        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |  MDBDC - THE MULTIOBJECTIVE DOUBLE BUNDLE METHOD FOR NONSMOOTH DC OPTIMIZATION   | | 
        !| |                                 (version 2)                                      | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                       by Kaisa Joki (last modified June 2020)                    | |
        !| |                                                                                  | |
        !| |      Features :                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |           * Possibility to use simple stepsize determination after               | |
        !| |             each 'main iteration'.                                               | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |           * During each round of 'main iteration' OpenMP can be utilized         | |
        !| |             to calculate subproblems in parallel. However if you DO NOT          | |
        !| |             WANT to use this feature then                                        | |
        !| |                1) in Makefile DELETE '-fopenmp' from line 5                      | |
        !| |                2) in mdbdc.f95 COMMENT lines 615-618                             | |
        !| |                                                                                  | |         
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |     The software is free for academic teaching and research purposes but I       | |
        !| |     ask you to refer the reference given below, if you use it.                   | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|    Utilizes new version of PLQDF1 by Ladislav Luksan as a quadratic solver.          |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !|                                                                                      |
        !|                                                                                      |
        !|   Codes include:                                                                     |
        !|                                                                                      |
        !|   tmdbdc.f95         - Main program for MDBDC (this file)                            |
        !|   constants.f95      - Double precision (also some parameters)                       |
        !|   bundle1.f95        - Bundle of DC component f_1                                    |
        !|   bundle2.f95        - Bundle of DC component f_2                                    |
        !|   functions.f95      - User-specified DC components f_1 and f_2 together with        |
        !|                        subgradients of DC components. Contains also user-specified   |
        !|                        initial values for parameters                                 |
        !|   mdbdc.f95          - MDBDC method                                                  |
        !|                                                                                      |
        !|   plqdf1.f           - Quadratic solver by Ladislav Luksan                           |
        !|                                                                                      |
        !|   Makefile           - Makefile                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|   To USE the software MODIFY   tmdbdc.f95   and   functions.f95   as needed          |
        !|                                                                                      |
        !|                                                                                      |
        !|   References:                                                                        |
        !|                                                                                      |
        !|   [1] Outi Montonen and Kaisa Joki: "Bundle-based descent method for nonsmooth       |
        !|       multiobjective DC optimization with inequality constraint". Journal of         |
        !|       Global Optimization, 72 (2018), 403–429.                                       |
        !|       https://doi.org/10.1007/s10898-018-0651-0                                      |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*       


      PROGRAM tmdbdc
      
         USE constants, ONLY : dp   ! double precision (i.e. accuracy)
         USE functions              ! INFORMATION from the USER
         USE bundle1                ! The BUNDLE of the DC component f_1
         USE bundle2                ! The BUNDLE of the DC component f_2
         USE mdbdc                  ! MDBDC method
         
        IMPLICIT NONE   
        
        ! 'user_n' is the number of variables in the problem (USER specifies this in MODULE functions.f95)
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: x_0           ! The starting point in the MDBDC method (specified by USER)
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: x_solution    ! The solution obtained to the problem
        
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f_solution    ! The objective function value at the solution 'x_solution'
        
        REAL(KIND=dp) :: time    ! The elapsed time
        
        INTEGER, DIMENSION(8) :: counter    ! Contains the values of different counteres:
                                            !   counter(1) = iter_counter:         the number of 'main iterations' executed
                                            !   counter(2) = subprob_counter:      the number of subproblems solved (in 'main iteration')
                                            !   counter(3) = f_counter:            the number of function values evaluated for the DC components (in 'main iteration')
                                            !   counter(4) = subgrad1_counter:     the number of subgradients calculated for f_1 (in 'main iteration')
                                            !   counter(5) = subgrad2_counter:     the number of subgradients calculated for f_2 (in 'main iteration')
                                            !--------------------------------------------------------------------------------------------------------------------------                 
                                            !   counter(6) = stop_cond_counter:    the number of times 'Clarke stationary algorithm' is executed 
                                            !   counter(7) = clarke_f_counter:     the number of function values evaluated for f (in 'Clarke stationary algorithms')
                                            !   counter(8) = clarke_sub_counter:   the number of subgradients caluculated for f (in 'Clarke stationary algorithms')
 
 
        INTEGER :: iprint                   ! Variable that specifies print option (specified by USER): 
                                            !   iprint = 0 : print is suppressed
                                            !   iprint = 1 : basic print of final result 
                                            !   iprint = -1: basic print of final result (without the solution vector)
                                            !   iprint = 2 : extended print of final result 
                                            !   iprint = -2: extended print of final result (without the solution vector)
                                            !   iprint = 3 : basic print of intermediate results and extended print of final results
                                            !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                            !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                            !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                            !   iprint = 5 : prints each step of the MDBDC algorithm (i.e. everything is printed step by step) (NOT recommended)
                                            !
                                            ! If 'iprint' <= -5 .OR. 'iprint' >= 6 then DEFAULT value 'iprint'=1 is used    
                                            
        
        INTEGER :: mit                      ! The maximum number of 'main iterations' (specified by USER).
                                            ! If 'mit' <=0 then DEFAULT value 'mit'=1000 is used

        INTEGER :: mrounds                  ! The maximum number of rounds during one 'main iteration' (specified by USER).
                                            ! If 'mrounds' <=0 then DEFAULT value 'mrounds'=5000 is used

        INTEGER :: mrounds_clarke           ! The maximum number of rounds during one 'Clarke stationarity' algorithm (specified by USER).
                                            ! If 'mrounds_clarke' <=0 then DEFAULT value 'mrounds_clarke'=5000 is used
                                            
        INTEGER :: termination              ! The reason for termination in MDBDC:
                                            !   1 - the stopping condition is satisfied (i.e. approximate Clarke stationarity)
                                            !   2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < epsilon)
                                            !   3 - the maximum number 'mrounds' of rounds is executed in one main iteration
                                            !   4 - the maximum number of 'main iterations' is executed 
                                            !   5 - the maximum number 'mrounds_clarke' of rounds is executed in one 'Clarke stationary' algorithm

        LOGICAL :: stepsize_used            ! If .TRUE. then simple stepsize determination is done in MDBDC after each 'main iteration' (specified by USER).
        LOGICAL :: scale_func               ! If .TRUE. then tha scaling procedure is performed at the beginning of MDBDC
        
        INTEGER :: i,j,k

        !--------------------------------------------------------------------------------------------------------------------------------
        !      *****   WHAT IS NEEDED TO DEFINE PROBLEM !    *****
        !--------------------------------------------------------------------------------------------------------------------------------
        INTEGER :: number_of_obj                           ! The number of objectives in the multiobjective problem
                                                            ! IF 'number_of_obj'=1 THEN we have a single objective case.
        INTEGER :: number_of_const                         ! The number of constraintss in the multiobjective problem
                                                            ! IF 'number_of_const'=0 THEN we have an unconstrained problem    
                                                                        
        INTEGER :: number_of_func                          ! The sum of objective and constraint functions                                                               

        !---------------------------------------------------------------------------------------------------------------------------------
        ! DC component lists where the first 'number_of_obj' places tell the DC components f1 for the objectives
        ! and the latter 'number_of_const' places tell the DC components g1 for the constraints         
        INTEGER, DIMENSION(:), ALLOCATABLE :: problem10      ! First is DC components f1 of objective functions used and then DC components g1 of constraint functions used
         
        INTEGER, DIMENSION(:), ALLOCATABLE :: problem20      ! DC components f2 and g2 of objective and constraint functions used 
        
        INTEGER :: startpoint                                ! Starting point use

        !--------------------------------------------------------------------------------------------------------------------------------
      
        INTEGER :: user_n                        ! The number of variables in the problem
        
        INTEGER :: dim_loop                      ! Defines the dimension(s) in the considered problem 
                                                 !    1 - the dimension is 2
                                                 !    2 - the dimension is 4
                                                 !    3 - the dimension is 10
                                                 !    4 - the dimensions are 10, 50, 100, 250 and 500
        
        INTEGER :: ind_start
        INTEGER :: ind_finish        
        
        CHARACTER(LEN=80) :: outfile1 

        outfile1 = 'results.txt'
        OPEN(46,file=outfile1)    
        
        !--------------------------------------------------------------------------------------------------------------------------------
 
       WRITE(46,*)  ' Problem ', ' user_n ', ' | ', ' obj ', ' obj2 ', ' obj3 ', ' | ',  &
                      & ' n_f_obj | ',& 
                      & ' n_sub | ', &
                     ! & ' n_sub2 | ', &
                      & ' n_f_escape_H ', ' n_sub_escape_H ', ' times_escape_H ', ' | ', &
                      & ' time ', ' n_iter ', ' termination ', ' | ', &
                      & ' total_f ', &
                      & ' total_sub | '

                      
        
      ! Different test problems are looked through
        DO k = 1, 21
        
            SELECT CASE(k)

              !-------------------------------------
              !           Problem   1
              !-------------------------------------           
               CASE(1)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                  
                ! The functions used 
                  problem10 = (/ 2,6 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 1
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1

              !-------------------------------------
              !           Problem   2
              !-------------------------------------           
               CASE(2)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 2
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1


              !-------------------------------------
              !           Problem   3
              !-------------------------------------           
               CASE(3)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 3
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1
                           
              !-------------------------------------
              !           Problem   4
              !-------------------------------------           
               CASE(4)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 6,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 4
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1

              !-------------------------------------
              !           Problem   5
              !-------------------------------------           
               CASE(5)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,9 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 5
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2


              !-------------------------------------
              !           Problem   6
              !-------------------------------------           
               CASE(6)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,13 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 6
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 3


              !-------------------------------------
              !           Problem   7
              !-------------------------------------           
               CASE(7)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 12,13 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 7
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 3


              !-------------------------------------
              !           Problem   8
              !-------------------------------------           
               CASE(8)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 8
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4
                  !dim_loop = 6

              !-------------------------------------
              !           Problem   9
              !-------------------------------------           
               CASE(9)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,12 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 9
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4

              !-------------------------------------
              !           Problem   10
              !-------------------------------------           
               CASE(10)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 14,15 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 10
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4
                  !dim_loop = 6

              !-------------------------------------
              !           Problem   11
              !-------------------------------------           
               CASE(11)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,7 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 11
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1

              !-------------------------------------
              !           Problem   12
              !-------------------------------------           
               CASE(12)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,4,9 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 12
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2

              !-------------------------------------
              !           Problem   13
              !-------------------------------------           
               CASE(13)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10,12 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 13
                  
                ! The dimension(s) of the problem looked through 
                 dim_loop = 4
                  !dim_loop = 3

              !-------------------------------------
              !           Problem   14
              !-------------------------------------           
               CASE(14)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10,16 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 14
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4
                          
              !-------------------------------------
              !           Problem   15
              !-------------------------------------           
               CASE(15)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 0
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,14,15 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 15
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4

              !-------------------------------------
              !           Problem   16
              !-------------------------------------           
               CASE(16)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,7,1 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 16
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1

              !-------------------------------------
              !           Problem   17
              !-------------------------------------           
               CASE(17)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,9,2 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 17
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2

              !-------------------------------------
              !           Problem   18
              !-------------------------------------           
               CASE(18)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 10,12,3 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 18
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4

              !-------------------------------------
              !           Problem   19
              !-------------------------------------           
               CASE(19)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,7,1 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 19
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1

              !-------------------------------------
              !           Problem   20
              !-------------------------------------           
               CASE(20)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 3,4,9,2 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 20
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 2

              !-------------------------------------
              !           Problem   21
              !-------------------------------------           
               CASE(21)
            
                ! The number of objectives 
                  number_of_obj = 3        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 4,10,16,3 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 21
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 4                  
                  
             !-------------------------------------
              !           Problem   22
              !-------------------------------------           
               CASE(22)
            
                ! The number of objectives 
                  number_of_obj = 2        
                ! The number of constraints               
                  number_of_const = 1
                ! The total number of objectives and constraints  
                  number_of_func = number_of_obj + number_of_const
                  
                  ALLOCATE(problem10(number_of_func),problem20(number_of_func))
                ! The functions used 
                  problem10 = (/ 2,6,4 /)
                  problem20 = problem10
                  
                ! The starting point   
                  startpoint = 22
                  
                ! The dimension(s) of the problem looked through 
                  dim_loop = 1                    
 
            END SELECT              
            
            ALLOCATE(f_solution(number_of_func))
            
            CALL allocate_prob_data(number_of_obj, number_of_const, number_of_func, &
                  & problem10, problem20, startpoint)

            SELECT CASE(dim_loop)   
 
                CASE(1) ! (n=2)
                    ind_start = 1
                    ind_finish = 1
                    
                CASE(2) ! (n=4)
                    ind_start = 2
                    ind_finish = 2
                    
                CASE(3) ! (n=10)
                    ind_start = 3
                    ind_finish = 3
                    
                CASE(4) ! (n= 10, 50, 100, 250, 500)
                    ind_start = 3
                    ind_finish = 7 
                    
                CASE(5) ! (n= 250)
                    ind_start = 3
                    ind_finish = 6

                CASE(6) ! (n= 10, 50, 100, 250)
                    ind_start = 3
                    ind_finish = 6

            END SELECT          
            
            DO j = ind_start, ind_finish
            
                SELECT CASE(j)
                   
                    CASE(1)
                       user_n = 2
                    
                    CASE (2)
                       user_n = 4
                    
                    CASE (3)
                       user_n = 10
                    
                    CASE (4)
                       user_n = 50
                    
                    CASE (5)
                       user_n = 100
                    
                    CASE (6)
                       user_n = 250
                    
                    CASE (7)
                       user_n = 500
                           
                END SELECT
                              
                CALL allocate_bundle_sizes(user_n)  
              
                ALLOCATE(x_0(user_n),x_solution(user_n))
          
                SELECT CASE(startpoint)    

                        
            ! Problem 1.  
            ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                      CASE(1)
                        x_0 = 1.0_dp
                        x_0(1) = -1.20_dp
             
            ! Problem 2.  
            ! x_0 = (/ -0.5_dp, 1.0_dp /)   
                      CASE(2)
                        x_0 = 1.0_dp
                        x_0(1) = -0.5_dp 
                        
            ! Problem 3.
            ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                      CASE(3)
                        x_0 = 1.0_dp
                        x_0(1) = -1.20_dp
                        
            ! Problem 4.
            ! x_0 = (/ -2.0_dp, 1.0_dp /)
                      CASE(4)
                        x_0 = 1.0_dp 
                        x_0(1) = -2.0_dp 
                        
            ! Problem 5.
            ! x_0= (/ 4.0_dp, 2.0_dp, 4.0_dp, 2.0_dp/)
                      CASE(5)
                        DO i = 1, user_n
                           IF ( (i == 1) .OR. (i == 3)) THEN
                              x_0(i) = 4.0_dp
                           ELSE
                              x_0(i) = 2.0_dp              
                           END IF
                        END DO              
                        
            ! Problem 6.
                      CASE(6)
                        DO i = 1, user_n
                           x_0(i)  = 0.1_dp * i
                        END DO    

            ! Problem 7.
                      CASE(7)
                        DO i = 1, user_n
                           x_0(i) = 2.0_dp * i 
                        END DO                  
                        
            ! Problem 8.
                      CASE(8)
                        DO i = 1, user_n
                           x_0(i)  = 0.1_dp * i
                        END DO              
                        
            ! Problem 9.
                      CASE(9)
                        DO i = 1, user_n
                           x_0(i) = 2.0_dp * i 
                        END DO          
                            
            ! Problem 10.
                      CASE(10)
                        x_0 = -1.0_dp
                        DO i = 1,user_n, 2
                          x_0(i) = 1.0_dp
                        END DO         

            ! Problem 11.
            ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                      CASE(11)
                        x_0 = 1.0_dp
                        x_0(1) = -1.20_dp   

            ! Problem 12.
            ! x_0 = (/ 1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp /)
                      CASE(12)
                        x_0 = 3.0_dp
                        x_0(1) = 1.0_dp
                        x_0(user_n) = 1.0_dp

            ! Problem 13.
                      CASE(13)
                        DO i = 1, user_n
                           x_0(i) = 0.1_dp * i 
                        END DO                  
                        
            ! Problem 14.
                      CASE(14)
                        DO i = 1, user_n
                           x_0(i) = 0.1_dp * i 
                        END DO                  
                        
            ! Problem 15.
                      CASE(15)
                        DO i = 1, user_n
                           x_0(i) = 0.1_dp * i 
                        END DO                  
                        
            ! Problem 16.  
            ! x_0 = (/ -0.5_dp, 1.0_dp /)   
                      CASE(16)
                        x_0 = 1.0_dp
                        x_0(1) = -0.5_dp    

            ! Problem 17.
            ! x_0= (/ 4.0_dp, 2.0_dp, 4.0_dp, 2.0_dp/)
                      CASE(17)
                        DO i = 1, user_n
                           IF ( (i == 1) .OR. (i == 3)) THEN
                              x_0(i) = 4.0_dp
                           ELSE
                              x_0(i) = 2.0_dp              
                           END IF
                        END DO              

            ! Problem 18.
                      CASE(18)
                        DO i = 1, user_n
                           x_0(i) = 2.0_dp * i 
                        END DO    

            ! Problem 19.
            ! x_0 = (/ -1.20_dp, 1.0_dp /)   
                      CASE(19)
                        x_0 = 1.0_dp
                        x_0(1) = -1.20_dp               

            ! Problem 20.
            ! x_0 = (/ 1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp /)
                      CASE(20)
                        x_0 = 3.0_dp
                        x_0(1) = 1.0_dp
                        x_0(user_n) = 1.0_dp
                        
            ! Problem 21.
                      CASE(21)
                        DO i = 1, user_n
                           x_0(i) = 0.1_dp * i 
                        END DO       
                        
            ! Problem 22.
                      CASE(22)
                        x_0(1) = 0.0_dp
                        x_0(2) = 2.0_dp

                        
                 END SELECT
            !--------------------------------------------------------------------------------------

  
               ! Some parameters of the method 
                  mrounds = 10000           ! maximum number of rounds during one 'main iteration'
                  mit = 10000               ! maximum number of 'main iterations'
                  mrounds_clarke = 10000    ! maximum number of rounds in Clarke stationary algorithm 
                  
                  iprint = 3                ! basic print of intermediate results and extended print of final results
                  
                  stepsize_used = .FALSE.   ! Simple stepsize determination is not used
                  !scale_func = .FALSE.     ! The scaling procedure is NOT performed
                  scale_func = .TRUE.       ! The scaling procedure is preformed
          

                IF (iprint >= 1) THEN 
                WRITE(*,*) '------------------------------------------------------------------'
                WRITE(*,*) '** START ** START ** START ** START ** START ** START ** START **'  
                WRITE(*,*) '------------------------------------------------------------------'         
                WRITE(*,*) ' '
                WRITE(*,*) ' ', 'Problem', k, 'user_n', user_n
                WRITE(*,*) ' '
                WRITE(*,*) ' ', 'Objectives:', problem10
                WRITE(*,*) ' '

                CALL bundle_algorithm( x_0, x_solution, f_solution, mit, mrounds, &
                            & mrounds_clarke, termination, counter,  &
                            & stepsize_used, iprint, scale_func, user_n, startpoint, & 
                            & number_of_obj, number_of_const, number_of_func, time)
                    

                WRITE(*,*) ' '
                WRITE(*,*) '------------------------------------------------------------------'
                WRITE(*,*) '** END ** END ** END ** END ** END ** END ** END ** END ** END **'  
                WRITE(*,*) '------------------------------------------------------------------'
                WRITE(*,*) ' '      
                END IF 


                IF (number_of_obj == 1) THEN 
                   WRITE(46,*)   k, user_n , ' | ', f_solution(1), ' . ', ' . ',  ' | ', &
                              & counter(3), ' | ', & 
                              & counter(4), ' | ', & 
                              !& counter(5), ' | ', &           
                              & counter(7), counter(8), counter(6), ' | ', &                                                    
                              & time, counter(1), termination , ' | ', &
                              & counter(3) + counter(7) , &
                              & Max(counter(4), counter(5)) + counter(8), ' | '
                  
                ELSE IF (number_of_obj == 2) THEN
                   WRITE(46,*)  k, user_n ,' | ', f_solution(1), f_solution(2), ' . ',  ' | ', &
                              & counter(3), ' | ', & 
                              & counter(4), ' | ', & 
                              !& counter(5), ' | ', &           
                              & counter(7), counter(8), counter(6), ' | ', &                                                    
                              & time, counter(1), termination , ' | ', &
                              & counter(3) + counter(7) , &
                              & Max(counter(4), counter(5)) + counter(8), ' | '                  
                ELSE IF (number_of_obj == 3) THEN
                   WRITE(46,*)  k, user_n ,' | ', f_solution(1), f_solution(2), f_solution(3), ' | ', &
                              & counter(3), ' | ', & 
                              & counter(4), ' | ', & 
                              !& counter(5), ' | ', &           
                              & counter(7), counter(8), counter(6), ' | ', &                                                    
                              & time, counter(1), termination , ' | ', &
                              & counter(3) + counter(7) , &
                              & Max(counter(4), counter(5)) + counter(8), ' | ' 
                              
                END IF              
                
                DEALLOCATE(x_solution, x_0)
                
            END DO   ! END: Different dimensions
      
            DEALLOCATE(problem10,problem20,f_solution)
            CALL deallocate_prob_data()
             
        END DO   ! END: Different test problems
        
        WRITE(46,*)
        WRITE(46,*) '------------------------------------------------------------------------------', &
                     & '------------------------------------------------------------------------------', &
                     & '------------------------------------------------------------------------------'
        
        
        CLOSE(46)  
          
      END PROGRAM tmdbdc






















