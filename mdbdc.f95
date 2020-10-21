      MODULE mdbdc
      
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                    MDBDC - THE MULTIOBJECTIVE DOUBLE BUNDLE METHOD               | |
        !| |                            FOR NONSMOOTH DC OPTIMIZATION                         | | 
        !| |                                 (version 2)                                      | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |      Features :                                                                  | |
        !| |                                                                                  | |
        !| |           * Possibility to use simple stepsize determination after               | |
        !| |             each 'main iteration'.                                               | |
        !| |                                                                                  | | 
        !| |           * During each round of 'main iteration' OpenMP can be utilized         | |
        !| |             to calculate subproblems in parallel. However if you DO NOT          | |
        !| |             WANT to use this feature then                                        | |
        !| |                1) in Makefile DELETE '-fopenmp' from line 5                      | |
        !| |                2) in mdbdc.f95 COMMENT lines 615-618                             | |
        !| |                                                                                  | |        
        !| |                                                                                  | |                
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|    Utilizes the new version of PLQDF1 by Ladislav Luksan as a quadratic solver.      |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
 

 
        USE omp_lib
        
        USE constants, ONLY : dp    ! Double precision (i.e. accuracy)
        USE bundle1                 ! The BUNDLE of the DC component f_1
        USE bundle2                 ! The BUNDLE of the DC component f_2
        USE functions               ! Contains INFORMATION from the USER
        
        IMPLICIT NONE    
        
        EXTERNAL PLQDF1             ! The QUADRATIC SOLVER by Ladislav Luksan
        
        CONTAINS
                
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                            CONTAINS SUBROUTINES:                                 | | 
        !| |                                                                                  | |
        !| |      bundle_algorithm    : The double bundle method for nonsmooth                | |
        !| |                            multiobjective DC optimization.                       | | 
        !| |      main_iteration      : The main iteration algorithm needed in the double     | |
        !| |                            bundle method for DC optimization.                    | |
        !| |      guaranteeing_clarke : The algorithm guaranteeing Clarke stationarity        | |       
        !| |                            for a solution obtained.                              | |
        !| |                                                                                  | |       
        !| |      quadratic_solver    : The solver for the quadratic norm minimization        | |       
        !| |                            problem. Needed in the Clarke stationary algorithm    | |       
        !| |      subproblem_solver   : The solver for subproblems in the search direction    | |
        !| |                            problem. Needed in the main iteration algorithm.      | |   
        !| |                                                                                  | |   
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*   
        
    
    
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |           THE DOUBLE BUNDLE ALGORITHM FOR NONSMOOTH DC OPTIMIZATION            |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************
        
           SUBROUTINE bundle_algorithm( x_0, x_solution, f_solution, mit,&
                            & mrounds, mrounds_clarke, termination, counter,  &
                            & stepsize_used, iprint, scale_func, user_n, startpoint, &
                            & number_of_obj, number_of_const, number_of_func, CPU)
            ! 
            ! Solves the unconstrained nonsmooth DC minimization problem
            !
            ! INPUT:  * 'x_0'            : A starting point
            !         * 'mit'            : The maximum number of 'main iterations'      
            !         * 'mrouds'         : The maximum number of rounds during one 'main iteration'
            !         * 'mrouds_clarke'  : The maximum number of rounds during one 'Clarke stationary' algorithm
            !         * 'stepsize_used'  : If .TRUE. then simple stepsize determination is used in the algorithm            
            !         * 'iprint'         : Specifies the print          
            !         * 'user_n'         : The dimension of the problem          
            !
            ! OUTPUT: * 'x_solution' : The solution obtained to the minimizatin problem
            !         * 'f_solution' : The objective function value at the solution 'x_solution'
            !         * 'termination': The cause of the termination in the algorithm
            !         * 'counter'    : Gives the values of different counters
            !
            ! NOTICE: * The dimension of vectors 'x_0' and 'x_solution' has to be 'user_n' defined by USER in MODULE functions.
            !         * The dimension of the vector 'counter' has to be 5.
            !         * 'mit', 'mrounds' and 'mrounds_clarke' have to be integers.
            !         * IF ('mit' <= 0) THEN DEFAULT value 1000 is used
            !         * IF ('mrounds' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * IF ('mrounds_clarke' <= 0 ) THEN DEFAULT value 5000 is used.
            !         * 'iprint' has to be -4, -3, -2, -1, 0, 1, 2, 3 or 4 (5). If it is NOT then DEFAULT value 1 is used. 
            !         * If 'iprint = 5' then everything is printed step by step (this is NOT recommended!!!)
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM USER ************************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_solution  ! the solution obtained to the problem
               
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: f_solution  ! the objective function value at the solution 'x_solution'
               
               REAL(KIND=dp), INTENT(OUT) :: CPU               ! the CPU time
               
               INTEGER, INTENT(INOUT) :: mit                  ! the maximum number of 'main iterations'
               INTEGER, INTENT(INOUT) :: mrounds              ! the maximum number of rounds during one 'main iteration'
               INTEGER, INTENT(INOUT) :: mrounds_clarke       ! the maximum number of rounds during one 'Clarke stationary' algorithm
               
               INTEGER, INTENT(OUT) :: termination        ! 1 - the stopping condition is satisfied (i.e. approximate Clarke stationarity)
                                                          ! 2 - the approximate stopping condition is satisfied (i.e. the step-length beta* < beta_0)
                                                          ! 3 - the maximum number 'mrounds' of rounds executed in one main iteration
                                                          ! 4 - the maximum number of 'main iterations' is executed  
                                                          ! 5 - the maximum number 'mrounds_clarke' of rounds executed in one 'Clarke stationary' alqorithm
                                                          
               INTEGER, DIMENSION(8), INTENT(OUT) :: counter  ! contains the values of different counteres: 
                                                              !   counter(1) = iter_counter         the number of 'main iterations' executed
                                                              !   counter(2) = subprob_counter      the number of subproblems solved
                                                              !   counter(3) = f_counter            the number of function values evaluated for DC component
                                                              !   counter(4) = subgrad1_counter     the number of subgradients calculated for f_1 
                                                              !   counter(5) = subgrad2_counter     the number of subgradients calculated for f_2 
                                                              !--------------------------------------------------------------------------------------------------------------------------                   
                                                              !   counter(6) = stop_cond_counter    the number of times 'Clarke stationary algorithm' is executed 
                                                              !   counter(7) = clarke_f_counter     the number of function values evaluated for f (in 'Clarke stationary algorithms')
                                                              !   counter(8) = clarke_sub_counter   the number of subgradients caluculated for f (in 'Clarke stationary algorithms')
                                                              
               LOGICAL, INTENT(IN) :: stepsize_used     ! .TRUE. if simple stepsize determination is used in the algorithm. Otherwise .FALSE.                                                             
               LOGICAL, INTENT(IN) :: scale_func        ! .TRUE. if scaling procedure is used in the algorithm. Otherwise .FALSE.                                                             
                        
               INTEGER, INTENT(INOUT) :: iprint ! variable that specifies print option:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result 
                                                !   iprint = -1: basic print of final result (without the solution vector)
                                                !   iprint = 2 : extended print of final result 
                                                !   iprint = -2: extended print of final result (without the solution vector)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results
                                                !   iprint = -3: basic print of intermediate results and extended print of final results (without the solution vector)
                                                !   iprint = 4 : extended print of intermediate results and extended print of final results 
                                                !   iprint = -4: extended print of intermediate results and extended print of final results (without the solution vectors)
                                                !   iprint = 5 : prints each step of the bundle algorithm (i.e. everything is printed step by step) 
               
               INTEGER, INTENT(IN) :: user_n        ! The dimension of the problem
               INTEGER, INTENT(IN) :: startpoint    ! The startpoint

               INTEGER, INTENT(IN) :: number_of_obj      ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_const    ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_func     ! The dimension of the problem
           
           !***************************** LOCAL VARIABLES ************************************  

               REAL(KIND=dp), DIMENSION(user_n)  :: x_0         ! the starting point
           
               TYPE(kimppu1) :: B1             ! The bundle B_1 for the DC component f_1
               TYPE(kimppu2) :: B2             ! The bundle B_2 for the DC component f_2

               REAL(KIND=dp), DIMENSION(user_n) :: x_current     ! the current iteration point (the dimension 'user_n' is the number of variables)
               REAL(KIND=dp), DIMENSION(user_n) :: x_new         ! the new iteration point (obtained from the previous 'main iteration')
               
               REAL(KIND=dp), DIMENSION(user_n) :: vect          ! a 'help' vector
               
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad1_k_all         ! the subgradients of the DC components f1 and g1 at x_k 
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad2_k_all         ! the subgradients of the DC components f2 and g2 at x_k 
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad1_new_all       ! the subgradients of the DC components f1 and g1 at x_new   
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad2_new_all       ! the subgradients of the DC components f2 and g2 at x_new   
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_k_all         ! the values of the DC components f1 and g1 at x_k      
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_k_all         ! the values of the DC components f2 and g2 at x_k  
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg_k_all          ! the values of the DC components f and g at x_k  
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_new_all       ! the values of the DC components f1 and g1 at x_new        
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_new_all       ! the values of the DC components f2 and g2 at x_new  
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: factor                ! the scaling factors for f_i and g_l              
               REAL(KIND=dp), DIMENSION(number_of_func) :: factor_old            ! the old scaling factors for f_i and g_l              
               INTEGER, DIMENSION(number_of_func) :: fac_ind            ! the scaling factor indices for f_i and g_l              
               INTEGER, DIMENSION(number_of_func) :: ero_ind            ! the difference between scaling factor indices for f_i and g_l              
               
               REAL(KIND=dp), DIMENSION(number_of_obj) :: f_0_all              ! the values of f_i at x_0   
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: AandB_k_all         ! the values of the components A_i and B_l at x_k                  
               REAL(KIND=dp), DIMENSION(number_of_func) :: AandB_new_all       ! the values of the components A_i and B_l at x_new                    
           
               REAL(KIND=dp), DIMENSION(user_n) :: grad1, grad2  ! subgradients (the dimenison 'user_n' is the length of subgradient)       
               REAL(KIND=dp), DIMENSION(user_n) :: dd            ! the search direction obtained from the 'main iteration': dd=x_new -x_current
                                                                 ! (the dimension 'user_n' is the length of subgradient)               
            
               REAL(KIND=dp) :: c1, c2, c3     ! Decrease parameter used in 'main iteration'
               REAL(KIND=dp) :: crit_tol, eps  ! 'crit_tol'=the stopping tolerance and 'eps'=the proximity measure 
               REAL(KIND=dp) :: m              ! the descent parameter used in 'main iteration'
               REAL(KIND=dp) :: m_extra        ! the extra descent parameter used in 'main iteration'
               REAL(KIND=dp) :: r_dec, r_inc   ! 'r_dec'=the decrease parameter and 'R_inc'=the increase parameter  
               
               REAL(KIND=dp) :: alpha          ! help variable when linearization error is calculated   
               REAL(KIND=dp) :: decrease       ! approximation of decrease
               
               REAL(KIND=dp) :: dec_tol        ! Decraese tolerance used in main iteration                  
               
               REAL(KIND=dp) :: step_tol       ! Step-length tolerance used in 'Clarke stationary' algorithm       
               REAL(KIND=dp) :: m_clarke       ! the descent parameter used in 'Clarke stationary' algorithm
           
               REAL(KIND=dp) :: t                           ! the value of the proximity parameter
               REAL(KIND=dp) :: t_min                       ! the lower bound for the proximity parameter
               REAL(KIND=dp) :: t_max                       ! the upper bound for the proximity parameter
               REAL(KIND=dp) :: norm1                       ! the norm of the subgradient for H_1 at the current iteration point 
               REAL(KIND=dp) :: max_norm                    ! the maximum norm of subgradient in B_2
               INTEGER :: i_t                               ! the index used in the update process of the proximity parameter
               
               REAL(KIND=dp) :: H1_current, H2_current      ! the value of H_1 and H_2 at the current solution  
               REAL(KIND=dp) :: H1_new, H2_new              ! the value of H_1 and H_2 at the new iteration point x_new (obtained from the previous main iteration)        
               REAL(KIND=dp) :: change                      ! 'change' = H(x_new) - H(x_current)  (i.e. the change in the objective function value)
               REAL(KIND=dp) :: change1                     ! 'change1' = H1(x_new) - H1(x_current)  (i.e. the change in the DC component H1)
               REAL(KIND=dp) :: change2                     ! 'change2' = H2(x_new) - H2(x_current)  (i.e. the change in the DC component H2)
               REAL(KIND=dp) :: norm                        ! the distance between two consecutive iteration points        
                   
               REAL(KIND=dp) :: start_time, finish_time                   ! start and finish CPU time
               REAL(KIND=dp) :: start_time_main_it, finish_time_main_it   ! start and finish CPU time in one 'main iteration'
               
               REAL(KIND=dp) :: elapsed_time                  ! elapsed 'clock' time in seconds
               INTEGER :: clock_start, clock_end, clock_rate  ! start and finish 'clock' time   

               INTEGER :: size_b1 , size_b2    ! The biggest possible size of the bundles B_1 and B_2 
               
               ! Parameter used in proximity parameter update procedure
               REAL(KIND=dp) :: u       
               REAL(KIND=dp) :: u_k
               REAL(KIND=dp) :: t_next
               REAL(KIND=dp) :: u_next
               REAL(KIND=dp) :: u_min
               LOGICAL :: extra_dec
                              
               INTEGER  :: iter_counter        ! the number of main iterations executed 
               INTEGER  :: subprob_counter     ! the number of subproblems solved during the execution of the bundle algorithm
               INTEGER  :: stop_cond_counter   ! the number of times 'Clarke stationary' algorithm is used during the algorithm
               INTEGER  :: f_counter           ! the number of function values evaluated for a DC component during 
                                               ! the execution of the bundle algorithm (same for the DC component f_1 and f_2) (without Clarke stationary algorithm)

               INTEGER :: clarke_f_counter     ! the number of function values evaluated for f in 'Clarke stationary' algorithm during the execution of the bundle algorithm
               INTEGER :: clarke_sub_counter   ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm during the execution of the bundle algorithm
               
               INTEGER  :: subgrad1_counter    ! the number of subgradients calculated for f_1 during the bundle algorithm
               INTEGER  :: subgrad2_counter    ! the number of subgradients calculated for f_2 during the bundle algorithm                 
               
               INTEGER :: help_mit_iter_counter    ! the number of iteration rounds executed in one 'main iteration'
               INTEGER :: help_subprob_counter     ! the number of subproblems solved in one 'main iteration'
               INTEGER :: help_f_counter           ! the number of function values evaluated for a DC component in one 'main iteration' (same for f_1 and f_2)
               INTEGER :: help_subgrad1_counter    ! the number of subgradients calculated for f_1 in one 'main iteration'
               INTEGER :: help_subgrad2_counter    ! the number of subgradients calculated for f_2 in one 'main iteration'
               INTEGER :: help_stop_cond_counter   ! the number of times 'Clarke stationarity' was tested during the 'main iteration'
               
               INTEGER :: help_clarke_f_counter     ! the number of function values evaluated for f_i in 'Clarke stationary' algorithm
               INTEGER :: help_clarke_sub_counter   ! the number of subgradients evaluated for f_i in 'Clarke stationary' algorithm
               
               
               INTEGER :: reason_for_stop       ! the reason for stop during the 'main iteration'
                                                ! 0 - a new iteration point found
                                                ! 1 - the stopping condition satisfied (i.e. approximate Clarke stationarity)
                                                ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < beta_0)
                                                ! 3 - the biggest possible number of rounds executed in the 'main iteration'
                                                ! 5 - the biggest possible number of rounds executed in the 'Clarke stationary' algorithm
               
               INTEGER :: max_threads           ! the maximum number of threads that can be used in parallellization
               INTEGER :: threads               ! the number of threads used in parallellization
               INTEGER :: max_sub_prob          ! the maximum number of subproblems solved at the same time
               
               INTEGER :: i, j         
               INTEGER :: ind
               INTEGER :: small_ind
               
               LOGICAL :: stop_alg      ! .TRUE. if the proximal bundle algorithm can be terminated                                             


               CHARACTER*30 outfi/'results.txt'/
               OPEN(40,file=outfi)             
               
               CALL cpu_time(start_time)                ! Start CPU timing     

               CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
               CALL SYSTEM_CLOCK(COUNT=clock_start)     ! Start timing 'clock' time            
               
                
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  

            
           
               !************ PARAMETER VALUES FROM USER ARE LOOKED THROUGH ********************
               !
               !     IF VALUES ARE SET WRONG THEN DEFAULT PARAMETER VALUES ARE USED
               !  
               
               ! ***** descent parameter 'm' *****
               IF ( (user_m<=0.0_dp) .OR. (user_m>=1.0_dp) ) THEN
                   m = 0.2_dp
               ELSE       
                   m = user_m 
               END IF
               
               ! ***** extra descent parameter 'm' *****
               IF ( (user_m_extra<=m) .OR. (user_m_extra>=1.0_dp) ) THEN
                   m_extra = 2.0_dp*m
                   IF (m_extra > 1.0_dp) THEN 
                       m_extra = 0.5_dp*(m+1.0_dp)
                   END IF
               ELSE       
                   m_extra = user_m_extra 
               END IF
               
               ! ***** decrease parameter 'c1' *****
               IF ( (user_c1<=0.0_dp) .OR. (user_c1>=1.0_dp) ) THEN
                   c1 = 0.5_dp
               ELSE       
                   c1 = user_c1 
               END IF       

               ! ***** decrease parameter 'c2' *****
               IF ( (user_c2<-0.0000000000001_dp) .OR. (user_c2>=1.0_dp) ) THEN
                 IF (number_of_obj==2) THEN
                   IF (user_n <= 10) THEN 
                     c2 = 0.4_dp
                   ELSE IF (user_n > 10 .AND. user_n < 100) THEN
                     c2 = 0.25_dp                 
                   ELSE IF (user_n>99 .AND. user_n < 200) THEN
                     c2 = 0.1_dp
                   ELSE IF (user_n > 199 .AND. user_n < 500) THEN    
                     c2 = 0.01_dp
                   ELSE IF (user_n > 499 ) THEN  
                     c2 = 0.001_dp
                   END IF       
                 ELSE 
                   IF (user_n < 100) THEN 
                     c2 = 0.5_dp              
                   ELSE IF (user_n>99 .AND. user_n < 200) THEN
                     c2 = MIN(0.1_dp*(number_of_obj - 1),0.5_dp)
                   ELSE IF (user_n > 199 .AND. user_n < 500) THEN    
                     c2 = MIN(0.01_dp*(number_of_obj - 1),0.5_dp)
                   ELSE IF (user_n > 499 ) THEN  
                     c2 = MIN(0.001_dp*(number_of_obj - 1),0.5_dp)                   
                   END IF   
                 END IF
               ELSE       
                   c2 = user_c2 
               END IF                  
               
               ! ***** decrease parameter 'c3' *****
               IF ( (user_c3<-0.0000000000000001_dp) .OR. (user_c3>=1.0_dp) ) THEN
                   c3 = 0.1_dp        
               ELSE       
                   c3 = user_c3 
               END IF   
                
              
               ! ***** stopping tolerance 'crit_tol' *****
               IF (user_crit_tol <= 0.0_dp .OR. user_crit_tol >= 1.0_dp ) THEN
                    crit_tol = (10.0_dp)**(-5)
               ELSE 
                   crit_tol = user_crit_tol
               END IF
              
               ! ***** proximity parameter 'eps' *****
               IF (user_eps <= 0.0_dp) THEN
                   eps = 0.00005_dp
               ELSE
                   eps = user_eps
               END IF
              
               ! ***** decrease parameter 'r_dec' *****
               IF ((user_r_dec <= 0.0_dp) .OR. (user_r_dec >= 1.0_dp) ) THEN
                   IF( user_n < 10) THEN 
                      r_dec = 0.75_dp
                   ELSE IF ( user_n >= 10 .AND. user_n <300 ) THEN 
                      r_dec = user_n/(user_n+5.0_dp)
                      !WRITE(*,*) r_dec
                   ELSE IF ( user_n >= 300 ) THEN 
                      r_dec = 0.99                    
                   END IF                     
               ELSE
                   r_dec = user_r_dec
               END IF

               ! ***** increase parameter 'r_inc' *****
               IF ( user_r_inc <= 1.0_dp ) THEN
                   r_inc = (10.0_dp)**10
               ELSE
                   r_inc = user_r_inc
               END IF

               ! ***** bundle B1 size *****
               IF ( user_size_b1 <= 0 ) THEN
                   size_b1 = MIN((user_n+5)*number_of_func,1000) 
               ELSE
                   size_b1 = user_size_b1
               END IF

               ! ***** bundle B2 size *****
               IF ( user_size_b2 <= 0 ) THEN
                   size_b2 = 3
               ELSE
                   size_b2 = user_size_b2
               END IF          
               
               !--------------------------------------------------------------------------
               
               ! ***** Descent parameter 'm_clarke' *****
               IF ((user_m_clarke<=0.0_dp) .OR. (user_m_clarke>=1.0_dp)) THEN
                   m_clarke = 0.01_dp
               ELSE
                   m_clarke = user_m_clarke
               END IF
               
               ! ***** Step-length tolerance 'step_tol' *****
               IF (user_step_tol <= 0.0_dp) THEN
                  step_tol = 0.0001_dp               
               ELSE
                   step_tol = user_step_tol
               END IF 
               
               ! ***** Decrease tolerance 'dec_tol' *****
               IF (user_dec_tol > 0.0_dp) THEN
                   IF (user_n <=100) THEN 
                     dec_tol = 0.0_dp
                   ELSE IF (user_n>100 .AND. user_n <250) THEN
                     dec_tol = -(10.0_dp)**(-5)
                   ELSE IF (user_n>=250 ) THEN
                     dec_tol = -(10.0_dp)**(-4)                  
                   END IF     
               ELSE
                   dec_tol = user_dec_tol
               END IF              
               
               !----------------------------------------------------------------------------
               ! ***** maximum number of main iterations 'mit' *****
               IF ( mit <= 0 ) THEN
                   mit = 1000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds' in the 'main iteration' *****
               IF ( mrounds <= 0 ) THEN
                   mrounds = 5000
               END IF   
               
               ! ***** maximum number of rounds 'mrounds' in the 'main iteration' *****
               IF ( mrounds_clarke <= 0 ) THEN
                   mrounds_clarke = 5000
               END IF                  
               
               ! ***** print option 'iprint' *****
               IF ( (iprint < -4 ) .OR. (iprint > 5) ) THEN   !executed if print value is wrong
                   iprint = 1
               END IF              
              !----------------------------------------------------------------------------
              
           !_______________________________________________________________________________
           !************************ STEP 0: PARAMETER INITIALIZATION *********************                
           
               factor = 1.0_dp                          ! the initial scaling factor
               x_current = x_0                          ! the current iteration point is the starting point x_0
 
               DO i = 1, number_of_obj                                      ! the values of DC components f1_i and f2_i at x_0
                  fg1_k_all(i) = f1(x_0,problem1(i),factor(i),user_n) 
                  fg2_k_all(i) = f2(x_0,problem2(i),factor(i),user_n) 
               END DO  
               DO i = number_of_obj+1, number_of_func                       ! the values of DC components g1_l and g2_l at x_0
                  fg1_k_all(i) = g1(x_0,problem1(i),factor(i),user_n) 
                  fg2_k_all(i) = g2(x_0,problem2(i),factor(i),user_n) 
               END DO              
               fg_k_all = fg1_k_all - fg2_k_all
                
               
               small_ind = 1
               DO i = 2, number_of_obj                      
                  IF (ABS(fg_k_all(small_ind))> ABS(fg_k_all(i))) THEN 
                    small_ind = i
                  END IF
               END DO       
               
               IF (scale_func) THEN
                DO i = 1, number_of_obj
                  IF ( ABS(fg_k_all(i)) < (10.0_dp)**(0)) THEN
                    fac_ind(i) = 0                      
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(1)) THEN
                    fac_ind(i) = 1  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(2)) THEN
                    fac_ind(i) = 2  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(3)) THEN
                    fac_ind(i) = 3  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(4)) THEN
                    fac_ind(i) = 4  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(5)) THEN
                    fac_ind(i) = 5  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(6)) THEN
                    fac_ind(i) = 6  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(7)) THEN
                    fac_ind(i) = 7  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(8)) THEN
                    fac_ind(i) = 8  
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(9)) THEN
                    fac_ind(i) = 9
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(10)) THEN
                    fac_ind(i) = 10 
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(11)) THEN
                    fac_ind(i) = 11 
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(12)) THEN
                    fac_ind(i) = 12 
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(13)) THEN
                    fac_ind(i) = 13
                  ELSE IF ( ABS(fg_k_all(i)) < (10.0_dp)**(14)) THEN
                    fac_ind(i) = 14                     
                  END IF    
                END DO

                DO i = 1, number_of_obj
                  ero_ind(i) = fac_ind(i)-fac_ind(small_ind)
                   IF(ero_ind(i) > 1) THEN
                     ero_ind(i) = ero_ind(i)-1
                   END IF
                   factor(i) = (10.0_dp)**(-ero_ind(i))                   
                END DO
               
               END IF

               DO i = 1, number_of_obj                       ! the values of DC components f1_i and f2_i at x_0
                  fg1_k_all(i) = f1(x_0,problem1(i),factor(i),user_n) 
                  fg2_k_all(i) = f2(x_0,problem2(i),factor(i),user_n)
                  f_0_all(i) = fg1_k_all(i) - fg2_k_all(i)   ! the values of the objectives f_i at x_0    
               END DO         
                              
               DO i = number_of_obj+1, number_of_func        ! the values of DC components g1_l and g2_l at x_0
                  fg1_k_all(i) = g1(x_0,problem1(i),factor(i),user_n) 
                  fg2_k_all(i) = g2(x_0,problem2(i),factor(i),user_n)
               END DO                  
               f_counter = 1                                 ! one function value evaluated for each DC component f1_i, f2_i, g1_l and g2_l             
               
               DO i = 1, number_of_obj                       ! the values of subgradients for DC components f1_i and f2_i at x_0
                   grad1 = subgradient_f1(x_0, problem1(i),factor(i),user_n)
                   grad2 = subgradient_f2(x_0, problem2(i),factor(i),user_n)
                   DO j = 1, user_n
                      grad1_k_all(i,j) = grad1(j)
                      grad2_k_all(i,j) = grad2(j)
                   END DO
               END DO  
               
               DO i = number_of_obj+1, number_of_func        ! the values of subgradients for DC components g1_l and g2_l at x_0
                   grad1 = subgradient_g1(x_0, problem1(i),factor(i),user_n)
                   grad2 = subgradient_g2(x_0, problem2(i),factor(i),user_n)
                   DO j = 1, user_n
                      grad1_k_all(i,j) = grad1(j)
                      grad2_k_all(i,j) = grad2(j)
                   END DO
               END DO              
               subgrad1_counter = 1                     ! one subgradient calculated for each DC component f1_i and g1_l
               subgrad2_counter = 1                     ! one subgradient calculated for each DC component f2_i and g2_l
               
               DO i = 1, number_of_func                  ! the values of components A_i and B_l for H_1 at x_0 
                  AandB_k_all(i) = AorB_i(fg1_k_all, fg2_k_all, fg1_k_all, fg2_k_all ,i)
               END DO
               
               H1_current = H1(AandB_k_all)              ! the value of the DC component H_1 at x_0
               H2_current = H2(fg2_k_all)                ! the value of the DC component H_2 at x_0
  
               ind = index_H1(AandB_k_all)                                        ! the index of the component A_i or B_l yielding the subgradient of H_1
               grad1 = subgradient_AorB_i(grad1_k_all, grad2_k_all, ind, user_n)  ! the subgradient of H_1 at x_0
               grad2 = subgradient_H2(grad2_k_all, user_n)                        ! the subgradient of H_2 at x_0
               
               ! The bundles B_1 and B_2 are initialized
               CALL init_bundle_b1(B1, size_b1, user_n)    
               CALL init_bundle_b2(B2, size_b2, user_n)
               
               ! The first bundle element is added into the bundle B_1 and B_2 (i.e. the one corresponding to the starting point)
               CALL add_first_element_b1(B1, grad1, ind)
               CALL add_first_element_b2(B2, grad2)
               
               DO j = 1, number_of_func
                  IF ( j/= ind) THEN 
                     grad1 = subgradient_AorB_i(grad1_k_all, grad2_k_all, j, user_n)        ! the subgradients of A_i and B_l at x_0
                     alpha = 0
                     CALL add_element_b1(B1, grad1, alpha, j)
                  END IF
               END DO
               
               iter_counter = 0             ! the number of 'main iterations' executed so far is zero
               subprob_counter = 0          ! the number of 'subproblems' solved so far is also zero
               stop_cond_counter = 0        ! the number of times 'Clarke stationary' algorithm is used is zero
               stop_alg = .FALSE.           ! we cannot stop the proximal bundle algorithm
               
               clarke_f_counter = 0         ! the initialization of counter for 'Clarke stationary' algorithm
               clarke_sub_counter = 0       ! the initialization of counter for 'Clarke stationary' algorithm
               
               
              ! -- Proximity parameter initialization
                vect = give_subgrad_b1(B1,0)                          ! the vector \bxi_1(x_k)
                norm1 = SQRT(DOT_PRODUCT(vect,vect))                  ! the value of the norm ||\bxi_1(x_k)||
                max_norm = max_norm_value(B2)                         ! the value of the maximum subgradient norm in the bundle B_2
                t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm ) ! the lower bound for t
                t_max = r_inc * t_min                                 ! the upper bound for t
                i_t = 0                                               ! index initialization
                
                IF (user_n < 250) THEN
                  t = 1.0_dp                                            ! the parameter t is selected  
                ELSE
                 IF (number_of_obj == 2) THEN
                    t = 2.0_dp
                 ELSE   
                    t = 10.0_dp
                 END IF 
                END IF

! --- --- --- Needed in OpenMP when we use PARALLELLIZATION --- --- ---   
               max_threads = omp_get_max_threads()
               max_sub_prob = give_max_size_b2(B2)+1
               threads = MIN(max_threads, max_sub_prob)
               CALL omp_set_num_threads(threads)  
! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
               
               IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 5) ) THEN        ! the basic/extended print of indermediate results     

                    WRITE(*,*) 'Main Iter:', 0, 'f(x):', fg1_k_all - fg2_k_all  
                   ! WRITE(*,*)  'x:', x_current  
                   ! WRITE(*,*)  '------------------------------------------------------------------------'  

               END IF
           !_______________________________________________________________________________
           !************************ STEP 0: END ******************************************    
           
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------  
           
            DO WHILE ( ( .NOT. stop_alg) .AND. (iter_counter < mit )  )  ! is repeated until the bundle algorithm can be terminated 
                                                                         
               !_______________________________________________________________________________
               !************************ STEP 1: MAIN ITERATION ******************************* 

                iter_counter = iter_counter + 1                     ! a new 'main iteration' is executed
                
                IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '         MAIN ITERATION', iter_counter, ' STARTS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
                END IF 

                CALL cpu_time(start_time_main_it)
                                
                CALL main_iteration( x_current, H1_current, H2_current,&
                            & fg1_k_all, fg2_k_all, x_new, H1_new, H2_new, &
                            & fg1_new_all, fg2_new_all, f_0_all, AandB_k_all, m_extra, &
                            & B1, B2, c1, c2, c3, t, t_min, t_max, i_t, decrease, & 
                            & crit_tol, eps, m, r_dec, r_inc, m_clarke, step_tol, mrounds,& 
                            & mrounds_clarke, iprint, reason_for_stop, help_mit_iter_counter,&
                            & help_subprob_counter, help_f_counter, help_subgrad1_counter, &
                            & help_subgrad2_counter, help_stop_cond_counter, &
                            & help_clarke_f_counter, help_clarke_sub_counter, &
                            & stepsize_used, extra_dec, factor, dec_tol, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)   

                IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '         MAIN ITERATION', iter_counter, ' ENDS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
                END IF 
                 
                
                CALL cpu_time(finish_time_main_it)          
                
                ! Different counters are updated
                subprob_counter = subprob_counter + help_subprob_counter
                f_counter = f_counter + help_f_counter
                stop_cond_counter = stop_cond_counter + help_stop_cond_counter
                subgrad1_counter = subgrad1_counter + help_subgrad1_counter
                subgrad2_counter = subgrad2_counter + help_subgrad2_counter 
                
                clarke_f_counter = clarke_f_counter + help_clarke_f_counter
                clarke_sub_counter = clarke_sub_counter + help_clarke_sub_counter

               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               IF (reason_for_stop == 0) THEN   ! in this case, a new iteration point is found in the 'main iteration'
               !_______________________________________________________________________________
               !************************ STEP 2: BUNDLE UPDATE*******************************

                 !----------------------------------------------------------                  
                 ! Update of the proximity parameter BEGINS
                 !----------------------------------------------------------

                    
                    u_k = 1.0_dp / t
                    u = u_k
                    u_min = 1.0_dp/t_max
                       
                    IF (extra_dec .AND. i_t>0 ) THEN     
                        u = (H1_new-H2_new-(H1_current-H2_current))
                        u = u / decrease 
                        u = 2.0_dp * u_k * (1-u)
                    ELSE IF ( i_t > 3 ) THEN  
                        u = u_k * 0.5_dp
                    END IF  
                     
                    u_next = MAX(u, 0.1_dp * u_k )   
                    u_next = MAX(u_next, u_min)
                  
                    i_t = MAX(i_t+1,1)
                  
                    IF (u_k /= u_next) THEN 
                       i_t  = 1
                    END IF
                  
                    t_next = 1.0_dp / u_next

                    IF (t_next >= t_max) THEN 
                        t = t_max             
                    ELSE IF (t_next < t_max .AND. t_next >= t_min) THEN 
                        t = t_next
                    ELSE IF (t_next < t_min) THEN  
                        t = t_min                         
                    END IF
 

                 !----------------------------------------------------------                  
                 ! Update of the proximity parameter ENDS 
                 !----------------------------------------------------------    

                    DO i = 1, number_of_func                ! the values of components A_i and B_l for H_1 at x_new 
                       AandB_new_all(i) = AorB_i(fg1_new_all, fg2_new_all, fg1_new_all, fg2_new_all, i)
                    END DO                  
                 
                    H1_new = H1(AandB_new_all)   ! a new function value for H_1
                    H2_new = H2(fg2_new_all)     ! a new function value for H_2
               
                    DO i = 1, number_of_obj                       ! the values of subgradients for DC components f1_i and f2_i at x_0
                      grad1 = subgradient_f1(x_new, problem1(i),factor(i),user_n)
                      grad2 = subgradient_f2(x_new, problem2(i),factor(i),user_n)
                      DO j = 1, user_n
                        grad1_new_all(i,j) = grad1(j)
                        grad2_new_all(i,j) = grad2(j)
                      END DO
                    END DO 

                    DO i = number_of_obj+1, number_of_func        ! the values of subgradients for DC components g1_l and g2_l at x_0
                      grad1 = subgradient_g1(x_new, problem1(i),factor(i),user_n)
                      grad2 = subgradient_g2(x_new, problem2(i),factor(i),user_n)
                      DO j = 1, user_n
                        grad1_new_all(i,j) = grad1(j)
                        grad2_new_all(i,j) = grad2(j)
                      END DO
                    END DO    
                    
                    subgrad1_counter = subgrad1_counter + 1                     ! one subgradient calculated for each DC component f1_i and g1_l
                    subgrad2_counter = subgrad2_counter + 1                     ! one subgradient calculated for each DC component f2_i and g2_l           
               
                    ! Subgradients of H_1 and H_2 at the new point x_new
                    ind = index_H1(AandB_new_all)                   
                    grad1 = subgradient_AorB_i(grad1_new_all, grad2_new_all, ind, user_n)
                    change1 = H1_new - H1_current
                    
                    grad2 = subgradient_H2(grad2_new_all, user_n)
                    change2 = H2_new - H2_current
                    
                    dd = x_new - x_current                              ! the search direction
                    change = H1_new - H2_new - (H1_current -H2_current) ! the change in the objective function value        
                
                    CALL update_b1(B1, grad1, ind, dd, AandB_new_all, AandB_k_all, &  ! bundle update for B_1
                                    & fg1_k_all, fg2_k_all, fg1_new_all, fg2_new_all, &
                                    & number_of_obj, number_of_func )                             
                    CALL update_b2(B2, grad2, dd, change2)                            ! bundle update for B_2
                    
                    
                     DO j = 1, number_of_func                  ! The rest of the new elements are added into the bundle B_1
                       IF ( j/= ind) THEN 
                         grad1 = subgradient_AorB_i(grad1_new_all, grad2_new_all, j, user_n)    ! the subgradient of A_i and B_l at x_new
                         alpha = 0.0_dp
                         CALL add_element_b1(B1, grad1, alpha, j)
                       END IF
                     END DO 
                    
                    
                    x_current = x_new                   ! update of the current iteration point
                    H1_current = H1_new                 ! update of the function value H_1
                    H2_current = H2_new                 ! update of the function value H_2
                    
                    fg1_k_all = fg1_new_all             ! update of the values f1_i and g1_l 
                    fg2_k_all = fg2_new_all             ! update of the values f2_i and g2_l
                    
                    fg_k_all = fg1_k_all -fg2_k_all     ! update of the values f_i and g_l
                    
                    AandB_k_all = AandB_new_all         ! update of the component values A_i and B_l
                    
                    norm = SQRT(DOT_PRODUCT(dd,dd))     ! Distance between two consecutive iteration points
                  
 
            !----------------------------------     
            ! Alteranative stopping condition utilizing the distance between consecutive iteration points
            !      IF ( change > - (10.0_dp)**(-6) ) THEN       ! The difference between objective function values in consecutive steps is small enough
            !      IF ( norm < (10.0_dp)**(-6) ) THEN           ! The difference between iteration points in consecutive steps is small enough
            !       stop_alg = .TRUE.
            !       reason_for_stop = 6
            !     END IF                                    
            !----------------------------------     


                    IF ( (ABS(iprint) >= 3) .AND. (ABS(iprint) <= 5) ) THEN        ! the basic/extended print of indermediate results     
                      IF (number_of_const == 0) THEN 
                        WRITE(*,*) 'Main Iter:', iter_counter, 'f(x):', fg_k_all, 't:', t    
                       ! WRITE(*,*) 'x:', x_current    
                       !WRITE(*,*)  '------------------------------------------------------------------------'  
                        
   
                      ELSE                  
                        WRITE(*,*) 'Main Iter:', iter_counter, 'f(x) and g(x):', fg1_k_all - fg2_k_all , 't:', t
                       ! WRITE(*,*) 'x:', x_current    
                       ! WRITE(*,*)  '------------------------------------------------------------------------'  
                        
                      END IF
                                    
                        IF (iprint > 3) THEN    ! IF iprint =  4 or 5 then the iteration point is printes
                            WRITE(*,*) 'The new iteration point:'
                            DO i = 1, user_n
                                WRITE(*,*) 'x(', i ,')*=',  x_current(i)
                            END DO
                            WRITE(*,*) ' '
                        END IF 
                        
                        IF ( (ABS(iprint) == 4 ) .OR. (ABS(iprint) == 5 ) ) THEN
                            WRITE(*,*) '--------------------------------------------------------------------------- '
                            WRITE(*,*) ' * * * * * DETAILIED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * ' 
                            WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                            WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                            WRITE(*,*) 'the number of function values calculated for DC components:', help_f_counter
                            WRITE(*,*) 'the number of subgradients calculated for DC components:', help_subgrad1_counter
                            WRITE(*,*) 'the number of times Clarke stationary algorithm was tested:',help_stop_cond_counter 
                            
                            IF (help_stop_cond_counter > 0) THEN
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) ' * * * * * * * * DETAILS ABOUT CLARKE STATIONARY ALGORITHM * * * * * * * *  '   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) 'The number of function values calculated for f', help_clarke_f_counter  
                                WRITE(*,*) 'The number of subgradients calculated for f', help_clarke_sub_counter   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                            END IF
                            
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                             
                            WRITE(*,*) '--------------------------------------------------------------------------- '   
                            WRITE(*,*) ' '                          
                        END IF
                        

                        
                    END IF  
                            
               !_______________________________________________________________________________
               !************************ STEP 2: END ****************************************       
               
               ELSE                     ! one of the stopping conditions is fulfilled
               
                    stop_alg = .TRUE.   ! the minimization algorithm can be STOPPED
                    
                    IF((ABS(iprint) >= 4) .AND. (ABS(iprint) <= 5)) THEN ! the extended print of indermediate results
                        WRITE(*,*) ' '
                        IF ( (ABS(iprint) == 4 ) .OR. (ABS(iprint) == 5 ) ) THEN
                            WRITE(*,*) '--------------------------------------------------------------------------- '
                            WRITE(*,*) ' * * * * * DETAILIED INFORMATION ABOUT THE MAIN ITERATION', iter_counter ,' * * * * * '
                            WRITE(*,*) 'the number of rounds needed:', help_mit_iter_counter
                            WRITE(*,*) 'the number of subproblems solved:', help_subprob_counter
                            WRITE(*,*) 'the number of function values calculated for DC components:', help_f_counter
                            WRITE(*,*) 'the number of subgradients calculated for DC components:', help_subgrad1_counter
                            WRITE(*,*) 'the number of times Clarke stationary algorithm was tested:',help_stop_cond_counter 
                            
                            IF (help_stop_cond_counter > 0) THEN    
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) ' * * * * * * * * DETAILS ABOUT CLARKE STATIONARY ALGORITHM * * * * * * * *  '   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                                WRITE(*,*) 'The number of function values calculated for f', help_clarke_f_counter  
                                WRITE(*,*) 'The number of subgradients calculated for f', help_clarke_sub_counter   
                                WRITE(*,*) '--------------------------------------------------------------------------- '   
                            END IF                          
                            
                            WRITE(*,*) '--------------------------------------------------------------------------- '       
                            WRITE(*,*) 'CPU time used:', finish_time_main_it-start_time_main_it                             
                            WRITE(*,*) '--------------------------------------------------------------------------- '
                        END IF      
                        WRITE(*,*) ' '
                        WRITE(*,*) 'During Main iteration round:', iter_counter
                        WRITE(*,*) 'Some of the stopping conditions is fulfilled and the algorithm is stopped.'                     
                        WRITE(*,*) ' '
                    END IF 
                    
               END IF
               
           END DO
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------         
           

           termination = reason_for_stop   ! the cause of the termination
           
           IF ( iter_counter >= mit ) THEN  ! the maximum number of 'main iterations' is executed and this causes the termination
               termination = 4 
           END IF

           factor_old = factor
           factor = 1.0_dp
           
           x_solution = x_current               ! the solution to the minimization problem

           DO i = 1, number_of_func
              f_solution(i) = f1(x_solution,problem1(i),factor(i),user_n)- & 
                  & f2(x_solution,problem2(i),factor(i),user_n) 
           END DO
           
           ! the values of different counters
           counter(1) = iter_counter 
           counter(2) = subprob_counter
           counter(3) = f_counter
           counter(4) = subgrad1_counter
           counter(5) = subgrad2_counter 
           counter(6) = stop_cond_counter
           counter(7) = clarke_f_counter
           counter(8) = clarke_sub_counter
            
           CALL cpu_time(finish_time)         ! Stop CPU timing
           CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing 'clock' time
           
          ! Calculate the elapsed 'clock' time in seconds:
          elapsed_time=(1.0_dp*clock_end-clock_start)/clock_rate        
          CPU = elapsed_time

           IF ( (ABS(iprint) == 1) ) THEN ! basic print of the final result
               WRITE(*,*) '---------------------------------------------------------------------------'        
               WRITE(*,*) '           * * * * * BASIC PRINT OF THE FINAL SOLUTION * * * * * '
               WRITE(*,*) '---------------------------------------------------------------------------'            
               WRITE(*,*) 'The cause of termination: ', termination
               IF (iprint > 0) THEN
                   WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i ,')*=',  x_solution(i)
                   END DO
               END IF
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '               
               WRITE(*,*) 'Objective functions:'
               DO i = 1, number_of_obj
                 WRITE(*,*) 'f(x)', problem1(i), ':', f1(x_solution,problem1(i),factor(i),user_n)- & 
                  & f2(x_solution,problem2(i),factor(i),user_n) 
               END DO                
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ' 
               IF (number_of_const > 0) THEN 
                 WRITE(*,*) 'Constraint functions:'
                 DO i = number_of_obj+1, number_of_func
                   WRITE(*,*) 'g(x)', problem1(i), ':', g1(x_solution,problem1(i),factor(i),user_n)- & 
                   & g2(x_solution,problem2(i),factor(i),user_n) 
                 END DO
               END IF    
               WRITE(*,*) '---------------------------------------------------------------------------' 
               WRITE(*,*) 'CPU time used:', finish_time-start_time                              
               WRITE(*,*) '---------------------------------------------------------------------------'
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'            
           END IF
           
           IF ((ABS(iprint) /= 1) .AND. (iprint /= 0)) THEN  ! extended print of the final result
               WRITE(*,*) '---------------------------------------------------------------------------'        
               WRITE(*,*) '        * * * * * EXTENDED PRINT OF THE FINAL SOLUTION * * * * * '
               WRITE(*,*) '---------------------------------------------------------------------------'            
               WRITE(*,*) 'The cause of termination: ', termination
               IF (iprint > 0) THEN
                   WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
                   DO i = 1, user_n
                       WRITE(*,*) 'x(', i ,')*=',  x_solution(i)
                   END DO
               END IF
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
               WRITE(*,*) 'Objective functions:'
               DO i = 1, number_of_obj
                 WRITE(*,*) 'f(x)', problem1(i), ':', f1(x_solution,problem1(i),factor(i),user_n) & 
                    & -f2(x_solution,problem2(i),factor(i),user_n) 
               END DO                
               WRITE(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ' 
               IF (number_of_const > 0) THEN 
                 WRITE(*,*) 'Constraint functions:'
                 DO i = number_of_obj+1, number_of_func
                   WRITE(*,*) 'g(x)', problem1(i), ':', g1(x_solution,problem1(i),factor(i),user_n) & 
                     & -g2(x_solution,problem2(i),factor(i),user_n) 
                 END DO
               END IF   
               WRITE(*,*) '---------------------------------------------------------------------------' 
               WRITE(*,*) 'scaling:', factor_old 
               WRITE(*,*) '---------------------------------------------------------------------------'                
               WRITE(*,*) 'The number of main iterations:', iter_counter
               WRITE(*,*) 'The number of subproblems solved:', subprob_counter             
               WRITE(*,*) 'The number of function values evaluated for DC components:', f_counter 
               WRITE(*,*) 'The number of subgradients calculated for DC components:', subgrad1_counter 
               WRITE(*,*) '---------------------------------------------------------------------------'
               WRITE(*,*) 'The number of times Clarke stationary algorithm (CSA) was used:', stop_cond_counter  
               WRITE(*,*) 'The number of function values computed for f in CSA:', clarke_f_counter
               WRITE(*,*) 'The number of subgradients calculated for f in CSA: ', clarke_sub_counter
               WRITE(*,*) '---------------------------------------------------------------------------'     
               WRITE(*,*) 'CPU time used:', finish_time-start_time                          
               WRITE(*,*) '--------------------------------------------------------------------------- '
               WRITE(*,*) 'Elapsed time:', elapsed_time
               WRITE(*,*) '---------------------------------------------------------------------------'                        
               !WRITE(*,*) 'total n_f: ', f_counter+clarke_f_counter                       
               !WRITE(*,*) 'total n_xi ', subgrad1_counter+clarke_sub_counter                       
               !WRITE(*,*) '---------------------------------------------------------------------------'   
           END IF 
      
      
      WRITE(40,*) 'The values of the objectives at final point:   '
               DO i = 1, number_of_obj
                 WRITE(40,*) 'f(x)', problem1(i), ':', f1(x_solution,problem1(i),factor(i),user_n) & 
                 & -f2(x_solution,problem2(i),factor(i),user_n) 
               END DO            
        WRITE(40,*) '-----------------------------------------------------------------------------&
         &---------------------------'                     
      WRITE(40,142) subprob_counter
  142  FORMAT(' The total number of subproblems solved:         ',i12)  
      WRITE(40,143) f_counter + clarke_f_counter
  143  FORMAT(' The total number of the objective evaluations:  ',i12)
      WRITE(40,144) subgrad1_counter + clarke_sub_counter
  144  FORMAT(' The total number of gradient evaluations for DC components:',i12)
        WRITE(40,*) '-----------------------------------------------------------------------------&
        &---------------------------'
      WRITE(40,145) f_counter 
  145  FORMAT(' The number of the objective evaluations in main iterations:  ',i12)
      WRITE(40,146) subgrad1_counter 
  146  FORMAT(' The number of gradient evaluations for DC components in main iterations:',i12)
        WRITE(40,*) '-----------------------------------------------------------------------------&
         &---------------------------'      
      WRITE(40,147) clarke_f_counter
  147  FORMAT(' The number of the objective evaluations in Clarke stationary algorithms:  ',i12)
      WRITE(40,148) clarke_sub_counter
  148  FORMAT(' The number of gradient evaluations for DC components in Clarke stationary algortihms:',i12)  
        WRITE(40,*) '-----------------------------------------------------------------------------&
        &---------------------------'       
      WRITE(40,149) finish_time-start_time
  149  FORMAT(' The CPU time:                                 ',f10.4)
      WRITE(40,150) elapsed_time
  150 FORMAT(' The elapsed time:                             ',f10.4)

           
           CLOSE(40)   
           
           END SUBROUTINE bundle_algorithm      
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|


        
        
        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                             THE MAIN ITERATION                                 |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE main_iteration( x_k, H1_k, H2_k, fg1_k_all, fg2_k_all,&
                            & x_new, H1_new, H2_new, fg1_new_all, fg2_new_all, f_0_all,&
                            & AB_k_all, m_extra, B1, B2, c1, c2, c3, t, t_min, t_max, i_t, &
                            & decrease, crit_tol, eps, m, r_dec, r_inc, m_clarke, step_tol, &
                            & mrounds, mrounds_clarke, iprint, reason_for_stop, iter_counter, &
                            & subprob_counter, fi_counter, subgrad1_counter, subgrad2_counter,&
                            & stop_cond_counter, clarke_f_counter, clarke_sub_counter, &
                            & stepsize_used, extra_dec, factor, dec_tol, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)                 
            !
            ! Executes the 'main iteration' algorithm. It is needed to find the search direction at the current iteration point. 
            !
            ! INPUT: * 'x_k'               : The current iteration point
            !        * 'H1_k' and 'H2_k'   : The values of the DC components H_1 and H_2 at the current iteration point
            !        * 'fg1_k_all'         : The values of the DC components f1_i and g1_l at the current iteration point
            !        * 'fg2_k_all'         : The values of the DC components f2_i and g2_l at the current iteration point
            !        * 'AB_k_all'          : The values of the components A_i and B_l at the current iteration point
            !        * 'f_val'             : The value used in descent condition
            !        * 'f_0_all'           : The values of the objective functions f_i at the starting point
            !        * 'factor'            : The values of sclaing factors of the objective functions f_i and constraints g_l
            !        * 'user_n'            : The dimension of the problem
            !
            !        * 'c1', 'c2' and 'c3' : The decrease parameters used in 'main iteration'
            !        * 'crit_tol' and 'eps': The stopping tolerance and the proximity measure
            !        * 'm'                 : The descent parameter
            !        * 'm_etxra'           : The extra descent parameter
            !        * 'r_dec' and 'r_inc' : The decrease and increase parameters   
            !        * 'm_clarke'          : The descent parameter used in 'Clarke stationary' algorithm 
            !        * 'step_tol'          : The step-length parameter used in 'Clarke stationary' algorithm
            !        * 'dec_tol'           : The tolerace for the approximate decrease obtained from the model
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'main iteration'
            !        * 'mrounds_clarke'    : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            !        * 'iprint'            : Specifies print
            !
            !        * 'stepsize_used'     : If .TRUE. then simple stepsize determination is used in the algorithm after each main iteration
            !
            !
            ! OUTPUT: * 'x_new'              : the new iteration point
            !         * 'H1_new' and 'H2_new': the new value of DC components H_1 and H_2
            !         * 'fg1_new_all'        : The values of the DC components f1_i and g1_l at the new iteration point
            !         * 'fg2_new_all'        : The values of the DC components f2_i and g2_l at the new iteration point          
            !         * 'reason_for_stop'    : Indicates the cause of the termination in the 'main iteration' 
            !         * 'decrease '          : approximation of descent in objective H
            !         * 'extra_dec'          : .TRUE. if the decrease in the objective satisfies the extra condtion
            !
            !         * 'iter_counter'       : the number of rounds needed in 'main iteration'
            !         * 'subprob_counter'    : the number of subproblems solved in 'main iteration'
            !         * 'fi_counter'         : the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2)
            !         * 'subgrad1_counter'   : the number of subgradients calculated for f_1 in 'main iteration'
            !         * 'subgrad2_counter'   : the number of subgradients calculated for f_2 in 'main iteration'
            !         * 'stop_cond_counter'  : the number of times approximate stopping condition is tested during the 'main iteration' 
            !         * 'clarke_f_counter'   : the number of function values evaluated for f in 'Clarke stationary' algorithm
            !         * 'clarke_sub_counter' : the number of subgradients evaluated for f in 'Clarke stationary' algorithm
            !
            ! INOUT: * 'B_1' and B_2'       : the bundles of the DC components f_1 and f_2
            !        * 't'                  : the value of the proximity parameter
            !        * 't_min' and 't_max'  : the lower and upper bounds of the proximity parameter 't' 
            !        * 'i_t'                : the index used in the update procedure of the proximity parameter 't' 
            !
            ! NOTICE: * The dimensions of vectors 'x_k' and 'x_new' has to be same. (i.e. the dimensio of 'x_k' and 'x_new' has to be
            !          'user_n' when SUBROUTINE main_iteration() is used in SUBROUTINE bundle_method.
            !         * The dimension of vectors 'fg1_k_all', 'fg2_k_all', 'fg1_new_all', 'fg2_new_all' and AB_k_all' has to be 'number_of_func'    
            !         * The dimention of vector 'f_0_all' has to be 'number_of_obj'         
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER *********************************  
            
               TYPE(kimppu1), INTENT(INOUT) :: B1                    ! the bundle B_1 for the DC component f_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                    ! the bundle B_2 for the DC component f_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k       ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fg1_k_all  ! the values of DC components f1_i at the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fg2_k_all  ! the values of DC components f2_i at the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new     ! the new iteration point if 'reason_for_stop'=0 in the 'main iteration'
               REAL(KIND=dp), INTENT(IN)  :: H1_k, H2_k               ! the value of H_1 and H_2 at x_k
               
               REAL(KIND=dp), INTENT(OUT) :: H1_new, H2_new  ! the value of H_1 and H_2 at x_new if 'reason_for_stop=0' in the 'main iteration'
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg1_new_all  ! the values of DC components f2_i at the new iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg2_new_all  ! the values of DC components f2_i at the new iteration point               
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: f_0_all         ! the values of DC components f_i at the starting point x_0
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AB_k_all        ! the values of components A_i at the starting point x_k
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: factor         ! the values of scaling factors of f_i and g_l
           
               REAL(KIND=dp), INTENT(IN) :: c1, c2, c3    ! the decrease parameters used in main iteration
               REAL(KIND=dp), INTENT(IN) :: crit_tol, eps ! crit_tol=stopping tolerance and eps=proximity measure (needed in stopping conditions)
               REAL(KIND=dp), INTENT(IN) :: m             ! descent parameter
               REAL(KIND=dp), INTENT(IN) :: m_extra       ! extra descent parameter
               REAL(KIND=dp), INTENT(IN) :: r_dec, r_inc  ! r_dec=decrease parameter and r_inc=increase parameter
               REAL(KIND=dp), INTENT(IN) :: m_clarke      ! The descent parameter used in 'Clarke stationary' algorithm 
               REAL(KIND=dp), INTENT(IN) :: step_tol      ! The step-length parameter used in 'Clarke stationary' algorithm
               REAL(KIND=dp), INTENT(IN) :: dec_tol       ! The tolerace for the approximate decrease obtained from the model
               
               REAL(KIND=dp), INTENT(OUT) :: decrease        ! approximation of descent in objective H
               
               REAL(KIND=dp), INTENT(INOUT) :: t             ! the proximity parameter t
               REAL(KIND=dp), INTENT(INOUT) :: t_min, t_max  ! the bounds for proximity parameter t            
               INTEGER, INTENT(INOUT) :: i_t                 ! the index used in update proxcedure of proximity parameter t            
               
               INTEGER, INTENT(IN) :: user_n                          ! the dimension of the problem'
               
               INTEGER, INTENT(IN) :: number_of_obj      ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_const    ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_func     ! The dimension of the problem            
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'main iteration'
                                                          ! If mrounds<=0, then DEFAULT value 500 is used (This is also done in the SUBROUTINE bundle_method()
               INTEGER :: mrounds_clarke                  ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm

                                                          
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (i.e. approximate Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < beta_0)
                                                          ! 3 - maximum number of rounds executed in 'main iteration'
                                                          ! 5 - maximum number of rounds executed in 'Clarke stationary' algorithm
                                                          
               INTEGER, INTENT(OUT) :: iter_counter        ! the number of rounds used in 'main iteration'
               INTEGER, INTENT(OUT) :: subprob_counter     ! the number of subproblems solved in 'main iteration'
               INTEGER, INTENT(OUT) :: fi_counter          ! the number of function values evaluated for a DC component in 'main iteration' (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad1_counter    ! the number of subgradients calculated for f_1 in 'main iteration'
               INTEGER, INTENT(OUT) :: subgrad2_counter    ! the number of subgradients calculated for f_2 in 'main iteration'  
               INTEGER, INTENT(OUT) :: stop_cond_counter   ! the number of times approximate stopping condition is tested during the 'main iteration'  

               INTEGER, INTENT(OUT) :: clarke_f_counter    ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               INTEGER, INTENT(OUT) :: clarke_sub_counter  ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm

               
               INTEGER, INTENT(IN) :: iprint    ! the variable that specifies print:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result (nothing is printed in main iteration)
                                                !   iprint = -1: basic print of final result without the solution vector (nothing is printed in main iteration)                                             
                                                !   iprint = 2 : extended print of final result (nothing is printed in main iteration)
                                                !   iprint = -2: extended print of final result without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results (nothing is printed in main iteration)
                                                !   iprint = -3: basic print of intermediate results and extended print of final results without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 4 : extended print of intermediate and final results (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = -4: extended print of intermediate and final results without the solution vector (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = 5 : everything is printed step by step (this is the only option which causes print in the 'main iteration') 
       
               LOGICAL, INTENT(OUT) :: extra_dec        ! If .TRUE. then the descent satisfies an extra condition
               LOGICAL, INTENT(IN) :: stepsize_used     ! If .TRUE. then simple stepsize determination is used in the algorithm after each main iteration
            !   
               

           !***************************** LOCAL VARIABLES ************************************
               
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad1_y_all            !  subgradients for DC components f1_i and g1_l at y
               REAL(KIND=dp), DIMENSION(number_of_func, user_n) :: grad2_y_all            !  subgradients for DC components f2_i and g2_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_y_all                   ! values of DC components f1_i and g1_l at y               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_y_all                   ! values of DC components f2_i and g2_l at y           
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_y_all                    ! values of components A_i at y 
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_koe_all                 ! values of DC components f1_i and g1_l at x_koe               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_koe_all                 ! values of DC components f2_i and g2_l at x_koe               
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_koe_all                  ! values of components A_i and B_l at x_koe               
               
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: d_t                    ! the search direction
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: y                      ! the new auxilary point
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: new_grad1, new_grad2   ! the subgradients of A_i or B_l and H_2 at y
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: vect                   ! 'help' vector        
               
               REAL(KIND=dp) :: norm1        ! ||\bxi_1(x_k)||
               REAL(KIND=dp) :: max_norm     ! the value of the maximum subgradient norm ||\bxi_{2,max}||
               
               REAL(KIND=dp) :: d_norm             ! the norm of the direction vector d_t
               REAL(KIND=dp) :: delta1, delta2     ! the predicted changes of H_1 and H_2, respectively
               REAL(KIND=dp) :: H1_y, H2_y         ! the function values of H_1 and H_2 at y
               REAL(KIND=dp) :: H_k                ! the function value of H at x_k
               REAL(KIND=dp) :: real_decrease      ! the real decrease of the objective function H
               REAL(KIND=dp) :: lin_err1, lin_err2 ! the linearization errors of A_i or B_l and H_2 (calculated using y)
               
               REAL(KIND=dp) :: erotus           ! the biggest difference obtained when function values at x_0 and at the new iteration point are compared (each f_i and g_l are taken into account)          
               REAL(KIND=dp) :: apu_erotus       ! the help variable used to obtain 'erotus'             
               REAL(KIND=dp) :: help             ! 'help' variable  

               INTEGER :: f_counter              ! help counter for the number of function values
               INTEGER :: sub_counter            ! help counter for the number of subgradients 
               INTEGER :: clarke_iter_counter    ! help counter for the number of subgradients 
               
               INTEGER :: s_counter              ! the number of subproblems solved during the current iteration round of the 'main iteration' algorithm
               INTEGER :: i, j                   ! 'help' variables
               INTEGER :: i_dec                  ! the number of bundle update steps made successively
                           
               LOGICAL :: was_b1_full            ! .TRUE. if during the previous round something was added to B_1 and B_1 was full before this insertion
               LOGICAL :: stop_main_it           ! .TRUE. if the current 'main iteration' can be stopped
               
               LOGICAL :: stop_askelpituus       ! .TRUE. if the current 'line_search' can be stopped
               REAL(KIND=dp) :: askelpituus      ! stepsize determined during 'line search' into descent direction
               REAL(KIND=dp) :: koe_H1           ! the value of f1 at the point tested in 'line search'
               REAL(KIND=dp) :: koe_H2           ! the value of f2 at the point tested in 'line search'
               REAL(KIND=dp), DIMENSION(give_n_b1(B1)) :: koe_y  ! a point tested during 'line search'
                       

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           ! <>  <>  <>  <>  <>  <>  <>  MAIN ITERATION STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           
           !__________________________________________________________________________________     
           !************************ STEP 0: STOPPING CONDITION ******************************
                
               iter_counter = 1                   ! the first round is executed
               subprob_counter = 0                ! the number of subproblems solved
               fi_counter = 0                     ! the number of function values evaluated
               subgrad1_counter = 0               ! the number of subgradients calculated for f_1
               subgrad2_counter = 0               ! the number of subgradients calculated for f_2
               stop_cond_counter = 0              ! the number of times approximate stopping condition is tested
               clarke_f_counter = 0               ! the number of function values evaluated for f in 'Clarke stationary' algorithm
               clarke_sub_counter = 0             ! the number of subgradients evaluated for f in 'Clarke stationary' algorithm
               
               i_dec = 0                          ! the number of bundle update steps made successively
                             
               stop_main_it = .FALSE.             ! We are NOT rady to stop
               
               vect = give_subgrad_b1(B1,0) - give_subgrad_b2(B2,0)       ! vector: ' \bxi_1(x_k) - \bxi_2(x_k) '

               start: IF ( SQRT(DOT_PRODUCT(vect,vect)) < (crit_tol)) THEN  ! ||\bxi_1(x_k) - \bxi_2(x_k)|| <  crit_tol 
               !**>>**>>**>> APPROXIMATE CRITICALITY OBTAINED <<**<<**<<**
                   
                   d_t = 1.0_dp
                   
                   CALL guaranteeing_clarke( x_k, H1_k, H2_k, fg1_k_all, fg2_k_all, &
                            & reason_for_stop, x_new , H1_new, H2_new, fg1_new_all, fg2_new_all,&
                            &  d_t, crit_tol, m_clarke, step_tol, iprint, mrounds_clarke, &
                            & clarke_iter_counter, f_counter, sub_counter, factor, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)
                    
                   clarke_f_counter = clarke_f_counter + f_counter                  
                   clarke_sub_counter = clarke_sub_counter + sub_counter                    
           
                   stop_cond_counter = stop_cond_counter + 1   ! The times that 'Clarke stationary' algorithm has been executed is increased by one
                   stop_main_it = .TRUE.                       ! the 'main iteration' can be stopped
                   extra_dec = .FALSE.
                   i_t = 0
                   

        
                   IF(iprint == 5) THEN  ! print everything        
                       WRITE(*,*) 'At STEP 0: the Clarke stationary algorithm is executed'
                       SELECT CASE(reason_for_Stop)
                         CASE(0)                              
                           WRITE(*,*) ' '
                           WRITE(*,*) 'A BETTER iteration point is found.'
                           WRITE(*,*) ' ' 
                         CASE(1)                              
                           WRITE(*,*) ' '
                           WRITE(*,*) 'The solution satisfies the CLARKE STATIONARITY CONDITION.'
                           WRITE(*,*) ' ' 
                         CASE(2)                              
                           WRITE(*,*) ' '
                           WRITE(*,*) 'The solution satisfies the APPROXIMATE STOPPING CONDITION.'
                           WRITE(*,*) ' ' 
                         CASE(5)                              
                           WRITE(*,*) ' '
                           WRITE(*,*) 'The maximum number of rounds is executed in Clarke stationary algorithm.'
                           WRITE(*,*) ' '                             
                       END SELECT                     
                   END IF

                   
           !__________________________________________________________________________________         
           !************************ STEP 0: END *********************************************         
               ELSE start   
               !_______________________________________________________________________________
               !************************ STEP 1: PARAMETER INITIALIZATION *********************
           
                   H_k = H1_k - H2_k                    ! the function value of H at x_k
                   
                   vect = give_subgrad_b1(B1,0)                           ! the vector \bxi_1(x_k)
                   norm1 = SQRT(DOT_PRODUCT(vect,vect))                   ! the value of the norm ||\bxi_1(x_k)||
                   max_norm = max_norm_value(B2)                          ! the value of the maximum subgradient norm in the bundle B_2
                   t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm )  ! the lower bound for t
                   t_max = r_inc * t_min                                  ! the upper bound for t     
                   
                   ! Checking that the parameter t is from the right interval
                   IF (t < t_min) THEN
                     t = t_min
                   END IF
                   IF (t > t_max) THEN
                     t = t_max 
                   END IF
                   
                   ! print everything   
                   IF(iprint == 5) THEN  
                      WRITE(*,*) ' '                   
                      WRITE(*,*) 'STEP 1: INITIALIZATION OF PARAMETERS'
                      WRITE(*,*) ' '                   
                      WRITE(*,*) 't=', t
                      WRITE(*,*) 't is from the interval [t_min, t_max]=[', t_min, ',', t_max, ']' 
                      WRITE(*,*) 'eps=', eps
                      WRITE(*,*) 'The maximum subgradient norm in B_2: ', max_norm, ' The norm |||\bxi_1(x_k)||: ', norm1 
                      WRITE(*,*) ' '
                      WRITE(*,*) '-------------------------------------------------------------------------------'  
                   END IF                 
            
                   was_b1_full = .FALSE.              ! B_1 was not full during the previous round since this is the first round of 'main iteration'
                       
                                   
                   IF ( mrounds <= 0 ) THEN           
                       mrounds = 5000                 ! the DEFAULT value for 'mrounds' is used
                   END IF
               !_______________________________________________________________________________
               !************************ STEP 1: END ****************************************** 
               
               END IF start

               
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------
           
               DO WHILE ( (.NOT. stop_main_it ) .AND. (iter_counter <= mrounds ))   ! is repeated until the 'main iteration' can be terminated         
               !______________________________________________________________________________
               !************************ STEP 2: SEARCH DIRECTION ****************************
                   
                     
                   CALL subproblem_solver(x_k, AB_k_all, give_n_b1(B1),&  ! subproblems are solved
                                        & give_size_b1(B1)+1, B1, B2, t, s_counter, number_of_func)

 
                   subprob_counter = subprob_counter + s_counter    ! the number of subproblems solved so far
                                   
                   CALL add_glob_index(B2)              ! calculates the index of the subproblem yielding the global solution        
                   d_t = give_solution(B2)              ! the global solution d_t is selected   
                   d_norm = SQRT(DOT_PRODUCT(d_t,d_t))  ! the norm ||d_t|| of the solution d_t is calculated      
                   
           
                   ! print everything                   
                   IF(iprint == 5) THEN
                      WRITE(*,*) ' '
                      WRITE(*,*) 'STEP 2: SEARCH DIRECTION' 
                      WRITE(*,*) ' '                  
                      WRITE(*,*) 'Round:' , iter_counter  , 'paramter t = ', t
                      WRITE(*,*) 'Search direction: [', d_t , ']' 
                      WRITE(*,*) 'The norm of search direction: ', d_norm                     
                      WRITE(*,*) 'The number of subproblems solved during this round: ', s_counter
                      WRITE(*,*) 'The number of subproblems solved till now: ', subprob_counter
                      WRITE(*,*) ' '                  
                      WRITE(*,*) '-------------------------------------------------------------------------------'
                   END IF 
                   
               !______________________________________________________________________________
               !************************ STEP 2: END *****************************************     
                   
                   
               !->->->->->-> EXECUTION OF BRANCH BEGINS (2 POSSIBLE BRANCHES) <-<-<-<-<-<-<-<-
                   branches: IF (d_norm < crit_tol .OR. (give_decrease(B2)-H1_k)>dec_tol  ) THEN                
               !->->->->->->->->->->->->->-> BRANCH 1 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-             
                   !__________________________________________________________________________
                   !******************** STEP 3: CLARKE STATIONARITY CHECK *******************
                       
                        CALL guaranteeing_clarke( x_k, H1_k, H2_k, fg1_k_all, fg2_k_all,&
                            & reason_for_stop, x_new, H1_new, H2_new, fg1_new_all, fg2_new_all,&
                            & d_t, crit_tol, m_clarke, step_tol, iprint, mrounds_clarke, &
                            & clarke_iter_counter, f_counter, sub_counter, factor, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)
                       
                        !Counters are updated
                        clarke_f_counter = clarke_f_counter + f_counter                 
                        clarke_sub_counter = clarke_sub_counter + sub_counter   
                       
                        stop_cond_counter = stop_cond_counter + 1   ! the update of the stopping condition counter
                        
                        stop_main_it = .TRUE.
                        extra_dec = .FALSE.
                        i_t = 0
                               
                           ! print everything
                           IF(iprint == 5) THEN
                              WRITE(*,*) ' '
                              WRITE(*,*) 'STEP 3: Clarke stationarity test' 
                              WRITE(*,*) ' '                  
                              WRITE(*,*) 'The number of the Clarke stationary algorithm' , stop_cond_counter
                              
                              SELECT CASE(reason_for_Stop)
                                CASE(0)                           
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'A BETTER iteration point is found.'
                                  WRITE(*,*) ' ' 
                                CASE(1)                           
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'The solution satisfies the CLARKE STATIONARITY CONDITION.'
                                  WRITE(*,*) ' ' 
                                CASE(2)                           
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'The solution satisfies the APPROXIMATE STOPPING CONDITION.'
                                  WRITE(*,*) ' ' 
                                CASE(5)                           
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'The maximum number of rounds is executed in algorithm.'
                                  WRITE(*,*) ' '                              
                              END SELECT 
                           END IF 

 

                   !__________________________________________________________________________
                   !******************** STEP 3: APPROXIMATE STOPPING CONDITION **************  
                
               !->->->->->->->->->->->->->-> BRANCH 1 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                        
                   ELSE branches
               !->->->->->->->->->->->->->-> BRANCH 2 BEGIN <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                   !__________________________________________________________________________
                   !******************** STEP 4: DESCENT TEST ********************************
       
                       y = x_k + d_t               ! a new auxilary point y
                       
                       DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at y
                          fg1_y_all(i) = f1(y,problem1(i),factor(i),user_n) 
                          fg2_y_all(i) = f2(y,problem2(i),factor(i),user_n) 
                       END DO 
                       
                       DO i = number_of_obj+1, number_of_func   ! the values of DC components g1_l and g2_l at y
                          fg1_y_all(i) = g1(y,problem1(i),factor(i),user_n) 
                          fg2_y_all(i) = g2(y,problem2(i),factor(i),user_n) 
                       END DO                      
                       fi_counter = fi_counter + 1              ! one more objective function value calculated for a DC component          
                                   
                       DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                          AB_y_all(i) = AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
                       END DO                  
                       
                       H1_y = H1(AB_y_all)    ! the value of H_1 at y
                       H2_y = H2(fg2_y_all)   ! the value of H_2 at y
                                   
                       real_decrease = ( H1_y - H2_y ) -  H_k  ! the real decrease in the value of H    
                  
                       
                   ! print everything
                   IF(iprint == 5) THEN
                       WRITE(*,*) ' '
                       WRITE(*,*) 'STEP 4: DESCENT TEST'    
                       WRITE(*,*) ' '                 
                       WRITE(*,*) 'A new AUXILIARY point: [', y, ']'                  
                       WRITE(*,*) ' '                 
                       WRITE(*,*) 'Real decrease:' , real_decrease, ' Predicted decrease: ', give_decrease(B2)            
                       WRITE(*,*) ' Is DECREASE enough: ', real_decrease, '<=', (m *( give_decrease(B2)-H1_k)), '?' 
                       WRITE(*,*) ' ' 
                   END IF
                       
                   
                      branch2: IF ( real_decrease <= (m * (give_decrease(B2)-H1_k)))THEN   ! is the descent in the function H enough?
                       !**>>**>>**>> NEW ITERATION POINT FOUND <<**<<**<<**  
                           
                           IF(stepsize_used) THEN 
                           !-**--**-  STEPSIZE DETERMINATION  -**--**-
                           IF (real_decrease < -0.1_dp) THEN
                              askelpituus = 1.0_dp
                              stop_askelpituus = .FALSE.
                              DO WHILE((.NOT. stop_askelpituus))
                                 koe_y = y + d_t
                                 
                                 DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at koe_y
                                   fg1_koe_all(i) = f1(koe_y,problem1(i),factor(i),user_n) 
                                   fg2_koe_all(i) = f2(koe_y,problem2(i),factor(i),user_n) 
                                 END DO  
                                 DO i = number_of_obj+1, number_of_func   ! the values of DC components g1_l and g2_l at koe_y
                                   fg1_koe_all(i) = g1(koe_y,problem1(i),factor(i),user_n) 
                                   fg2_koe_all(i) = g2(koe_y,problem2(i),factor(i),user_n) 
                                 END DO                                  
                                 fi_counter = fi_counter + 1              ! one more objective function value calculated for a DC component          
                                   
                                 DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                                   AB_koe_all(i) = AorB_i(fg1_koe_all, fg2_koe_all, fg1_k_all, fg2_k_all ,i)
                                 END DO                
                       
                                 koe_H1 = H1(AB_koe_all)    ! the value of H_1 at y
                                 koe_H2 = H2(fg2_koe_all)   ! the value of H_2 at y

                                 IF ( (koe_H1 - koe_H2 - H_k) <= ( 0.7_dp * m * (give_decrease(B2)-H1_k))) THEN  !the case where the stepsize can be increased
                                    askelpituus = askelpituus + 1.0_dp
                                    y = koe_y
                                    H1_y = koe_H1
                                    H2_y = koe_H2   
                                    fg1_y_all = fg1_koe_all
                                    fg2_y_all = fg2_koe_all
                                 ELSE
                                    stop_askelpituus = .TRUE.  
                                 END IF
                              END DO
                           END IF 
                           !-**--**- STEPSIZE DETERMINATION END -**--**-
                           END IF
                           
                           x_new = y                 ! the new iteration point is the current auxilary point
                           H1_new = H1_y             ! the value of H_1 at x_new
                           H2_new = H2_y             ! the value of H_2 at x_new
                           fg1_new_all = fg1_y_all   ! the values of DC components f1_i and g1_l at x_new
                           fg2_new_all = fg2_y_all   ! the values of DC components f2_i and g2_l at x_new
                           
                           stop_main_it = .TRUE.   ! the 'main iteration' can be stopped
                           reason_for_stop = 0     ! the reason for stopping is the new iteration point
                           
                           decrease  = give_decrease(B2)-H1_k
                           
                           IF ( real_decrease <= ( m_extra * decrease )) THEN
                                extra_dec = .TRUE.
                           ELSE
                                extra_dec = .FALSE.
                           END IF
                           
                           ! print everything
                           IF(iprint == 5) THEN
                              WRITE(*,*) 'A NEW ITERATION POINT is found: x_new= [', x_new, ']'  
                              WRITE(*,*) 'A NEW OBJECTIVE VALUE: f_new=', H1_new - H2_new                 
                              WRITE(*,*) 'The old point: [', x_k,']'              
                              WRITE(*,*) 'The old objective value:', H_k              
                              WRITE(*,*) ' '                  
                           END IF
                           
                           
                   !__________________________________________________________________________
                   !******************** STEP 4: END *****************************************         
                   
                       ELSE branch2 
                       !----------- SOME UPDATES IN VARIABLES BEGINS -------------------------  

                          ! print everything                       
                          IF (iprint == 5) THEN
                              WRITE(*,*) '-----------------------------------------------------------------------------------'
                          END IF                       
                           
                       
                   !__________________________________________________________________________     
                   !******************** STEP 5: BUNDLE UPDATE *******************************
                          
                          erotus = fg1_y_all(1)-fg2_y_all(1)-f_0_all(1)             !the biggest difference between the new auxiliary point y and x_0 when we compare function values of f_i and g_l 
                          DO i =2, number_of_obj
                              apu_erotus = fg1_y_all(i)-fg2_y_all(i)-f_0_all(i)
                              IF (apu_erotus > erotus) THEN
                                 erotus = apu_erotus
                              END IF
                          END DO
                          
                          update: IF (( erotus>0 ) .AND. & 
                                                      & (d_norm > eps)) THEN 
                          !----------- PARAMETER t REDUCTION --------------------------------
                              
                              t = t - c1 * ( t - t_min )   ! the parameter t is updated
                              i_dec = 0                    ! the numebr of bundle update steps made successively
                              i_t = -1                     ! null step is done in main iteration
                              
                          ! print everything
                          IF(iprint == 5) THEN
                              WRITE(*,*) ' '
                              WRITE(*,*) 'STEP 5a: PARAMETER t UPDATE'  
                              WRITE(*,*) ' '                  
                              WRITE(*,*) 'Updated parameter t=', t                
                              WRITE(*,*) '-----------------------------------------------------------------------------------'
                          END IF                  

                               
                          ELSE update
                          !----------- BUNDLE AND PARAMETER UPDATE --------------------------
                              
                              DO i = 1, number_of_obj                       ! the values of subgradients for DC components f1_i and f2_i at y
                                new_grad1 = subgradient_f1(y, problem1(i),factor(i),user_n)
                                new_grad2 = subgradient_f2(y, problem2(i),factor(i),user_n)
                                DO j = 1, user_n
                                   grad1_y_all(i,j) = new_grad1(j)
                                   grad2_y_all(i,j) = new_grad2(j)
                                END DO
                              END DO      
                              DO i = number_of_obj+1, number_of_func        ! the values of subgradients for DC components g1_l and g2_l at y
                                new_grad1 = subgradient_g1(y, problem1(i),factor(i),user_n)
                                new_grad2 = subgradient_g2(y, problem2(i),factor(i),user_n)
                                DO j = 1, user_n
                                   grad1_y_all(i,j) = new_grad1(j)
                                   grad2_y_all(i,j) = new_grad2(j)
                                END DO
                              END DO                              
                              subgrad1_counter = subgrad1_counter + 1          ! the subgradient counter is updated
                              subgrad2_counter = subgrad2_counter + 1          ! the subgradient counter is updated
                              
                              DO j = 1, number_of_func

                                new_grad1 = subgradient_AorB_i(grad1_y_all, grad2_y_all, j, user_n)    ! a subgradient of A_i or B_l at y

                                lin_err1 = AB_k_all(j) - AB_y_all(j) + &
                                           & DOT_PRODUCT(d_t,new_grad1)                        ! a new linearization error for A_i or B_l                       
                              
                               !->>>>>>>>>>>> BUNDLE B_1 UPDATE BEGINS <<<<<<<<<<<<<<<<<<<<<<<<<<- 
                                 bundle1: IF (is_full_b1(B1)) THEN   
                                 ! bundle B_1 is full and a new element is added to the bundle B_1                                
                                     CALL add_element_b1(B1, new_grad1, lin_err1, j)
                                     was_b1_full = .TRUE.
                                
                                 ELSE bundle1
                                 ! bundle B_1 is NOT full and a new element is added to the bundle B_1

                                     CALL add_element_b1(B1, new_grad1, lin_err1, j)
                                     was_b1_full = .FALSE.

                                 END IF bundle1                          
                               !->>>>>>>>>>>> BUNDLE B_1 UPDATE ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<-
                            
                              END DO
                              
                              ! print everything 
                              IF(iprint == 5) THEN
                                  WRITE(*,*) ' '
                                  WRITE(*,*) 'STEP 5b: BUNDLE UPDATE'   
                                  WRITE(*,*) ' '                  
                                  WRITE(*,*) 'delta2: ', delta2, ' delta1: ', delta1, ' glob solution index in B2: ', i 
                                  WRITE(*,*) ' ' 
                                  WRITE(*,*) 'A new bundle element is inserted into the bundle B_1. '   
                              END IF
                              
                                  
                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE BEGINS <<<<<<<<<<<-             

                              !bundle2: IF( delta2 >= 0 .AND. (give_max_size_b2(B2)>0) ) THEN   ! bundle B_2 is updated (happens only when delta2 >= 0 and bundle size is larger than 1) (SOME OTHER INSERTION RULES ALSO POSSIBLE!)
                              bundle2: IF( (give_max_size_b2(B2)>0) ) THEN   ! bundle B_2 is updated (happens only when delta2 >= 0 and bundle size is larger than 1) (SOME OTHER INSERTION RULES ALSO POSSIBLE!)

                                  new_grad2 = subgradient_H2(grad2_y_all, user_n)              ! a subgradient of f_2 at y  
                                  lin_err2 = H2_k - H2_y + DOT_PRODUCT(d_t,new_grad2)  ! a new linearization error for f_2  
                                
                                  CALL add_element_b2(B2, new_grad2, lin_err2)         ! a new element is inserted into the bundle B_2
                                  !In the algorithm we never overwrite the 'bundle element' yielding the previous global solution of the search direction problem.

                              
                                  ! print everything                                  
                                  IF(iprint == 5) THEN
                                      WRITE(*,*) 'A new bundle element is inserted into the bundle B_2'
                                      WRITE(*,*) ' '                  
                                      WRITE(*,*) ' subgrad2 = [', new_grad2,']'               
                                      WRITE(*,*) ' '                  
                                      WRITE(*,*) '--------------------&
                                        &------------------------------&
                                        &--------------------------------'
                                  END IF 

                                  
                   !__________________________________________________________________________     
                   !******************** STEP 5: END *****************************************
                   
                   !__________________________________________________________________________     
                   !******************** STEP 6: PARAMETER UPDATE ****************************
                   
                                  help = SQRT(DOT_PRODUCT(new_grad2,new_grad2)) ! the norm of the new subgradient of f_2 (subgradient is calculated 
                                                                                ! at the new auxilary point y)
                                  IF(help > max_norm) THEN 
                                      max_norm = help            ! the updated value of the maximum subgradient norm in the bundle B_2
                                      t_min = (0.5_dp * r_dec * eps) / ( norm1 + max_norm ) ! the updated value of the lower bound t_min
                                  END IF    
                              
                                  
                                  ! print everything
                                  IF(iprint == 5) THEN
                                      WRITE(*,*) '---------------------&
                                      &----------------------------------------------------&
                                      &----------'
                                      WRITE(*,*) ' '
                                      WRITE(*,*) 'STEP 6: PARAMETER UPDATE'           
                                      WRITE(*,*) ' '                        
                                  END IF 
  
                                  
                   !__________________________________________________________________________     
                   !******************** STEP 6: END *****************************************
                   
                            END IF bundle2
                            
 
                              ! --> the update of the proximity parameter BEGINS <--  
                                  i_t = -1                         ! Null step is done in main iteration
                                  
                                  IF (i_dec <= 50) THEN
                                      t = t - c2 * ( t - t_min )    ! the parameter t is decreased
                                  ELSE IF (i_dec >= 51 ) THEN
                                      t = t - c3 * ( t - t_min ) ! the parameter t is decreased     
                                      !WRITE(*,*) 'taalla'                                    
                                  END IF 
                                  
                                  i_dec = i_dec + 1                ! one more bundle update step executed

                              ! --> the update of the proximity parameter ENDS <--  
                              

                            ! print everything
                            IF(iprint == 5) THEN
                                WRITE(*,*) '-----------------------------------------------------------------------------------'    
                            END IF                          

                            !->>>>>>>>>>>> BUNDLE B_2 AND PARAMETER UPDATE ENDS <<<<<<<<<<<<<-
                                                              
                          END IF update        
                        !------ PARAMETER t REDUCTION & BUNDLE AND PARAMETER UPDATE ENDS -----
                        
                          iter_counter = iter_counter + 1  ! update of the iteration counter                        
                              
                       END IF branch2   
                     !--------- NEW ITERATION POINT & SOME UPDATES IN VARIABLES ENDS --------- 
                     
               !->->->->->->->->->->->->->-> BRANCH 2 END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-   
                   END IF branches
               !->->->->->->-> EXECUTION OF BRANCH ENDS  (2 BRANCHES) <-<-<-<-<-<-<-<-<-<-<-<-
              END DO 
              
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <>  REPEATABLE PART ENDS  <>  <>  <>  <>  <>  <>  <>  <> 
           !----------------------------------------------------------------------------------           
              
              IF ( (.NOT. stop_main_it ) ) THEN 
                  reason_for_stop = 3            ! the maximum number of rounds have been executed 
                  iter_counter = iter_counter -1 
              END IF        
              
           END SUBROUTINE main_iteration
        !.......................................................................................           
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>   
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
        


  
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |
        !  |                      Guaranteeing Clarke stationarity                          |  |
        !  |                                                                                |  |
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************        
           
            SUBROUTINE guaranteeing_clarke( x_k, H1_k, H2_k, fg1_k_all, fg2_k_all, &
                            & reason_for_stop, x_new , H1_new, H2_new, fg1_new_all, fg2_new_all,&
                            & d_cause, crit_tol, m, step_tol, iprint, mrounds, &
                            & iter_counter, f_counter, subgrad_counter, factor, user_n, & 
                            & number_of_obj, number_of_const, number_of_func)
                                                
            !
            ! Executes the 'Clarke stationary' algorithm. It is needed to guarantee Clarke stationarity at the current iteration point. 
            ! If the current point is not Clarke stationary then the algorithm generates a descent direction.
            !
            ! INPUT: * 'x_k'               : The current iteratioin point
            !        * 'H1_k' and 'H2_k'   : The values of DC components H_1 and H_2 at the current iteration point
            !  
            !        * 'fg1_k_all' and 'fg2_k_all'   : The values of DC components f1_i, f2_i, g1_l and g2_l at the current iteration point
            !
            !        * 'factor'            : The values of scaling factors for f_i and g_l
            !        * 'd_cause'           : The search direction which caused the execution of the clarke stationary algorithm
            !        * 'crit_tol'          : The stopping tolerance
            !        * 'm'                 : The descent parameter
            !        * 'step_tol'          : The step-length tolerance  
            !
            !        * 'mrounds'           : The maximum number of possible rounds in the 'Clarke stationary' algorithm
            !        * 'iprint'            : Specifies print
            
            !        * 'user_n'            : The dimension of the problem
            !
            !
            ! OUTPUT: * 'x_new'                 : Tthe new iteration point (if it is found, i.e., 'reason_for_stop' = 0)
            !         * 'H1_new' and 'H2_new'   : The values of DC components H_1 and H_2 at the new iteration point
            !         * 'fg1_new_all' and 'fg2_new_all': the new value of DC components f1_i, f2_i, g1_l and g2_l (if new iteration point is found, i.e., 'reason_for_stop' = 0)
            !
            !         * 'reason_for_stop'       : Indicates the cause of the termination in the 'Clarke stationarity' algorithm  
            !
            !         * 'iter_counter'          : the number of rounds needed in 'Clarke stationary' algorithm
            !         * 'f_counter'             : the number of function values evaluated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !         * 'subgrad_counter'       : the number of subgradients calculated for f in 'Clarke stationary' algorithm (same for f_1 and f_2)
            !
            ! NOTICE: * The dimensions of vectors 'x_k', 'x_new' and 'd_cause' has to be same. (i.e. the dimension of 'x_k', 'x_new' and 'd_cause' has to be
            !           'user_n' when SUBROUTINE guaranteeing_clarke() is used in SUBROUTINE main_iteration().
            !         * The dimension of vectors 'fg1_k_all', 'fg2_k_all', 'fg1_new_all' and 'fg2_new_all' has to be 'number_of_func'
            !
            !***********************************************************************************
               IMPLICIT NONE
            !**************************** NEEDED FROM THE USER ********************************* 
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x_k      ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: d_cause  ! the search direction
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_new    ! the new iteration point if 'reason_for_stop'=0 in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg1_k_all    ! the values of DC components f1_i and g1_l at x_k
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: fg2_k_all    ! the values of DC components f2_i and g2_l at x_k
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg1_new_all ! the values of DC components f1_i and g1_l at x_new
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: fg2_new_all ! the values of DC components f2_i and g2_l at x_new

               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: factor       ! the values of scaling factors for f_i and g_l           
               
               REAL(KIND=dp), INTENT(IN)  :: H1_k, H2_k      ! the value of H_1 and H_2 at x_k
               REAL(KIND=dp), INTENT(OUT) :: H1_new, H2_new  ! the value of H_1 and H_2 at x_new if 'reason_for_stop=0' in the 'Clarke stationary' algorithm
               
               REAL(KIND=dp), INTENT(IN) :: crit_tol      ! 'crit_tol'=stopping tolerance (needed in stopping condition)
               REAL(KIND=dp), INTENT(IN) :: m             ! a descent parameter
               REAL(KIND=dp), INTENT(IN) :: step_tol      ! a step-length tolerance
               
               INTEGER, INTENT(IN) :: user_n              ! the dimension of the problem
               
               INTEGER, INTENT(IN) :: number_of_obj      ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_const    ! The dimension of the problem
               INTEGER, INTENT(IN) :: number_of_func     ! The dimension of the problem            
               
               INTEGER :: mrounds                         ! the maximum number of possible rounds in the 'Clarke statinonary' algorithm
 
               INTEGER, INTENT(OUT) :: reason_for_stop    ! 0 - a new iteration point found
                                                          ! 1 - stopping condition satisfied (approximate Clarke stationarity)
                                                          ! 2 - approximate stopping condition satisfied (i.e. the step-length beta* < beta_0)
                                                          ! 5 - the maximum number of rounds executed in 'Clarke stationary' algorithm
                                                          
               INTEGER, INTENT(OUT) :: iter_counter       ! the number of rounds used in 'Clarke stationary' algorithm
               INTEGER, INTENT(OUT) :: f_counter          ! the number of function values evaluated for H in 'Clarke stationary' algorithm (same for f_1 and f_2)
               INTEGER, INTENT(OUT) :: subgrad_counter    ! the number of subgradients calculated for H in 'Clarke stationary' algorithm 
               
               INTEGER, INTENT(IN) :: iprint    ! the variable that specifies print:
                                                !   iprint = 0 : print is suppressed
                                                !   iprint = 1 : basic print of final result (nothing is printed in main iteration)
                                                !   iprint = -1: basic print of final result without the solution vector (nothing is printed in main iteration)                                             
                                                !   iprint = 2 : extended print of final result (nothing is printed in main iteration)
                                                !   iprint = -2: extended print of final result without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 3 : basic print of intermediate results and extended print of final results (nothing is printed in main iteration)
                                                !   iprint = -3: basic print of intermediate results and extended print of final results without the solution vector (nothing is printed in main iteration)
                                                !   iprint = 4 : extended print of intermediate and final results (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = -4: extended print of intermediate and final results without the solution vector (something is printed about the 'main iteration' (this print is done in the bundle algorithm))               
                                                !   iprint = 5 : everything is printed step by step (this is the only option which causes print in the 'main iteration') 
               

           !***************************** LOCAL VARIABLES ************************************
           
               TYPE(kimppu1) :: B                                ! the bundle B of H containing subgradients from \partial H(x)
               
               REAL(KIND=dp), DIMENSION(user_n) :: d                      ! the search direction
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad               ! the subgradient of H at x
               REAL(KIND=dp), DIMENSION(user_n) :: new_grad1, new_grad2   ! the subgradients of H_1 and H_2 at x_d
               REAL(KIND=dp), DIMENSION(user_n) :: u_k                    ! the vector yielding minimum norm for quadratic norm minimization problem
               REAL(KIND=dp), DIMENSION(user_n) :: y                      ! the new auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: test_y                 ! the test auxiliary point            
               REAL(KIND=dp), DIMENSION(user_n) :: lower_y                ! the 'hepl' auxiliary point            
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_y_all    ! the values of DC components f1_i and g1_l at y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_y_all    ! the values of DC components f2_i and g2_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_y_all     ! the values of components A_i and B_l at y
               
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_test_all     ! the values of DC components f1_i and g1_l at test_y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_test_all     ! the values of DC components f2_i and g2_l at test_y   
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_test_all      ! the values of components A_i and B_l at test_y              
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg1_lower_all    ! the values of DC components f1_i and g1_l at lower_y
               REAL(KIND=dp), DIMENSION(number_of_func) :: fg2_lower_all    ! the values of DC components f2_i and g2_l at lower_y  
               REAL(KIND=dp), DIMENSION(number_of_func) :: AB_lower_all     ! the values of components A_i and B_l at lower_y              

               REAL(KIND=dp), DIMENSION(number_of_func,user_n) :: fg1_grad_all    ! the subgradients of DC components f1_i and g1_l
               REAL(KIND=dp), DIMENSION(number_of_func,user_n) :: fg2_grad_all    ! the subgradients of DC components f2_i and g2_l
               
               REAL(KIND=dp) :: H_k              ! the value of the objective function H at the current iteration point 
               REAL(KIND=dp) :: obj_u            ! the optimal value of the quadratin norm minimization problem 
               REAL(KIND=dp) :: norm_u           ! the norm for the vector u_k  
               REAL(KIND=dp) :: eps              ! stepsize used in subgradient calculation 
               REAL(KIND=dp) :: div_eps          ! 1.0_dp / eps     
               REAL(KIND=dp) :: real_decrease    ! the real decrease in the objective function f
               REAL(KIND=dp) :: H1_y, H2_y       ! the function values of H_1 and H_2 at y
               REAL(KIND=dp) :: test_H1, test_H2 ! the function values of H_1 and H_2 at test_y
               REAL(KIND=dp) :: lower_H1, lower_H2 ! the function values of H_1 and H_2 at lower_y
               REAL(KIND=dp) :: direc_der        ! directional derivative of H 
               REAL(KIND=dp) :: stepsize         ! stepsize 
               REAL(KIND=dp) :: upper, lower     ! bounds for stepsize 
               REAL(KIND=dp) :: descent_app      ! approximation of descent in the objective function value
                              
               INTEGER :: size_b                 ! the biggest possible bundle size for B
               INTEGER :: N                      ! the current bundle size for B with the current element
               
               INTEGER :: i, j, l 
               
               LOGICAL :: ylos               ! .TRUE. if 'Clarke stationary' algorithm can be stopped
               LOGICAL :: stop_alg               ! .TRUE. if 'Clarke stationary' algorithm can be stopped
               LOGICAL :: stop_step_det          ! .TRUE. if the stepsize determination can be stopped
               

           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM STARTS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             
              IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '              CLARKE STATIONARY ALGORITHM STARTS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
              END IF             
             
           !_______________________________________________________________________________
           !************************ STEP 0: INITIALIZATION *******************************              
                        
            ! the size of the bundle B
               IF (user_size < 2) THEN
                   size_b = (user_n + 5)*2
               ELSE
                   size_b = user_size
               END IF
               
               eps = (10.0_dp)**(-5)     ! Stepsize used when subgradient is calculated
               div_eps = 1.0_dp / eps
           
            ! Initialization of iteration counters
               iter_counter = 1                       ! The first iteration round is executed
               f_counter = 0                          ! The number of function values evaluated for f
               subgrad_counter = 0                    ! The numeber of subgradients evaluated for f
               l = 0
               
               H_k = H1_k - H2_k
            
               d =  -d_cause / SQRT(DOT_PRODUCT(d_cause, d_cause))             
               
            ! The bundles B is initialized
               CALL init_bundle_b1(B, user_size, user_n) 
            ! Nothing is stored into bundle B      
               N = 0                        
               
               stop_alg = .FALSE.
               
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 0: INITIALIATION '
                   WRITE(*,*)''
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 0: END **********************************************      
           
               
               y = x_k + (eps) * d 
               
               DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at y
                   fg1_y_all(i) = f1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = f2(y,problem2(i),factor(i),user_n) 
               END DO  
               
               DO i = number_of_obj+1, number_of_func   ! the values of DC components g1_l and g2_l at y
                   fg1_y_all(i) = g1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = g2(y,problem2(i),factor(i),user_n) 
               END DO              
               f_counter = f_counter + 1                ! one more objective function value calculated for DC components          
                               
               DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                   AB_y_all(i) = AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
               END DO     

  
           !----------------------------------------------------------------------------------   
           ! <>  <>  <>  <>  <>  <>  <> REPEATABLE PART BEGINS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------         
           DO WHILE ( (.NOT. stop_alg) .AND. (iter_counter <= mrounds))                  
           !_______________________________________________________________________________
           !************************ STEP 1: NEW SUBGRADIENT ******************************          
               
               IF (iter_counter > 100) THEN 
                   eps = (10.0_dp)**(-6)
               END IF
              
                   
               DO i = 1, number_of_obj                       ! the values of subgradients for DC components f1_i and f2_i at x_0
                  new_grad1 = subgradient_f1(y, problem1(i),factor(i),user_n)
                  new_grad2 = subgradient_f2(y, problem2(i),factor(i),user_n)
                  DO j = 1, user_n
                     fg1_grad_all(i,j) = new_grad1(j)
                     fg2_grad_all(i,j) = new_grad2(j)
                  END DO
               END DO   
               
               DO i = number_of_obj+1,number_of_func          ! the values of subgradients for DC components g1_l and g2_l at x_0
                  new_grad1 = subgradient_g1(y, problem1(i),factor(i),user_n)
                  new_grad2 = subgradient_g2(y, problem2(i),factor(i),user_n)
                  DO j = 1, user_n
                     fg1_grad_all(i,j) = new_grad1(j)
                     fg2_grad_all(i,j) = new_grad2(j)
                  END DO
               END DO                
               subgrad_counter = subgrad_counter+1            ! one more subgradient calculated for DC components 

               new_grad1 = subgradient_H1(fg1_grad_all, fg2_grad_all, AB_y_all, user_n) 
               new_grad2 = subgradient_H2(fg2_grad_all, user_n)
               
               ! new subgradient for H at x 
               new_grad = new_grad1 - new_grad2
               
               IF (N == 0) THEN 
                  CALL add_first_element_b1(B, new_grad,0)  
                  N = N + 1
               ELSE
                  CALL add_element_b1(B, new_grad, 0.0_dp,0)
               END IF 
           
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 1: NEW SUBGRADIENT '
                   WRITE(*,*)''
                   WRITE(*,*)' Iteration round: ', iter_counter
                   WRITE(*,*)' New subgardient = [', new_grad, ']' 

               END IF
               
           !_______________________________________________________________________________    
           !******************** STEP 1: END **********************************************            
           
           !_______________________________________________________________________________
           !************************ STEP 2: CALRKE STATIONARITY **************************            
           
               N = give_size_b1(B) + 1      ! the current bundle size
               
               CALL quadratic_solver(x_k, user_n, N, B, 1.0_dp, u_k, obj_u) 
               ! The solution u_k is has already negative sign. Thus d = u_k/norm_u
                               
               norm_u = SQRT( DOT_PRODUCT(u_k,u_k))

               
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 2: CLARKE STATIONARITY '
                   WRITE(*,*)''
                   WRITE(*,*)' The value of minimum norm = ', norm_u
                   WRITE(*,*)' The solution to quadratic norm minimization problem  = [',u_k, ']'               
               END IF              
               
               IF (norm_u < crit_tol) THEN
                  
                  reason_for_stop = 1 
                  stop_alg = .TRUE.
                  
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)''
                   WRITE(*,*)' CLARKE STATIONARITY ACHIEVED'
                   WRITE(*,*)''                
               END IF                     
                  
           !_______________________________________________________________________________    
           !******************** STEP 2: END **********************************************     
               
               ELSE            
           !_______________________________________________________________________________
           !************************ STEP 3: SERACH DIRECTION *****************************            
               
               d =  u_k / norm_u
               y = x_k + eps * d           
               
               DO i = 1, number_of_obj                  ! the values of DC components f1_i anf f2_i at y
                   fg1_y_all(i) = f1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = f2(y,problem2(i),factor(i),user_n) 
               END DO      
               
               DO i = number_of_obj+1,number_of_func    ! the values of DC components g1_l and g2_l at y
                   fg1_y_all(i) = g1(y,problem1(i),factor(i),user_n) 
                   fg2_y_all(i) = g2(y,problem2(i),factor(i),user_n) 
               END DO                  
               f_counter = f_counter + 1                ! one more objective function value calculated for DC components          
                               
               DO i = 1, number_of_func                 ! the values of components A_i and B_l for H_1 at y
                   AB_y_all(i) = AorB_i(fg1_y_all, fg2_y_all, fg1_k_all, fg2_k_all ,i)
               END DO   
               
               H1_y = H1(AB_y_all) 
               H2_y = H2(fg2_y_all) 
               
               real_decrease = H1_y - H2_y - H_k

               direc_der = div_eps * real_decrease

            
               IF(iprint==5) THEN  !prints everything
                   WRITE(*,*)'-------------------------------------------------------------------------------'             
                   WRITE(*,*)''
                   WRITE(*,*)' STEP 3: SEARCH DIRECTION '
                   WRITE(*,*)''               
                   WRITE(*,*)' The search direction d = [', d, ']' 
                   WRITE(*,*)''                
               END IF           
               
               IF (direc_der <= (-m * norm_u)) THEN
                     reason_for_stop = 10                 ! we will execute step-length determination at Step 4 
                     stop_alg = .TRUE.
                     
                     IF(iprint==5) THEN  !prints everything
                         WRITE(*,*)''
                         WRITE(*,*)' SUFFICIENT DESCENT IS OBTAINED.'
                         WRITE(*,*)''                  
                     END IF                      
                     
               END IF
               
               iter_counter = iter_counter +1 
               
               END IF 
               
           !_______________________________________________________________________________    
           !******************** STEP 3: END ********************************************** 
            END DO 
           !----------------------------------------------------------------------------------   
           !  <>  <>  <>  <>  <>  <>  <> REPEATABLE PART ENDS <>  <>  <>  <>  <>  <>  <>  <>
           !----------------------------------------------------------------------------------           
            
           !_______________________________________________________________________________
           !************************ STEP 4: STEP-LENGTH **********************************            
               IF ( (reason_for_stop > 1) .AND. stop_alg ) THEN 
           
               !-**--**- STEPSIZE DETERMINATION START -**--**-         
                  upper = 1.0_dp
                  lower = user_step_tol
                  stepsize = user_step_tol
                  stop_step_det = .FALSE.
                  ylos = .FALSE.
                
                  lower_y = x_k + stepsize * d
                  DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at lower_y
                      fg1_lower_all(i) = f1(lower_y,problem1(i),factor(i),user_n) 
                      fg2_lower_all(i) = f2(lower_y,problem2(i),factor(i),user_n) 
                  END DO  
                  DO i = number_of_obj+1,number_of_func    ! the values of DC components g1_l and g2_l at lower_y
                      fg1_lower_all(i) = g1(lower_y,problem1(i),factor(i),user_n) 
                      fg2_lower_all(i) = g2(lower_y,problem2(i),factor(i),user_n) 
                  END DO                      
                  f_counter = f_counter+1             
                               
                  DO i = 1, number_of_func              ! the values of components A_i and B_l for H_1 at lower_y
                      AB_lower_all(i) = AorB_i(fg1_lower_all, fg2_lower_all,  & 
                                              & fg1_k_all, fg2_k_all ,i)
                  END DO                          
                      
                 lower_H1 = H1(AB_lower_all)     ! the value of H_1 at test_y
                 lower_H2 = H2(fg2_lower_all)    ! the value of H_2 at test_y 

                 descent_app = -m * stepsize * norm_u       
                 
                 IF ((lower_H1 - lower_H2 - H_k) > descent_app ) THEN 
                    stop_step_det = .TRUE.
                    stepsize = -1.0_dp                
                 END IF
                 l = 1  
     
                  
                 DO WHILE((.NOT. stop_step_det))
                  
                      IF (l == 1) THEN
                        stepsize = 100.0_dp * stepsize
                      ELSE IF (l == 2) THEN 
                        IF (upper < 0.02_dp) THEN 
                          stepsize = 0.1_dp * stepsize
                          !stepsize = 0.5_dp*(upper +lower)
                        ELSE
                          stepsize = 10.0_dp*stepsize
                        END IF
                      ELSE IF (l > 2) THEN
                        stepsize = 0.5_dp*(upper +lower)
                      END IF 
                      
                      test_y = x_k + stepsize * d
                      DO i = 1, number_of_obj                  ! the values of DC components f1_i and f2_i at test_y
                         fg1_test_all(i) = f1(test_y,problem1(i),factor(i),user_n) 
                         fg2_test_all(i) = f2(test_y,problem2(i),factor(i),user_n) 
                      END DO  
                      DO i = number_of_obj+1,number_of_func    ! the values of DC components g1_l and g2_l at test_y
                         fg1_test_all(i) = g1(test_y,problem1(i),factor(i),user_n) 
                         fg2_test_all(i) = g2(test_y,problem2(i),factor(i),user_n) 
                      END DO                      
                      f_counter = f_counter + 1             ! one more objective function value calculated for a DC component          
                               
                      DO i = 1, number_of_func              ! the values of components A_i and B_l for H_1 at test_y
                        AB_test_all(i) = AorB_i(fg1_test_all, fg2_test_all,  & 
                                              & fg1_k_all, fg2_k_all ,i)
                      END DO                          
                      
                      test_H1 = H1(AB_test_all)     ! the value of H_1 at test_y
                      test_H2 = H2(fg2_test_all)    ! the value of H_2 at test_y 

                      IF (l < 7) THEN
                        descent_app = -m * stepsize * norm_u  
                        IF ( (test_H1 - test_H2 - H_k) < descent_app) THEN
                      
                         fg1_lower_all = fg1_test_all
                         fg2_lower_all = fg2_test_all
                         lower_H1 = test_H1
                         lower_H2 = test_H2
                         lower_y = test_y
                         
                         l = l + 1 
                         lower = stepsize
                        ELSE 
                         upper = stepsize
                         l = l+1
                        END IF
                      ELSE  

                         stop_step_det = .TRUE. 
                         y = lower_y
                         H1_y = lower_H1
                         H2_y = lower_H2
                         fg1_y_all = fg1_lower_all
                         fg2_y_all = fg2_lower_all

                      END IF

                      
                  END DO
               !-**--**- STEPSIZE DETERMINATION END -**--**-
               
                  IF(iprint==5) THEN  !prints everything
                     WRITE(*,*)'-------------------------------------------------------------------------------'               
                     WRITE(*,*)''
                     WRITE(*,*)' STEP 4: STEP-LENGTH '
                     WRITE(*,*)''                 
                     WRITE(*,*)' The step-lenght t =', stepsize
                  END IF                   
               
       
                  IF (stepsize >= step_tol) THEN    ! the step-length is big enough
                 
                     x_new = y                 ! the new iteration point is the auxilary point y
                     H1_new = H1_y             ! the value of H_1 at x_new
                     H2_new = H2_y             ! the value of H_2 at x_new
                     fg1_new_all = fg1_y_all   ! the values of f1_i and g1_l at x_new
                     fg2_new_all = fg2_y_all   ! the values of f2_i and g2_l at x_new
                     
                     reason_for_stop = 0     ! the reason for stopping is the new iteration point   
                     
                     IF(iprint==5) THEN  !prints everything
                        WRITE(*,*)' New itaration point is found and x_new = [', x_new, ']' 
                        WRITE(*,*)' New objective function value f_new = ', H1_new -H2_new                 
                        WRITE(*,*)''                   
                     END IF                      
                     
                  ELSE
                     reason_for_stop = 2
                     
                     IF(iprint==5) THEN  !prints everything
                        WRITE(*,*)' APPROXIMATE STOPPING CONDITION HOLD at the current iteration point.'               
                        WRITE(*,*)''                   
                     END IF
                     
                  END IF 
               
                  iter_counter = iter_counter -1 
                  
               END IF
           !_______________________________________________________________________________    
           !******************** STEP 4: END **********************************************             
           
              IF ( (.NOT. stop_alg ) ) THEN 
                  reason_for_stop = 5            ! the maximum number of rounds have been executed 
                  
                  IF(iprint==5) THEN  !prints everything
                     WRITE(*,*)'--------------------------------------&
                     &-----------------------------------------'                       
                     WRITE(*,*)''                                     
                     WRITE(*,*)' THE MAXIMUM NUMBER OF ROUNDS HAS BEEN EXECUTED.'              
                     WRITE(*,*)''                  
                  END IF                  
              END IF           
           
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
           !   <>  <>  <>  <>  <>  <>  <>  ALGORITHM ENDS  <>  <>  <>  <>  <>  <>  <>  <>
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         

              IF(iprint==5) THEN  !prints everything
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) '                    CLARKE STATIONARY ALGORITHM ENDS'
                    WRITE(*,*) '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
                    WRITE(*,*) ' '
              END IF               
           
           END SUBROUTINE guaranteeing_clarke
        !.......................................................................................   
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END **  
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|    


        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                          QUADRATIC NORM SOLVER                                 |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE quadratic_solver(x, NF, NA, B, t, d, obj) 
               ! Solves the quadratic norm problem in the 'Clarke stationary' algorithm.
               ! Calls PLQDF1 by Ladislav Luksan
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B                   ! the bundle B
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! the current iteration point
               REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: d       ! the vector giving minimum norm
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! the proximity parameter t
               REAL(KIND=dp), INTENT(OUT) :: obj                   ! the optimal value of the quadratic norm problem
               
               INTEGER, INTENT(IN)  :: NF    ! the number of variables
               INTEGER, INTENT(IN)  :: NA    ! the bundle size of B is 'give_size_b1(B) + 1' 
                
            !*************************** LOCAL VARIABLES ************************************
            
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  

               REAL(KIND=dp) :: u       ! 'help' variable
               
               INTEGER :: j, k, l
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               !INTEGER :: IDECF=10                     ! IDECF=10 diagonal matrix
               !INTEGER :: MFP=2                        ! MFP=2 optimum feasible point
               INTEGER :: KBC=0, KBF=0                  ! KBC=0 no linear constraints; KBF=0 no simple bounds
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b           ! subgradient matrix of B


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               grad_m_b = grad_matrix(B)       ! subgradient matrix of B 
               
               ! NO linearization errors
               AF = 0.0_dp      
  
               ! matrix containing subgradients   
               DO j = 1, NA
                   k = (j-1)*NF
                   DO l = 1, NF
                      AG(k+l) = grad_m_b(k+l) 
                   END DO
               END DO   
                        
               !Calls PLQDF1 by Ladislav Luksan
               CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                           & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,2,KBF, &
                           & KBC,10,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                       
                           & UMAX,GMAX,N,ITERQ)

               ! the vector giving minimum norm               
               DO j = 1, NF
                   d(j) = S(j) 
               END DO
              
              ! the value of minimum norm
               obj = 0.5_dp * DOT_PRODUCT(d, d)         
    
                   
               
           END SUBROUTINE quadratic_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|        
        
        


        
          
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
        ! START ** START ** START ** START ** START ** START ** START ** START ** START ** START 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|         
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>
        !***************************************************************************************
        !  ----------------------------------------------------------------------------------  |
        !  |                                                                                |  |        
        !  |                           SUBPROBLEM SOLVER                                    |  |
        !  |                                                                                |  |        
        !  ----------------------------------------------------------------------------------  |
        !***************************************************************************************           

            SUBROUTINE subproblem_solver(x, AB_k_all, NF, NA, B1, B2, t, subprob_counter, number_of_func) 
               ! Solves the subproblems in the original search direction problem.
               ! Calls PLQDF1 by Ladislav Luksan
               IMPLICIT NONE               
            !**************************** NEEDED FROM USER ***********************************             
 
               TYPE(kimppu1), INTENT(IN)    :: B1                  ! bundle B_1
               TYPE(kimppu2), INTENT(INOUT) :: B2                  ! bundle B_2
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x        ! current iteration point
               
               REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: AB_k_all  ! the values of components A_i and B_l at the current iteration point
               
               REAL(KIND=dp), INTENT(IN) :: t                      ! proximity parameter t
               
               INTEGER, INTENT(IN)  :: NF    ! number of variables
               INTEGER, INTENT(IN)  :: NA    ! bundle size of B1  
               
               INTEGER, INTENT(IN)  :: number_of_func    ! bundle size of B1         
               
               INTEGER, INTENT(OUT) :: subprob_counter    ! number of subproblems solved    
                       
              
            !*************************** LOCAL VARIABLES ************************************
            
               
               ! .. Parameters used in PLQDF1 ..
               REAL(KIND=dp), PARAMETER :: ETA0  = 1.0E-15_dp  
               REAL(KIND=dp), PARAMETER :: ETA2  = 1.0E-12_dp
               REAL(KIND=dp), PARAMETER :: ETA9  = 1.0E+60_dp
               REAL(KIND=dp), PARAMETER :: EPS7  = 1.0E-14_dp
               REAL(KIND=dp), PARAMETER :: EPS9  = 1.0E-12_dp 
               
               ! .. Scalar Arguments in PLQDF1 ..               
               REAL(KIND=dp) :: GMAX   ! output of PLQDF1: maximum absolute value of a partial derivative
               REAL(KIND=dp) :: UMAX   ! output of PLQDF1: maximum absolute value of a negative lagrange multiplier
               REAL(KIND=dp) :: XNORM  ! output of PLQDF1: value of linearized minimax function = delta_1 + delta_2  
               
               REAL(KIND=dp) :: alpha_b2     ! a linearization error of B_2
               REAL(KIND=dp) :: u, obj       ! 'help' variable
               
               REAL(KIND=dp) :: a
               INTEGER :: i, j, k, l
               
               INTEGER, PARAMETER :: NC = 0             ! number of constraints is zero
               !INTEGER :: IDECF=10                     ! IDECF=10 diagonal matrix;
               !INTEGER :: MFP=2                        ! MFP=2 optimum feasible point
               INTEGER :: KBC=0, KBF=0                  ! KBC=0 no linear constraints; KBF=0 no simple bounds
               INTEGER :: ITERQ, N                      ! output values of PLQDF1: ITERQ=type of feasible point 
                                                        ! N=dimension of manifold defined by active constraints 
               ! .. Array Arguments in PLQDF1..
               INTEGER, DIMENSION(NF) :: IX     ! vector containing types of bounds
               INTEGER, DIMENSION(NA) :: IA     ! vector containing types of deviations 
               INTEGER, DIMENSION(NC) :: IC     ! vector containing types of constraints. NOT significant because NC=0
               INTEGER, DIMENSION(NF+1) :: IAA  ! Output of PLQDF1: vector containing indicies of active functions
               
               REAL(KIND=dp), DIMENSION(NF) :: XL, XU          ! lower and upper bounds for x. NOT significant because variables are unbounded
               REAL(KIND=dp), DIMENSION(NF) :: H               ! diagonal matrix (1/t)*I
               REAL(KIND=dp), DIMENSION(NC) :: CF, CL, CU      ! NOT significant since NC=0 (are related to constraints)
               REAL(KIND=dp), DIMENSION(NA) ::  AF             ! vector of bundle function values (-alpha)
               REAL(KIND=dp), DIMENSION(NA) :: AFD             ! Output of PLQDF1: vector containing increments of the approximated functions
               REAL(KIND=dp), DIMENSION(NF+1) :: AZ            ! Output of PLQDF1: vector of Lagrange multipliers
               REAL(KIND=dp), DIMENSION(NF+1) :: S             ! Output of PLQDF1: direction vector
               REAL(KIND=dp), DIMENSION(NF) :: direction       ! Actual direction vector used (notice dimension)
               REAL(KIND=dp), DIMENSION(NF+1) :: G             ! Output of PLQDF1: gradient of the Lagrangian function
               REAL(KIND=dp), DIMENSION(NF*NA) :: AG           ! matrix whose columns are bundle subgradients 
               REAL(KIND=dp), DIMENSION((NF+1)*(NF+2)/2) :: AR ! Output of PLQDF1: triangular decomposition of kernel of the orthogonal projection
               REAL(KIND=dp), DIMENSION(NF*NC) :: CG           ! NOT significant since NC=0. matrix whose columns are normals of the linear constraints            
               
               ! .. Some other varibles ..    
               REAL(KIND=dp), DIMENSION(NF*NA) :: grad_m_b1           ! subgradient matrix of B_1
               REAL(KIND=dp), DIMENSION(NA) :: alpha_m_b1             ! linearization error matrix of B_1
               REAL(KIND=dp), DIMENSION(NF) :: grad_b2                ! a subgradient of B_2


           !************************** SUBPROBLEM SOLVER STARTS *********************************   
           
           !****************************** INITIALIZATIONS **************************************
           
           
               IX = 0             ! types of bounds: 0 - unbounded variables
               IA = 2             ! types of deviations
               
               u = 1.0_dp / t
               H = u              ! diagonal matrix
               
               grad_m_b1 = grad_matrix(B1)                     ! subgradient matrix of B_1
               alpha_m_b1 = lin_err_and_f1_matrix(B1,AB_k_all,number_of_func) ! linearization error matrix of B_1 with component values A_i(x_h) and B_l(x_h)
               
               subprob_counter = give_size_b2(B2) 
              

               !$OMP PARALLEL DO PRIVATE(grad_b2,alpha_b2,direction,a,obj) & 
               !$OMP FIRSTPRIVATE(NA,X,IX,XL,XU,AFD,IA,IAA) &
               !$OMP FIRSTPRIVATE(AR,AZ,CF,IC,CL,CU,CG,G,H,S,KBF) &
               !$OMP FIRSTPRIVATE(KBC,XNORM,UMAX,GMAX,N,ITERQ) &
               !$OMP PRIVATE(i, AG, AF ) &
               !$OMP SHARED(grad_m_b1,alpha_m_b1,B2)               
                       
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-
                    subproblems1: DO i = 0, subprob_counter   ! each subproblem is looked through
                    

                        grad_b2  = give_subgrad_b2(B2, i)  ! subgradient of B_2 in the subproblem i
                        alpha_b2 = give_linerr_b2(B2, i)   ! linearization error of B_2 in the subproblem i

                        
                        !$OMP CRITICAL  
                        DO j = 1, NA
                           k = (j-1)*NF
                           DO l = 1, NF
                               AG(k+l) = grad_m_b1(k+l) - grad_b2(l)
                           END DO
                           AF(j) = - alpha_m_b1(j) + alpha_b2 
                        END DO
                        !$OMP END CRITICAL  
                        
                        !Calls PLQDF1 by Ladislav Luksan
                        CALL PLQDF1(NF,NA,NC,X,IX,XL,XU,AF,AFD,IA,IAA, &
                              & AG,AR,AZ,CF,IC,CL,CU,CG,G,H,S,2,KBF, &
                              & KBC,10,ETA0,ETA2,ETA9,EPS7,EPS9,XNORM, &                        
                              & UMAX,GMAX,N,ITERQ)

                              
                        DO j = 1, NF
                           direction(j) = S(j) 
                        END DO
              
                        a = DOT_PRODUCT(direction, direction)           

                        obj =   XNORM + (a * u) / 2
                        !$OMP CRITICAL
                        CALL add_solution(B2, i , direction, XNORM, obj )   
                        !$OMP END CRITICAL
                    
                    END DO subproblems1
               !->->->->->->->->->->->->->-> EACH SUBPROBLEM SOLVED END <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-                      
               !$OPM END PARALLEL DO

               subprob_counter = subprob_counter + 1 
               
               
           END SUBROUTINE subproblem_solver
        !.......................................................................................
        ! <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>        
        ! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  
        !| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |        
        !** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** END ** 
        !|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

        






      END MODULE mdbdc