PROGRAM ADVECT_FIBERS
  !! Main program
  !! Given an inital configuration, the fiber system is advected in time
  USE quadrature
  USE MatrixAssembly
  USE write_Matlab
  USE coeffs_and_vels
  USE Int_analytical
  IMPLICIT NONE

  !****************************************************************************!
  !                                                                            !
  ! Variable Declarations                                                      !
  !                                                                            !
  !****************************************************************************!

  !me:  The number of individual fibers
  INTEGER M
  !me:  @todo The number of terms sued in the force expansion? Not absolutely
  !     sure what the really means yet
  INTEGER N
  !me:  The total number of quadrature points, and thus also the number of
  !     legrende polynominals.
  INTEGER LQ
  !me:  Number of intervals used in combination with the gaussian quadrature
  !     The integral to be estimated is split into this many subintervals and
  !     the gaussian quadrature is used for each interval INDEPENDENTLY.
  !     @todo Is there a name for this quadrature scheme?
  INTEGER NoQI
  !me:  @todo ???  
  INTEGER nocc

  !me:  REAL*8 is a double precision floating point type i.e. double in other 
  !     languages the *8 comes from the number of bytes (here 8) that are used 
  !     to store the number

  !me:  The slenderness parameter modelling how "thin"/slender the fiber is
  !     epsilon = a / 2*L, where L is the half length of the fiber, i.e. 2*L is
  !     the total fiber length
  !     @todo what is a?
  REAL*8 eeps
  !me:  @todo Probably the total elapsed time
  REAL*8 t
  !me:  The timestep
  REAL*8 dt
  !me:  @todo ???
  REAL*8 tvecMod
  !me:  The tolerance parameter used for GMRES
  REAL*8 tol

  !me:  REAL*8,ALLOCATABLE,DIMENSION(:,:):: defines a 2-dimensional array of
  !     type REAL*8 (i.e. double). ALLOCATABLE allows the size of the array
  !     to be set (aka allocate) later dynamically

  !me:  The position vectors of the fiber centers
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::XcVecs
  !me:  The orientation of the fiber, i.e. the unit tangent vector  
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::tVecs
  !me:  The linear velocity vectors of the fibers
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::VelVecs
  !me:  The rotational velocity vectors of the fibers
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::RotVecs

  !me:  The points used for integral estimation using gaussian quadrature along
  !     subintervals
  REAL*8,ALLOCATABLE,DIMENSION(:)::pv
  !me:  The weights used for integral estimation using gaussian quadrature along
  !     subintervals
  REAL*8,ALLOCATABLE,DIMENSION(:)::wv

  !me:  The legendre polynomials evaluated at the the quadrature points
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::LvecMat

  !me:  @todo The external forces acting on all fibers, e.g. gravity? 
  REAL*8,ALLOCATABLE,DIMENSION(:)::ExtForce
  !me:  @todo ???
  REAL*8,ALLOCATABLE,DIMENSION(:)::coeffvec
  !me:  @todo Probably used as an initial guess for some kind of integration 
  !     scheme
  REAL*8,ALLOCATABLE,DIMENSION(:)::init_guess

  !me:  A temporary variable holding the length of a single orientation vector.
  !     @todo Why does this need to be a vector? A simple REAL*8 should be
  !     enough.
  REAL*8,ALLOCATABLE,DIMENSION(:)::tmod0
  !me:  A temporary variable holding the length of the orientation vectors.
  REAL*8,ALLOCATABLE,DIMENSION(:)::tmod

  !me:  @todo ???
  INTEGER one
  !me:  @todo ???
  INTEGER two
  !me:  @todo ???
  INTEGER three
  !me:  @todo ???
  INTEGER old
  !me:  @todo ???
  INTEGER new
  !me:  The maximum number of iterations used for GMRES
  INTEGER max_iters
  !me:  The restart parameter used for GMRES
  INTEGER restart
  !me:  @todo ???
  INTEGER label
  !me:  The interval after which the current state is saved to the output file
  INTEGER save_ival
  !me:  The number of timesteps to simulate
  INTEGER no_ts
  !me:  @todo The number of saves in each file? Maybe this outputs to multiple
  !     files in case one file gets to large?
  INTEGER no_saves_in_file
  !me:  @todo ???
  INTEGER nos
  !me:  @todo ???
  INTEGER label2
  !me:  Flag indicating whether to solve the integral (@todo ref eq)
  !     analytically (=1) or numerically (=0)
  INTEGER int_an
  !me:  Flag indicating whether to use a direct solver (@todo what does direct
  !     mean?) (=1) or to use GMRES (=0)
  INTEGER dir_solve
  !me:  @todo ???
  INTEGER i
  !me:  @todo ???
  INTEGER j
  !me:  @todo ???
  INTEGER ind
  !me:  @todo ???
  INTEGER tt
  !me:  @todo ???
  INTEGER pp
  !me:  @todo ???
  INTEGER count_rate
  !me:  @todo ???
  INTEGER count_max
  !me:  @todo ???
  INTEGER count1
  !me:  @todo ???
  INTEGER count2

  !me:  @todo ???
  CHARACTER (LEN=10)::labstr
  !me:  @todo ???
  CHARACTER (LEN=10)::filenostr
  !me:  @todo ???
  CHARACTER (LEN=10)::labstr2

  !me:  @todo ???
  CHARACTER (LEN=25)::abbr
  !me:  @todo ???
  CHARACTER (LEN=25)::abbr2
  !me:  @todo ???
  CHARACTER (LEN=25)::repfile
  !me:  @todo ???
  CHARACTER (LEN=25)::filename
  !me:  @todo ???
  CHARACTER (LEN=25)::filename2
  !me:  @todo ???
  CHARACTER (LEN=25)::filename3
  !me:  @todo ???
  CHARACTER (LEN=25)::abbr3
  !me:  @todo ???
  CHARACTER (LEN=25)::filein
  !me:  @todo ???
  CHARACTER (LEN=25)::parfname  

  !me:  @todo ???
  REAL*8 lim
  !me:  @todo ???
  REAL*8,ALLOCATABLE,DIMENSION(:)::distVec
  !me:  @todo ???
  REAL*8 CPU_p


  !****************************************************************************!
  !                                                                            !
  ! Read Parameters from Standard In                                           !
  !                                                                            !
  ! This also allows to simply pipe in a configuration file                    !
  !                                                                            !
  !****************************************************************************!

  PRINT *,"Give label for run: "
  READ *,label
  PRINT *,"Give label for indatafile: "
  READ *,label2
  CALL int2str(labstr,label)
  CALL int2str(labstr2,label2)
  abbr="XcTData"//trim(labstr)//"_p"
  abbr2="VelRotData"//trim(labstr)//"_p"
  abbr3="Coeffvec"//trim(labstr)//"_p"
  !M=4;
  filein="XcT_init"//trim(labstr2)//".in"
  PRINT *,"Reading M and initial fiber positions from ",filein
  OPEN(33,file=trim(filein));
  READ (33,*) M
  PRINT *,"M= ",M,"."
  PRINT *,"Give no of terms in force expansion (N):"
  READ *,N 
  PRINT *,"Give eeps:"
  READ *,eeps 
  PRINT *,"Give dt:"
  READ *,dt 
  PRINT *,"Give number of time steps:"
  READ *,no_ts 
  PRINT *,"Give number of time steps between each save:"
  READ *,save_ival 
  PRINT *,"Give number of saves in each file:"
  READ *,no_saves_in_file
  PRINT*, "Give number of quad. intervals:"
  READ*,NoQI
  PRINT*, "Use analytical integration (1/0):"
  READ*,int_an
  PRINT*, "Use direct solver (1/0):"
  READ*,dir_solve
  IF (dir_solve==0) THEN
     PRINT*, "Give restart paramter for GMRES:"
     READ*,restart
     PRINT*, "Give maximum number of itrerations for GMRES:"
     READ*,max_iters
     PRINT*, "Give tolerance for GMRES:"
     READ*,tol
  END IF
  !!NoQI=8

  !me:  Initalize number of quadrature points. For a 3rd order gaussian
  !     quadrature this is 3 times the number of subintervals.
  LQ=3*NoQI
  nocc=3*N*M

  !****************************************************************************!
  !                                                                            !
  ! Allocate matrices and vectors                                              !
  !                                                                            !
  !****************************************************************************!

  !me:  The fiber matrices are all 3*Mx3. The 2-dimension is simply used as a
  !     kind of triple buffer for the timestepping
  !me:  Allocate the positions matrix.
  ALLOCATE(XcVecs(3*M,3))
  !me:  Allocate the orientation matrix.
  ALLOCATE(tVecs(3*M,3))
  !me:  Allocate the linear velocity matrix.
  ALLOCATE(VelVecs(3*M,3))
  !me:  Allocate the rotational velocity matrix.
  ALLOCATE(RotVecs(3*M,3))

  !me:  Allocate external force vector with 3 components per fiber.
  !     So probably the first 1-3 components are the external force vector
  !     acting on fiber 0 and components 4-6 act on fiber 1 at so on...
  !
  !     Why is the external force not the same for all particles? Is this maybe
  !     NOT a uniform external force but instead the 'external' force acting on
  !     the fiber as a result of all the other fibers?
  !     Answer it really appears to be simply gravity and for now at least IS
  !     the same for all fibers, as it is initalized only once and the same for
  !     each fiber. However that might be subject to change in the future. For
  !     this is an easy way to save some memory.   
  ALLOCATE(ExtForce(3*M))
  !me:  Allocate the coefficient vector with 3*N components per fiber.
  !     N here is the configurated number of terms used in the force expansion.
  !     @todo Still have to understand what exactly that means.
  ALLOCATE(coeffvec(3*N*M))
  !me:  Allocate the initial guess for the coefficients used for solver.
  ALLOCATE(init_guess(3*N*M))

  !me:  Allocate a temporary vector holding the length of each orientation
  !     vectors as a scalar
  ALLOCATE(tmod0(M))
  !me:  Allocate a temporary vector holding the length of each orientation
  !     vectors as in vector form, i.e. 3-components for each vector.
  !     @todo Why do we need a 3-component vector here per orientation vector?
  !     Does that make normalization easier/faster? One length per orientation
  !     vector should be enough? Potential memory saving opportunity
  ALLOCATE(tmod(3*M))

  !me:  Allocate the points and weights for the gaussian quadrature for all
  !     subintervals 
  ALLOCATE(pv(LQ),wv(LQ))
  !me:  Allocate legendre polynomial matrix where each row represents a
  !     quadrature point and each column entry the corresponding legendre
  !     polynominal evaluated at that point.
  ALLOCATE(LvecMat(LQ,N))
  !me:  @todo ???
  ALLOCATE(distVec(M*(M-1)))

  !****************************************************************************!
  !                                                                            !
  ! Initializing data matrix and vectors                                       !
  !                                                                            !
  !****************************************************************************!

  !me:  @todo these all appear to be some variables used for indexing however
  !     how exaclty that's used later is not yet clear
  one=1;
  two=2;
  three=3;
  old=1; 
  new=2;

  !me:  Fortran by default uses 1-based indexing.
  !me:  Initalizes the external force for all fibers. This sets the third
  !     component (i.e. z) to -1 with probably models gravity. Be aware of the
  !     index calculation 0-based vs 1-based!
  !     Result: (0,0,-1),(0,0,-1),...
  !me:  Fortran allows arrays (even multidimensional arrays) to be initalized
  !     with a single assignment
  ExtForce=0.0d0;
  !ExtForce(3:3:3*M)=-1.0d0
  DO j=1,M
     ind=(j-1)*3
     ExtForce(ind+3)=-1.0d0
  END DO

  !me:  Intializes the positons and @todo tVecs with data from the specified
  !     input file.
  XcVecs=0.0d0
  tVecs=0.0d0
  DO j=1,M
     ind=(j-1)*3
     READ(33,*) XcVecs(ind+1,three),XcVecs(ind+2,three),XcVecs(ind+3,three),tVecs(ind+1,three),tVecs(ind+2,three),tVecs(ind+3,three)
  END DO
  CLOSE(33)
  PRINT *,"Read initial data from file. "
  
  !! Set initial guess for iterative solver to zero
  DO j=1,3*M*N
     init_guess(j)=0.0d0
  END DO
  
  
  !! Normalize the tVecs
  DO i=1,M
     ind=(i-1)*3
     
     tmod0(i)=sqrt(tVecs(ind+1,three)**2+tVecs(ind+2,three)**2+tVecs(ind+3,three)**2)
     tmod(ind+1)=tmod0(i)
     tmod(ind+2)=tmod0(i)
     tmod(ind+3)=tmod0(i)
  END DO
  tVecs(:,three)=tVecs(:,three)/tmod;
  
  
  
  
  !!Computes the quadrature weights and 
  !! point to be used for numerical integration
  !! based on a three point Gauss quadrature. Section 4.2 in
  !! AKT-KG 2006
  CALL quad_pts_and_wts(NoQI,pv,wv) 
  LvecMat=compute_LVecs(N,LQ,pv);
  
  !=================================================================
  PRINT *,"=========== advect_fibers =============="
  PRINT *,"label= ", label,  "."
  PRINT *,"M= ", M,  ", N= " ,N, "."
  PRINT *,"eeps= ", eeps,  ", dt= " ,dt, "."
  PRINT *,"no_ts= ", no_ts,  ", save_ival= " ,save_ival, "."
  !=================================================================
  repfile="rep_"//trim(labstr)//".txt"
  PRINT *,"Writing report to file ",trim(repfile),"."
  PRINT*, int_an
  PRINT*, NoQI
  CALL write_Rep_file(repfile,abbr,label,M,N,NoQI,eeps,&
       dt,int_an,no_ts,save_ival,no_saves_in_file);
  
  parfname="pars"//trim(labstr)//"_.m"
  CALL write_Pars(parfname,label,M,N,eeps,&
       dt,no_ts,save_ival,no_saves_in_file);
  
  t=0.0d0;
  !!Need to save to file. 
  
  pp=1
  filename=trim(abbr)//"1"
  filename2=trim(abbr2)//"1"
  filename3=trim(abbr3)//"1"
  
  OPEN(10,file=trim(filename));
  OPEN(30,file=trim(filename2));
  OPEN(70,file=trim(filename3));
  PRINT *,"Saving to file ",trim(filename)
  PRINT *,"Saving Vel and Rotvel to file ",trim(filename2)
  DO i=1,M
     ind=(i-1)*3
     WRITE(10,'(2F24.16)') XcVecs(ind+1,three)
     WRITE(10,'(2F24.16)') XcVecs(ind+2,three)
     WRITE(10,'(2F24.16)') XcVecs(ind+3,three)
     WRITE(10,*) ' '
     WRITE(10,'(2F24.16)') tVecs(ind+1,three)
     WRITE(10,'(2F24.16)') tVecs(ind+2,three)
     WRITE(10,'(2F24.16)') tVecs(ind+3,three)
     WRITE(10,*) ' '
     WRITE(30,'(2F24.16)') 0.0;
     WRITE(30,'(2F24.16)') 0.0;
     WRITE(30,'(2F24.16)') 0.0;
     WRITE(30,*) ' '
     WRITE(30,'(2F24.16)') 0.0;
     WRITE(30,'(2F24.16)') 0.0;
     WRITE(30,'(2F24.16)') 0.0;
     WRITE(30,*) ' '
     
  END DO
  nos=1
  
  lim=5.0d0;
  
  
  !!================================================================================
  !!Starting the time-stepping. 
  !!===============================================================================
  
  CALL SYSTEM_CLOCK(count1,count_rate,count_max);

  PRINT *,"=========== Starting the time-stepping ==================="
  DO tt=1,no_ts
     IF (mod(tt,5)==0) THEN
        PRINT *,"time step no ",tt
     END IF
     !! Compute the coefficient in the force expansion using Legendre polynomials
     !! Use either direct solver for the linear system (dir_solve==1) or GMRES
     
     IF (dir_solve==1) THEN  !! Use direct solver
        IF (int_an==1) THEN  !! Use analytical integration for inner integral 
                             !! in Eq. (23) and numerical for outer
           coeffvec=solve_coeffs_an(eeps,M,N,XcVecs(:,three),tVecs(:,three),ExtForce,LQ,pv,wv,LvecMat);
           
        ELSE  !! Use only numerical quadrature
           coeffvec=solve_coeffs(eeps,M,N,XcVecs(:,three),tVecs(:,three),ExtForce,LQ,pv,wv,LvecMat);
           
        END IF
     ELSE  !! Use GMRES
        IF (int_an==1) THEN
           coeffvec=solve_coeffs_an_iter(eeps,M,N,XcVecs(:,three),tVecs(:,three),ExtForce,LQ,pv,wv,LvecMat,&
                restart,max_iters,tol,init_guess);
         ELSE
            
            coeffvec=solve_coeffs_iter(eeps,M,N,XcVecs(:,three),tVecs(:,three),ExtForce,LQ,pv,wv,LvecMat,&
                 restart,max_iters,tol,init_guess)
         END IF
      END IF
      
      init_guess=coeffvec;
      old=new; 
      new=3-old; 
      
      !! Compute velocitites i.e. the right hand sides in Eqns. (24) (25)
      IF (int_an==1) THEN
         CALL compute_velocities_an(eeps,M,N,LQ,pv,wv,LvecMat,&
              XcVecs(:,three),tVecs(:,three),ExtForce,coeffvec,&
              VelVecs(:,new),RotVecs(:,new))
      ELSE
         CALL compute_velocities(eeps,M,N,LQ,pv,wv,LvecMat,&
              XcVecs(:,three),tVecs(:,three),ExtForce,coeffvec,&
              VelVecs(:,new),RotVecs(:,new))
      END IF
      !Result returned in VelVecs and RotVecs
      one=two;
      two=three;
      three=6-(one+two);

      !me:  The triple buffer ordering for the first 4 iterations.
      !     During initialization only column 3 is filled.
      !     For the first loop a simple euler step is taken.
      !     @todo What is the reasoning behind choosing this particular
      !     integration scheme?
      !     
      !     original_index -> new_index => used_columns -> target_column
      !     1,2,3 -> 2,3,1 => 3     -> 1
      !     2,3,1 -> 3,1,2 => 1,3   -> 2
      !     3,1,2 -> 1,2,3 => 2,1   -> 3
      !     1,2,3 -> 2,3,1 => 3,2   -> 1
      
      !! Update postition and orientation by solving (24) and (25) in time
      IF (tt==1) THEN
         !!first time step. 
         XcVecs(:,three)=XcVecs(:,two)+dt*VelVecs(:,new);
         tVecs(:,three)=tVecs(:,two)+dt*RotVecs(:,new);
         PRINT *,"One first order ts."
         
      ELSE 
         XcVecs(:,three)=4.0d0/3.0d0*XcVecs(:,two)-1.0d0/3.0d0*XcVecs(:,one)+&
              2.0d0/3.0d0*dt*(2.0d0*VelVecs(:,new)-VelVecs(:,old));
         tVecs(:,three)=4.0d0/3.0d0*tVecs(:,two)-1.0d0/3.0d0*tVecs(:,one)+&
              2.0d0/3.0d0*dt*(2.0d0*RotVecs(:,new)-RotVecs(:,old));
      END IF
      
      !! Normalize, make orientation vectors have length one
      DO i=1,M
         ind=(i-1)*3
         
         tmod0(i)=sqrt(tVecs(ind+1,three)**2+tVecs(ind+2,three)**2+tVecs(ind+3,three)**2)
         tmod(ind+1)=tmod0(i)
         tmod(ind+2)=tmod0(i)
         tmod(ind+3)=tmod0(i)
      END DO
      
      tmod0=sqrt(tVecs(1:3:3*M,three)**2+tVecs(2:3:3*M,three)**2+tVecs(3:3:3*M,three)**2);
      
      tVecs(:,three)=tVecs(:,three)/tmod;
      
      
      t=t+dt;
      
      IF (mod(tt,save_ival)==0) THEN
         !Time to save
         nos=nos+1;
         PRINT *,"========================================"
         PRINT *,"Time step no ",tt,"."
         PRINT *,"========================================"
         PRINT *,"Saving to file ",trim(filename),", nos=",nos
         PRINT *,"Saving to file ",trim(filename2),", nos=",nos
         DO i=1,M
            ind=(i-1)*3
            WRITE(10,'(2F24.16)') XcVecs(ind+1,three)
            WRITE(10,'(2F24.16)') XcVecs(ind+2,three)
            WRITE(10,'(2F24.16)') XcVecs(ind+3,three)
            WRITE(10,*) ' '
            WRITE(10,'(2F24.16)') tVecs(ind+1,three)
            WRITE(10,'(2F24.16)') tVecs(ind+2,three)
            WRITE(10,'(2F24.16)') tVecs(ind+3,three)
            WRITE(10,*) ' '
            WRITE(30,'(2F24.16)') VelVecs(ind+1,new)
            WRITE(30,'(2F24.16)') VelVecs(ind+2,new)
            WRITE(30,'(2F24.16)') VelVecs(ind+3,new)
            WRITE(30,*) ' '
            WRITE(30,'(2F24.16)') RotVecs(ind+1,new)
            WRITE(30,'(2F24.16)') RotVecs(ind+2,new)
            WRITE(30,'(2F24.16)') RotVecs(ind+3,new)
            WRITE(30,*) ' '
         END DO
         DO i=1,3*M*N
            WRITE(70,'(2F24.16)') coeffvec(i)
         END DO
         WRITE(70,*) ' '
      END IF
     IF (mod(tt,save_ival*no_saves_in_file)==0 .AND. tt<no_ts) THEN
        pp=pp+1;
        CALL int2str(filenostr,pp)
        
        filename=trim(abbr)//filenostr
        filename2=trim(abbr2)//filenostr
        CLOSE(10)
        CLOSE(30)
        CLOSE(70)
        OPEN(10,file=trim(filename));
        OPEN(30,file=trim(filename2));
        OPEN(70,file=trim(filename3));
        PRINT *,"Opened new file: ",filename, "and", filename2
     END IF
     
  END DO
  
  
  
  PRINT *,"DONE!"
  PRINT *,"M= ",M,", N= ",N
  CALL SYSTEM_CLOCK(count2,count_rate,count_max);
  CPU_p=real(count2-count1)/count_rate

  PRINT *,"Taking ", no_ts," time steps took ",CPU_p," seconds." 
END PROGRAM ADVECT_FIBERS

