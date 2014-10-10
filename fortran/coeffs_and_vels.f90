MODULE coeffs_and_vels
USE quadrature
USE Int_analytical
USE MatrixAssembly

IMPLICIT NONE

CONTAINS

FUNCTION solve_coeffs(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat)
  !! Sets up and solve for the coefficients in the Legendre expansion.
  !! Basically setting up and solving Eq. (20). This function uses 
  !! only numerical integration and solve the system by a direct solver. 
  !! Returns the coefficients, a_m^n in Eq. (19)


  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce

  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M*N)::solve_coeffs

  REAL*8,DIMENSION(3*M*N,3*M*N):: AMat
  REAL*8,DIMENSION(3*M*N):: Brhs
  INTEGER info,nocc,j,i
  REAL*8,DIMENSION(3*M*N)::IPIV
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p
  nocc=3*M*N

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  !! Assembly of the coefficient matrix for solving Aa=f Eq. (20)
  CALL assemble_matrix(eeps,M,N,XcVecs,tVecs,LQ,pv,wv,LvecMat,AMat);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling the matrix took ",CPU_p," seconds."

!  OPEN(10,file="AMat.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (AMat(i,j),j=1,3*M*N)
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  !! Assembly of the right hand side for solving Aa=f Eq. (20)
  CALL assemble_rhs(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,Brhs);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling the right hand side took ",CPU_p," seconds."
  
!  OPEN(10,file="BVec.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (Brhs(i))
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL dgesv(nocc,1,AMat,nocc,IPIV,Brhs,nocc,info)
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Solving system using direct solver took ",CPU_p," seconds."
 
  solve_coeffs=Brhs
  OPEN(10,file="XVec.out");
  DO i=1,3*M*N
    WRITE(10,'(*(F16.8))') (solve_coeffs(i))
  END DO
  CLOSE(10)

  !DO j=1,M
  !   PRINT *,"coeff, element ",(j-1)*3*N+6,":",solve_coeffs((j-1)*3*N+6)
  !END DO

END FUNCTION solve_coeffs

FUNCTION solve_coeffs_an(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat)

  !! Sets up and solve for the coefficients in the Legendre expansion.
  !! Basically setting up and solving Eq. (20). This function uses 
  !! analytical and numerical integration and solve the system by a direct solver. 

  !me:  Unused
  INTEGER count_rate,count_max,count1,count2
  !me:  Unused
  REAL*8 CPU_p

  !me:  The slenderness parameter
  REAL*8,INTENT(IN)::eeps
  !me:  The number of fibers to be simulated
  INTEGER, INTENT(IN)::M
  !me:  @todo The number of terms (e.g. 5) sued in the force expansion? Not absolutely
  !     sure what the really means yet
  INTEGER, INTENT(IN)::N
  !me:  The total number of quadrature points
  INTEGER, INTENT(IN)::LQ
  !me:  The position vectors of the fiber centers
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs
  !me:  The orientation of the fiber, i.e. the unit tangent vector  
  REAL*8,DIMENSION(3*M),INTENT(IN)::tVecs
  !me:  The external forces acting on each fiber
  REAL*8,DIMENSION(3*M),INTENT(IN)::ExtForce

  !pv,wv,and Lvec_m should be of length LQ. 
  !me:  The legendre polynominals evaluated at each point
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  !me:  The points and weights of the piecewise gaussian quadrature
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  !me:  The final solutions for the coefficients
  REAL*8,DIMENSION(3*M*N)::solve_coeffs_an

  !me:  The final matrix A of the linear system of equations in Ax=b
  REAL*8,DIMENSION(3*M*N,3*M*N):: AMat
 
  !me:  The final right hand side b of the linear of equations in Ax=b
  REAL*8,DIMENSION(3*M*N):: Brhs

  !me:  Unused
  INTEGER info,j,i
  !me:  @todo
  INTEGER nocc

  !me:  Unused
  REAL*8,DIMENSION(3*M*N)::IPIV
  
  !me:  @todo
  nocc=3*M*N
  
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL assemble_matrix_an(eeps,M,N,XcVecs,tVecs,LQ,pv,wv,LvecMat,AMat);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling the matrix took ",CPU_p," seconds."

!  OPEN(10,file="AMat.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (AMat(i,j),j=1,3*M*N)
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL assemble_rhs_an(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,Brhs);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling the right hand side took ",CPU_p," seconds."

!  OPEN(10,file="BVec.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (Brhs(i))
!  END DO
!  CLOSE(10)
  
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL dgesv(nocc,1,AMat,nocc,IPIV,Brhs,nocc,info)
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Solving system using direct solver took ",CPU_p," seconds."
  
!  solve_coeffs_an=Brhs
!  OPEN(10,file="XVec.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (solve_coeffs_an(i))
!  END DO
!  CLOSE(10)

  

END FUNCTION solve_coeffs_an

FUNCTION solve_coeffs_iter(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,&
     restart,max_iters,tol,init_guess)

  !! Sets up and solve for the coefficients in the Legendre expansion.
  !! Basically setting up and solving Eq. (20). This function uses 
  !! only numerical integration and solve the system by GMRES. 
  !! Returns the coefficients, a_m^n in Eq. (19)

  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce

  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,INTENT(IN)::tol
  INTEGER,INTENT(IN)::restart,max_iters
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::init_guess
  REAL*8,DIMENSION(3*M*N)::solve_coeffs_iter

  REAL*8,DIMENSION(3*M*N,3*M*N):: AMat
  REAL*8,DIMENSION(3*M*N):: Brhs
  REAL*8,DIMENSION(3*M*N):: ccvec

  INTEGER LUinfo,nocc,j,i
  REAL*8,DIMENSION(3*M*N)::IPIV

  INTEGER lwork,done,iter

  INTEGER revcom,colx,coly,colz,nbscal
  INTEGER,DIMENSION(5)::irc
  INTEGER,DIMENSION(8)::icntl
  INTEGER,DIMENSION(3)::info

  !PARAMETER (lwork=restart**2+restart*(3*M*N+5)+5*(3*M*N)+1)
  !REAL*8,DIMENSION(lwork):: work
  !Added one to minimum
  REAL*8,DIMENSION(restart**2+restart*(3*M*N+5)+5*(3*M*N)+2):: work
  REAL*8,DIMENSION(5):: cntl
  REAL*8,DIMENSION(2):: rinfo

  INTEGER matvec,precondLeft,precondRight,dotProd
  PARAMETER (matvec=1,precondLeft=2,precondRight=3,dotProd=4)

  REAL*8 ZERO,ONE
  PARAMETER(ZERO=0.0d0,ONE=1.0d0)

  INTEGER count_rate,count_max,count1,count2,count3,count4
  REAL*8 CPU_p

  nocc=3*M*N


  CALL SYSTEM_CLOCK(count1,count_rate,count_max);
  
  CALL assemble_matrix(eeps,M,N,XcVecs,tVecs,LQ,pv,wv,LvecMat,AMat);
  CALL SYSTEM_CLOCK(count2,count_rate,count_max);
  CPU_p=real(count2-count1)/count_rate
  PRINT *,"Assemble matrix took ",CPU_p

!  OPEN(10,file="AMat.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (AMat(i,j),j=1,3*M*N)
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1,count_rate,count_max);
  CALL assemble_rhs(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,Brhs);
  CALL SYSTEM_CLOCK(count2,count_rate,count_max);
  CPU_p=real(count2-count1)/count_rate
  PRINT *,"Assemble rhs took ",CPU_p

!  OPEN(10,file="BVec.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (Brhs(i))
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1,count_rate,count_max);
  !PRINT *,"restart= ",restart
  lwork=restart**2+restart*(3*M*N+5)+5*(3*M*N)+2
  work=0.0d0
  work(1:nocc)=init_guess
  work(nocc+1:2*nocc)=Brhs
  !PRINT *,"lwork= ",lwork

  !Initialize control parameters to default values. 
  call init_dgmres(icntl,cntl)
  !Tune some parameters. 
  !Tolerance
  cntl(1)=tol
  !Do not print optimal calue for lwork. 
  icntl(2)=0
  !Output in file fort.20
  icntl(3)=20
  !No preconditioning
  icntl(4)=0

  !orthogonalization choice. 
  !0/1/2/3 yields MGS/IMGS/CGS/ICGS
  !Started with 3.
  icntl(5)=0

  !Initial guess supplied if icntl(6) set to 1. 
  icntl(6)=1
  !max no of iters. 
  icntl(7)=max_iters

  done=0
  iter=1
  DO WHILE (done==0) 
     CALL drive_dgmres(nocc,nocc,restart,lwork,work,irc,icntl,cntl,info,rinfo)
     revcom=irc(1)
     colx=irc(2)
     coly=irc(3)
     colz=irc(4)
     nbscal=irc(5)
     !PRINT *,"iter= ",iter,", revcom= ",revcom

     IF (revcom == matvec) THEN 
        !Perform the matrix vector product
        !work(colz)<--A*work(colx)
        !PRINT *,"Computing matrix vector product"
        !CALL SYSTEM_CLOCK(count3,count_rate,count_max);
        CALL dgemv('N',nocc,nocc,ONE,AMat,nocc,work(colx),1,&
             ZERO,work(colz),1)
        !CALL SYSTEM_CLOCK(count4,count_rate,count_max);
        !CPU_p=real(count4-count3)/count_rate
        !PRINT *,"Computing matrix-vecor product ",CPU_p
        !!alpha*A*x+beta*y, alpha set to 1 (arg no 4), beta to 0 (arg no 9). 
        !!The 1 after work(colx) is increment in that vector. 
        !!The 1 after work(colz) is increment in that vector. 
        iter=iter+1
     ELSEIF (revcom == dotProd) THEN 
        !PRINT *,"Computing dot product"
        !work(colz)<--work(colx) work(coly)
        CALL dgemv('C',nocc,nbscal,ONE,work(colx),nocc,&
             work(coly),1,ZERO,work(colz),1)
        
     ELSE
        !PRINT *,"NO MORE matrix products. revcom= ",revcom
        !PRINT *,"No of matrix-vec products: ",iter
        done=1
     END IF
  END DO
  PRINT*,"info", info(1)
  IF (info(1)<0) THEN 
     PRINT *,"GMRES returned with info(1)= ",info(1)
     PRINT*,"With a backward error of: ", rinfo(2)
     IF (info(1)==-3) THEN 
        PRINT *,"Workspace too small. lwork= ",lwork,", need ",info(2)
     END IF
     ELSE IF (info(1)==0) THEN
     PRINT*,"GMRES conv. in no. of iter.: ", info(2)
     PRINT*,"With a backward error of: ", rinfo(2)
  END IF
  IF (info(1)==-4) THEN 
     PRINT*,"Solve with direct solver"
     CALL dgesv(3*M*N,1,AMat,3*M*N,IPIV,Brhs,3*N*M,info)
     
     solve_coeffs_iter=Brhs
  ELSE
     ccvec=work(1:nocc)
     solve_coeffs_iter=ccvec
  END IF

  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Solving system using GMRES took ",CPU_p," seconds."


  OPEN(10,file="XVec.out");
  DO i=1,3*M*N
    WRITE(10,'(*(F16.8))') (solve_coeffs_iter(i))
  END DO
  CLOSE(10)


  

 


END FUNCTION solve_coeffs_iter


FUNCTION solve_coeffs_an_iter(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,&
     restart,max_iters,tol,init_guess)
  !! Sets up and solve for the coefficients in the Legendre expansion.
  !! Basically setting up and solving Eq. (20). This function uses 
  !! analytical and  numerical integration and solve the system by GMRES. 
  !! Returns the coefficients, a_m^n in Eq. (19)


  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce

  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,INTENT(IN)::tol
  INTEGER,INTENT(IN)::restart,max_iters
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::init_guess
  REAL*8,DIMENSION(3*M*N)::solve_coeffs_an_iter

  REAL*8,DIMENSION(3*M*N,3*M*N):: AMat
  REAL*8,DIMENSION(3*M*N):: Brhs
  REAL*8,DIMENSION(3*M*N):: ccvec

  INTEGER LUinfo,nocc,j,i
  REAL*8,DIMENSION(3*M*N)::IPIV

  INTEGER lwork,done,iter

  INTEGER revcom,colx,coly,colz,nbscal
  INTEGER,DIMENSION(5)::irc
  INTEGER,DIMENSION(8)::icntl
  INTEGER,DIMENSION(3)::info

  !PARAMETER (lwork=restart**2+restart*(3*M*N+5)+5*(3*M*N)+1)
  !REAL*8,DIMENSION(lwork):: work
  !Added one to minimum
  REAL*8,DIMENSION(restart**2+restart*(3*M*N+5)+5*(3*M*N)+2):: work
  REAL*8,DIMENSION(5):: cntl
  REAL*8,DIMENSION(2):: rinfo

  INTEGER matvec,precondLeft,precondRight,dotProd
  PARAMETER (matvec=1,precondLeft=2,precondRight=3,dotProd=4)

  REAL*8 ZERO,ONE
  PARAMETER(ZERO=0.0d0,ONE=1.0d0)

  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p

  nocc=3*M*N
 
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL assemble_matrix_an(eeps,M,N,XcVecs,tVecs,LQ,pv,wv,LvecMat,AMat);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling the matrix took ",CPU_p," seconds."

!  OPEN(10,file="AMat.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (AMat(i,j),j=1,3*M*N)
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL assemble_rhs_an(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,Brhs);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling the right hand side took ",CPU_p," seconds."
 
!  OPEN(10,file="BVec.out");
!  DO i=1,3*M*N
!    WRITE(10,'(*(F16.8))') (Brhs(i))
!  END DO
!  CLOSE(10)

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  !PRINT *,"restart= ",restart
  lwork=restart**2+restart*(3*M*N+5)+5*(3*M*N)+2
  work=0.0d0
  work(1:nocc)=init_guess
  work(nocc+1:2*nocc)=Brhs
  !PRINT *,"lwork= ",lwork

  !Initialize control parameters to default values. 
  call init_dgmres(icntl,cntl)
  !Tune some parameters. 
  !Tolerance
  cntl(1)=tol
  !Do not print optimal calue for lwork. 
  icntl(2)=0
  !Output in file fort.20
  icntl(3)=20
  !No preconditioning
  icntl(4)=0

  !orthogonalization choice. 
  !0/1/2/3 yields MGS/IMGS/CGS/ICGS
  !Started with 3.
  icntl(5)=0

  !Initial guess supplied if icntl(6) set to 1. 
  icntl(6)=1
  !max no of iters. 
  icntl(7)=max_iters

  done=0
  iter=1
  DO WHILE (done==0) 
     CALL drive_dgmres(nocc,nocc,restart,lwork,work,irc,icntl,cntl,info,rinfo)
     revcom=irc(1)
     colx=irc(2)
     coly=irc(3)
     colz=irc(4)
     nbscal=irc(5)
     !PRINT *,"iter= ",iter,", revcom= ",revcom

     IF (revcom == matvec) THEN 
        !Perform the matrix vector product
        !work(colz)<--A*work(colx)
        !PRINT *,"Computing matrix vector product"
        CALL dgemv('N',nocc,nocc,ONE,AMat,nocc,work(colx),1,&
             ZERO,work(colz),1)
        !!alpha*A*x+beta*y, alpha set to 1 (arg no 4), beta to 0 (arg no 9). 
        !!The 1 after work(colx) is increment in that vector. 
        !!The 1 after work(colz) is increment in that vector. 
        iter=iter+1
     ELSEIF (revcom == dotProd) THEN 
        !PRINT *,"Computing dot product"
        !work(colz)<--work(colx) work(coly)
        CALL dgemv('C',nocc,nbscal,ONE,work(colx),nocc,&
             work(coly),1,ZERO,work(colz),1)
        
     ELSE
        !PRINT *,"NO MORE matrix products. revcom= ",revcom
        !PRINT *,"No of matrix-vec products: ",iter
        done=1
     END IF
  END DO
  PRINT*,"info", info(1)
  IF (info(1)<0) THEN 
     PRINT *,"GMRES returned with info(1)= ",info(1)
     PRINT*,"With a backward error of: ", rinfo(2)
     IF (info(1)==-3) THEN 
        PRINT *,"Workspace too small. lwork= ",lwork,", need ",info(2)
     END IF
  ELSE IF (info(1)==0) THEN
     PRINT*,"GMRES conv. in no. of iter.: ", info(2)
     PRINT*,"With a backward error of: ", rinfo(2)
  END IF
  IF (info(1)==-4) THEN 
     PRINT*,"Solve with direct solver"
     CALL dgesv(3*M*N,1,AMat,3*M*N,IPIV,Brhs,3*N*M,info)
     
     solve_coeffs_an_iter=Brhs
  ELSE
     ccvec=work(1:nocc)
     solve_coeffs_an_iter=ccvec
  END IF

  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Solving system using GMRES took ",CPU_p," seconds."

  OPEN(10,file="XVec.out");
  DO i=1,3*M*N
    WRITE(10,'(*(F16.8))') (solve_coeffs_an_iter(i))
  END DO
  CLOSE(10)

END FUNCTION solve_coeffs_an_iter










SUBROUTINE compute_velocities(eeps,M,N,LQ,pv,wv,LvecMat,&
                XcVecs,tVecs,ExtForce,coeffvec,VelVecs,RotVecs)
  !! Compute the right hand side in Eqns. (24) and (25)
  !! using only numerical integration.

  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::coeffvec
  REAL*8,DIMENSION(3*M),INTENT(INOUT)::VelVecs,RotVecs

  REAL*8::c,d,e,sh,mubar
  INTEGER:: filno,fno,ind,ind2
  REAL*8,DIMENSION(LQ,M):: Fmatx,Fmaty,Fmatz
  REAL*8,DIMENSION(LQ):: LvecZero
  REAL*8::ta_dot_Fa
  REAL*8,DIMENSION(3):: xc,ta,xcb,tb,IF_a0,IF_a1
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p

  INTEGER i
  
  !PRINT *,"IN compute_velocities"
  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  mubar=d;
  sh=0.0d0;
  LvecZero=1.0d0

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL assemble_forces(M,N,LQ,LvecMat,ExtForce,coeffvec,&
             Fmatx,Fmaty,Fmatz);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling forces took ",CPU_p," seconds."
  
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  !PRINT *,"Assembled forces"
  !!From fil b to fil a. 
  DO filno=1,M
    !PRINT *,"filno= ",filno

    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    ta_dot_Fa=sum(ta*ExtForce(ind+1:ind+3));

    !!Local contribution to velocity. 
    VelVecs(ind+1:ind+3)=0.5d0*((d+e)*ExtForce(ind+1:ind+3)+(d-e)*ta_dot_Fa*ta);
    VelVecs(ind+1)=VelVecs(ind+1)+sh*mubar*xc(2);
                                          !!Y-coord of xc for this fil. 
    RotVecs(ind+1:ind+3)=-sh*mubar*ta(2)*ta(1)*ta;
    RotVecs(ind+1)=RotVecs(ind+1)+sh*mubar*ta(2);

    DO fno=1,M
      IF (fno/=filno) THEN
        !!Add IF_a0 and IF_a1 contributions to the velocities. 
        !!From fil fno to fil filno.  
        ind2=(fno-1)*3;
        xcb=XcVecs(ind2+1:ind2+3);
        tb=tVecs(ind2+1:ind2+3);
        
        
        
        IF_a0=TH_f_GQ(xcb,tb,xc,ta,eeps,Fmatx(:,fno),Fmaty(:,fno),Fmatz(:,fno),&
              LQ,pv,wv,LvecZero);
        IF_a1=TH_f_GQ(xcb,tb,xc,ta,eeps,Fmatx(:,fno),Fmaty(:,fno),Fmatz(:,fno),&
             LQ,pv,wv,LvecMat(:,1));
        
        
        !!Adding on non-local contributions to velocities. 
        VelVecs(ind+1:ind+3)=VelVecs(ind+1:ind+3)+0.5d0*IF_a0;
        RotVecs(ind+1:ind+3)=RotVecs(ind+1:ind+3)+1.5d0*(IF_a1-sum(ta*IF_a1)*ta);

      END IF
    END DO
  END DO

  VelVecs=VelVecs/mubar;
  RotVecs=RotVecs/mubar;

  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Updating velocities took ",CPU_p," seconds."

  OPEN(10,file="TRANSVel.out");
  DO i=1,3*M
    WRITE(10,'(*(F16.8))') (VelVecs(i))
  END DO
  CLOSE(10)
  OPEN(10,file="ROTVel.out");
  DO i=1,3*M
    WRITE(10,'(*(F16.8))') (RotVecs(i))
  END DO
  CLOSE(10)
  
END SUBROUTINE compute_velocities


SUBROUTINE compute_velocities_an(eeps,M,N,LQ,pv,wv,LvecMat,&
                XcVecs,tVecs,ExtForce,coeffvec,VelVecs,RotVecs)

  !! Compute the right hand side in Eqns. (24) and (25)
  !! using analytical and  numerical integration.
  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::coeffvec
  REAL*8,DIMENSION(3*M),INTENT(INOUT)::VelVecs,RotVecs

  REAL*8::c,d,e,sh,mubar
  INTEGER:: filno,fno,ind,ind2
  REAL*8,DIMENSION(LQ,M):: Fmatx,Fmaty,Fmatz
  REAL*8,DIMENSION(LQ):: LvecZero
  REAL*8::ta_dot_Fa
  REAL*8,DIMENSION(3):: xc,ta,xcb,tb,IF_a0,IF_a1
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p
  INTEGER i


  
  !PRINT *,"IN compute_velocities"
  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  mubar=d;
  sh=0.0d0;
  LvecZero=1.0d0
  
  !! Do not need to do this when using analytical integration
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  CALL assemble_forces(M,N,LQ,LvecMat,ExtForce,coeffvec,&
             Fmatx,Fmaty,Fmatz);
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Assembling forces took ",CPU_p," seconds."

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  !!From fil b to fil a. 
  DO filno=1,M
    

    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    ta_dot_Fa=sum(ta*ExtForce(ind+1:ind+3));

    !!Local contribution to velocity. 
    VelVecs(ind+1:ind+3)=0.5d0*((d+e)*ExtForce(ind+1:ind+3)+(d-e)*ta_dot_Fa*ta);
    VelVecs(ind+1)=VelVecs(ind+1)+sh*mubar*xc(2);
                                          !!Y-coord of xc for this fil. 
    RotVecs(ind+1:ind+3)=-sh*mubar*ta(2)*ta(1)*ta;
    RotVecs(ind+1)=RotVecs(ind+1)+sh*mubar*ta(2);

    DO fno=1,M
      IF (fno/=filno) THEN
        !!Add IF_a0 and IF_a1 contributions to the velocities. 
        !!From fil fno to fil filno.  
        ind2=(fno-1)*3;
        xcb=XcVecs(ind2+1:ind2+3);
        tb=tVecs(ind2+1:ind2+3);
        
        IF_a0=TH_f_GQ_kg(M,N,xcb,tb,xc,ta,eeps,ExtForce,coeffvec,LQ,pv,wv,& 
             LvecZero,fno); !! added here
        
        IF_a1=TH_f_GQ_kg(M,N,xcb,tb,xc,ta,eeps,ExtForce,coeffvec,LQ,pv,wv,&
                       LvecMat(:,1),fno); !! added here
        
                 
        !!Adding on non-local contributions to velocities. 
        VelVecs(ind+1:ind+3)=VelVecs(ind+1:ind+3)+0.5d0*IF_a0;
        RotVecs(ind+1:ind+3)=RotVecs(ind+1:ind+3)+1.5d0*(IF_a1-sum(ta*IF_a1)*ta);
      END IF
    END DO
  END DO

  VelVecs=VelVecs/mubar;
  RotVecs=RotVecs/mubar;
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  CPU_p = real(count2-count1)/count_rate
  PRINT *,"Updating velocities took ",CPU_p," seconds."

!  OPEN(10,file="TRANSVel.out");
!  DO i=1,3*M
!    WRITE(10,'(*(F16.8))') (VelVecs(i))
!  END DO
!  CLOSE(10)
!  OPEN(10,file="ROTVel.out");
!  DO i=1,3*M
!    WRITE(10,'(*(F16.8))') (RotVecs(i))
!  END DO
!  CLOSE(10)

END SUBROUTINE compute_velocities_an



END MODULE coeffs_and_vels

