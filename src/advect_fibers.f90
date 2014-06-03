PROGRAM ADVECT_FIBERS
  !! Main program
  !! Given an inital configuration, the fiber system is advected in time
  USE quadrature
  USE MatrixAssembly
  USE write_Matlab
  USE coeffs_and_vels
  USE Int_analytical
  IMPLICIT NONE
  
  
  INTEGER M,N,LQ,NoQI,nocc
  REAL *8 eeps,t,dt,tvecMod,tol
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::XcVecs,tVecs,VelVecs,RotVecs
  REAL*8,ALLOCATABLE,DIMENSION(:)::pv,wv
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::LvecMat
  REAL*8,ALLOCATABLE,DIMENSION(:)::ExtForce,coeffvec,init_guess
  REAL*8,ALLOCATABLE,DIMENSION(:)::tmod0,tmod
  INTEGER one,two,three,old,new,max_iters,restart
  INTEGER label,save_ival,no_ts,no_saves_in_file,nos,label2
  INTEGER int_an, dir_solve
  INTEGER i,j,ind,tt,pp
  INTEGER count_rate,count_max,count1,count2
  CHARACTER (LEN=10):: labstr,filenostr,labstr2
  CHARACTER (LEN=25):: abbr,abbr2,repfile,filename, filename2,filename3,abbr3
  CHARACTER (LEN=25):: filein,parfname  
  REAL*8 lim
  REAL*8,ALLOCATABLE,DIMENSION(:)::distVec
  REAL*8 CPU_p
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
  LQ=3*NoQI
  nocc=3*N*M
  ALLOCATE(XcVecs(3*M,3),tVecs(3*M,3),VelVecs(3*M,3),RotVecs(3*M,3))
  ALLOCATE(ExtForce(3*M),coeffvec(3*N*M),init_guess(3*N*M))
  ALLOCATE(tmod0(M),tmod(3*M))
  ALLOCATE(pv(LQ),wv(LQ))
  ALLOCATE(LvecMat(LQ,N))
  ALLOCATE(distVec(M*(M-1)))

  one=1;
  two=2;
  three=3;
  old=1; 
  new=2;

  ExtForce=0.0d0;
  !ExtForce(3:3:3*M)=-1.0d0
  DO j=1,M
     ind=(j-1)*3
     ExtForce(ind+3)=-1.0d0
  END DO
  
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

