MODULE MatrixAssembly
USE quadrature
USE Int_analytical
IMPLICIT NONE

CONTAINS

FUNCTION compute_LVecs(N,LQ,pv);
  !! Computes and store the Legendre polynomials.
  !! Evaluated in the quadrature points. 

  !me:  see: http://en.wikipedia.org/wiki/Legendre_polynomials
  !     This also contains a listing of the formulas up to n = 10
  !     However we currently only allow up to n = 8 subintervals

  !me:  The number of subintervals
  INTEGER, INTENT(IN)::N
  !me:  The total number of quadrature points. For 3rd order gaussian quadrature
  !     this should always be 3 * N
  INTEGER, INTENT(IN)::LQ

  !me:  The points previously obtained by applying 3rd order gaussian quadrature
  !     to each interval
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv

  !me:  The output matrix containing the legendre polynomials evaluated at each
  !     quadrature point. So for each quadrature point we calculate each
  !     legendre polynomial up to the number of subintervals. The results is a
  !     matrix where each row represents a points and each column entry
  !     represents a legendre polynomial evaluated at that point.
  REAL*8,DIMENSION(LQ,N)::compute_LVecs

  !me: Temporary variable for iteration
  INTEGER kk

  DO kk=1,N
    IF (kk==1) THEN
      compute_Lvecs(:,kk)=pv;  
    ELSEIF (kk==2) THEN
      compute_Lvecs(:,kk)=0.5d0*(3.0d0*pv**2-1.0d0);  
    ELSEIF (kk==3) THEN
      compute_Lvecs(:,kk)=0.5d0*(5.0d0*pv**3-3.0d0*pv);  
    ELSEIF (kk==4) THEN
      compute_Lvecs(:,kk)=0.125d0*(35.0d0*pv**4-30.0d0*pv**2+3.0d0);
    ELSEIF (kk==5) THEN
      compute_Lvecs(:,kk)=0.125d0*(63.0d0*pv**5-70.0d0*pv**3+15.0d0*pv);
    ELSEIF (kk==6) THEN
       compute_Lvecs(:,kk)=0.0625d0*(231.0d0*pv**6-315.0d0*pv**4+105.0d0*pv**2-5.0d0);
    ELSEIF (kk==7) THEN
       compute_Lvecs(:,kk)=0.0625d0*(429.0d0*pv**7-693.0d0*pv**5+315.0d0*pv**3-35.0d0*pv);
    ELSEIF (kk==8) THEN
       compute_Lvecs(:,kk)=6435.0d0/128.0d0*pv**8-3003.0d0/32.0d0*pv**6+3465.0d0/64.0d0*pv**4-315.0d0/32.0d0*pv**2+35.0d0/128.0d0;
    ELSE
      PRINT *,"In compute_Lvecs: choice of N too large."
    END IF
  END DO
END FUNCTION compute_LVecs

!me:  UNUSED
SUBROUTINE compute_LvecMat(N,LQ,pv,LvecMat);
  !! Computes and store the Legendre polynomials.
  !! Evaluated in the quadrature points. 
  INTEGER, INTENT(IN):: N,LQ
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv
  REAL*8,DIMENSION(LQ,N),INTENT(OUT)::LvecMat

  INTEGER kk

  DO kk=1,N
    IF (kk==1) THEN
      LvecMat(:,kk)=pv;  
    ELSEIF (kk==2) THEN
      LvecMat(:,kk)=0.5d0*(3.0d0*pv**2-1.0d0);  
    ELSEIF (kk==3) THEN
      LvecMat(:,kk)=0.5d0*(5.0d0*pv**3-3.0d0*pv);  
    ELSEIF (kk==4) THEN
      LvecMat(:,kk)=0.125d0*(35.0d0*pv**4-30.0d0*pv**2+3.0d0);
    ELSEIF (kk==5) THEN
      LvecMat(:,kk)=0.125d0*(63.0d0*pv**5-70.0d0*pv**3+15.0d0*pv);
    ELSEIF (kk==6) THEN
       LvecMat(:,kk)=0.0625d0*(231.0d0*pv**6-315.0d0*pv**4+105.0d0*pv**2-5.0d0);
    ELSEIF (kk==7) THEN
       LvecMat(:,kk)=0.0625d0*(429.0d0*pv**7-693.0d0*pv**5+315.0d0*pv**3-35.0d0*pv);
    ELSEIF (kk==8) THEN
       LvecMat(:,kk)=6435.0d0/128.0d0*pv**8-3003.0d0/32.0d0*pv**6+3465.0d0/64.0d0*pv**4-315.0d0/32.0d0*pv**2+35.0d0/128.0d0;
    ELSE
      PRINT *,"In compute_LvecMat: choice of N too large."
    END IF
  END DO
END SUBROUTINE compute_LVecMat

SUBROUTINE assemble_force(N,LQ,LvecMat,ExtForce,ax_coeffs,ay_coeffs,az_coeffs,&
  fx,fy,fz);
  
  !! Given the coefficient, a_m^n, in Eq. (19). This subroutin compute
  !! the force f_m
  !!N: no of terms in force expansion 0..N.
  !LvecMat(:,k) contains Lvec no k, in LQ quadrature pts. 
  !Not Lvec no 0 (which is constant equal to 1).

  !ExtForce: 3 components, x,y, and z of external force. 
  !ax_coeffs: N x-coefficients, etc.  Vector of length k. 

  INTEGER, INTENT(IN):: N,LQ
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(N),INTENT(IN)::ax_coeffs,ay_coeffs,az_coeffs
  REAL*8,DIMENSION(3),INTENT(IN)::ExtForce

  REAL*8,DIMENSION(LQ),INTENT(OUT)::fx,fy,fz

  INTEGER kk
  fx=0.5d0*ExtForce(1)
  fy=0.5d0*ExtForce(2)
  fz=0.5d0*ExtForce(3)

  DO kk=1,N
    fx=fx+ax_coeffs(kk)*LvecMat(:,kk);
    fy=fy+ay_coeffs(kk)*LvecMat(:,kk);
    fz=fz+az_coeffs(kk)*LvecMat(:,kk);
  END DO
END SUBROUTINE assemble_force

SUBROUTINE assemble_forces(M,N,LQ,LvecMat,ExtForce,coeffvec,&
             Fmatx,Fmaty,Fmatz);
  !! Given the coefficient, a_m^n, in Eq. (19). This subroutin compute
  !! the force f_m for all fibers, m=1..M.
  !!N: no of terms in force expansion 0..N.
  !LvecMat(:,k) contains Lvec no k, in LQ quadrature pts. 
  !Not Lvec no 0 (which is constant equal to 1).

  !ExtForce: 3M components, x,y, and z of external force for each fiber. 
  !coeffvec: ax,ay,az: First for fiber no 1, for force term no 1..N, 
  !                    then for fiber no 2...

  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::coeffvec
  REAL*8,DIMENSION(3*M),INTENT(IN)::ExtForce

  REAL*8,DIMENSION(LQ,M),INTENT(OUT)::Fmatx,Fmaty,Fmatz
  INTEGER kk,i,ind

  DO i=1,M
    ind=(i-1)*3;
    Fmatx(:,i)=0.5d0*ExtForce(ind+1)
    Fmaty(:,i)=0.5d0*ExtForce(ind+2)
    Fmatz(:,i)=0.5d0*ExtForce(ind+3)

    DO kk=1,N
      ind=(i-1)*3*N+(kk-1)*3;
      Fmatx(:,i)=Fmatx(:,i)+coeffvec(ind+1)*LvecMat(:,kk);
      Fmaty(:,i)=Fmaty(:,i)+coeffvec(ind+2)*LvecMat(:,kk);
      Fmatz(:,i)=Fmatz(:,i)+coeffvec(ind+3)*LvecMat(:,kk);
    END DO
  END DO

END SUBROUTINE assemble_forces

FUNCTION mat_prod(xvec,eeps,M,N,XcVecs,tVecs,LQ,LvecMat,pv,wv);
  !!Computing A*xvec where A is the matrix for the linear 
  !!system for the force coefficients. Never assembling A.                             
  !pv, wv: Quad pts and weights. 
  !LvecMat(:,k) contains Lvecs from 1 to N, def in the quad pts.  

  !!sh set to 1: include shear flow.
  !!sh set to 0: do not include. 
  
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::xvec
  REAL*8,INTENT(IN)::eeps
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs

  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M*N):: mat_prod

  REAL*8::c,d,e,cc,D1,gammak,Ek
  REAL*8,DIMENSION(LQ,M):: Fmatx,Fmaty,Fmatz
  REAL*8,DIMENSION(3*M):: ExtForceZero
  REAL*8,DIMENSION(N,M):: ax_coeffs,ay_coeffs,az_coeffs
  

  REAL*8,DIMENSION(3):: xc,ta,xbar,IF_a,aac,xcb,tb

  REAL*8,DIMENSION(N):: lambdavec,Ekvec,gammak_vec
  REAL*8,DIMENSION(LQ,3):: Gmat
  INTEGER:: i,j,ii,ind,filno,fno,kk

  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  cc=1.0d0;
  D1=0.75d0/(d-2.0d0*cc);

  !Coeffs 0 to N, but computing for 1 to N.
  lambdavec(1)=2;
  DO i=2,N
    lambdavec(i)=lambdavec(i-1)+2.0d0/i;
  END DO
  Ekvec=(d-e-cc*lambdavec)/2.0d0/(d-cc*lambdavec);

  !!Forces stored columnwise for each fiber 1..M.
  !!Should not include ExtForce in assembly, this is only matrix, not 
  !!incl rhs. 
  ExtForceZero=0.0d0;

  !!CONTINUE FROM HERE.
  CALL assemble_forces(M,N,LQ,LvecMat,ExtForceZero,xvec,&
           Fmatx,Fmaty,Fmatz);

  mat_prod=xvec;
  !!Identiy times xvec. Now, add the rest. 
  DO filno=1,M
    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    !xc=XcMat(:,filno);
    !ta=PaMat(:,filno);
    DO fno=1,M
      IF (fno /= filno) THEN
        !!From fil fno to fil filno. 
        !disp(['filno = ' num2str(filno) ', fno= ' num2str(fno) '.']);
        ind=(fno-1)*3;
        xcb=XcVecs(ind+1:ind+3);
        tb=tVecs(ind+1:ind+3);
        DO i=1,LQ
          xbar=xc+pv(i)*ta;  
          Gmat(i,:)=G_f_GQ(xcb,tb,xbar,eeps,&
                     LQ,Fmatx(:,fno),Fmaty(:,fno),Fmatz(:,fno),pv,wv);  

        END DO
  
        DO i=1,3
          IF_a(i)=sum(wv*(Gmat(:,i)*LvecMat(:,1)));
        END DO

        !!Contr to x,y,z comp for first coefficient.
        aac=D1*sum(ta*IF_a)*ta;
        !!aac=-D1*sum(ta*IF_a)*ta+1.5d0*tau_c_ta;
        ind=(filno-1)*3*N;
        mat_prod(ind+1:ind+3)=mat_prod(ind+1:ind+3)+aac;

        DO kk=2,N
          Ek=Ekvec(kk);
          gammak=0.5d0*(2.0d0*kk+1.0d0)/(d+e-cc*lambdavec(kk));
          !!gammak=gammak_vec(kk);
          !!From fil fno to fil filno. 
          !!Int with L_k(s). 
          !!Inner integral over forces computed further up in loop.
 	  !!Now, integrating this towards Lvec_k
          DO i=1,3
            IF_a(i)=sum(wv*(Gmat(:,i)*LvecMat(:,kk)));
          END DO
          aac=gammak*(IF_a-Ek*sum(ta*IF_a)*ta);
          ind=(filno-1)*3*N+(kk-1)*3;
          mat_prod(ind+1:ind+3)=mat_prod(ind+1:ind+3)+aac;
        END DO;  !!for kk. 
      END IF; !!if fno!=filno
    END DO; !!for fno=1:M
  END DO; !!for filno=1:M


END FUNCTION mat_prod

SUBROUTINE assemble_matrix(eeps,M,N,XcVecs,tVecs,LQ,pv,wv,LvecMat,AMat);

  !! Assembles the coefficent matrix in Eq. (20).
  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs

  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M*N,3*M*N),INTENT(OUT):: AMat

  REAL*8::c,d,e,cc,D1,gammak,Ek
  REAL*8,DIMENSION(LQ,M):: Fmatx,Fmaty,Fmatz
  REAL*8,DIMENSION(LQ):: ExtForceZero

  REAL*8,DIMENSION(3):: xc,ta,xcb,tb,xbar,IF_a,aac

  REAL*8,DIMENSION(N):: lambdavec,Ekvec,gammak_vec
  REAL*8,DIMENSION(LQ,6):: Gmat
  REAL*8,DIMENSION(6):: thvec
  INTEGER:: i,j,ii,ind,filno,fno,kk,l,nocc,jj
  INTEGER:: rowno,p
  REAL*8 Q1,Q2,Q3
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p
  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  cc=1.0d0;
  D1=0.75d0/(d-2.0d0*cc);

  !Coeffs 0 to M, but computing for 1 to M.
  lambdavec(1)=2;
  DO i=2,N
    lambdavec(i)=lambdavec(i-1)+2.0d0/i;
  END DO
  Ekvec=(d-e-cc*lambdavec)/2.0d0/(d-cc*lambdavec);

  nocc=3*M*N;
  AMat=0.0d0;
  DO i=1,nocc
    AMat(i,i)=1.0d0
  END DO
  !!Order: a_x^1,a_y^1,a_z^1,a_x^2...a_x^N, a_y^N...etc for first fiber.
  !!Then the same for second. 

  !! Loop over all fibers 
  !! For fiber filno compute the interaction with all other fibers fno=1..M
  DO filno=1,M
    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    
    DO fno=1,M
      IF (fno /= filno) THEN !! The fiber does not interact with itself
       
        ind=(fno-1)*3;
        xcb=XcVecs(ind+1:ind+3);
        tb=tVecs(ind+1:ind+3);
        DO l=1,N  !! Loop over all N
          
          kk=1; 
          rowno=(filno-1)*3*N+1;
          
          
          DO i=1,LQ
            xbar=xc+pv(i)*ta;  
            Gmat(i,:)=G_compute_GQ(xcb,tb,xbar,eeps,LQ,pv,wv,LvecMat(:,l));
            !! Inner integral over G*P_l(s') in Eq. (23)
         END DO
          
          
          DO i=1,6
            thvec(i)=sum(wv*Gmat(:,i)*LvecMat(:,1));  !! Theta in Eq. (23)
          END DO  
          
          Q1=thvec(1)*ta(1)+thvec(4)*ta(2)+thvec(5)*ta(3);
          Q2=thvec(4)*ta(1)+thvec(2)*ta(2)+thvec(6)*ta(3);
          Q3=thvec(5)*ta(1)+thvec(6)*ta(2)+thvec(3)*ta(3);
          
          !!Will add up as we loop over fno, and over l=1..N. 
          p=(fno-1)*N*3+3*(l-1)+1;
          
          AMat(rowno,p)=D1*ta(1)*Q1;
          AMat(rowno,p+1)=D1*ta(1)*Q2;
          AMat(rowno,p+2)=D1*ta(1)*Q3;
  
          !!Second row for a_y^1.
          AMat(rowno+1,p)=D1*ta(2)*Q1;
          AMat(rowno+1,p+1)=D1*ta(2)*Q2;
          AMat(rowno+1,p+2)=D1*ta(2)*Q3;
  
          !!Second row for a_z^1.
          AMat(rowno+2,p)=D1*ta(3)*Q1;
          AMat(rowno+2,p+1)=D1*ta(3)*Q2;
          AMat(rowno+2,p+2)=D1*ta(3)*Q3;

          !!For higher k, same formula, so now we can loop. 
          DO kk=2,N
            rowno=(filno-1)*3*N+3*(kk-1)+1;
            Ek=Ekvec(kk);
            gammak=0.5d0*(2.0d0*kk+1.0d0)/(d+e-cc*lambdavec(kk));
            !!From fno to filno. 
            DO i=1,6
              thvec(i)=sum(wv*Gmat(:,i)*LvecMat(:,kk));
            END DO  
            Q1=thvec(1)*ta(1)+thvec(4)*ta(2)+thvec(5)*ta(3);
            Q2=thvec(4)*ta(1)+thvec(2)*ta(2)+thvec(6)*ta(3);
            Q3=thvec(5)*ta(1)+thvec(6)*ta(2)+thvec(3)*ta(3);
  
            !!Will add up as we loop over fno, and over l=1..N. 
            !!Take the absolute value of each coefficient before adding. 
            !!To the row of a^k_x for fiber filno:
            p=(fno-1)*N*3+3*(l-1)+1;
            
            !!First row for a_x^k.
            AMat(rowno,p)=gammak*(thvec(1)-Ek*ta(1)*Q1);
            AMat(rowno,p+1)=gammak*(thvec(4)-Ek*ta(1)*Q2);
            AMat(rowno,p+2)=gammak*(thvec(5)-Ek*ta(1)*Q3);
            !!Second row for a_y^k.
            AMat(rowno+1,p)=gammak*(thvec(4)-Ek*ta(2)*Q1);
            AMat(rowno+1,p+1)=gammak*(thvec(2)-Ek*ta(2)*Q2);
            AMat(rowno+1,p+2)=gammak*(thvec(6)-Ek*ta(2)*Q3);
            !!Third row for a_z^k.
            AMat(rowno+2,p)=gammak*(thvec(5)-Ek*ta(3)*Q1);
            AMat(rowno+2,p+1)=gammak*(thvec(6)-Ek*ta(3)*Q2);
            AMat(rowno+2,p+2)=gammak*(thvec(3)-Ek*ta(3)*Q3);
          END DO !!for kk=2:N.
        END DO;  !!for l=1..N
      END IF;  !!if (fno~=filno)
    END DO;  !!For fno=1:M
  END DO; !!For filno=1..M. 

 

END SUBROUTINE assemble_matrix


SUBROUTINE assemble_matrix_an(eeps,M,N,XcVecs,tVecs,LQ,pv,wv,LvecMat,AMat);
  !! Assembles the matrix using analytical and numerical integration.
  !! Coefficient in the matrix given by Eq. (20).

  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs

  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M*N,3*M*N),INTENT(OUT):: AMat

  REAL*8::c,d,e,cc,D1,gammak,Ek
  REAL*8,DIMENSION(LQ,M):: Fmatx,Fmaty,Fmatz
  REAL*8,DIMENSION(LQ):: ExtForceZero

  REAL*8,DIMENSION(3):: xc,ta,xcb,tb,xbar,IF_a,aac

  REAL*8,DIMENSION(N):: lambdavec,Ekvec,gammak_vec
  REAL*8,DIMENSION(LQ,6):: Gmat
  REAL*8,DIMENSION(6):: thvec
  INTEGER:: i,j,ii,ind,filno,fno,kk,l,nocc
  INTEGER:: rowno,p
  REAL*8 Q1,Q2,Q3
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p
  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  cc=1.0d0;
  D1=0.75d0/(d-2.0d0*cc);

  !Coeffs 0 to M, but computing for 1 to M.
  lambdavec(1)=2;
  DO i=2,N
    lambdavec(i)=lambdavec(i-1)+2.0d0/i;
  END DO
  Ekvec=(d-e-cc*lambdavec)/2.0d0/(d-cc*lambdavec);

  nocc=3*M*N;
  AMat=0.0d0;
  DO i=1,nocc
    AMat(i,i)=1.0d0
  END DO
  
  !!Order: a_x^1,a_y^1,a_z^1,a_x^2...a_x^N, a_y^N...etc for first fiber.
  !!Then the same for second. 

  !! Loop over all fibers 
  !! For fiber filno compute the interaction with all other fibers fno=1..M
  DO filno=1,M
    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    
    DO fno=1,M
      IF (fno /= filno) THEN
        
        ind=(fno-1)*3;
        xcb=XcVecs(ind+1:ind+3);
        tb=tVecs(ind+1:ind+3);
        DO l=1,N
          
          !!k=1, different formula. For higer k's: loop. 
          kk=1; 
          rowno=(filno-1)*3*N+1;
          !!From fno to filno. 
          !Gmat=zeros(6,LQ);  
          
          DO i=1,LQ
             xbar=xc+pv(i)*ta;  
             
             Gmat(i,:)=G_compute_GQ_kg(xcb,tb,xbar,eeps,l);
            
             
          END DO
          
          

          DO i=1,6
            thvec(i)=sum(wv*Gmat(:,i)*LvecMat(:,1));
          END DO  

          Q1=thvec(1)*ta(1)+thvec(4)*ta(2)+thvec(5)*ta(3);
          Q2=thvec(4)*ta(1)+thvec(2)*ta(2)+thvec(6)*ta(3);
          Q3=thvec(5)*ta(1)+thvec(6)*ta(2)+thvec(3)*ta(3);
          
          !!Will add up as we loop over fno, and over l=1..N. 
          p=(fno-1)*N*3+3*(l-1)+1;
          
          !!First row for a_x^1.
          AMat(rowno,p)=D1*ta(1)*Q1;
          AMat(rowno,p+1)=D1*ta(1)*Q2;
          AMat(rowno,p+2)=D1*ta(1)*Q3;
  
          !!Second row for a_y^1.
          AMat(rowno+1,p)=D1*ta(2)*Q1;
          AMat(rowno+1,p+1)=D1*ta(2)*Q2;
          AMat(rowno+1,p+2)=D1*ta(2)*Q3;
  
          !!Second row for a_z^1.
          AMat(rowno+2,p)=D1*ta(3)*Q1;
          AMat(rowno+2,p+1)=D1*ta(3)*Q2;
          AMat(rowno+2,p+2)=D1*ta(3)*Q3;
  
          !!For higher k, same formula, so now we can loop. 
          DO kk=2,N
            rowno=(filno-1)*3*N+3*(kk-1)+1;
            Ek=Ekvec(kk);
            gammak=0.5d0*(2.0d0*kk+1.0d0)/(d+e-cc*lambdavec(kk));
            !!From fno to filno. 
            DO i=1,6
              thvec(i)=sum(wv*Gmat(:,i)*LvecMat(:,kk));
            END DO  
            Q1=thvec(1)*ta(1)+thvec(4)*ta(2)+thvec(5)*ta(3);
            Q2=thvec(4)*ta(1)+thvec(2)*ta(2)+thvec(6)*ta(3);
            Q3=thvec(5)*ta(1)+thvec(6)*ta(2)+thvec(3)*ta(3);
  
            !!Will add up as we loop over fno, and over l=1..N. 
            !!Take the absolute value of each coefficient before adding. 
            !!To the row of a^k_x for fiber filno:
            p=(fno-1)*N*3+3*(l-1)+1;
            
            !!First row for a_x^k.
            AMat(rowno,p)=gammak*(thvec(1)-Ek*ta(1)*Q1);
            AMat(rowno,p+1)=gammak*(thvec(4)-Ek*ta(1)*Q2);
            AMat(rowno,p+2)=gammak*(thvec(5)-Ek*ta(1)*Q3);
            !!Second row for a_y^k.
            AMat(rowno+1,p)=gammak*(thvec(4)-Ek*ta(2)*Q1);
            AMat(rowno+1,p+1)=gammak*(thvec(2)-Ek*ta(2)*Q2);
            AMat(rowno+1,p+2)=gammak*(thvec(6)-Ek*ta(2)*Q3);
            !!Third row for a_z^k.
            AMat(rowno+2,p)=gammak*(thvec(5)-Ek*ta(3)*Q1);
            AMat(rowno+2,p+1)=gammak*(thvec(6)-Ek*ta(3)*Q2);
            AMat(rowno+2,p+2)=gammak*(thvec(3)-Ek*ta(3)*Q3);
          END DO !!for kk=2:N.
        END DO;  !!for l=1..N
      END IF;  !!if (fno~=filno)
    END DO;  !!For fno=1:M
  END DO; !!For filno=1..M. 

  
END SUBROUTINE assemble_matrix_an


SUBROUTINE assemble_rhs(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,Brhs);
  !! Assemles the RHS in Eq. (20)
  REAL*8::ffvec_x,ffvec_y,ffvec_z !! added this since need scalar force!!!!!!
  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce
  
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M*N),INTENT(OUT):: Brhs

  REAL*8::c,d,e,cc,D1,gammak,Ek
  REAL*8,DIMENSION(LQ):: Fvec_x,Fvec_y,Fvec_z

  REAL*8,DIMENSION(3):: xc,ta,xbar,contr,xcb,tb

  REAL*8,DIMENSION(N):: lambdavec,Ekvec
  REAL*8,DIMENSION(LQ,3):: Gmat
  INTEGER:: i,j,ii,ind,filno,fno,kk,nocc
  INTEGER:: rowno,p
  REAL*8 Q1,Q2,Q3,ta_dot_ac
  INTEGER count_rate,count_max,count1,count2
  REAL*8 CPU_p
  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  cc=1.0d0;
  D1=0.75d0/(d-2.0d0*cc);
  
  !Coeffs 0 to N, but computing for 1 to N.
  lambdavec(1)=2;
  DO i=2,N
    lambdavec(i)=lambdavec(i-1)+2.0d0/i;
  END DO;
  Ekvec=(d-e-cc*lambdavec)/2.0d0/(d-cc*lambdavec);

  nocc=3*M*N;
  !!Order: a_x^1,a_y^1,a_z^1,a_x^2,etc for first fib.
  !!Then the same for second. 

  Brhs=0.0d0
  DO filno=1,M
    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    
    DO fno=1,M
      IF (fno /= filno) THEN
        !!From fno to filno. 
  
        !!kk=1, different formula. For higer kk's: loop. 
        kk=1; 
        rowno=(filno-1)*3*N+1;
        ind=(fno-1)*3;
        Fvec_x=0.5d0*ExtForce(ind+1)
        Fvec_y=0.5d0*ExtForce(ind+2)
        Fvec_z=0.5d0*ExtForce(ind+3)
        
       
        xcb=XcVecs(ind+1:ind+3);
        tb=tVecs(ind+1:ind+3);
        
        DO i=1,LQ
          xbar=xc+pv(i)*ta;  
          
          
          Gmat(i,:)=G_f_GQ(xcb,tb,xbar,eeps,&
               LQ,Fvec_x,Fvec_y,Fvec_z,pv,wv)
             
          
          END DO
          
        DO i=1,3
          contr(i)=sum(wv*Gmat(:,i)*LvecMat(:,1));
        END DO
  
        ta_dot_ac=ta(1)*contr(1)+ta(2)*contr(2)+ta(3)*contr(3);
        Brhs(rowno)=Brhs(rowno)-D1*ta_dot_ac*ta(1);
        Brhs(rowno+1)=Brhs(rowno+1)-D1*ta_dot_ac*ta(2);
        Brhs(rowno+2)=Brhs(rowno+2)-D1*ta_dot_ac*ta(3);
        
  
        !!For higher k, same formula, so now we can loop. 
        DO kk=2,N
          rowno=(filno-1)*3*N+3*(kk-1)+1;
          Ek=Ekvec(kk);
          gammak=0.5d0*(2.0d0*kk+1.0d0)/(d+e-cc*lambdavec(kk));
          DO i=1,3
            contr(i)=sum(wv*Gmat(:,i)*LvecMat(:,kk));
          END DO
          ta_dot_ac=ta(1)*contr(1)+ta(2)*contr(2)+ta(3)*contr(3);
  
          
          Brhs(rowno)=Brhs(rowno)-gammak*(contr(1)-Ek*ta_dot_ac*ta(1));
          Brhs(rowno+1)=Brhs(rowno+1)-gammak*(contr(2)-Ek*ta_dot_ac*ta(2));
          Brhs(rowno+2)=Brhs(rowno+2)-gammak*(contr(3)-Ek*ta_dot_ac*ta(3));
          !disp(['For the third of these rows, adding ' num2str(-gammak*(contr(3)-Ek*ta_dot_ac*ta(3))) '.']);
  
        END DO !!for kk=2:N.
      END IF;  !!if (fno~=filno)
    END DO;  !!For fno=1:M
  END DO; !!For filno=1..M. 

END SUBROUTINE assemble_rhs

SUBROUTINE assemble_rhs_an(eeps,M,N,XcVecs,tVecs,ExtForce,LQ,pv,wv,LvecMat,Brhs);
  !M no of fils. N: no of terms in force exp. 
  !Total length of vector: N*3*M.
  REAL*8::ffvec_x,ffvec_y,ffvec_z !! added this since need scalar force!!!!!!
  REAL*8,INTENT(IN)::eeps
  INTEGER, INTENT(IN):: M,N,LQ
  REAL*8,DIMENSION(3*M),INTENT(IN)::XcVecs,tVecs,ExtForce
  
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ,N),INTENT(IN)::LvecMat
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  REAL*8,DIMENSION(3*M*N),INTENT(OUT):: Brhs

  REAL*8::c,d,e,cc,D1,gammak,Ek
  REAL*8,DIMENSION(LQ):: Fvec_x,Fvec_y,Fvec_z

  REAL*8,DIMENSION(3):: xc,ta,xbar,contr,xcb,tb

  REAL*8,DIMENSION(N):: lambdavec,Ekvec
  REAL*8,DIMENSION(LQ,3):: Gmat
  INTEGER:: i,j,ii,ind,filno,fno,kk,nocc
  INTEGER:: rowno,p
  REAL*8 Q1,Q2,Q3,ta_dot_ac

  c=log(eeps**2*exp(1.0d0));
  d=-c;
  e=2.0d0;
  cc=1.0d0;
  D1=0.75d0/(d-2.0d0*cc);
  
  !Coeffs 0 to N, but computing for 1 to N.
  lambdavec(1)=2;
  DO i=2,N
    lambdavec(i)=lambdavec(i-1)+2.0d0/i;
  END DO;
  Ekvec=(d-e-cc*lambdavec)/2.0d0/(d-cc*lambdavec);

  nocc=3*M*N;
  !!Order: a_x^1,a_y^1,a_z^1,a_x^2,etc for first fil.
  !!Then the same for second. 

  Brhs=0.0d0
  DO filno=1,M
    ind=(filno-1)*3;
    xc=XcVecs(ind+1:ind+3);
    ta=tVecs(ind+1:ind+3);
    
    DO fno=1,M
      IF (fno /= filno) THEN
        !!From fno to filno. 
  
        !!kk=1, different formula. For higer kk's: loop. 
        kk=1; 
        rowno=(filno-1)*3*N+1;
        ind=(fno-1)*3;
        Fvec_x=0.5d0*ExtForce(ind+1)
        Fvec_y=0.5d0*ExtForce(ind+2)
        Fvec_z=0.5d0*ExtForce(ind+3)
        
        
        ffvec_x=0.0d0*ExtForce(ind+1);
        ffvec_y=0.0d0*ExtForce(ind+2);
        ffvec_z=0.5d0*ExtForce(ind+3);
       
        xcb=XcVecs(ind+1:ind+3);
        tb=tVecs(ind+1:ind+3);
        DO i=1,LQ
          xbar=xc+pv(i)*ta;  
          
         
          Gmat(i,:)=G_f_GQ_kg(0,xcb,tb,xbar,&
               ffvec_x,ffvec_y,ffvec_z,eeps);
          
          
          
          END DO
  
        DO i=1,3
          contr(i)=sum(wv*Gmat(:,i)*LvecMat(:,1));
        END DO
  
        ta_dot_ac=ta(1)*contr(1)+ta(2)*contr(2)+ta(3)*contr(3);
        Brhs(rowno)=Brhs(rowno)-D1*ta_dot_ac*ta(1);
        Brhs(rowno+1)=Brhs(rowno+1)-D1*ta_dot_ac*ta(2);
        Brhs(rowno+2)=Brhs(rowno+2)-D1*ta_dot_ac*ta(3);
        
  
        !!For higher k, same formula, so now we can loop. 
        DO kk=2,N
          rowno=(filno-1)*3*N+3*(kk-1)+1;
          Ek=Ekvec(kk);
          gammak=0.5d0*(2.0d0*kk+1.0d0)/(d+e-cc*lambdavec(kk));
          DO i=1,3
            contr(i)=sum(wv*Gmat(:,i)*LvecMat(:,kk));
          END DO
          ta_dot_ac=ta(1)*contr(1)+ta(2)*contr(2)+ta(3)*contr(3);
  
          
          Brhs(rowno)=Brhs(rowno)-gammak*(contr(1)-Ek*ta_dot_ac*ta(1));
          Brhs(rowno+1)=Brhs(rowno+1)-gammak*(contr(2)-Ek*ta_dot_ac*ta(2));
          Brhs(rowno+2)=Brhs(rowno+2)-gammak*(contr(3)-Ek*ta_dot_ac*ta(3));
          
  
        END DO !!for kk=2:N.
      END IF;  !!if (fno~=filno)
    END DO;  !!For fno=1:M
  END DO; !!For filno=1..M. 

END SUBROUTINE assemble_rhs_an


END MODULE MatrixAssembly

