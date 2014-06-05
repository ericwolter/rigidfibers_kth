MODULE Int_analytical
IMPLICIT NONE
CONTAINS

SUBROUTINE Analytical_int(xb,pb,xbar,N,L11,L12,L13,L14,L22,L23,L24)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::N
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,DIMENSION(3)::R_0
  REAL*8::clim
  !REAL*8,DIMENSION(N+3)::Analytical_int,I,J,S
  REAL*8,DIMENSION(30)::I,J,S
  REAL*8,DIMENSION(N+1),INTENT(OUT)::L11,L12,L13,L14,L22,L23,L24
  !!REAL*8::Analytical_int
  REAL*8::b,c,se,sb,ue,ub,te,tb,d,nn
  INTEGER::ii
  I=0.0d0;
  J=0.0d0;
  S=0.0d0;
  clim=10;
  R_0=(xbar-xb);
  b=-2.0d0*SUM(pb*R_0);
  c=SUM(R_0**2.0d0);
  se=1.0d0;
  sb=-1.0d0;
  ue=sqrt(1.0d0+b+c);
  ub=sqrt(1.0d0-b+c);
  te=1.0d0+0.5d0*b;
  tb=-1.0d0+0.5d0*b;
  d=c-0.25d0*b*b;
  
  
  I(1) = LOG(ABS(2.0d0*se+b+2.0d0*ue))-LOG(ABS(2.0d0*sb + b + 2.0d0*ub));
  I(2) = ue-ub-1/2.0d0*b*I(1);
  IF (c<clim) THEN
     nn=2.0d0; 
     DO ii=3,N+3
        
        nn=nn+1.0d0;
        I(ii) = -se**(ii-2)*ue/(1.0d0-nn) + sb**(ii-2)*ub/(1.0d0-nn) &
             - (0.5d0-(nn-1.0d0))*b/(1.0d0-nn)*I(ii-1) &
             + (nn-2.0d0)*c/(1.0d0-nn)*I(ii-2);
     END DO
  ELSE
     I(30)=0.0d0;
     I(29)=0.0d0;
     nn=27.0d0;
     DO ii=28,3,-1
        
        I(ii) =(nn+2)/((nn+1)*c)*(se**(nn+1)*ue/(nn+2) - sb**(nn+1)*ub/(nn+2) + &
             (1-2*(nn+2))/(2*(nn+2))*b*I(ii+1)- I(ii+2));
        nn=nn-1.0d0;
     END DO
  END IF
 
 
  
  IF (d<1e-14) THEN
     
     J(1) = -1.0d0/(2.0d0*(se+b/2.0d0)**2.0d0) + 1.0d0/(2.0d0*(sb+b/2.0d0)**2.0d0);
     
  ELSE
     J(1) = te/(d*ue)-tb/(d*ub);
     
  END IF
  J(2) =  -1.0d0/ue+1.0d0/ub - b/2.0d0*J(1);

  IF (c<clim) THEN
     DO ii=3,N+3
        J(ii) = I(ii-2) - b*J(ii-1)-c*J(ii-2);
        
     END DO
  ELSE
     J(30) = 0;
     J(29)= 0;
     DO ii=28,3,-1
        J(ii) = 1/c*(I(ii)-b*J(ii+1)-J(ii+2));
     END DO
  END IF
  

 
  
  
  IF (d<1e-14) THEN 
     S(1) = -4.0d0/(2.0d0*se+b)**4.0d0 + 4.0d0/(2.0d0*sb+b)**4.0d0;
     S(2) = -8.0d0/(3.0d0*(2.0d0*se+b)**3.0d0)+8.0d0/(3.0d0*(2.0d0*sb+b)**3.0d0)-b/2.0d0*S(1);
  
  ELSE
     S(1) = (2.0d0*se+b)/(6.0d0*d*ue**3.0d0)-(2.0d0*sb+b)/(6.0d0*d*ub**3.0d0)+2.0d0/(3.0d0*d)*J(1);
     S(2) = -(b*se+2.0d0*c)/(6.0d0*d*ue**3.0d0)+(b*sb+2.0d0*c)/(6.0d0*d*ub**3.0d0)-b/(3.0d0*d)*J(1);
  END IF
  IF (c<clim) THEN
     DO ii=3,N+3
        S(ii) = J(ii-2) - b*S(ii-1)-c*S(ii-2);
     END DO
  ELSE
  S(30)=0;
  S(29)=0;
  DO  ii=28,3,-1
     S(ii)=1/c*(J(ii)-b*S(ii+1)-S(ii+2));
  END DO
END IF
 
 

DO ii=1,N+1
  IF (ii==1 .OR. ii==2) THEN
     L11(ii) = I(ii);
     L12(ii) = J(ii);
     L13(ii) = J(ii+1);
     L14(ii) = J(ii+2);
     L22(ii) = S(ii);
     L23(ii) = S(ii+1);
     L24(ii) = S(ii+2);
    
  ELSEIF (ii==3) THEN
     L11(ii) = 0.5*(3.0d0*I(ii)-I(ii-2));
     L12(ii) = 0.5*(3.0d0*J(ii)-J(ii-2));
     L13(ii) = 0.5*(3.0d0*J(ii+1)-J(ii-1));
     L14(ii) = 0.5*(3.0d0*J(ii+2)-J(ii));
     L22(ii) = 0.5*(3.0d0*S(ii)-S(ii-2));
     L23(ii) = 0.5*(3.0d0*S(ii+1)-S(ii-1));
     L24(ii) = 0.5*(3.0d0*S(ii+2)-S(ii));
     
  ELSEIF (ii==4) THEN
     L11(ii) = 0.5d0*(5.0d0*I(ii)-3.0d0*I(ii-2));
     L12(ii) = 0.5d0*(5.0d0*J(ii)-3.0d0*J(ii-2));
     L13(ii) = 0.5d0*(5.0d0*J(ii+1)-3.0d0*J(ii-1));
     L14(ii) = 0.5d0*(5.0d0*J(ii+2)-3.0d0*J(ii));
     L22(ii) = 0.5d0*(5.0d0*S(ii)-3.0d0*S(ii-2));
     L23(ii) = 0.5d0*(5.0d0*S(ii+1)-3.0d0*S(ii-1));
     L24(ii) = 0.5d0*(5.0d0*S(ii+2)-3.0d0*S(ii));
     
  ELSEIF (ii==5) THEN
     L11(ii) = 0.125d0*(35.0d0*I(ii)-30.0d0*I(ii-2)+3.0d0*I(ii-4));
     L12(ii) = 0.125d0*(35.0d0*J(ii)-30.0d0*J(ii-2)+3.0d0*J(ii-4));
     L13(ii) = 0.125d0*(35.0d0*J(ii+1)-30.0d0*J(ii-1)+3.0d0*J(ii-3));
     L14(ii) = 0.125d0*(35.0d0*J(ii+2)-30.0d0*J(ii)+3.0d0*J(ii-2));
     L22(ii) = 0.125d0*(35.0d0*S(ii)-30.0d0*S(ii-2)+3.0d0*S(ii-4));
     L23(ii) = 0.125d0*(35.0d0*S(ii+1)-30.0d0*S(ii-1)+3.0d0*S(ii-3));
     L24(ii) = 0.125d0*(35.0d0*S(ii+2)-30.0d0*S(ii)+3.0d0*S(ii-2));
  ELSEIF (ii==6) THEN
     L11(ii) = 0.125d0*(63.0d0*I(ii)-70.0d0*I(ii-2)+15.0d0*I(ii-4));
     L12(ii) = 0.125d0*(63.0d0*J(ii)-70.0d0*J(ii-2)+15.0d0*J(ii-4));
     L13(ii) = 0.125d0*(63.0d0*J(ii+1)-70.0d0*J(ii-1)+15.0d0*J(ii-3));
     L14(ii) = 0.125d0*(63.0d0*J(ii+2)-70.0d0*J(ii)+15.0d0*J(ii-2));
     L22(ii) = 0.125d0*(63.0d0*S(ii)-70.0d0*S(ii-2)+15.0d0*S(ii-4));
     L23(ii) = 0.125d0*(63.0d0*S(ii+1)-70.0d0*S(ii-1)+15.0d0*S(ii-3));
     L24(ii) = 0.125d0*(63.0d0*S(ii+2)-70.0d0*S(ii)+15.0d0*S(ii-2));
  ELSEIF (ii==7) THEN
     L11(ii) = 0.0625d0*(231.0d0*I(ii)-315.0d0*I(ii-2)+105.0d0*I(ii-4)-5.0d0*I(ii-6));
     L12(ii) = 0.0625d0*(231.0d0*J(ii)-315.0d0*J(ii-2)+105.0d0*J(ii-4)-5.0d0*J(ii-6));
     L13(ii) = 0.0625d0*(231.0d0*J(ii+1)-315.0d0*J(ii-1)+105.0d0*J(ii-3)-5.0d0*J(ii-5));
     L14(ii) = 0.0625d0*(231.0d0*J(ii+2)-315.0d0*J(ii)+105.0d0*J(ii-2)-5.0d0*J(ii-4));
     L22(ii) = 0.0625d0*(231.0d0*S(ii)-315.0d0*S(ii-2)+105.0d0*S(ii-4)-5.0d0*S(ii-6));
     L23(ii) = 0.0625d0*(231.0d0*S(ii+1)-315.0d0*S(ii-1)+105.0d0*S(ii-3)-5.0d0*S(ii-5));
     L24(ii) = 0.0625d0*(231.0d0*S(ii+2)-315.0d0*S(ii)+105.0d0*S(ii-2)-5.0d0*S(ii-4)); 
  ELSEIF (ii==8) THEN
     L11(ii) = 0.0625d0*(429.0d0*I(ii)-693.0d0*I(ii-2)+315.0d0*I(ii-4)-35.0d0*I(ii-6));
     L12(ii) = 0.0625d0*(429.0d0*J(ii)-693.0d0*J(ii-2)+315.0d0*J(ii-4)-35.0d0*J(ii-6));
     L13(ii) = 0.0625d0*(429.0d0*J(ii+1)-693.0d0*J(ii-1)+315.0d0*J(ii-3)-35.0d0*J(ii-5));
     L14(ii) = 0.0625d0*(429.0d0*J(ii+2)-693.0d0*J(ii)+315.0d0*J(ii-2)-35.0d0*J(ii-4));
     L22(ii) = 0.0625d0*(429.0d0*S(ii)-693.0d0*S(ii-2)+315.0d0*S(ii-4)-35.0d0*S(ii-6));
     L23(ii) = 0.0625d0*(429.0d0*S(ii+1)-693.0d0*S(ii-1)+315.0d0*S(ii-3)-35.0d0*S(ii-5));
     L24(ii) = 0.0625d0*(429.0d0*S(ii+2)-693.0d0*S(ii)+315.0d0*S(ii-2)-35.0d0*S(ii-4));
  ELSEIF (ii==9) THEN
     L11(ii) = 1.0d0/128.0d0*(6435.0d0*I(ii)-12012.0d0*I(ii-2)+6930.0d0*I(ii-4)-1260.0d0*I(ii-6)+35.0d0*I(ii-8));
     L12(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii)-12012.0d0*J(ii-2)+6930.0d0*J(ii-4)-1260.0d0*J(ii-6)+35.0d0*J(ii-8));
     L13(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii+1)-12012.0d0*J(ii-1)+6930.0d0*J(ii-3)-1260.0d0*J(ii-5)+35.0d0*J(ii-7));
     L14(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii+2)-12012.0d0*J(ii)+6930.0d0*J(ii-2)-1260.0d0*J(ii-4)+35.0d0*J(ii-6));
     L22(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii)-12012.0d0*S(ii-2)+6930.0d0*S(ii-4)-1260.0d0*S(ii-6)+35.0d0*S(ii-8));
     L23(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii+1)-12012.0d0*S(ii-1)+6930.0d0*S(ii-3)-1260.0d0*S(ii-5)+35.0d0*S(ii-7));
     L24(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii+2)-12012.0d0*S(ii)+6930.0d0*S(ii-2)-1260.0d0*S(ii-4)+35.0d0*S(ii-6));

     ELSEIF (ii>9) THEN 
     PRINT*,"K  too large"
  END IF

END DO




END SUBROUTINE Analytical_int


SUBROUTINE Analytical_intII(xb,pb,xbar,N,L11,L12,L13,L14,L22,L23,L24)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::N
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,DIMENSION(3)::R_0
  !REAL*8,DIMENSION(N+3)::Analytical_int,I,J,S
  REAL*8,DIMENSION(N+3)::I,J,S
  REAL*8,DIMENSION(N+1),INTENT(OUT)::L11,L12,L13,L14,L22,L23,L24
  !!REAL*8::Analytical_int
  REAL*8::b,c,se,sb,ue,ub,te,tb,d
  INTEGER::ii
  
  R_0=(xbar-xb);
  b=-2.0d0*SUM(pb*R_0);
  c=SUM(R_0**2);
  se=1.0d0;
  sb=-1.0d0;
  ue=1.0d0+b+c;
  ub=1.0d0-b+c;
  te=1.0d0+0.5d0*b;
  tb=-1.0d0+0.5d0*b;
  d=c-0.25d0*b*b;
  
  IF (d<1e-14) THEN
     
     I(1) = LOG(ABS(se+b/2.0d0))-LOG(ABS(sb+b/2.0d0));
     DO ii=2,N+3
        I(ii) = 1.0d0/(ii-1)*(se**(ii-1)-sb**(ii-1) - b/2.0d0*(ii-1)*I(ii-1));
     END DO
  ELSE
     
     I(1) = LOG(ABS(2.0d0*se+b+2.0d0*SQRT(ue)))-LOG(ABS(2.0d0*sb + b + 2.0d0*SQRT(ub))); 
     
     I(2) = SQRT(ue)-SQRT(ub) - b/2.0d0*I(1);
     I(3) = se*SQRT(ue)/2.0d0-sb*SQRT(ub)/2.0d0 - 3.0d0*b/4.0d0*I(2) - c/2*I(1);
     
     DO ii=4,N+3
        I(ii) = -se**(ii-2)*SQRT(ue)/(1-ii) + sb**(ii-2)*SQRT(ub)/(1-ii) - (0.5d0-(ii-1))*b/(1-ii)*I(ii-1) + (ii-2)*c/(1-ii)*I(ii-2);
     END DO
  END IF
!!$  DO ii=1,N+3
!!$     PRINT*,"I(",ii,") = ", I(ii)
!!$  END DO
  IF (d<1e-14) THEN
    
   J(1) = -1.0d0/(2.0d0*((se+b/2.0d0)**2)) + 1.0d0/(2.0d0*((sb+b/2.0d0)**2));
   
   J(2) = b/(2.0d0*(2.0d0*(se+b/2.0d0)**2)) - 1.0d0/(se+b/2.0d0)-b/(2.0d0*(2.0d0*(sb+b/2.0d0)**2))+1.0d0/(sb+b/2.0d0);
   
   J(3) = LOG(ABS(se+b/2.0d0))-LOG(ABS(sb+b/2.0d0))+b/(se+b/2.0d0)-b/(sb+b/2.0d0)-b**2/(8.0d0*(se+b/2.0d0)**2)
   J(3)=J(3)+b**2/(8.0d0*(sb+b/2.0d0)**2);

   
   DO ii=4,N+3
      J(ii) = -1.0d0/(2.0d0-(ii-1))*(se**(ii-1)/(se+b/2.0d0)**2 - sb**(ii-1)/(sb+b/2.0d0)**2 -(ii-1)*b/2.0d0*J(ii-1));
   END DO
   
ELSE
   J(1) = te/(d*SQRT(ue))-tb/(d*SQRT(ub));
  
   J(2) = -1.0d0/SQRT(ue)+1.0d0/SQRT(ub) - b/2.0d0*J(1);
   
   DO ii=3,N+3
      J(ii) = I(ii-2) - b*J(ii-1)-c*J(ii-2);
       
   END DO
END IF

IF (d<1e-14) THEN 
   
   S(1) = -1.0d0/(4.0d0*(se+b/2.0d0)**4) + 1.0d0/(4.0d0*(sb+b/2.0d0)**4);
   S(2) = -1.0d0/(3.0d0*(se+b/2.0d0)**3)+1.0d0/(3.0d0*(sb+b/2.0d0)**3) 
   S(2) = S(2) +b/2.0d0*1.0d0/(4.0d0*(se+b/2.0d0)**4)-b/2.0d0*1.0d0/(4.0d0*(sb+b/2.0d0)**4);
ELSE
  
   S(1) = 2.0d0*te/(6.0d0*d*ue**(1.5d0)) -2.0d0*tb/(6.0d0*d*ub**(1.5d0))+2.0d0/(3.0d0*d)*J(1);
   S(2) = -(b*(te-b/2.0d0)+2.0d0*c)/(6.0d0*d*ue**(1.5d0)) +(b*(tb-b/2.0d0)+2.0d0*c)/(6.0d0*d*ub**(1.5d0)) -b/3.0d0/d*J(1);
END IF

DO ii=3,N+3
   S(ii) = J(ii-2) - b*S(ii-1)-c*S(ii-2);
END DO


DO ii=1,N+3
   PRINT *,"I(",ii,") = ", I(ii)
END DO

DO ii=1,N+3
   PRINT *,"J(",ii,") = ", J(ii)
END DO
DO ii=1,N+3
   PRINT *,"S(",ii,") = ", S(ii)
END DO

DO ii=1,N+1
  IF (ii==1 .OR. ii==2) THEN
     L11(ii) = I(ii);
     L12(ii) = J(ii);
     L13(ii) = J(ii+1);
     L14(ii) = J(ii+2);
     L22(ii) = S(ii);
     L23(ii) = S(ii+1);
     L24(ii) = S(ii+2);
    
  ELSEIF (ii==3) THEN
     L11(ii) = 0.5*(3.0d0*I(ii)-I(ii-2));
     L12(ii) = 0.5*(3.0d0*J(ii)-J(ii-2));
     L13(ii) = 0.5*(3.0d0*J(ii+1)-J(ii-1));
     L14(ii) = 0.5*(3.0d0*J(ii+2)-J(ii));
     L22(ii) = 0.5*(3.0d0*S(ii)-S(ii-2));
     L23(ii) = 0.5*(3.0d0*S(ii+1)-S(ii-1));
     L24(ii) = 0.5*(3.0d0*S(ii+2)-S(ii));
     
  ELSEIF (ii==4) THEN
     L11(ii) = 0.5d0*(5.0d0*I(ii)-3.0d0*I(ii-2));
     L12(ii) = 0.5d0*(5.0d0*J(ii)-3.0d0*J(ii-2));
     L13(ii) = 0.5d0*(5.0d0*J(ii+1)-3.0d0*J(ii-1));
     L14(ii) = 0.5d0*(5.0d0*J(ii+2)-3.0d0*J(ii));
     L22(ii) = 0.5d0*(5.0d0*S(ii)-3.0d0*S(ii-2));
     L23(ii) = 0.5d0*(5.0d0*S(ii+1)-3.0d0*S(ii-1));
     L24(ii) = 0.5d0*(5.0d0*S(ii+2)-3.0d0*S(ii));
     
  ELSEIF (ii==5) THEN
    L11(ii) = 0.125d0*(35.0d0*I(ii)-30.0d0*I(ii-2)+3.0d0*I(ii-4));
    L12(ii) = 0.125d0*(35.0d0*J(ii)-30.0d0*J(ii-2)+3.0d0*J(ii-4));
    L13(ii) = 0.125d0*(35.0d0*J(ii+1)-30.0d0*J(ii-1)+3.0d0*J(ii-3));
    L14(ii) = 0.125d0*(35.0d0*J(ii+2)-30.0d0*J(ii)+3.0d0*J(ii-2));
    L22(ii) = 0.125d0*(35.0d0*S(ii)-30.0d0*S(ii-2)+3.0d0*S(ii-4));
    L23(ii) = 0.125d0*(35.0d0*S(ii+1)-30.0d0*S(ii-1)+3.0d0*S(ii-3));
    L24(ii) = 0.125d0*(35.0d0*S(ii+2)-30.0d0*S(ii)+3.0d0*S(ii-2));
  ELSEIF (ii==6) THEN
     L11(ii) = 0.125d0*(63.0d0*I(ii)-70.0d0*I(ii-2)+15.0d0*I(ii-4));
     L12(ii) = 0.125d0*(63.0d0*J(ii)-70.0d0*J(ii-2)+15.0d0*J(ii-4));
     L13(ii) = 0.125d0*(63.0d0*J(ii+1)-70.0d0*J(ii-1)+15.0d0*J(ii-3));
     L14(ii) = 0.125d0*(63.0d0*J(ii+2)-70.0d0*J(ii)+15.0d0*J(ii-2));
     L22(ii) = 0.125d0*(63.0d0*S(ii)-70.0d0*S(ii-2)+15.0d0*S(ii-4));
     L23(ii) = 0.125d0*(63.0d0*S(ii+1)-70.0d0*S(ii-1)+15.0d0*S(ii-3));
     L24(ii) = 0.125d0*(63.0d0*S(ii+2)-70.0d0*S(ii)+15.0d0*S(ii-2));
  ELSEIF (ii==7) THEN
     L11(ii) = 0.0625d0*(231.0d0*I(ii)-315.0d0*I(ii-2)+105.0d0*I(ii-4)-5.0d0*I(ii-6));
     L12(ii) = 0.0625d0*(231.0d0*J(ii)-315.0d0*J(ii-2)+105.0d0*J(ii-4)-5.0d0*J(ii-6));
     L13(ii) = 0.0625d0*(231.0d0*J(ii+1)-315.0d0*J(ii-1)+105.0d0*J(ii-3)-5.0d0*J(ii-5));
     L14(ii) = 0.0625d0*(231.0d0*J(ii+2)-315.0d0*J(ii)+105.0d0*J(ii-2)-5.0d0*J(ii-4));
     L22(ii) = 0.0625d0*(231.0d0*S(ii)-315.0d0*S(ii-2)+105.0d0*S(ii-4)-5.0d0*S(ii-6));
     L23(ii) = 0.0625d0*(231.0d0*S(ii+1)-315.0d0*S(ii-1)+105.0d0*S(ii-3)-5.0d0*S(ii-5));
     L24(ii) = 0.0625d0*(231.0d0*S(ii+2)-315.0d0*S(ii)+105.0d0*S(ii-2)-5.0d0*S(ii-4)); 
  ELSEIF (ii==8) THEN
     L11(ii) = 0.0625d0*(429.0d0*I(ii)-693.0d0*I(ii-2)+315.0d0*I(ii-4)-35.0d0*I(ii-6));
     L12(ii) = 0.0625d0*(429.0d0*J(ii)-693.0d0*J(ii-2)+315.0d0*J(ii-4)-35.0d0*J(ii-6));
     L13(ii) = 0.0625d0*(429.0d0*J(ii+1)-693.0d0*J(ii-1)+315.0d0*J(ii-3)-35.0d0*J(ii-5));
     L14(ii) = 0.0625d0*(429.0d0*J(ii+2)-693.0d0*J(ii)+315.0d0*J(ii-2)-35.0d0*J(ii-4));
     L22(ii) = 0.0625d0*(429.0d0*S(ii)-693.0d0*S(ii-2)+315.0d0*S(ii-4)-35.0d0*S(ii-6));
     L23(ii) = 0.0625d0*(429.0d0*S(ii+1)-693.0d0*S(ii-1)+315.0d0*S(ii-3)-35.0d0*S(ii-5));
     L24(ii) = 0.0625d0*(429.0d0*S(ii+2)-693.0d0*S(ii)+315.0d0*S(ii-2)-35.0d0*S(ii-4));
  ELSEIF (ii==9) THEN
     L11(ii) = 1.0d0/128.0d0*(6435.0d0*I(ii)-12012.0d0*I(ii-2)+6930.0d0*I(ii-4)-1260.0d0*I(ii-6)+35.0d0*I(ii-8));
     L12(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii)-12012.0d0*J(ii-2)+6930.0d0*J(ii-4)-1260.0d0*J(ii-6)+35.0d0*J(ii-8));
     L13(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii+1)-12012.0d0*J(ii-1)+6930.0d0*J(ii-3)-1260.0d0*J(ii-5)+35.0d0*J(ii-7));
     L14(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii+2)-12012.0d0*J(ii)+6930.0d0*J(ii-2)-1260.0d0*J(ii-4)+35.0d0*J(ii-6));
     L22(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii)-12012.0d0*S(ii-2)+6930.0d0*S(ii-4)-1260.0d0*S(ii-6)+35.0d0*S(ii-8));
     L23(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii+1)-12012.0d0*S(ii-1)+6930.0d0*S(ii-3)-1260.0d0*S(ii-5)+35.0d0*S(ii-7));
     L24(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii+2)-12012.0d0*S(ii)+6930.0d0*S(ii-2)-1260.0d0*S(ii-4)+35.0d0*S(ii-6));

     ELSEIF (ii>9) THEN 
     PRINT*,"K  too large"
  END IF
END DO

END SUBROUTINE Analytical_intII


SUBROUTINE Analytical_intIII(xb,pb,xbar,N,L11,L12,L13,L14,L22,L23,L24)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::N
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,DIMENSION(3)::R_0
  !REAL*8,DIMENSION(N+3)::Analytical_int,I,J,S
  REAL*8,DIMENSION(N+3)::I,J,S
  REAL*8,DIMENSION(N+1),INTENT(OUT)::L11,L12,L13,L14,L22,L23,L24
  !!REAL*8::Analytical_int
  REAL*8::b,c,se,sb,ue,ub,te,tb,d,nn
  INTEGER::ii
  
  R_0=(xbar-xb);
  b=-2.0d0*SUM(pb*R_0);
  c=SUM(R_0**2.0d0);
  se=1.0d0;
  sb=-1.0d0;
  ue=sqrt(1.0d0+b+c);
  ub=sqrt(1.0d0-b+c);
  te=1.0d0+0.5d0*b;
  tb=-1.0d0+0.5d0*b;
  d=c-0.25d0*b*b;
  
  
  I(1) = LOG(ABS(2.0d0*se+b+2.0d0*ue))-LOG(ABS(2.0d0*sb + b + 2.0d0*ub));
  I(2) = ue-ub-b/2.0d0*I(1);
  
  nn=2.0d0; 
  DO ii=3,N+3
     nn=nn+1.0d0;
     !I(ii) = -se**(ii-2)*ue/(1-ii) + sb**(ii-2)*ub/(1-ii) - (0.5d0-(ii-1))*b/(1-ii)*I(ii-1) + (ii-2)*c/(1-ii)*I(ii-2);
     
     I(ii) = -se**(ii-2)*ue/(1.0d0-nn) + sb**(ii-2)*ub/(1.0d0-nn) &
          - (0.5d0-(nn-1.0d0))*b/(1.0d0-nn)*I(ii-1) &
          + (nn-2.0d0)*c/(1.0d0-nn)*I(ii-2);
  END DO
!!$  DO ii=1,N+3
!!$     PRINT*,"I_new(",ii,") = ", I(ii)
!!$  END DO
!!$  DO ii=1,N+3
!!$     PRINT*,"I", I(ii)
!!$  END DO
  
  IF (d<1e-14) THEN
     
     J(1) = -1.0d0/(2.0d0*(se+b/2.0d0)**2.0d0) &
          + 1.0d0/(2.0d0*(sb+b/2.0d0)**2.0d0);
     
  ELSE
     J(1) = te/(d*ue)-tb/(d*ub);
     
  END IF
  J(2) =  -1.0d0/ue+1.0d0/ub - b/2.0d0*J(1);
 
  DO ii=3,N+3
     J(ii) = I(ii-2) - b*J(ii-1)-c*J(ii-2);
     
  END DO
!!$  DO ii=1,N+3
!!$     PRINT*,"J", J(ii)
!!$  END DO
  
  
  
  IF (d<1e-14) THEN 
     S(1) = -4.0d0/(2.0d0*se+b)**4.0d0 + 4.0d0/(2.0d0*sb+b)**4.0d0;
     S(2) = -8.0d0/(3.0d0*(2.0d0*se+b)**3.0d0) &
     +8.0d0/(3.0d0*(2.0d0*sb+b)**3.0d0)-b/2.0d0*S(1);
  
  ELSE
     S(1) = (2.0d0*se+b)/(6.0d0*d*ue**3.0d0)&
          -(2.0d0*sb+b)/(6.0d0*d*ub**3.0d0)+2.0d0/(3.0d0*d)*J(1);
     S(2) = -(b*se+2.0d0*c)/(6.0d0*d*ue**3.0d0)&
          +(b*sb+2.0d0*c)/(6.0d0*d*ub**3.0d0)-b/(3.0d0*d)*J(1);
  END IF

DO ii=3,N+3
   S(ii) = J(ii-2) - b*S(ii-1)-c*S(ii-2);
END DO

!!$DO ii=1,N+3
!!$   PRINT *,"Inew(",ii,") = ", I(ii)
!!$END DO
!!$
!!$DO ii=1,N+3
!!$   PRINT *,"Jnew(",ii,") = ", J(ii)
!!$END DO
!!$DO ii=1,N+3
!!$   PRINT *,"Snew(",ii,") = ", S(ii)
!!$END DO


DO ii=1,N+1
  IF (ii==1 .OR. ii==2) THEN
     L11(ii) = I(ii);
     L12(ii) = J(ii);
     L13(ii) = J(ii+1);
     L14(ii) = J(ii+2);
     L22(ii) = S(ii);
     L23(ii) = S(ii+1);
     L24(ii) = S(ii+2);
    
  ELSEIF (ii==3) THEN
     L11(ii) = 0.5*(3.0d0*I(ii)-I(ii-2));
     L12(ii) = 0.5*(3.0d0*J(ii)-J(ii-2));
     L13(ii) = 0.5*(3.0d0*J(ii+1)-J(ii-1));
     L14(ii) = 0.5*(3.0d0*J(ii+2)-J(ii));
     L22(ii) = 0.5*(3.0d0*S(ii)-S(ii-2));
     L23(ii) = 0.5*(3.0d0*S(ii+1)-S(ii-1));
     L24(ii) = 0.5*(3.0d0*S(ii+2)-S(ii));
     
  ELSEIF (ii==4) THEN
     L11(ii) = 0.5d0*(5.0d0*I(ii)-3.0d0*I(ii-2));
     L12(ii) = 0.5d0*(5.0d0*J(ii)-3.0d0*J(ii-2));
     L13(ii) = 0.5d0*(5.0d0*J(ii+1)-3.0d0*J(ii-1));
     L14(ii) = 0.5d0*(5.0d0*J(ii+2)-3.0d0*J(ii));
     L22(ii) = 0.5d0*(5.0d0*S(ii)-3.0d0*S(ii-2));
     L23(ii) = 0.5d0*(5.0d0*S(ii+1)-3.0d0*S(ii-1));
     L24(ii) = 0.5d0*(5.0d0*S(ii+2)-3.0d0*S(ii));
     
  ELSEIF (ii==5) THEN
    L11(ii) = 0.125d0*(35.0d0*I(ii)-30.0d0*I(ii-2)+3.0d0*I(ii-4));
    L12(ii) = 0.125d0*(35.0d0*J(ii)-30.0d0*J(ii-2)+3.0d0*J(ii-4));
    L13(ii) = 0.125d0*(35.0d0*J(ii+1)-30.0d0*J(ii-1)+3.0d0*J(ii-3));
    L14(ii) = 0.125d0*(35.0d0*J(ii+2)-30.0d0*J(ii)+3.0d0*J(ii-2));
    L22(ii) = 0.125d0*(35.0d0*S(ii)-30.0d0*S(ii-2)+3.0d0*S(ii-4));
    L23(ii) = 0.125d0*(35.0d0*S(ii+1)-30.0d0*S(ii-1)+3.0d0*S(ii-3));
    L24(ii) = 0.125d0*(35.0d0*S(ii+2)-30.0d0*S(ii)+3.0d0*S(ii-2));
  ELSEIF (ii==6) THEN
     L11(ii) = 0.125d0*(63.0d0*I(ii)-70.0d0*I(ii-2)+15.0d0*I(ii-4));
     L12(ii) = 0.125d0*(63.0d0*J(ii)-70.0d0*J(ii-2)+15.0d0*J(ii-4));
     L13(ii) = 0.125d0*(63.0d0*J(ii+1)-70.0d0*J(ii-1)+15.0d0*J(ii-3));
     L14(ii) = 0.125d0*(63.0d0*J(ii+2)-70.0d0*J(ii)+15.0d0*J(ii-2));
     L22(ii) = 0.125d0*(63.0d0*S(ii)-70.0d0*S(ii-2)+15.0d0*S(ii-4));
     L23(ii) = 0.125d0*(63.0d0*S(ii+1)-70.0d0*S(ii-1)+15.0d0*S(ii-3));
     L24(ii) = 0.125d0*(63.0d0*S(ii+2)-70.0d0*S(ii)+15.0d0*S(ii-2));
  ELSEIF (ii==7) THEN
     L11(ii) = 0.0625d0*(231.0d0*I(ii)-315.0d0*I(ii-2)&
          +105.0d0*I(ii-4)-5.0d0*I(ii-6));
     L12(ii) = 0.0625d0*(231.0d0*J(ii)-315.0d0*J(ii-2)&
          +105.0d0*J(ii-4)-5.0d0*J(ii-6));
     L13(ii) = 0.0625d0*(231.0d0*J(ii+1)-315.0d0*J(ii-1)&
          +105.0d0*J(ii-3)-5.0d0*J(ii-5));
     L14(ii) = 0.0625d0*(231.0d0*J(ii+2)-315.0d0*J(ii)&
          +105.0d0*J(ii-2)-5.0d0*J(ii-4));
     L22(ii) = 0.0625d0*(231.0d0*S(ii)-315.0d0*S(ii-2)&
          +105.0d0*S(ii-4)-5.0d0*S(ii-6));
     L23(ii) = 0.0625d0*(231.0d0*S(ii+1)-315.0d0*S(ii-1)&
          +105.0d0*S(ii-3)-5.0d0*S(ii-5));
     L24(ii) = 0.0625d0*(231.0d0*S(ii+2)-315.0d0*S(ii)&
          +105.0d0*S(ii-2)-5.0d0*S(ii-4)); 
  ELSEIF (ii==8) THEN
     L11(ii) = 0.0625d0*(429.0d0*I(ii)-693.0d0*I(ii-2)&
     +315.0d0*I(ii-4)-35.0d0*I(ii-6));
     L12(ii) = 0.0625d0*(429.0d0*J(ii)-693.0d0*J(ii-2)&
          +315.0d0*J(ii-4)-35.0d0*J(ii-6));
     L13(ii) = 0.0625d0*(429.0d0*J(ii+1)-693.0d0*J(ii-1)&
          +315.0d0*J(ii-3)-35.0d0*J(ii-5));
     L14(ii) = 0.0625d0*(429.0d0*J(ii+2)-693.0d0*J(ii)&
          +315.0d0*J(ii-2)-35.0d0*J(ii-4));
     L22(ii) = 0.0625d0*(429.0d0*S(ii)-693.0d0*S(ii-2)&
          +315.0d0*S(ii-4)-35.0d0*S(ii-6));
     L23(ii) = 0.0625d0*(429.0d0*S(ii+1)-693.0d0*S(ii-1)&
          +315.0d0*S(ii-3)-35.0d0*S(ii-5));
     L24(ii) = 0.0625d0*(429.0d0*S(ii+2)-693.0d0*S(ii)&
          +315.0d0*S(ii-2)-35.0d0*S(ii-4));
  ELSEIF (ii==9) THEN
     L11(ii) = 1.0d0/128.0d0*(6435.0d0*I(ii)-12012.0d0*I(ii-2)&
          +6930.0d0*I(ii-4)-1260.0d0*I(ii-6)+35.0d0*I(ii-8));
     L12(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii)-12012.0d0*J(ii-2)&
          +6930.0d0*J(ii-4)-1260.0d0*J(ii-6)+35.0d0*J(ii-8));
     L13(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii+1)-12012.0d0*J(ii-1)&
          +6930.0d0*J(ii-3)-1260.0d0*J(ii-5)+35.0d0*J(ii-7));
     L14(ii) = 1.0d0/128.0d0*(6435.0d0*J(ii+2)-12012.0d0*J(ii)&
          +6930.0d0*J(ii-2)-1260.0d0*J(ii-4)+35.0d0*J(ii-6));
     L22(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii)-12012.0d0*S(ii-2)&
          +6930.0d0*S(ii-4)-1260.0d0*S(ii-6)+35.0d0*S(ii-8));
     L23(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii+1)-12012.0d0*S(ii-1)&
          +6930.0d0*S(ii-3)-1260.0d0*S(ii-5)+35.0d0*S(ii-7));
     L24(ii) = 1.0d0/128.0d0*(6435.0d0*S(ii+2)-12012.0d0*S(ii)&
          +6930.0d0*S(ii-2)-1260.0d0*S(ii-4)+35.0d0*S(ii-6));

     ELSEIF (ii>9) THEN 
     PRINT*,"K  too large"
  END IF
END DO

END SUBROUTINE Analytical_intIII





FUNCTION G_compute_GQ_kg(xb,pb,xbar,eeps,N);

  !!Contr. from filament b to point xbar. 
  !!G is integral over filament b, with kernel 
  !!multiplied by L_m(s).
  !!Result are 6 values stored in Gvec:
  !!G11,G22,G33,G12,G13,G23.
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,INTENT(IN)::eeps
  REAL*8,DIMENSION(3)::R_0

  INTEGER, INTENT(IN):: N
  REAL*8,DIMENSION(N+1)::L11,L12,L13,L14,L22,L23,L24
  !!REAL*8,DIMENSION(N+1)::L11n,L12n,L13n,L14n,L22n,L23n,L24n
  REAL*8,DIMENSION(6):: G_compute_GQ_kg
  REAL*8,DIMENSION(3,3):: R_00,Identity,R_0_pb, pb_R_0,Gvec_tmp,pb_pb
  INTEGER i,j
  !Gvec=zeros(6,1);
  R_0=xbar-xb;
  DO i=1,3
     DO j=1,3
        R_00(i,j) = R_0(i)*R_0(j);
        pb_R_0(i,j) = pb(i)*R_0(j);
        R_0_pb(j,i) = R_0(j)*pb(i);
        pb_pb(i,j)=pb(i)*pb(j);
     END DO
  END DO
  Identity=RESHAPE((/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/), (/3, 3/))

  
  CALL Analytical_int(xb,pb,xbar,N,L11,L12,L13,L14,L22,L23,L24);
  
  
  Gvec_tmp =L11(N+1)*Identity + R_00*L12(N+1) - (pb_R_0+R_0_pb)*L13(N+1)+pb_pb*L14(N+1) &
  +2.0d0*eeps**2*(Identity*L12(N+1)-3.0d0*R_00*L22(N+1)+3.0d0*(pb_R_0+R_0_pb)*L23(N+1)-3.0d0*pb_pb*L24(N+1));
  
 
  
 
  G_compute_GQ_kg(1)=Gvec_tmp(1,1);
  G_compute_GQ_kg(2)=Gvec_tmp(2,2);
  G_compute_GQ_kg(3)=Gvec_tmp(3,3);
  G_compute_GQ_kg(4)=Gvec_tmp(1,2);
  G_compute_GQ_kg(5)=Gvec_tmp(1,3);
  G_compute_GQ_kg(6)=Gvec_tmp(2,3);

  

END FUNCTION G_compute_GQ_kg




FUNCTION G_compute_GQ_f_kg(M,N,xb,pb,xbar,eeps,ExtForce, coeffvec,fno);
  !!Contr. from filament b to point xbar. 
  !!G is integral over filament b, with kernel 
  !!multiplied by L_m(s).
  !!Result are 6 values stored in Gvec:
  !!G11,G22,G33,G12,G13,G23.
  INTEGER, INTENT(IN):: M,N,fno
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,INTENT(IN)::eeps
  REAL*8,DIMENSION(3)::R_0
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::coeffvec
  REAL*8,DIMENSION(3*M),INTENT(IN)::ExtForce
  
  REAL*8,DIMENSION(N+1)::L11,L12,L13,L14,L22,L23,L24
  REAL*8,DIMENSION(6,N+1):: Gvec
  REAL*8,DIMENSION(3,3):: R_00,Identity,R_0_pb, pb_R_0,Gvec_tmp,pb_pb
  REAL*8,DIMENSION(3)::G_compute_GQ_f_kg
  INTEGER i,j,ind
  
  R_0=xbar-xb;
  DO i=1,3
     DO j=1,3
        R_00(i,j) = R_0(i)*R_0(j);
        pb_R_0(i,j) = pb(i)*R_0(j);
        R_0_pb(j,i) = R_0(j)*pb(i);
        pb_pb(i,j)=pb(i)*pb(j);
     END DO
  END DO
  Identity=RESHAPE((/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/), (/3, 3/))
  
  CALL Analytical_int(xb,pb,xbar,N,L11,L12,L13,L14,L22,L23,L24);
  
  

  
  
  DO i=1,N+1
     Gvec_tmp =L11(i)*Identity + R_00*L12(i) - (pb_R_0+R_0_pb)*L13(i)+pb_pb*L14(i) &
    +2.0d0*eeps**2*(Identity*L12(i)-3*R_00*L22(i)+3.0d0*(pb_R_0+R_0_pb)*L23(i)-3.0d0*pb_pb*L24(i));
    
     Gvec(1,i)=Gvec_tmp(1,1);
     Gvec(2,i)=Gvec_tmp(2,2);
     Gvec(3,i)=Gvec_tmp(3,3);
     Gvec(4,i)=Gvec_tmp(1,2);
     Gvec(5,i)=Gvec_tmp(1,3);
     Gvec(6,i)=Gvec_tmp(2,3);
     
    

  END DO
  


ind=(fno-1)*3
  
  
G_compute_GQ_f_kg(1) = 0.5d0*ExtForce(ind+1)*Gvec(1,1) + 0.5d0*ExtForce(ind+2)*Gvec(4,1)&
     + 0.5d0*ExtForce(ind+3)*Gvec(5,1); 
G_compute_GQ_f_kg(2) = 0.5d0*ExtForce(ind+1)*Gvec(4,1) + 0.5d0*ExtForce(ind+2)*Gvec(2,1)&
     + 0.5d0*ExtForce(ind+3)*Gvec(6,1);
G_compute_GQ_f_kg(3) = 0.5d0*ExtForce(ind+1)*Gvec(5,1) + 0.5d0*ExtForce(ind+2)*Gvec(6,1)&
     + 0.5d0*ExtForce(ind+3)*Gvec(3,1);
 


DO i=2,N+1
   ind=(fno-1)*3*N+(i-2)*3;
   
   G_compute_GQ_f_kg(1) = G_compute_GQ_f_kg(1) + coeffvec(ind+1)*Gvec(1,i) + coeffvec(ind+2)*Gvec(4,i) + coeffvec(ind+3)*Gvec(5,i);
   
  G_compute_GQ_f_kg(2) = G_compute_GQ_f_kg(2) + coeffvec(ind+1)*Gvec(4,i) +coeffvec(ind+2)*Gvec(2,i) + coeffvec(ind+3)*Gvec(6,i);
   
   G_compute_GQ_f_kg(3) = G_compute_GQ_f_kg(3) + coeffvec(ind+1)*Gvec(5,i) + coeffvec(ind+2)*Gvec(6,i) + coeffvec(ind+3)*Gvec(3,i);
   
END DO



END FUNCTION G_compute_GQ_f_kg




FUNCTION G_f_GQ_kg(N,xb,pb,xbar,fvec_x,fvec_y,fvec_z,eeps);
  !!Assuming all vectors in are column vectors. 
  !!Result is column vector as well. 
  
  !!Contr. from filament b to point xbar. 
  !!G is integral over filament b, with kernel 
  !!multiplied fvec.
  !!pv,wv are quad pts and weights.
  !!Result are 3 values stored in Gvec:
  !!KF_x,KF_y,KF_z.
  INTEGER, INTENT(IN)::N
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,INTENT(IN)::eeps
  REAL*8,DIMENSION(3)::R_0
  REAL*8,INTENT(IN)::fvec_x,fvec_y,fvec_z
  REAL*8,DIMENSION(N+1)::L11,L12,L13,L14,L22,L23,L24
  REAL*8,DIMENSION(3)::G_f_GQ_kg
  REAL*8,DIMENSION(3,3):: R_00,Identity,R_0_pb, pb_R_0,Gvec_tmp,pb_pb
  
  INTEGER i,j,ind
  
  
  R_0=(xbar-xb);
  
  DO i=1,3
     DO j=1,3
        R_00(i,j) = R_0(i)*R_0(j);
        pb_R_0(i,j) = pb(i)*R_0(j);
        R_0_pb(j,i) = R_0(j)*pb(i);
        pb_pb(i,j)=pb(i)*pb(j);
     END DO
  END DO
  Identity=RESHAPE((/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/), (/3, 3/))
  
  CALL Analytical_int(xb,pb,xbar,0,L11,L12,L13,L14,L22,L23,L24);
  
  
 
   Gvec_tmp =L11(1)*Identity + R_00*L12(1) - (pb_R_0+R_0_pb)*L13(1)+pb_pb*L14(1) &
 + 2.0d0*eeps**2*(Identity*L12(1)-3.0d0*R_00*L22(1)+3.0d0*(pb_R_0+ R_0_pb)*L23(1)-3.0d0*pb_pb*L24(1));
   
   G_f_GQ_kg(1)=Gvec_tmp(1,1)*fvec_x+Gvec_tmp(1,2)*fvec_y+Gvec_tmp(1,3)*fvec_z; 
   G_f_GQ_kg(2)=Gvec_tmp(1,2)*fvec_x+Gvec_tmp(2,2)*fvec_y+Gvec_tmp(2,3)*fvec_z;
   G_f_GQ_kg(3)=Gvec_tmp(1,3)*fvec_x+Gvec_tmp(2,3)*fvec_y+Gvec_tmp(3,3)*fvec_z;



  
  END FUNCTION G_f_GQ_kg 



FUNCTION TH_f_GQ_kg(M,N,xb,pb,xa,pa,eeps,ExtForce,coeffvec,LQ,pv,wv,Lvec,fno);

  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xa,pa
  REAL*8,INTENT(IN)::eeps
  REAL*8,DIMENSION(3*M*N),INTENT(IN)::coeffvec
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv,Lvec
  INTEGER, INTENT(IN):: M,N,LQ,fno
  INTEGER i
  REAL*8,DIMENSION(3):: TH_f_GQ_kg
  REAL*8,DIMENSION(3*M),INTENT(IN)::ExtForce
  REAL*8,DIMENSION(3)::xbar
  REAL*8,DIMENSION(LQ,3)::Gmat
  !!Contr. from filament b to filament a. 
  !!First G is integral over filament b, with kernel 
  !!multiplied by fvec.
  !!This yield result that depends on what point of filament a
  !!we are at. Next step is to integrate this result multiplied 
  !!by L_k(s) over filament a to get a final result. 
  !!Result is 3 components in thvec, TH_f_GQ

  !!Quad pts in pv, weights in wv. 
    
  !DO i=1,LQ
  DO i=1,LQ
     xbar=xa+pv(i)*pa;  
     Gmat(i,:)=G_compute_GQ_f_kg(M,N,xb,pb,xbar,eeps,ExtForce,coeffvec,fno);
     
  END DO
  
  DO i=1,3
     TH_f_GQ_kg(i)=sum(wv*Gmat(:,i)*Lvec);
  END DO
END FUNCTION TH_f_GQ_kg






END MODULE Int_analytical