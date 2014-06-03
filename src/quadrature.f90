MODULE quadrature
IMPLICIT NONE

CONTAINS

SUBROUTINE quad_pts_and_wts(NoQI,pv,wv)
!!points and weights for 3rd order Gauss Quad with NoQI ivals 
!!on -1 to 1. 

  INTEGER, INTENT(IN):: NoQI
  REAL*8,DIMENSION(3*NoQI),INTENT(OUT)::pv,wv

  REAL*8,DIMENSION(3):: pv0,wv0
  REAL*8 a,iv
  INTEGER i,l0

  l0=3
  !pv0=[ 0.0 sqrt(15.0)/5.0]';
  pv0(1)=-sqrt(15.0d0)/5.0d0
  pv0(2)=0.0d0
  pv0(3)=sqrt(15.0d0)/5.0d0

  !!wv0=1.0/9.0*[5.0 8.0 5.0]';
  wv0(1)=5.0d0/9.0d0
  wv0(2)=8.0d0/9.0d0
  wv0(3)=5.0d0/9.0d0

  a=-1.0d0;
  iv=2.0d0/NoQI;
  
  DO i=1,NoQI
    pv((i-1)*l0+1:i*l0)=(2.0d0*a+iv+pv0*iv)/2.0d0;
    wv((i-1)*l0+1:i*l0)=wv0/NoQI;
    a=a+iv;
 END DO
!!$ DO i=1,(3*NoQI-1)/2
!!$    pv(i+(3*NoQI-1)/2+1)=-pv((3*NoQI-1)/2-i+1);
!!$ END DO
!!$ pv((3*NoQI-1)/2+1)=0.0d0;
!!$ DO i=1,3*NoQI
!!$    PRINT*,"pv",pv(i)
!!$ END DO
!!$ DO i=1,3*NoQI
!!$    PRINT*,"Wv",wv(i)
!!$ END DO
 
END SUBROUTINE quad_pts_and_wts



FUNCTION G_compute_GQ(xb,pb,xbar,eeps,LQ,pv,wv,Lvec_m);
  
  !!Contr. from filament b to point xbar. 
  !!G is integral over filament b, with kernel 
  !!multiplied by L_m(s).
  !!Result are 6 values stored in Gvec:
  !!G11,G22,G33,G12,G13,G23.
    
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,INTENT(IN)::eeps

  INTEGER, INTENT(IN):: LQ
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv,Lvec_m

  REAL*8,DIMENSION(6):: G_compute_GQ

  REAL*8,DIMENSION(LQ)::Rvec_x,Rvec_y,Rvec_z,Rmod
  REAL*8,DIMENSION(LQ)::K11,K22,K33,K12,K13,K23
  INTEGER i

  DO i=1,LQ
    Rvec_x(i)=xbar(1)-(xb(1)+pb(1)*pv(i));  
    Rvec_y(i)=xbar(2)-(xb(2)+pb(2)*pv(i));  
    Rvec_z(i)=xbar(3)-(xb(3)+pb(3)*pv(i));  
  END DO
  Rmod=sqrt(Rvec_x**2+Rvec_y**2+Rvec_z**2);
  
  !Rvec=xbar'*ones(size(pv))-(xb'*ones(size(pv))+pb'*pv);
  !Rmod=sqrt(sum(Rvec.^2));
  !K11=1.0./Rmod+Rvec(1,:).^2./Rmod.^3;
  !K22=1.0./Rmod+Rvec(2,:).^2./Rmod.^3;
  !K33=1.0./Rmod+Rvec(3,:).^2./Rmod.^3;
  !K12=Rvec(1,:).*Rvec(2,:)./Rmod.^3;
  !K13=Rvec(1,:).*Rvec(3,:)./Rmod.^3;
  !K23=Rvec(2,:).*Rvec(3,:)./Rmod.^3;

  K11=1.0d0/Rmod+Rvec_x**2/Rmod**3+&
             2.0d0*eeps**2*(1.0d0/Rmod**3-3.0d0*Rvec_x**2/Rmod**5);
  K22=1.0d0/Rmod+Rvec_y**2/Rmod**3+&
             2.0d0*eeps**2*(1.0d0/Rmod**3-3.0d0*Rvec_y**2/Rmod**5);
  K33=1.0d0/Rmod+Rvec_z**2/Rmod**3+&
             2.0d0*eeps**2*(1.0d0/Rmod**3-3.0d0*Rvec_z**2/Rmod**5);
  K12=Rvec_x*Rvec_y/Rmod**3-&
                    6.0d0*eeps**2*Rvec_x*Rvec_y/Rmod**5;
  K13=Rvec_x*Rvec_z/Rmod**3-&
                    6.0d0*eeps**2*Rvec_x*Rvec_z/Rmod**5;
  K23=Rvec_y*Rvec_z/Rmod**3-&
                    6.0d0*eeps**2*Rvec_y*Rvec_z/Rmod**5;
  
  G_compute_GQ(1)=sum(wv*K11*Lvec_m);
  G_compute_GQ(2)=sum(wv*K22*Lvec_m);
  G_compute_GQ(3)=sum(wv*K33*Lvec_m);
  G_compute_GQ(4)=sum(wv*K12*Lvec_m);
  G_compute_GQ(5)=sum(wv*K13*Lvec_m);
  G_compute_GQ(6)=sum(wv*K23*Lvec_m);
!!$  PRINT*,"G_compute_GQ(1)",G_compute_GQ(1)
!!$  PRINT*,"G_compute_GQ(2)",G_compute_GQ(2)
!!$  PRINT*,"G_compute_GQ(3)",G_compute_GQ(3)
!!$  PRINT*,"G_compute_GQ(4)",G_compute_GQ(4)
!!$  PRINT*,"G_compute_GQ(5)",G_compute_GQ(5)
!!$  PRINT*,"G_compute_GQ(6)",G_compute_GQ(6)
!!$STOP  

END FUNCTION G_COMPUTE_GQ

FUNCTION  G_f_GQ(xb,pb,xbar,eeps,LQ,fvec_x,fvec_y,fvec_z,pv,wv);
  !!Contr. from filament b to point xbar. 
  !!G is integral over filament b, with kernel 
  !!multiplied fvec.
  !!pv,wv are quad pts and weights.
  !!Result are 3 values stored in Gvec:
  !!KF_x,KF_y,KF_z.
    
  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xbar
  REAL*8,INTENT(IN)::eeps

  INTEGER, INTENT(IN):: LQ
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ),INTENT(IN)::fvec_x,fvec_y,fvec_z
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv

  INTEGER i
  REAL*8,DIMENSION(3):: G_f_GQ

  REAL*8,DIMENSION(LQ)::Rvec_x,Rvec_y,Rvec_z,Rmod
  REAL*8,DIMENSION(LQ)::K11,K22,K33,K12,K13,K23

  DO i=1,LQ
    Rvec_x(i)=xbar(1)-(xb(1)+pb(1)*pv(i));  
    Rvec_y(i)=xbar(2)-(xb(2)+pb(2)*pv(i));  
    Rvec_z(i)=xbar(3)-(xb(3)+pb(3)*pv(i));  
  END DO
  Rmod=sqrt(Rvec_x**2+Rvec_y**2+Rvec_z**2);

  !Rvec=xbar'*ones(size(pv))-(xb'*ones(size(pv))+pb'*pv);
  !Rmod=sqrt(sum(Rvec.^2));
  !K11=1.0d0./Rmod+Rvec(1,:).^2./Rmod.^3;
  !K22=1.0d0./Rmod+Rvec(2,:).^2./Rmod.^3;
  !K33=1.0d0./Rmod+Rvec(3,:).^2./Rmod.^3;
  !K12=Rvec(1,:).*Rvec(2,:)./Rmod.^3;
  !K13=Rvec(1,:).*Rvec(3,:)./Rmod.^3;
  !K23=Rvec(2,:).*Rvec(3,:)./Rmod.^3;

  K11=1.0d0/Rmod+Rvec_x**2/Rmod**3+&
             2.0d0*eeps**2*(1.0d0/Rmod**3-3.0d0*Rvec_x**2/Rmod**5);
  K22=1.0d0/Rmod+Rvec_y**2/Rmod**3+&
             2.0d0*eeps**2*(1.0d0/Rmod**3-3.0d0*Rvec_y**2/Rmod**5);
  K33=1.0d0/Rmod+Rvec_z**2/Rmod**3+&
             2.0d0*eeps**2*(1.0d0/Rmod**3-3.0d0*Rvec_z**2/Rmod**5);
  K12=Rvec_x*Rvec_y/Rmod**3-&
                    6.0d0*eeps**2*Rvec_x*Rvec_y/Rmod**5;
  K13=Rvec_x*Rvec_z/Rmod**3-&
                    6.0d0*eeps**2*Rvec_x*Rvec_z/Rmod**5;
  K23=Rvec_y*Rvec_z/Rmod**3-&
                    6.0d0*eeps**2*Rvec_y*Rvec_z/Rmod**5;

  G_f_GQ(1)=sum(wv*(K11*fvec_x+K12*fvec_y+K13*fvec_z));
  G_f_GQ(2)=sum(wv*(K12*fvec_x+K22*fvec_y+K23*fvec_z));
  G_f_GQ(3)=sum(wv*(K13*fvec_x+K23*fvec_y+K33*fvec_z));
  
END FUNCTION G_f_GQ

FUNCTION TH_f_GQ(xb,pb,xa,pa,eeps,fvec_x,fvec_y,fvec_z,LQ,pv,wv,Lvec);

  REAL*8,DIMENSION(3),INTENT(IN)::xb,pb,xa,pa
  REAL*8,INTENT(IN)::eeps

  INTEGER, INTENT(IN):: LQ
  !pv,wv,and Lvec_m should be of length LQ. 
  REAL*8,DIMENSION(LQ),INTENT(IN)::fvec_x,fvec_y,fvec_z
  REAL*8,DIMENSION(LQ),INTENT(IN)::pv,wv,Lvec

  INTEGER i
  REAL*8,DIMENSION(3):: TH_f_GQ

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
    
  DO i=1,LQ
     xbar=xa+pv(i)*pa;  
     Gmat(i,:)=G_f_GQ(xb,pb,xbar,eeps,LQ,fvec_x,fvec_y,fvec_z,pv,wv);
     
  END DO
  
  DO i=1,3
     TH_f_GQ(i)=sum(wv*Gmat(:,i)*Lvec);
  END DO
END FUNCTION TH_f_GQ

END MODULE quadrature

