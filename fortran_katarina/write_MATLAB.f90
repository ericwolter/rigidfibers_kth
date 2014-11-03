MODULE write_MATLAB
IMPLICIT NONE
  
CONTAINS
  
!!============================================================
  SUBROUTINE int2str(str,num)
    !use constants
    implicit none
    
    CHARACTER(len=10) :: str
    INTEGER :: num,nstp,ipos,ndigits,idigit
    str = ''
    nstp = ABS(num)
    IF (num == 0) THEN
       str = CHAR(ICHAR('0'))
    ELSEIF (num <= -1) THEN
       ndigits = FLOOR(LOG10(REAL(ABS(num))))+2
       DO ipos = ndigits,2,-1
          idigit = MOD(nstp,10)
          str(ipos:ipos) = CHAR(ICHAR('0')+idigit)
          nstp = nstp/10
       END DO
       str(1:1) = '-'
    ELSEIF (num >= 1) THEN
       ndigits = FLOOR(LOG10(REAL(ABS(num))))+1
       DO ipos = ndigits,1,-1
          idigit = MOD(nstp,10)
          str(ipos:ipos) = CHAR(ICHAR('0')+idigit)
          nstp = nstp/10
       END DO
    END IF
  END SUBROUTINE int2str
  
  SUBROUTINE save_array_to_file(filename,arr)
    REAL*8, DIMENSION(:),INTENT(IN)::arr
    CHARACTER (LEN=20)::filename
    
    INTEGER j,M
    M=size(arr);
    
    OPEN(10,file=trim(filename));
    DO j=1,M
       WRITE(10,'(2F24.16)') arr(j)
    END DO
    CLOSE(10)
    
  END SUBROUTINE save_array_to_file
  
  SUBROUTINE save2D_array_to_file(filename,arr)
    REAL*8, DIMENSION(:,:),INTENT(IN)::arr
    CHARACTER (LEN=20)::filename
    
    INTEGER i,j,M1,M2
    M1=size(arr,1);
    M2=size(arr,2);
    
    OPEN(10,file=trim(filename));
    DO j=1,M2
       DO i=1,M1
          WRITE(10,'(2F24.16)') arr(i,j)
       END DO
    END DO
    CLOSE(10)
    
  END SUBROUTINE save2D_array_to_file
    
  SUBROUTINE write_Rep_file(repfile,abbr,label,M,N,NoQI,eeps,&
                            dt,int_an,no_ts,save_ival,no_saves_in_file)

    CHARACTER (len=25),INTENT(IN)::repfile,abbr
    INTEGER,INTENT(IN)::label,M,N,no_ts,save_ival,no_saves_in_file,NoQI,int_an
    REAL*8,INTENT(IN):: eeps,dt
    PRINT*, "int_an= ",int_an
    PRINT*, "NoQI ",NoQI
    OPEN(10,file=trim(repfile));
    WRITE (10,'(A40)')"========================================"
    WRITE (10,'(A30)')repfile
    WRITE (10,*) "Run with label ",label
    WRITE (10,'(A40)')"========================================"

    WRITE (10,'(A40)')"========================================"
    WRITE (10,'(A25,A15,A4)')"Results written to files ",abbr,&
         "+pX"
    WRITE (10,'(A30)')"Parameters for the run: "
    WRITE (10,'(A20,I6)') "No of fibers:",M
    WRITE (10,'(A33,I4)') "No of terms in force expansion:",N
    WRITE (10,'(A30,I4)') "No of quad. intervals:",NoQI
    WRITE (10,'(A10,E10.4)') "eeps:",eeps
    WRITE (10,'(A10,F10.4)') "dt:",dt
    WRITE (10,'(A10,I6)') "no_ts:",no_ts
    WRITE (10,'(A33,I4)') "Use analytical integration (1/0):", int_an
    WRITE (10,'(A10,I4)') "save_ival:",save_ival
    WRITE (10,'(A26,I4)') "No of saves in each file:",no_saves_in_file
    WRITE (10,'(A40)')"========================================"

    CLOSE(10)
    
  END SUBROUTINE write_Rep_file
!!=================================================================
  SUBROUTINE write_Pars(parsfile,label,M,N,eeps,&
                        dt,no_ts,save_ival,no_saves_in_file)
    
    CHARACTER (len=25),INTENT(IN)::parsfile
    INTEGER,INTENT(IN)::label,M,N,no_ts,save_ival,no_saves_in_file
    REAL*8,INTENT(IN):: eeps,dt

    CHARACTER (len=10)::fstr
    CHARACTER (len=10)::labstr
    CALL int2str(labstr,label);
    
    OPEN(10,file=trim(parsfile));

    WRITE (10,'(A8,A2,I4,A2)') "M"//trim(labstr),"=",M,";"
    WRITE (10,'(A8,A2,I4,A2)') "N"//trim(labstr),"=",N,";"
    WRITE (10,'(A8,A2,F20.16,A3)') "eeps"//trim(labstr),"=",eeps,";"
    WRITE (10,'(A8,A2,F20.16,A3)') "dt"//trim(labstr),"=",dt,";"
    WRITE (10,'(A8,A2,I6,A2)') "no_ts"//trim(labstr),"=",no_ts,";"
    WRITE (10,'(A12,A2,I4,A2)') "save_ival"//trim(labstr),"=",save_ival,";"
    WRITE (10,'(A12,A2,I4,A2)') "nos_if"//trim(labstr),"=",no_saves_in_file,";"
    
  END SUBROUTINE write_Pars
!!=================================================================




!!================================================================
END MODULE write_MATLAB

