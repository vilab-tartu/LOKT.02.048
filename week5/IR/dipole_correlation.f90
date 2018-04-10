!
!
IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(0.0D0)
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: correlation
REAL(KIND=dp) :: integral,omega,Pi,timestep
REAL(KIND=dp), DIMENSION(3,1000000) :: dipder
INTEGER :: N,I,J,Nmax,Fmax
CHARACTER(LEN=100) :: line,filename

Pi=4.0D0*ATAN(1.0D0)

READ(5,*) filename
READ(5,*) timestep

OPEN(10,FILE=filename)
N=0
DO 
 READ(10,'(A100)',END=999) line
 IF (INDEX(line,' DIPOLE [Non Periodic] DERIVATIVE(A.U.)|').NE.0) THEN
    N=N+1
    READ(line(45:),*) dipder(:,N)
 ENDIF
ENDDO
999 CONTINUE

CLOSE(10)

! use only the first 10% of the data for the correlation function, as the rest is not statistically meaningful
! and our quadratic algorithm becomes too slow
Nmax=N/10

print *, Nmax

ALLOCATE(correlation(0:Nmax))
correlation=0.0_dp

DO I=1,N-Nmax
 DO J=I,I+Nmax
    correlation(J-I)=correlation(J-I)+DOT_PRODUCT(dipder(:,I),dipder(:,J))
 ENDDO
ENDDO

DO I=0,Nmax
   correlation(I)=correlation(I)/(REAL(N-I,kind=dp)*REAL(N,kind=dp))
ENDDO

OPEN(UNIT=10,FILE="dip_dip_correlation.time")
write(10,*) "# correlation in the time domain, first column time in fs"
DO I=-Nmax,Nmax
   write(10,*) I*timestep,correlation(ABS(I))/correlation(0)
ENDDO
CLOSE(10)

OPEN(UNIT=10,FILE="dip_dip_correlation.freq")
write(10,*) "# correlation in the frequency domain [cm^-1]"
!Fmax up to 4000 cm^-1
Fmax=N*((4000D0/(2*Pi))*(2.0D0*Pi*timestep*1.0D-15*29979245800.0))
DO I=0,Fmax
   omega=(2.0D0*Pi*I)/N
   integral=0.0_dp
   DO J=0,Nmax
      integral=integral+cos(omega*J)*correlation(J)
   ENDDO
   ! frequency in cm^-1
   write(10,*) omega/(2.0D0*Pi*timestep*1.0D-15*29979245800.0),integral
ENDDO
CLOSE(10)

END
