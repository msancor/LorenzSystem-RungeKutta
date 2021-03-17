PROGRAM lorenz
REAL*8::a,b,h,f,g,u,sigma,beta,rho
REAL*8::x(0:7000),y(0:7000),z(0:7000),t(0:7000),k0(0:7000),k1(0:7000)
REAL*8::k2(0:7000),k3(0:7000),l0(0:7000),l1(0:7000),l2(0:7000),l3(0:7000)
REAL*8::m0(0:7000),m1(0:7000),m2(0:7000),m3(0:7000),m(0:7000)
INTEGER::k,n

!We define our parameters and initial conditions

h=0.01
t(0)=0.0
sigma=10.0
beta=8.0/3.0


WRITE(*,*) "Insert the interval [a,b]"
READ(*,*) a,b

WRITE(*,*) "Insert the initial conditions for x, y and z:"
READ(*,*) x(0),y(0),z(0)

WRITE(*,*) "Insert the value for Rayleigh's number:"
READ(*,*) rho

n=INT((b-a)/h)

!We perform the 4th Order Runge-Kutta's method

DO k=0,n
k0(k)=h*f(x(k),y(k),sigma)
l0(k)=h*g(x(k),y(k),z(k),rho)
m0(k)=h*u(x(k),y(k),z(k),beta)
k1(k)=h*f(x(k)+(0.5)*k0(k),y(k)+(0.5)*l0(k),sigma)
l1(k)=h*g(x(k)+(0.5)*k0(k),y(k)+(0.5)*l0(k),z(k)+(0.5)*m0(k),rho)
m1(k)=h*u(x(k)+(0.5)*k0(k),y(k)+(0.5)*l0(k),z(k)+(0.5)*m0(k),beta)
k2(k)=h*f(x(k)+(0.5)*k1(k),y(k)+(0.5)*l1(k),sigma)
l2(k)=h*g(x(k)+(0.5)*k1(k),y(k)+(0.5)*l1(k),z(k)+(0.5)*m1(k),rho)
m2(k)=h*u(x(k)+(0.5)*k1(k),y(k)+(0.5)*l1(k),z(k)+(0.5)*m1(k),beta)
k3(k)=h*f(x(k)+k2(k),y(k)+l2(k),sigma)
l3(k)=h*g(x(k)+k2(k),y(k)+l2(k),z(k)+m2(k),rho)
m3(k)=h*u(x(k)+k2(k),y(k)+l2(k),z(k)+m2(k),beta)
x(k+1)=x(k)+(k0(k)+2*k1(k)+2*k2(k)+k3(k))/(6.0)
y(k+1)=y(k)+(l0(k)+2*l1(k)+2*l2(k)+l3(k))/(6.0)
z(k+1)=z(k)+(m0(k)+2*m1(k)+2*m2(k)+m3(k))/(6.0)
t(k+1)=t(k)+h
END DO

!We create and delete a datafile to save the maximum values of the trajectories

OPEN(1,FILE="maximos.txt",STATUS="new")
DO j=1,n
IF (z(j).ge.z(j+1) .AND. z(j).ge.z(j-1)) THEN
WRITE(1,*) z(j)
END IF
END DO
CLOSE(1)

OPEN(2,FILE="maximos.txt",status="old")
DO k=1,97
READ(2,*) m(k)
END DO
CLOSE(2,STATUS="DELETE")

!We save the data to graph the Lorenz map

OPEN(3,FILE="map.txt",STATUS="new")
DO j=1,73
WRITE(3,*) m(j),m(j+1)
END DO
CLOSE(3)

!We save the time series of the x, y and z trajectories and the coordinates of the 2D and 3D phase space

OPEN(4,FILE="atractor3D.txt",STATUS="new")
OPEN(5,FILE="atractor2D.txt",STATUS="new")
OPEN(6,FILE="x(t).txt",STATUS="new")
OPEN(7,FILE="y(t).txt",STATUS="new")
OPEN(8,FILE="z(t).txt",STATUS="new")
OPEN(9,FILE="x.txt",STATUS="new")
DO k=0,n
WRITE(4,*) x(k),y(k),z(k)
WRITE(5,*) x(k),z(k)
WRITE(6,*) t(k),x(k)
WRITE(7,*) t(k),y(k)
WRITE(8,*) t(k),z(k)
WRITE(9,*) x(k)
END DO
CLOSE(4)
CLOSE(5)
CLOSE(6)
CLOSE(7)
CLOSE(8)
CLOSE(9)







END PROGRAM


!We define the Lorenz equations

DOUBLE PRECISION FUNCTION f(x,y,sigma)
REAL*8::x,y,sigma
f=sigma*(y-x)
END FUNCTION

DOUBLE PRECISION FUNCTION g(x,y,z,rho)
REAL*8::x,y,z,rho
g=x*(rho-z)-y
END FUNCTION

DOUBLE PRECISION FUNCTION u(x,y,z,beta)
REAL*8::x,y,z,beta
u=x*y-beta*z
END FUNCTION


