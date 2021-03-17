PROGRAM lorenz3
REAL*8::a,b,h,f,g,u,sigma,beta
REAL*8::x(0:100000),y(0:100000),z(0:100000),t(0:100000),k0(0:100000),k1(0:100000)
REAL*8::k2(0:100000),k3(0:100000),l0(0:100000),l1(0:100000),l2(0:100000),l3(0:100000)
REAL*8::m0(0:100000),m1(0:100000),m2(0:100000),m3(0:100000),rho(0:1000)
INTEGER::k,n,i,l,j

!We define our parameters and initial conditions

h=0.001
t(0)=0.0
sigma=10.0
beta=8.0/3.0
a=0.0
b=100.0
x(0)=0.0
y(0)=1.0
z(0)=0.0
rho(0)=0

n=INT((b-a)/h)
l=INT((250)/0.25)

!We perform the 4th order Runge-Kutta method for every value of rho and save the bifurcation data

OPEN(1,FILE="bifurcacionesx.txt",STATUS="new")
OPEN(2,FILE="bifurcacionesy.txt",STATUS="new")
OPEN(3,FILE="bifurcacionesz.txt",STATUS="new")
DO i=0,l
DO k=0,n
k0(k)=h*f(x(k),y(k),sigma)
l0(k)=h*g(x(k),y(k),z(k),rho(i))
m0(k)=h*u(x(k),y(k),z(k),beta)
k1(k)=h*f(x(k)+(0.5)*k0(k),y(k)+(0.5)*l0(k),sigma)
l1(k)=h*g(x(k)+(0.5)*k0(k),y(k)+(0.5)*l0(k),z(k)+(0.5)*m0(k),rho(i))
m1(k)=h*u(x(k)+(0.5)*k0(k),y(k)+(0.5)*l0(k),z(k)+(0.5)*m0(k),beta)
k2(k)=h*f(x(k)+(0.5)*k1(k),y(k)+(0.5)*l1(k),sigma)
l2(k)=h*g(x(k)+(0.5)*k1(k),y(k)+(0.5)*l1(k),z(k)+(0.5)*m1(k),rho(i))
m2(k)=h*u(x(k)+(0.5)*k1(k),y(k)+(0.5)*l1(k),z(k)+(0.5)*m1(k),beta)
k3(k)=h*f(x(k)+k2(k),y(k)+l2(k),sigma)
l3(k)=h*g(x(k)+k2(k),y(k)+l2(k),z(k)+m2(k),rho(i))
m3(k)=h*u(x(k)+k2(k),y(k)+l2(k),z(k)+m2(k),beta)
x(k+1)=x(k)+(k0(k)+2*k1(k)+2*k2(k)+k3(k))/(6.0)
y(k+1)=y(k)+(l0(k)+2*l1(k)+2*l2(k)+l3(k))/(6.0)
z(k+1)=z(k)+(m0(k)+2*m1(k)+2*m2(k)+m3(k))/(6.0)
t(k+1)=t(k)+h
END DO
DO j=1,n
IF (x(j).le.x(j+1) .AND. x(j).le.x(j-1)) THEN
WRITE(1,*) rho(i),x(j)
ELSE IF (x(j).ge.x(j+1) .AND. x(j).ge.x(j-1)) THEN
WRITE(1,*) rho(i),x(j)
END IF

IF (y(j).le.y(j+1) .AND. y(j).le.y(j-1)) THEN
WRITE(2,*) rho(i),y(j)
ELSE IF (y(j).ge.y(j+1) .AND. y(j).ge.y(j-1)) THEN
WRITE(2,*) rho(i),y(j)
END IF

IF (z(j).le.z(j+1) .AND. z(j).le.z(j-1)) THEN
WRITE(3,*) rho(i),z(j)
ELSE IF (z(j).ge.z(j+1) .AND. z(j).ge.z(j-1)) THEN
WRITE(3,*) rho(i),z(j)
END IF
END DO

rho(i+1)=rho(i)+0.25
x(0)=x(n+1)
y(0)=y(n+1)
z(0)=z(n+1)
END DO
CLOSE(1)
CLOSE(2)
CLOSE(3)
END PROGRAM

!We define the Lorenz System



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



