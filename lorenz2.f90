PROGRAM lorenz2
REAL*8::a,b,h,f,g,u,sigma,beta,delta,m,ord,sumat,sumad,cuad,prod,rho
REAL*8::x(0:7000),y(0:7000),z(0:7000),t(0:7000),k0(0:7000),k1(0:7000)
REAL*8::k2(0:7000),k3(0:7000),l0(0:7000),l1(0:7000),l2(0:7000),l3(0:7000),x1(0:7000)
REAL*8::m0(0:7000),m1(0:7000),m2(0:7000),m3(0:7000),x0(0:7000),y0(0:7000),z0(0:7000)
REAL*8::d(0:7000),rect(0:7000),func(0:7000)
INTEGER::k,n

!We define our parameters and initial conditions

h=0.01
t(0)=0.0
sigma=10.0
beta=8.0/3.0
delta=0.000000000001


WRITE(*,*) "Insert the interval [a,b]:"
READ(*,*) a,b

WRITE(*,*) "Insert the initial conditions for x, y and z:"
READ(*,*) x0(0),y0(0),z0(0)

WRITE(*,*) "Inserte the value for Rayleigh's number:"
READ(*,*) rho

n=INT((b-a)/h)

x(0)=x0(0)+delta
y(0)=y0(0)+delta
z(0)=z0(0)+delta

!We perform the 4th Order Runge Kutta's Method on a perturbation of the initial conditions

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

!We obtain the Lyapunov's coefficient using a Least Squares fit on a log-log scale

OPEN(1,FILE="x.txt",status="old")
DO k=0,n
READ(1,*) x1(k)
END DO
CLOSE(1,STATUS="DELETE")


DO k=0,n
d(k)=ABS(x1(k)-x(k))
func(K)=LOG(d(k))
END DO

sumat=0.0
sumad=0.0
cuad=0.0
prod=0.0

DO k=1500,4500
sumat=sumat+t(k)
sumad=sumad+func(k)
cuad=cuad+(t(k))**2
prod=prod+(t(k))*(func(k))
END DO

m=(((3000)*prod)-(sumat*sumad))/(((3000)*cuad)-(sumat**2))
ord=(sumad-(sumat*m))/(3000)

WRITE(*,*) "Lyapunov's coefficient is: ",m

!We obtain the data of the displacement between the trajectories given a perturbation and it's linear fit in a log-log scale
DO k=0,n
rect(k)=m*t(k)+ord
END DO

OPEN(2,FILE="semilogdelta.txt",STATUS="new")
OPEN(3,FILE="recta.txt",STATUS="new")
DO k=0,n
WRITE(2,*) t(k),LOG(d(k))
WRITE(3,*) t(k),rect(k)
END DO
CLOSE(2)
CLOSE(3)

!We obtain the time series of the x, y and z trajectories and the 2D and 3D coordinates of the phase space


OPEN(4,FILE="atractor3Dd.txt",STATUS="new")
OPEN(5,FILE="atractor2Dd.txt",STATUS="new")
OPEN(6,FILE="x(t)d.txt",STATUS="new")
OPEN(7,FILE="y(t)d.txt",STATUS="new")
OPEN(8,FILE="z(t)d.txt",STATUS="new")
DO k=0,n
WRITE(4,*) x(k),y(k),z(k)
WRITE(5,*) x(k),z(k)
WRITE(6,*) t(k),x(k)
WRITE(7,*) t(k),y(k)
WRITE(8,*) t(k),z(k)
END DO
CLOSE(4)
CLOSE(5)
CLOSE(6)
CLOSE(7)
CLOSE(8)


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
