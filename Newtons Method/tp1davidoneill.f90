!-- Module containing some real constants --!
module constants
implicit none
real :: g, R, S, P, Tc, T, l, n, m, omega_sqrd, gam, theta0
end module

!-- Module acting as a function interface --!
module IFunc
implicit none
real :: x, f, df
end module

!-- Program entry point --!
program
implicit none
external :: SinFunc, SinApprox, ExpFunc, TanhFunc, Landau
real :: Newton
write(*,*) "Newton's Method : Zero Finder"
write(*,*)
write(*,*)
write(*,*) "6a) Finding zero for function sin(x) with x0 = 0.1 :"
Newton(SinFunc, 0.1, .0001)
write(*,*) "6b.1) Finding zero for function sin_approx(x) with x0 = 0.1 :"
Newton(SinApprox, 0.1, .0001)
write(*,*) "	  Finding zero for function sin(x) with x0 = 1.55 :"
Newton(SinFunc, 1.55, .0001)
write(*,*) "6c) Finding zero for function exp(x) with x0 = 0.1 :"
Newton(ExpFunc, 0.1, .0001)
write(*,*) "6d) Finding zero for function TanhFunc(x) with x0 = 0.1 :"
Newton(TanhFunc, 1.0, .0001)
write(*,*) "7) Finding zero for function Landau(x) with x0 = 0.1 :"
Newton(Landau, 0.1, .0001)
end program

!-- Iteratively estimates the root of a supplied function. Max iteration = 10,000 --!
real function Newton(func, x0, eps)
implicit none
integer :: i, iMax
real :: root, f, df, xNext, xCurrent, epsCurrent, x0, eps
logical :: converges
i = 0
iMax = 10000
f = 0
df = 0
xNext = 0
root = 0
xCurrent = x0
epsCurrent = 100
do while(abs(epsCurrent) > eps .and. i <= iMax)
	call func(xCurrent, f, df)
	if(f .ne. df) then
		xNext = xCurrent - f / df
		epsCurrent = xNext - xCurrentroot 
		root = xNext
		xCurrent = xNext
	endif
	else if(f <= eps) then
		cycle
	endif
	i = i + 1
enddo
converges = i .le. iMax
if(converges) then
	write(*,*) root
else
	write(*,*) "No root found - does not converge!"
endif
Newton = root
end function

!-- A subroutine that calculates sin(x) and sin'(x) [cos(x)] --!
subroutine SinFunc(x, f, df)
use IFunc
implicit none
f = sin(x)
df = cos(x)
end subroutine

!-- A subroutine that calculates sin(x) and estimates sin'(x) --!
subroutine SinApprox(x, f, df)
use IFunc
implicit none
step = 1e-4
f = sin(x)
df = (sin(x+step/2.0) - sin(x-step/2.0)) / step
end subroutine

!-- A subroutine that calculates exp(x) exp'(x) [exp(x)]--!
subroutine ExpFunc(x, f, df)
use IFunc
implicit none
f = exp(x)
df = exp(x)
end subroutine

!-- A subroutine that calculates (1/2) - tanh(x-1) and f'(x) --!
subroutine TanhFunc(x, f, df)
use IFunc
implicit none
f = (1.0 / 2.0) - tanh(x - 1)
df = -1.0 * (1/(cosh(1-x)))**2.0
end subroutine

!-- A subroutine that calculates f(x) and f'(x) for question 7 --!
subroutine Landau(x, f, df)
use IFunc
use constants
implicit none
g = 9.81
R = 8.32
S = 1e-4
P = 200.0
T = 300.0
l=0.1
n = 8e-7
Tc = (m * g * l * theta0**2.0) / (2.0 * n * R)
m = 6e-3
omega_sqrd = g / l
gam = 2.0 * n * R / (m * l**2.0)
theta0 = 1.0
f = (sin(x) / x) * (theta0**2.0 - x**2.0) - (gam * Tc / (2.0 * omega_sqrd))
df = (theta0**2.0 - x**2.0) * ((x * cos(x) - sin(x)) / x**2.0) - 2.0 * sin(x)
end subroutine