#-----------------------------------------------------------------------
# "Raw AA Library v.1.0" is a C/Matlab draft library, which may be   
# useful for simple computations involving calculus subroutines, 
# 1D heath transport (Finite Differences), 1D free surface flows (FD), 
# 1D pollutant dispersion (FD), 0D and 1D climate models (FD).
# Raw AA Library v.1.0 Copyright 2000,2003,2004,2008 Andrea Amicarelli
# email contact: Andrea.Amicarelli@gmail.com
#-----------------------------------------------------------------------
# This file is part of Raw AA Library v.1.0. 
# Raw AA Library v.1.0 is free software: you can redistribute it and/or 
# modify it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation, either version 3 of the 
# License, or (at your option) any later version.
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public 
# License along with this library. If not, see 
# <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
# Program unit name: Newton
# Program unit description: Newton’s bisection method for the function 
# 2x^3-4x+1. Year: 2000.
#-----------------------------------------------------------------------
#include<stdio.h>
#include<math.h>


double f(double x)
{double a, b, m, ERR;
double fa, fb, fm;

do{printf("A: ") ;
scanf("%f", &a);
printf("B: ");
scanf("%f", &b);
fa=2*a*a*a-4*a+1;
fb=2*b*b*b-4*b+1;
}
while(fa*fb>0);
do{
m=(a+b)/2;
fm=2*m*m*m-4*m+1;
if(fm!=0)
{if(fa*fm<0)
b=m;
else a=m;
}}
while(fabs(fm)>ERR);
printf("Newton dice: %f", m);
return  2*x*x*x-4*x+1      ;
}
