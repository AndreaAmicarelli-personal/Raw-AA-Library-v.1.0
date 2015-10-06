%-----------------------------------------------------------------------
% "Raw AA Library" is a C++/Matlab draft library, which may be useful  
% for simple computations involving calculus subroutines, 
% 1D heath transport (Finite Differences), 1D free surface flows (FD), 
% 1D pollutant dispersion (FD), 0D and 1D climate models (FD).
% Raw AA Library Copyright 2000,2003,2004,2008 Andrea Amicarelli
% email contact: Andrea.Amicarelli@gmail.com
%-----------------------------------------------------------------------
% This file is part of Raw AA Library. 
% Raw AA Library is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% You should have received a copy of the GNU Lesser General Public 
% License along with this library. If not, see 
% <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------
% Program unit name: North_75
% Program unit description: modello di North_75. (non ho utilizato il 
% controllo ricorsivo sull'albedo, ma qui l'errore è minimo ed equivale 
% a dire che la T_lim per il cambio di albedo è -12°C invece di -10°C 
% (non ho avuto tempo per l'iterazione). Year: 2008.
%-----------------------------------------------------------------------
% dati
c=90;
pigreco=3.14159;
phi=zeros(c+1);
for i=1:(c+1),
   phi(i)=(i-1)*(90/c);
end;
x_N=zeros(c+1);
for i=1:(c+1),
   x_N(i)=sin(phi(i)*(pigreco)/180);
end;
dphi=(90/c)*(pigreco)/180;
dx=zeros(c);
for i=1:c,
   dx(i)=(x_N(i+1))-(x_N(i));
end;
alb=zeros(c);
alb_min=0.32;
alb_max=0.62;
i_alb=ceil(72/90*c)+1;
for i=1:i_alb,
   alb(i)=alb_min;
end;
for i=(i_alb+1):c,
   alb(i)=alb_max;
end;
I_sun=334.4;
S=zeros(c);
for i=1:c,
   S(i)=1-0.482*((3*((x_N(i)+x_N(i+1))/2)^2-1)/2);
end;
% A=201.4;
A=211.2;
D=0.31;
%B=1.45;
B=1.55;

%metodo di fattorizzazione LU di Thomas (per sistemi lineari con matrice dei coefficienti tridiagonale:
% definisco le costanti del sistema lineare 
d(1)=1;
for i=2:(c-1),
    d(i)=-(((1-x_N(i)^2)/(dx(i)^2))+(B/D)+((1-x_N(i+1)^2)/(dx(i)^2)));
end;
d(c)=-(((1-x_N(c)^2)/(dx(c)^2))+(B/D));
s=zeros(c-1);
s(1)=-1
for i=2:(c-1),
    s(i)=(1-x_N(i+1)^2)/(dx(i)^2);
end;
a=zeros(c-1);
for i=1:(c-1);
    a(i)=(1-x_N(i+1)^2)/(dx(i+1)^2);
end;
b=zeros(c);
b(1)=0;
for i=2:(c-1),
    b(i)=-S(i)*I_sun*(1-alb(i))/D+A/D;
end;
b(c)=-S(c)*I_sun*(1-alb(c))/D+A/D;

% algoritmo di Thomas: definisco le costanti delle matrici L e U, fattori
% della matrice M tridiagonale dei cofficienti del sistema lineare
% (fattorizzazione LU)
u=zeros(c);
v=zeros(c-1);
u(1)=d(1);
for i=1:(c-1),
    v(i)=s(i);
    alfa(i)=a(i)/u(i);
    u(i+1)=d(i+1)-alfa(i)*v(i);
end;
v(c)=s(c);
alfa(c)=a(c)/u(c);

% algoritmo di Thomas: trovo il vettore Y=UX
x=zeros(c);
y=zeros(c);
y(1)=b(1);
for i=2:c,
    y(i)=b(i)-alfa(i-1)*y(i-1);
end;

% algoritmo di Thomas: trovo il vettore X(temperatura): MX=B
x(c)=y(c)/u(c);
for i=(c-1):-1:1,
    x(i)=(y(i)-v(i)*x(i+1))/u(i);
end;
end;
end;
