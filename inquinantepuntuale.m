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
% Program unit name: inquinantepuntuale
% Program unit description: 1D Solution to the point source pollutant 
% diffusion via Thomas’ algorithm. Year: 2004.
%-----------------------------------------------------------------------
l=input('lunghezza:');
K=input('inverso della diffusività turbolenta:');
D=input('scrivi 1:');
r=input('calore generato per unità di lunghezza per unità di tempo:');
C=input('scrivi 1:');
A=input('concentrazione di picco iniziale:');
t=input('numero di passi temporali (escluso il primo):');
z=input('indice spaziale di immissione: ');
h=0.15;
%metodo di Thomas, vettori per la fattorizzazione LU della matrice tridiagonale
c=l/h-1;
d=zeros(c);
for i=1:c,
    d(i)=2+2*K;
end;
s=zeros(c-1);
for i=1:(c-1),
    s(i)=-1;
end;
a=zeros(c-1);
for i=1:(c-1);
    a(i)=-1;
end;
alfa=zeros(c-1);
u=zeros(c);
v=zeros(c-1);
u(1)=d(1);
for i=1:(c-1),
    v(i)=s(i);
end;
for i=1:(c-1),
    alfa(i)=a(i)/u(i);
    u(i+1)=d(i+1)-alfa(i)*v(i);
end;
%metodo di Thomas, risultati, le T sono le x, prima iterazione
b=zeros(c);
b=zeros(c);
b(z-1)=A+K;
b(z)=A+K;
b(z+1)=A;
%for i=1:c,
 %   b(i)=A*sin((pi)*(i-1)*h/l)+(2*K-2)*A*sin((pi)*i*h/l)+A*sin(pi*(i+1)*h/l)+K*r/D/C;
 %end;
x=zeros(c+2,t+1);
%condizioni iniziali
for i=1:(c+1),
    x(i,1)=0;
end;
x(z,1)=A;
%for i=1:(c+1),
 %   x(i,1)=A*sin((pi)*(i*h-h)/l);
 %end;
%condizioni al contorno
for j=2:t+1,
    x(1,j)=0;
    x(c+2,j)=0;
end;
%risultati I iterazione
y=zeros(c);
y(1)=b(1);
for i=2:c,
    y(i)=b(i)-alfa(i-1)*y(i-1);
end;
x(c+1,2)=y(c)/u(c);
for i=c:-1:2,
    x(i,2)=(y(i-1)-v(i-1)*x(i+1,2))/u(i-1);
end;
%metodo di Thomas, risultati, iterazioni successive
for j=3:t+1,
    for i=1:c,
        b(i)=x(i,j-1)+(2*K-2)*x(i+1,j-1)+x(i+2,j-1)+K*r/D/C;
    end;
    y(1)=b(1);
    for i=2:c,
        y(i)=b(i)-alfa(i-1)*y(i-1);
    end;
    x(c+1,j)=y(c)/u(c);
    for i=c:-1:2,
        x(i,j)=(y(i-1)-v(i-1)*x(i+1,j))/u(i-1);
    end;
end
end;
end

