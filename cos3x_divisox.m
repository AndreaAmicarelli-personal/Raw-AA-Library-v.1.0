%-----------------------------------------------------------------------
% "Raw AA Library v.1.0" is a C/Matlab draft library, which may be   
% useful for simple computations involving calculus subroutines, 
% 1D heath transport (Finite Differences), 1D free surface flows (FD), 
% 1D pollutant dispersion (FD), 0D and 1D climate models (FD).
% Raw AA Library v.1.0 Copyright 2000,2003,2004,2008 Andrea Amicarelli
% email contact: Andrea.Amicarelli@gmail.com
%-----------------------------------------------------------------------
% This file is part of Raw AA Library v.1.0. 
% Raw AA Library v.1.0 is free software: you can redistribute it and/or 
% modify it under the terms of the GNU Lesser General Public License as 
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
% Program unit name: cos3x_divisox
% Program unit description: Computation of the interpolating polynome 
% for cos(3x)/x Language: C++. Year: 2003.
%-----------------------------------------------------------------------
nnodi=input('Inserisci il numero di nodi:');
%intervallo di interpolazione
a=1;
b=5;
%calcolo dei nodi
xnodi=linspace(a,b,nnodi);
fnodi=cos(3*xnodi)./xnodi;
%calcolo del polinomio nodale
x=linspace(a,b);
n=lenght(x);
for 1=1:n
    pnod(i)=prod(x(i)-xnodi);
end;

%calcolo dei polinomi di base di Lagrange
plag=zeros(n,nnodi);
for k=1:nnodi,
    for i:1:n
        plag(i,k)=prod(x(i)-xnodi(1:k-1))*prod(x(i)-xnodi(k+1:nnodi))/...
            prod(xnodi(k)-xnodi(1:k-1))*prod((xnodi(k)-xnodi(k+1:nnodi)));
    end;
end;
 %calcolo del polinomio interpolatore e grafici
 pn=fnodi*plag;
 f=cos(3*x)./x:plot(x,f,'k',x,pn,'r',xnodi,fnodi,'k*');
 pause;%tra un grafico e l'altro
 plot(x,f-pn,'b');
end.
