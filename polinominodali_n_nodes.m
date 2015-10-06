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
% Program unit name: polinominodali_n_nodes
% Program unit description: Computation of the nodal polynomials 
% (Chebychev). polinomi nodali in[-1,1] secondo nodi equispaziati e 
% di Chebychev. Year: 2004.
%-----------------------------------------------------------------------
n=input('inserisci n: ');
inodi=0:n;
xeq=-1+2*inodi/n;
xCh=cos((2*inodi+1)*pi/((n+1)*2));
x=linspace(-1,1);
for i=1:100,
    polnodeq(i)=prod(x(i)-xeq);
end
for i=1:100,
    polnodCh=prod(x(i)-xCh);
end
plot(x,polnodeq,'r',x,polnodCh,'b',x,zeros(1,100),'k');
end.
