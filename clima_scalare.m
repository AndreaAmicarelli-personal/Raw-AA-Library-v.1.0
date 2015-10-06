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
% Program unit name: clima_scalare
% Program unit description: modello climatologico scalare (Cushman_97), 
% con studio di sensibilità sul coefficiente di assorbimento 
% dell'atmosfera alle onde lunghe (effetto serra) e verifica delle 3 
% condizioni di equilibrio al variare della temperatura. Year: 2008.
%-----------------------------------------------------------------------
% dati
aAs=0.18;
IT=344;
tAs=0.49;
aTs=0.96;
Q=113;
sigma0=5.67e-8;
eps_max=1;
aTs_g=0.42;
aAl=zeros(11);
ET=zeros(11);
EA=zeros(11);
TT=zeros(11);
TA=zeros(11);

% analisi di sensibilità per quantificazione effetto serra; la soluzione
% principale è la (6), con aAl=0.95 (equilibrio attuale)
for i=1:11,
    aAl(i)=0.89+0.01*i
    ET(i)=((0.64*aAs+tAs*aTs)*IT-0.36*Q)/(1-0.64*aAl(i));
    EA(i)=aAl(i)*ET(i)+aAs*IT+Q;
    TT(i)=((ET(i)^0.25)/(sigma0^0.25))-273.16;
    TA(i)=((EA(i)^0.25)/((aAl(i)*sigma0)^0.25))-273.16;    
end;

% equilibrio ere glaciali 
ET_g=((0.64*aAs+tAs*aTs_g)*IT-0.36*Q)/(1-0.64*aAl(6));
EA_g=aAl(6)*ET_g+aAs*IT+Q;
TT_g=((ET_g^0.25)/(sigma0^0.25))-273.16;
TA_g=((EA_g^0.25)/((aAl(6)*sigma0)^0.25))-273.16;

% equilibrio instabile
a=((aTs-aTs_g)*tAs*IT)/(30*(sigma0^0.25)*(1-0.64*aAl(6)));
b=(0.64*aAs*IT+tAs*(aTs_g-(aTs-aTs_g)*253/30)*IT-0.36*Q)/(1-0.64*aAl(6));
ET_old=300;
ET_new=300;
eps=1000;
j=0;
while ((eps>eps_max)&(j<100)),
    ET_old=(ET_new+ET_old)/2;
    j=j+1;
    ET_new=(-(ET_old^0.25)*a-b+ET_old)/(0.25*a*(ET_old^(-0.75))-1)+ET_old;
    eps=(ET_new^0.25)*a+b-ET_new;
end;
if j>99,
    d=input(' errore: dopo 99 iterazioni la condizione per eps non è verificata ');
   else 
      ET_i=ET_new;
      TT_i=((ET_i^0.25)/(sigma0^0.25))-273.16;
      EA_i=aAl(6)*ET_i+aAs*IT+Q;
      TA_i=((EA_i^0.25)/((aAl(6)*sigma0)^0.25))-273.16;
end;

% grafici effetto serra
plot(aAl,ET,'k',aAl,EA,'g',aAl,TT,'r',aAl,TA,'b');
xlabel('aAl')
title('effetto serra')
legend('ET','EA','TT','TA')
colorbar;
save ('effetto_serra.out');
end;
end;
end;
