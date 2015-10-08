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
% Program unit name: canaleveloceagricolo
% Program unit description: metodo delle caratteristiche per il moto 
% vario nei corsi d'acqua a pelo libero con corrente veloce, sezione 
% uniforme (a scelta tra triangolare, rettangolare e trapezoidale), 
% regolato da un serbatoio di monte e da una portata a monte che varia 
% tra 0 e Q(tfi) linearmente e poi rimane costante, con flussi esterni 
% uniformi, con canale inizialmente asciutto. It is equivalent to  
% "canalelentoindustriale" with supercritical free surface flows and 
% channel initially dry. Thomas’ algorithm. Year: 2004.
%-----------------------------------------------------------------------
sez=input('Scrivi 1 per triangolare, 2 per rettangolare, 3 per trapezoidale: ');
if sez==1,
   a=input('pendenza pareti: ');
   elseif sez==2,
          b=input('base: ');
      else 
          b=input('base: ');
          a=input('pendenza pareti: ');
end;
L=input('lunghezza alveo: ');
ds=input('discretizzazione spaziale: ');
ts=input('tempo della simulazione: ');
dt=input('discretizzazione temporale: ');
n=input('coefficiente di Manning: ');
tf=input('tempo di fine aumento (multiplo di dt): ');
tq=input('tempo di inizio prelievi (multiplo di dt): ');
qag=input('prelievi per unità di lunghezza (segno meno): ');
ps=L/ds+1;
pt=ts/dt+1;
Y=zeros(ps,pt);
U=zeros(ps,pt);
g=9.81;
ifo=zeros(ps);
teta=dt/ds;
si=zeros(1,ps);
for i=1:ps
    si(i)=(i-1)*ds;
end;
ti=zeros(pt);
for i=1:pt,
    ti(i)=(i-1)*dt;
end;
Yaus=zeros(1);
%condizioni iniziali e pendenza del fondo
Yuni=input('Y uniforme: ');
for i=1:ps,
    Y(i,1)=Yuni;
end;
ifo(1)=input('pendenza del fondo: ');
for i=1:ps,
    ifo(i)=ifo(1);
end;
% condizioni al contorno
q=zeros(pt);
for i=1:(tq/dt+1),
    q(i)=0;
end;
for i=(tq/dt+1):pt,
    q(i)=qag;
end;
% calcolo le celerità relative, le aree, le larghezze del pelo libero
% ed i contorni bagnati in ogni sezione per le condizioni iniziali
c=zeros(ps,pt);
A=zeros(ps,pt);
B=zeros(ps,pt);
for j=1:ps,
    if sez==1,
       A(j,1)=(Y(j,1)^2)/a;
       B(j,1)=2*Y(j,1)/a;
       C(j,1)=2*(B(j,1)/2/cos(atan(a)));
       elseif sez==2,
              A(j,1)=b*Y(j,1);
              B(j,1)=b;
              C(j,1)=b+2*Y(j,1);
              else 
                  A(j,1)=b*Y(j,1)+(Y(j,1)^2)/a;
                  B(j,1)=2*A(j,1)/Y(j,1)-b;
                  C(j,1)=b+2*((B(j,1)-b)/2/cos(atan(a)));
    end;
    c(j,1)=sqrt(g*Y(j,1));
    J(j,1)=(((U(j,1))^2)*(n^2))/((A(j,1)/C(j,1))^(4/3));
end;
for i=1:ps,
    U(i,1)=((A(i,1)/(C(i,1)))^(2/3))/n*sqrt(ifo(1));
    Quni=U(1,1)*A(1,1);
end;
Yser=Yuni+(U(1,1)^2)/2/g;
Q=zeros(1,pt);
tfi=tf/dt;
Q(tfi+1)=input('portata affluente a regime: ');
for i=1:(tfi),
    Q(i)=(i-1)/(tfi)*(Q(tfi+1)-Quni)+Quni;
end;
for i=(tfi+2):pt,
    Q(i)=Q(tfi+1);
end;
% iterazione del procedimento risolutivo
for i=2:pt,
    for j=1:ps,
        %calcolo le soluzioni "centrali" e a valle
        if j<ps,
        UR=zeros(1);
        US=zeros(1);
        cR=zeros(1);
        cS=zeros(1);
        YR=zeros(1);
        YS=zeros(1);
        JR=zeros(1);
        JS=zeros(1);
        AR=zeros(1);
        AS=zeros(1);
        UR=(U(j+1,i-1)+teta*(-U(j+1,i-1)*c(j,i-1)+c(j+1,i-1)*U(j,i-1)))/(1+teta*(U(j+1,i-1)-U(j,i-1)+c(j+1,i-1)-c(j,i-1)));
        cR=(c(j+1,i-1)-teta*UR*(c(j+1,i-1)-c(j,i-1)))/(1+teta*(c(j+1,i-1)-c(j,i-1)));
        YR=Y(j+1,i-1)-teta*(UR+cR)*(Y(j+1,i-1)-Y(j,i-1));
        JR=J(j+1,i-1)-teta*(UR+cR)*(J(j+1,i-1)-J(j,i-1));
        AR=A(j+1,i-1)-teta*(UR+cR)*(A(j+1,i-1)-A(j,i-1));
        US=(U(j+1,i-1)+teta*(U(j+1,i-1)*c(j,i-1)-c(j+1,i-1)*U(j,i-1)))/(1+teta*(U(j+1,i-1)-U(j,i-1)-c(j+1,i-1)+c(j,i-1)));
        cS=(c(j+1,i-1)-teta*US*(c(j+1,i-1)-c(j,i-1)))/(1-teta*(c(j+1,i-1)-c(j,i-1)));
        YS=Y(j+1,i-1)-teta*(US-cS)*(Y(j+1,i-1)-Y(j,i-1));
        JS=J(j+1,i-1)-teta*(US-cS)*(J(j+1,i-1)-J(j,i-1));
        AS=A(j+1,i-1)-teta*(US-cS)*(A(j+1,i-1)-A(j,i-1));
        Y(j+1,i)=(1/(cR+cS))*(YS*cR+YR*cS+cR*cS*((UR-US)/g-dt*(JR-JS)-(q(i)*dt/g)*((UR-cR)/AR-(US+cS)/AS)));
        U(j+1,i)=UR-g*((Y(j+1,i)-YR)/cR)-g*dt*(JR-ifo(j+1))-q(i)*dt/AR*(UR-cR);
        end;
    end;
    % soluzioni "al contorno" a monte
    AAm=zeros(1);
    BBm=zeros(1);
    Yaus=Y(1,i-1);
    if sez==1,
       AAm=5*(Yaus^4)-4*Yser*(Yaus^3);
       BBm=-(Yaus^5)+Yser*(Yaus^4)-(Q(i)^2)*(a^2)/2/g+Yaus*AAm;
         elseif sez==2,
                AAm=3*(Yaus^2)-2*Yaus*Yser;
                BBm=-(Yaus^3)+(Yaus^2)*Yser-(Q(i)^2)/(b^2)/2/g+Yaus*AAm;
                else 
                     AAm=3*(b^2)*(Yaus^2)+5*(Yaus^4)/(a^2)+8*b*(Yaus^3)/a-2*(b^2)*Yser*Yaus-4*(Yaus^3)*Yser/(a^2)-6*(Yaus^2)*b*Yser/a;
                     BBm=-((b^2)*(Yaus^2)+(Yaus^4)/(a^2)+2*b*(Yaus^3)/a)*(Yaus-Yser)-(Q(i)^2)/2/g+AAm*Yaus;
      end;
      Y(1,i)=BBm/AAm;
      diff=zeros(1);
      diff=abs(Y(1,i)-Yaus);
      while diff>0.0001,
            Yaus=Y(1,i);
            if sez==1,
               AAm=5*(Yaus^4)-4*Yser*(Yaus^3);
               BBm=-(Yaus^5)+Yser*(Yaus^4)-(Q(i)^2)*(a^2)/2/g+Yaus*AAm;
               elseif sez==2,
                      AAm=3*(Yaus^2)-2*Yaus*Yser;
                      BBm=-(Yaus^3)+(Yaus^2)*Yser-(Q(i)^2)/(b^2)/2/g+Yaus*AAm;
                else 
                      AAm=3*(b^2)*(Yaus^2)+5*(Yaus^4)/(a^2)+8*b*(Yaus^3)/a-2*(b^2)*Yser*Yaus-4*(Yaus^3)*Yser/(a^2)-6*(Yaus^2)*b*Yser/a;
                      BBm=-((b^2)*(Yaus^2)+(Yaus^4)/(a^2)+2*b*(Yaus^3)/a)*(Yaus-Yser)-(Q(i)^2)/2/g+AAm*Yaus;
                  end;
       Y(1,i)=BBm/AAm;
      diff=abs(Y(1,i)-Yaus);
      end;
    % ricavo aree, larghezze del pelo libero, contorni bagnati, celerità
    % relative, velocità di monte e cadenti energetiche per ogni sezione,
    % ad ogni iterazione
    for j=1:ps,    
        if sez==1,
           A(j,i)=(Y(j,i)^2)/a;
           B(j,i)=2*Y(j,i)/a;
           C(j,i)=2*(B(j,i)/2/cos(atan(a)));
           elseif sez==2,
                  A(j,i)=b*Y(j,i);
                  B(j,i)=b;
                  C(j,i)=b+2*Y(j,i);
                  else 
                       A(j,i)=b*Y(j,i)+(Y(j,i)^2)/a;
                       B(j,i)=2*A(j,i)/Y(j,i)-b;
                       C(j,i)=b+2*((B(j,i)-b)/2/cos(atan(a)));
        end;
     c(j,i)=sqrt(g*Y(j,i));
     end;
     U(1,i)=Q(i)/A(1,i);
     for j=1:ps,
         J(j,i)=(((U(j,i))^2)*(n^2))/((A(j,i)/C(j,i))^(4/3));
     end;
end;
% verifica della discretizzazione
errdis=0;
for i=1:ps,
    for j=1:pt,
        if ((errdis==0)&(dt>(ds/c(i,j)*sqrt(1-dt/2*(g*J(i,j)/U(i,j)+q(i)/A(i,j)))))),
           errdis=1;
           d=input(' errore: dt è troppo elevato (scrivi 1) ');
        end;
    end;
end;
% verifica del numero di Froud
errFr=0;
Fr=zeros(ps,pt);
for i=1:ps,
    for j=1:pt,
        Fr(i,j)=U(i,j)/sqrt(g*Y(i,j));
        if ((Fr<1)&(errFr==0)),
           errFr=1;
           d=input(' errore: la corrente è lenta in qualche sezione (scrivi 1) ');
        end;
    end;
end;
% trovo le quote del pelo libero e le velocità estreme
iYmax=1;
iYmin=1;
jYmax=1;
jYmin=1;
iUmax=1;
iUmin=1;
jUmax=1;
jUmin=1;
for i=1:ps,
    for j=1:pt,
        if Y(i,j)>Y(iYmax,jYmax),
           iYmax=i;
           jYmax=j;
       end;
       if  Y(i,j)<Y(iYmin,jYmin),
           iYmin=i;
           jYmin=j;
       end;
       if  U(i,j)>U(iUmax,jUmax),
           iUmax=i;
           jUmax=j;
       end;
       if  U(i,j)<U(iUmin,jUmin),
           iUmin=i;
           jUmin=j;
       end;
   end;
end;
% portate di monte e di valle
Qris=zeros(ps,2);
for i=1:pt,
    for j=1:ps
    Qris(i,1)=U(1,i)*A(1,i);
    Qris(i,2)=U(ps,i)*A(ps,i);
end;
end;
% altezze del pelo libero e del fondo sul fondo della sezione iniziale
fon=zeros(ps);
for i=1:ps,
    fon(i)=si(i)*ifo(1);
end;
h=zeros(ps,pt);
for i=1:ps,
    for j=1:pt,
        h(i,j)=Y(i,j)-fon(i);
    end;
end;
% portate di monte e di valle
Qris=zeros(ps,pt);
for i=1:pt,
    for j=1:ps,
            Qris(i,j)=U(j,i)*A(j,i);
    Qris(i,j)=U(j,i)*A(j,i);
end;
FrYmax=Fr(iYmax,jYmax);
FrYmin=Fr(iYmin,jYmin);
FrUmax=Fr(iUmax,jUmax);
FrUmin=Fr(iUmin,jUmin);
cYmax=c(iYmax,jYmax);
cYmin=c(iYmin,jYmin);
cUmax=c(iUmax,jUmax);
cUmin=c(iUmin,jUmin);
AYmax=A(iYmax,jYmax);
AYmin=A(iYmin,jYmin);
AUmax=A(iUmax,jUmax);
AUmin=A(iUmin,jUmin);
BYmax=B(iYmax,jYmax);
BYmin=B(iYmin,jYmin);
BUmax=B(iUmax,jUmax);
BUmin=B(iUmin,jUmin);
CYmax=C(iYmax,jYmax);
CYmin=C(iYmin,jYmin);
CUmax=C(iUmax,jUmax);
CUmin=C(iUmin,jUmin);
JYmax=J(iYmax,jYmax);
JYmin=J(iYmin,jYmin);
JUmax=J(iUmax,jUmax);
JUmin=J(iUmin,jUmin);
Ymax=Y(iYmax,jYmax);
Ymin=Y(iYmin,jYmin);
Umax=U(iUmax,jUmax);
Umin=U(iUmin,jUmin);
Yfin=zeros(pt,1);
UYmax=U(iYmax,jYmax);
UYmin=U(iYmin,jYmin);
YUmax=Y(iUmax,jUmax);
YUmin=Y(iUmin,jUmin);
for i=1:pt,
    Yfin(i,1)=Y(ps,i);
end;
end
