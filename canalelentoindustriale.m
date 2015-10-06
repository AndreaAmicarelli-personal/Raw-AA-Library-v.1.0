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
% Program unit name: canalelentoindustriale
% Program unit description: metodo delle caratteristiche per il moto 
% vario nei corsi d'acqua a pelo libero con corrente lenta, sezione 
% uniforme (a scelta tra triangolare, rettangolare e trapezoidale), 
% regolato da un serbatoio di monte e da una portata a valle che cresce 
% da 0mc/s a Q(tfi) in tf secondi e poi rimane costante con condizioni 
% iniziali statiche; la pendenza del canale è uniforme, i flussi esterni
% costanti ed uniformi. In altre parole si tratta di canalelento con 
% portata a valle che cresce da 0mc/s a Q(tfi) in tf secondi e poi 
% rimane costante con condizioni iniziali statiche. Thomas’ algorithm.
% Year: 2004.
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
q=input('flussi esterni: ');
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
Y(1,1)=input('altezza di monte: ');
Yser=Y(1,1);
ifo(1)=input('pendenza del canale:');
for i=2:ps,
    ifo(i)=ifo(1);
end;
for i=2:ps
    Y(i,1)=Y(1,1)+ifo(1)*((i-1)*ds);
end;
for i=1:ps,
    U(i,1)=0;
end;
% condizioni al contorno
Q=zeros(1,pt);
tf=input('tempo di fine aumento (multiplo della discretizzazione temporale): ');
tfi=tf/dt;
Q(tfi+1)=input('portata a regime: ');
for i=1:(tfi),
    Q(i)=(i-1)/(tfi)*Q(tfi+1);
end;
for i=(tfi+2):pt,
    Q(i)=Q(tfi+1);
end;
% calcolo le celerità relative, le aree, i contorni bagnati, le larghezze
% del pelo libero e le cadenti energetiche in ogni sezione per le 
% condizioni iniziali
c=zeros(ps,pt);
A=zeros(ps,pt);
B=zeros(ps,pt);
J=zeros(ps,pt);
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
% iterazione del procedimento risolutivo
for i=2:pt;
    for j=1:ps;
        %calcolo le soluzioni "centrali"
        if j<ps-1,
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
         US=(U(j+1,i-1)-teta*(U(j+1,i-1)*c(j+2,i-1)-c(j+1,i-1)*U(j+2,i-1)))/(1-teta*(U(j+1,i-1)-U(j+2,i-1)-c(j+1,i-1)+c(j+2,i-1)));
         cS=(c(j+1,i-1)+teta*US*(c(j+1,i-1)-c(j+2,i-1)))/(1+teta*(c(j+1,i-1)-c(j+2,i-1)));
         YS=Y(j+1,i-1)+teta*(US-cS)*(Y(j+1,i-1)-Y(j+2,i-1));
         JS=J(j+1,i-1)-teta*(US-cS)*(J(j+2,i-1)-J(j+1,i-1));
         AS=A(j+1,i-1)-teta*(US-cS)*(A(j+2,i-1)-A(j+1,i-1));
         Y(j+1,i)=(1/(cR+cS))*(YS*cR+YR*cS+cR*cS*((UR-US)/g-dt*(JR-JS)-(q*dt/g)*((UR-cR)/AR-(US+cS)/AS)));
         U(j+1,i)=UR-g*((Y(j+1,i)-YR)/cR)-g*dt*(JR-ifo(j+1))-q*dt/AR*(UR-cR);
        end;
    end;
    % soluzioni "al contorno" a monte
    USmon=zeros(1);
    USmon=(U(1,i-1)-teta*(U(1,i-1)*c(2,i-1)-c(1,i-1)*U(2,i-1)))/(1-teta*(U(1,i-1)-U(2,i-1)-c(1,i-1)+c(2,i-1)));
    cSmon=(c(1,i-1)+teta*USmon*(c(1,i-1)-c(2,i-1)))/(1+teta*(c(1,i-1)-c(2,i-1)));
    YSmon=Y(1,i-1)+teta*(USmon-cS)*(Y(1,i-1)-Y(2,i-1));
    JSmon=J(1,i-1)-teta*(USmon-cSmon)*(J(2,i-1)-J(1,i-1));
    ASmon=A(1,i-1)-teta*(USmon-cSmon)*(A(2,i-1)-A(1,i-1));
     % metodo di Newton
        AAm=zeros(1);
        BBm=zeros(1);
        Uaus=U(1,i-1);
        AAm=1+Uaus/cSmon;
        BBm=-Uaus+USmon+g/cSmon*(Y(1,1)-(Uaus^2)/2/g-YSmon)-g*(JSmon-ifo(1))*dt-q*(USmon+cSmon)/ASmon*dt+AAm*Uaus;
        U(1,i)=BBm/AAm;
        diffm=zeros(1);
        diffm=abs(U(1,i)-Uaus);
        while diffm>0.0001,
              Uaus=U(1,i);
              AAm=1+Uaus/cSmon;
              BBm=-Uaus+USmon+g/cSmon*(Y(1,1)-(Uaus^2)/2/g-YSmon)-g*(JSmon-ifo(1))*dt-q*(USmon+cSmon)/ASmon*dt+AAm*Uaus;
              U(1,i)=BBm/AAm;
              diffm=abs(U(1,i)-Uaus);
        end;
        Y(1,i)=Y(1,1)-(U(1,i)^2)/2/g;
    % soluzioni "al contorno" a valle
    URval=(U(ps,i-1)+teta*(-U(ps,i-1)*c(ps-1,i-1)+c(ps,i-1)*U(ps-1,i-1)))/(1+teta*(U(ps,i-1)-U(ps-1,i-1)+c(ps,i-1)-c(ps-1,i-1)));
    cRval=(c(ps,i-1)-teta*URval*(c(ps,i-1)-c(ps-1,i-1)))/(1+teta*(c(ps,i-1)-c(ps-1,i-1))); 
    YRval=Y(ps,i-1)-teta*(URval+cRval)*(Y(ps,i-1)-Y(ps-1,i-1));
    JRval=J(ps,i-1)-teta*(URval+cRval)*(J(ps,i-1)-J(ps-1,i-1));
    ARval=A(ps,i-1)-teta*(URval+cRval)*(A(ps,i-1)-A(ps-1,i-1));
      % metodo di Newton
      AA=zeros(1);
      BB=zeros(1);
      Yaus=Y(ps,i-1);
      if sez==1,
         AA=-2*URval*Yaus+g*3*(Yaus^2)/cRval-2*g*YRval*Yaus/cRval+2*g*(JRval-ifo(ps))*dt*Yaus+2*q*(URval-cRval)/ARval*dt*Yaus;
         BB=-Q(i)*a+URval*(Yaus^2)-g/cRval*(Yaus-YRval)*(Yaus^2)-g*(Yaus^2)*(JRval-ifo(ps))*dt-q*(URval-cRval)/ARval*dt*(Yaus^2)+AA*Yaus;
         elseif sez==2,
                AA=-URval+2*g/cRval*Yaus-g/cRval*YRval+g*(JRval-ifo(ps))*dt+q*(URval-cRval)/ARval*dt;
                BB=-Q(i)/b+URval*Yaus-g/cRval*(Yaus^2)+g/cRval*YRval*Yaus-g*Yaus*(JRval-ifo(ps))*dt-q*(URval-cRval)/ARval*dt*Yaus+AA*Yaus;
                else 
                     AA=-URval*b-2*URval/a*Yaus+2*g/cRval*b*Yaus+3*g/cRval/a*(Yaus^2)-g/cRval*YRval*b-2*g*YRval/cRval/a*Yaus+g*(JRval-ifo(ps))*dt*b+2*g*(JRval-ifo(ps))*dt/a*Yaus+q*(URval-cRval)/ARval*dt*(b+2/a*Yaus);
                     BB=-Q(i)+URval*(b*Yaus+(Yaus^2)/a)-g/cRval*(Yaus-YRval)*(b*Yaus+(Yaus^2)/a)-g*(JRval-ifo(ps))*dt*(b*Yaus+(Yaus^2)/a)-q*(URval-cRval)/ARval*dt*(b*Yaus+(Yaus^2)/a)+AA*Yaus;
      end;
      Y(ps,i)=BB/AA;
      diff=zeros(1);
      diff=abs(Y(ps,i)-Yaus);
      while diff>0.0001,
            Yaus=Y(ps,i);
            if sez==1,
               AA=-2*URval*Yaus+g*3*(Yaus^2)/cRval-2*g*YRval*Yaus/cRval+2*g*(JRval-ifo(ps))*dt*Yaus+2*q*(URval-cRval)/ARval*dt*Yaus;
               BB=-Q(i)*a+URval*(Yaus^2)-g/cRval*(Yaus-YRval)*(Yaus^2)-g*(Yaus^2)*(JRval-ifo(ps))*dt-q*(URval-cRval)/ARval*dt*(Yaus^2)+AA*Yaus;
               elseif sez==2,
                      AA=-URval+2*g/cRval*Yaus-g/cRval*YRval+g*(JRval-ifo(ps))*dt+q*(URval-cRval)/ARval*dt;
                      BB=-Q(i)/b+URval*Yaus-g/cRval*(Yaus^2)+g/cRval*YRval*Yaus-g*Yaus*(JRval-ifo(ps))*dt-q*(URval-cRval)/ARval*dt*Yaus+AA*Yaus;
                      else 
                          AA=-URval*b-2*URval/a*Yaus+2*g/cRval*b*Yaus+3*g/cRval/a*(Yaus^2)-g/cRval*YRval*b-2*g*YRval/cRval/a*Yaus+g*(JRval-ifo(ps))*dt*b+2*g*(JRval-ifo(ps))*dt/a*Yaus+q*(URval-cRval)/ARval*dt*(b+2/a*Yaus);
                          BB=-Q(i)+URval*(b*Yaus+(Yaus^2)/a)-g/cRval*(Yaus-YRval)*(b*Yaus+(Yaus^2)/a)-g*(JRval-ifo(ps))*dt*(b*Yaus+(Yaus^2)/a)-q*(URval-cRval)/ARval*dt*(b*Yaus+(Yaus^2)/a)+AA*Yaus;
            end;   
      Y(ps,i)=BB/AA;
      diff=abs(Y(ps,i)-Yaus);
      end;
    % ricavo aree, larghezze del pelo libero, contorni bagnati, celerità
    % relative velocità di valle e cadenti energetiche per ogni sezione per
    % ogni iterazione
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
     U(ps,i)=Q(i)/A(ps,i);
     for j=1:ps,
         J(j,i)=(((U(j,i))^2)*(n^2))/((A(j,i)/C(j,i))^(4/3));
     end;
end;
% verifica della discretizzazione
errdis=0;
for i=1:ps,
    for j=2:pt,
        if ((errdis==0)&(dt>(ds/c(i,j)*sqrt(1-dt/2*(g*J(i,j)/U(i,j)+q/A(i,j)))))),
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
        if ((Fr>1)&(errFr==0)),
           errFr=1;
           d=input(' errore: la corrente è veloce in qualche sezione (scrivi 1) ');
        end;
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
for i=1:pt,
    Yfin(i,1)=Y(ps,i);
end;
end
