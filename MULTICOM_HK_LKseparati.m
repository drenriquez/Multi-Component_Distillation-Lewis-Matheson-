clc;
clear;

% DATI GENERICI

Rrif=1.5; %rapporto di riflusso
P=15; %Pressione totale P in atm:
N=6; %numero di componenti
Rid=0.00000000001; %Riduttore per i bilanci di massa  0.00000000001
Rid_dist=0.71;
%_________________________________________________________________________________________________________________________________________

% DATI FEED

% Fwt_tot=153189.2014; % portata totale feed in kg/h
Fwt_tot=176681;
XFwt=[0.003045959,0.004027002,0.202643965,0.222913649,0.251446959,0.315922466];% vettore delle frazioni in peso del feed

PM=[44.1,58.12,92.53908485,101.2788346,115.2837797,142.28]; %pesi molare in kg dei composti del feed

Fwt=Fwt_tot*XFwt;  % vettore delle portate in peso dei componenti del feed
Fmoli=Fwt./PM;  % vettore delle portate molari dei componenti del feed
Fmoli_tot=sum(Fmoli); % portata molare totale feed in Kmoli/h
XFmoli=Fmoli/Fmoli_tot; %vettore delle frazioni in moli del feed

Tc=[369.8,425.2,507.4,540.2,568.8,617.6]; % in °k
Pc=[42.5,38,29.7,27.4,24.8,21.1]; % in BAR
w=[171,0.197,0.316,0.344,0.42,0.502]; % fattore acentrico


%__________________________________________________________________________________________________________________________________________________________________________
  
% LIGHT KEY E HEAVY KEY

%fissare il componente chiave HK e LK, indicare il numero della colonna nel vettore F
HK=5;
LK=3;

%frazione in peso dell HK nel distillato:     0.02
DwtHK=0.04;%0.04
%frazione in peso dell LK nel bottom:         0.18
BwtLK=0.0065;%0.0065

%_____________________________________________________________________________________________________________________________________________________________

%DATI TERMICI PER I BILANCI ENERGETICI

% Ln (P°)= ANTA - ANTB/(T + ANTC), P°=[mmHg], T=[°k]
ant=[15.726	1872.46	-25.16
15.6782	2154.9	-34.42
15.8366	2697.55	-48.78
15.8737	2911.32	-56.51
15.9426	3120.29	-63.63
16.0114	3456.8	-78.67];
 
% temperatura di riferimento per l' ENTALPIA MOLARE, in °K.
Trif=400;
 
% Entalpa molare in Kj/mol
Hliq=[-102.5531,-126.1528,-176.5692,-199.1423,-222.0128,-265.5372];
Hgas=[-96.92885,-115.7603,-153.9882,-174.8683,-193.6997,-229.3598];

%CpmediL in kj/mol k
CpmedioL=[0.113426364,0.180620909,0.349664545,0.333766364,0.341537273,0.402278182];
CpmedioG=[0.106716818,0.142386364,0.218328182,0.269465455,0.295111818,0.313349091];

%__________________________________________________________________________________________________________________________________________________________________
% __________________________________________________________________________________________________________________________________________________________________

%BILANCI DI MATERIA

Dsum=Rid*(sum(Fwt(HK+1:N)));
Dsot=Rid*(sum(Fwt(1:LK-1)));

B1_=[sum(Fwt(1:LK))+Rid_dist*sum(Fwt(LK+1:HK-1));sum(Fwt(HK:N))+((1-Rid_dist)*sum(Fwt(LK+1:HK-1)))]; % vettore dei termini noti del sistema di equazioni del bilancio
B_=[B1_(1)+Dsum-Dsot B1_(2)-Dsum+Dsot]; % vettore dei termini noti del sistemi di equazione del bilanci 
    
A=[1-DwtHK   BwtLK
   DwtHK     1-BwtLK];  %matrice dei termini noti del sistema di equazione del bilancio   

%Sintassi matriciale del metodo di gauss della risoluzione di sistemi lineari
C=(inv(A))';
X=B_*C;  %vettore delle incognite D e B

Dwt_tot=X(1); %portata mssica distillato
Bwt_tot=X(2); %portata massica bottom


%____________________________________________________________________________________________________________________________________________________________________

%CARATTERIZZAZIONE DEL BOTTOM E DISTILLATO

%costruzione del vettore delle portate massiche dei componenti nel distillato
Dwt(1:LK-1)=Fwt(1:LK-1)-Rid*Fwt(1:LK-1);
Dwt(LK)=Fwt(LK)-Bwt_tot*BwtLK;
Dwt(LK+1:HK-1)=Rid_dist*sum(Fwt(LK+1:HK-1));
Dwt(HK)=Dwt_tot*DwtHK; 
Dwt(HK+1:N)=Rid*Fwt(HK+1:N);

%costruzione del vettore delle portate massiche dei componenti nel bottom 
Bwt(1:LK-1)=Rid*Fwt(1:LK-1);
Bwt(LK)=Bwt_tot*BwtLK;
Bwt(LK+1:HK-1)=(1-Rid_dist)*sum(Fwt(LK+1:HK-1));
Bwt(HK)=Fwt(HK)-Dwt_tot*DwtHK;
Bwt(HK+1:N)=Fwt(HK+1:N)-Rid*Fwt(HK+1:N);

%Vettori delle portate molari dei componenti nel distillato e nel bottom
Dmoli=Dwt./PM;



% _______________________________________________________________________________________________________________-
% _______________________________________________________________________________________________________________-
% _______________________________________________________________________________________________________________-

% MACRO LOOP CONVERGENZA LEWIS-MATHESON


for z=1:2

% Dmoli=[10.5806641208979 10.6128604712561 271.170837256415 9.31231511877650 2.59435093103905e-07 6.39573382032817e-08];
Bmoli=Fmoli-Dmoli;
%Portata molare totale nel bottom e nel distillato
Dmoli_tot=sum(Dmoli);
Bmoli_tot=sum(Bmoli);

% Vettori delle frazioni molari dei componenti nel bottom e nel distillato

XDmoli=Dmoli/Dmoli_tot;
XBmoli=Bmoli/Bmoli_tot;

YDmoli=XDmoli; %Condensatore parziale --> tutto D è in fase gas (e conosciamo tutto di DY)

%_______________________________________________________________________________________________________________-

%Vettore delle temperature di ebollizione dei componenti alla pressione operativa della colonna
for i=1:N
[Teb(i)]=ANTOINEINV(P*760,ant(i,1),ant(i,2),ant(i,3));
end


%__________________________________________________________________________________________________________________

% TEMPERATURA FEED

%Teb ribollitore 
Tmax=Teb(N);
Tmin=Teb(1);
% Tip ---> P°---> Kvalue=P°/P ---> Yd=Xd*Kvalue ---> Verifica: sommatoria(Yd)=1

for j=1:50
Tfeed=(Tmax+Tmin)/2;

    for i=1:N
    [Po(i)]=ANTOINE(Tfeed,ant(i,1),ant(i,2),ant(i,3));
%      [Po(i)]=PENG_ROBINSON(Tfeed,w(i),Tc(i),Pc(i));
    end
    
    Kvalue=Po/(P*760);
    YFmoli=XFmoli.*Kvalue;
    sommatoriaY=sum(YFmoli);
    
    if sommatoriaY>1
        Tmax=Tfeed;
    else Tmin=Tfeed;
    end
end




%__________________________________________________________________________________________________________________

% CONDENSATORE

% Dew point cal. Per stabilire la temperatura del condensatore e la composizione di Lo 
% Tip ---> P°---> Kvalue=P°/P ---> Xd=Yd/Kvalue ---> Verifica: sommatoria(Xd)=1
%Temperatura di ebollizione alla pressione operativa minima e massima.
%Servono per usare il metodo di bisezione
 Tmax=Teb(N);
 Tmin=Teb(1); 
 
for j=1:50
Tcond=(Tmax+Tmin)/2;

    for i=1:N
    [Po(i)]=ANTOINE(Tcond,ant(i,1),ant(i,2),ant(i,3));
%     [Po(i)]=PENG_ROBINSON(Tcond,w(i),Tc(i),Pc(i));
    end
    P1(j,:)=Po;
    Kvalue=Po/(P*760);
    XLomoli=YDmoli./Kvalue;
    sommatoriaX=sum(XLomoli);
    
    if sommatoriaX<1
        Tmax=Tcond;
    else Tmin=Tcond;
    end
end
      
%bilancio di materia dei componenti sul condensatore per caratterizzare la corrente V1

Lomoli_tot=Rrif*Dmoli_tot; %portata totale in moli di Lo
Lomoli=XLomoli*Lomoli_tot; %vettore delle portate molari dei componenti in Lo (bilancio sui componenti al condensatore)
Vmoli=Lomoli+Dmoli; %vettore delle portate in moli dei componenti nella corrente V1
YVmoli=Vmoli/sum(Vmoli); %vettore delle frazioni molari dei componenti della corrente V1

%Dew point calc per trovare la temperatura del primo piatto
Tmax=Teb(N);
Tmin=Teb(1);

%Dew point cal per calcolare la temperatura del primo piatto e la composizione di L1 
% Tip ---> P°---> Kvalue=P°/P ---> Xd=Yd/Kvalue ---> Verifica: sommatoria(Xd)=1

for j=1:50
T=(Tmax+Tmin)/2;

    for i=1:N
    [Po(i)]=ANTOINE(T,ant(i,1),ant(i,2),ant(i,3));
%     [Po(i)]=PENG_ROBINSON(T,w(i),Tc(i),Pc(i));
    end
    
    Kvalue=Po/(P*760);
    XLmoli=YVmoli./Kvalue;
    sommatoriaX=sum(XLmoli);
    
    if sommatoriaX<1
        Tmax=T;
    else Tmin=T;
    end
end

% Bilancio di energia del condensatore
HV=sum(Vmoli.*(Hgas+CpmedioG*(T-Trif))); % portata totala di calore V1
HDv=sum(Dmoli.*(Hgas+CpmedioG*(Tcond-Trif)));%portata totale di calore di Dv
hLo=sum(Lomoli.*(Hliq+CpmedioL*(Tcond-Trif)));%portata totale di di Lo

Qc=HV-HDv-hLo; %Bilancio di energia al condensatore per avere HV1 e Qc

%_________________________________________________________________________________________________________________________________

%RIBOLLITORE

%Teb ribollitore 
Tmax=Teb(N);
Tmin=Teb(1);


% Tip ---> P°---> Kvalue=P°/P ---> Yd=Xd*Kvalue ---> Verifica: sommatoria(Yd)=1

for j=1:50
Trib=(Tmax+Tmin)/2;

    for i=1:N
    [Po(i)]=ANTOINE(Trib,ant(i,1),ant(i,2),ant(i,3));
%     [Po(i)]=PENG_ROBINSON(Trib,w(i),Tc(i),Pc(i));
    end
    
    Kvalue=Po/(P*760);
    YVomoli=XBmoli.*Kvalue;
    sommatoriaY=sum(YVomoli);
    
    if sommatoriaY>1
        Tmax=Trib;
    else Tmin=Trib;
    end
end

% Bilancio di energia totale per ricavare  Qr

HDv=sum(Dmoli.*(Hgas+CpmedioG*(Tcond-Trif)));%portata totale di calore di Dv
hF=sum(Fmoli.*(Hliq+CpmedioL*(Tfeed-Trif)));%portata totale di calore del Feed
hB=sum(Bmoli.*(Hliq+CpmedioL*(Trib-Trif)));%portata totale di calore di B

Qr=hB+HDv-hF+Qc; %Bilancio di energia al condensatore per avere HV1 e Qc
%_________________________________________________________________________________________________________________________________

% LOOP SEZIONE D ARRICCHIMENTO

 %(il loop parte dal secondo piatto,  perchè il primo è stato già %caratterizzato)
 
Arresto=XFmoli(LK)/XFmoli(HK);% Valore che ferma il ciclo del loop sezione arricchimento (per trovare il piatto del feeed)
Tolleranza=0.15; %range di tolleranza per il valore di arresto

Lmoli(1,:)=zeros(1,N); %Inizializzazione matrice costituita dalle portate molari in fase liquida (colonne) nei vari piatti (righe)
XLmoli(1,:)=XLmoli;% Inizializzazione matrice costituita dalle frazioni molari in fase liquida(colonne)nei vari piatti(righe)
Vmoli(1,:)=Vmoli; % Inizializzazione matrice costituita dalle portate molari in fase vapore (colonne) nei vari piatti (righe)
YVmoli(1,:)=YVmoli;% Inizializzazione matrice costituita dalle frazioni molari in fase vapore(colonne)nei vari piatti(righe)

T(1)=T;% Inizializzazione vettore delle temperture dei vari piatti
k=1 ;

%while k<3
    
    

 while and(abs((XLmoli(k,LK)/XLmoli(k,HK))-Arresto)>Tolleranza,k<30)

     rapportoLK_HK(k)=XLmoli(k,LK)/XLmoli(k,HK);
     
     Lmoli_totmax=2*Lomoli_tot;
     Lmoli_totmin=0;
     
     for i=1:50
         Lmoli_tot(k)=(Lmoli_totmax+Lmoli_totmin)/2;
         Lmoli(k,:)=XLmoli(k,:)*Lmoli_tot(k);
         Vmoli(k+1,:)=Lmoli(k,:)+Dmoli;
         Vmoli_tot(k+1)=sum(Vmoli(k+1,:));
         YVmoli(k+1,:)=Vmoli(k+1,:)/Vmoli_tot(k+1);
           
       % Dew point cal. Per stabilire la temperatura e la composizione di L
       % Tip ---> P°---> Kvalue=P°/P ---> Xd=Yd/Kvalue ---> Verifica: sommatoria(Xd)=1

             Tmax=Teb(N);
             Tmin=Teb(1);

                for j=1:50
                T(k+1)=(Tmax+Tmin)/2;

                    for i=1:N
                   [Po(i)]=ANTOINE(T(k+1),ant(i,1),ant(i,2),ant(i,3));
%                     [Po(i)]=PENG_ROBINSON(T(k+1),w(i),Tc(i),Pc(i));
                    end
              
                Kvalue=Po/(P*760);
                XLmoli(k+1,:)=YVmoli(k+1,:)./Kvalue;
                sommatoriaX=sum(XLmoli(k+1,:));
    
                if sommatoriaX<1
                Tmax=T(k+1);
                 else Tmin=T(k+1);
                end
                end
         
             % Bilancio di energia che ingloba il piatto k-esimo ed il condensatore
            HVk=sum(Vmoli(k+1,:).*(Hgas+CpmedioG*(T(k+1)-Trif))); % portata totala di calore V1
            HDv=sum(Dmoli.*(Hgas+CpmedioG*(Tcond-Trif)));%portata totale di calore di Dv
            hL=sum(Lmoli(k,:).*(Hliq+CpmedioL*(T(k)-Trif)));%portata totale di di Lo

            Bilancio(k)=HVk-hL-HDv-Qc;
            if Bilancio(k)>0
                Lmoli_totmax=Lmoli_tot(k)
            else  Lmoli_totmin=Lmoli_tot(k)
            end
         
         
         
         
     end
      
     
   k=k+1;
   PIATTI_TESTA=k;
 end  

 %_________________________________________________________________________________________________________________________________

% LOOP SEZIONE DI ESAURIMENTO
Tolleranza=0.12

Lbmoli(1,:)=Bmoli; %Inizializzazione matrice costituita dalle portate molari in fase liquida (colonne) nei vari piatti (righe)
XLbmoli(1,:)=XBmoli;% Inizializzazione matrice costituita dalle frazioni molari in fase liquida(colonne)nei vari piatti(righe)
Vbmoli(1,:)=zeros(1,N); % Inizializzazione matrice costituita dalle portate molari in fase vapore (colonne) nei vari piatti (righe)
YVbmoli(1,:)=YVomoli;% Inizializzazione matrice costituita dalle frazioni molari in fase vapore(colonne)nei vari piatti(righe)

Tb(1)=Trib;% Inizializzazione vettore delle temperture dei vari piatti
k=1 ;
%while k<6
 while and(abs((YVbmoli(k,LK)/YVbmoli(k,HK))-Arresto)>Tolleranza,k<30)
     rapportobLK_HK(k)=YVbmoli(k,LK)/YVbmoli(k,HK);
     rapportobXLK_HK(k)=XLbmoli(k,LK)/XLbmoli(k,HK);
  
     
     Vbmoli_totmax=2*Bmoli_tot;
     Vbmoli_totmin=0;

         for i=1:50
         Vbmoli_tot(k)=(Vbmoli_totmax+Vbmoli_totmin)/2;
         Vbmoli(k,:)=YVbmoli(k,:)*Vbmoli_tot(k);
         Lbmoli(k+1,:)=Vbmoli(k,:)+Bmoli;
         Lbmoli_tot(k+1)=sum(Lbmoli(k+1,:));
         XLbmoli(k+1,:)=Lbmoli(k+1,:)/Lbmoli_tot(k+1);
           
       % bubble point cal. Per stabilire la temperatura e la composizione di V
       % Tip ---> P°---> Kvalue=P°/P --->Yd=Xd*Kvalue ---> Verifica: sommatoria(Yd)=1

             Tmax=Teb(N);
             Tmin=Teb(1);

                  for j=1:50
                  Tb(k+1)=(Tmax+Tmin)/2;

                    for i=1:N
                   [Po(i)]=ANTOINE(Tb(k+1),ant(i,1),ant(i,2),ant(i,3));
%                     [Po(i)]=PENG_ROBINSON(T(k+1),w(i),Tc(i),Pc(i));
                    end
              
                Kvalue=Po/(P*760);
                YVbmoli(k+1,:)=XLbmoli(k+1,:).*Kvalue;
                sommatoriaY=sum(YVbmoli(k+1,:));
    
                if sommatoriaY>1
                Tmax=Tb(k+1);
                 else Tmin=Tb(k+1);
                end
                end


            % Bilancio di energia che ingloba il piatto k-esimo ed il condensatore
            HVbk=sum(Vbmoli(k,:).*(Hgas+CpmedioG*(Tb(k)-Trif))); % portata totala di calore V1
            hB=sum(Bmoli.*(Hliq+CpmedioL*(Trib-Trif)));%portata totale di calore di B
            hLb=sum(Lbmoli(k+1,:).*(Hliq+CpmedioL*(Tb(k+1)-Trif)));%portata totale di di Lo

            Bilanciob(k)=Qr+hLb-hB-HVbk
            if Bilanciob(k)<0
                Vbmoli_totmax=Vbmoli_tot(k)
            else  Vbmoli_totmin=Vbmoli_tot(k)
            end
 
         
     end
      
     
   k=k+1  
   PIATTI_CODA=k-1;
end   
%CRITERIO DI CONVERGENZA

%vettore errore e=f(i)+l_testa(i)-l_coda(i)
epsilon=Fmoli+Lmoli(PIATTI_TESTA-1,:)-Lbmoli(PIATTI_CODA+1,:);
alfa=Lbmoli(PIATTI_CODA+1,:)./(Fmoli+Lmoli(PIATTI_TESTA-1,:));

Dmoli_calc=Dmoli.*alfa;
Bmoli_calc=Bmoli./alfa;
%SE(dk<=bk and dk+1_calc<=f) allora dk+1=dk+1_calc; else dk+1=f-bk+1_calc)
for i=1:N
    if and(Dmoli(i)<=Bmoli(i),Dmoli_calc(i)<=Fmoli(i))
        Dmoli_new(i)=Dmoli_calc(i);
    else
        Dmoli_new(i)=Fmoli(i)-Bmoli_calc(i);
    end
end
DmoliB(z,:)=Dmoli;%serve per tener traccia del vettore D nelle varie iterazioni z
Dmoli=Dmoli_new;
end
  DmoliB(z+1,:)=Dmoli;

%_________________________________________________________________________________________________________________________________
% PLOT




Tb1=Tb(PIATTI_CODA+1:-1:1);
Tpiatti=T;
Tpiatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA+1)=Tb1;

PIATTI=[PIATTI_TESTA+PIATTI_CODA+1:-1:1];
Tpiatti2=Tpiatti(2:PIATTI);%senza rib


XLK_piatticoda=XLbmoli((PIATTI_CODA+1:-1:1),LK)';
XLK_piatti=XLmoli(:,LK)';
XLK_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA+1)=XLK_piatticoda;
XLK_piatti2=XLK_piatti(1:PIATTI_TESTA+PIATTI_CODA);

XHK_piatticoda=XLbmoli((PIATTI_CODA+1:-1:1),HK)';
XHK_piatti=XLmoli(:,HK)';
XHK_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA+1)=XHK_piatticoda;
XHK_piatti2=XHK_piatti(1:PIATTI_TESTA+PIATTI_CODA);

Y1_piatticoda=YVbmoli((PIATTI_CODA+1:-1:1),1)';
Y1_piatti=YVmoli(:,1)';
Y1_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA+1)=Y1_piatticoda;
Y1_piatti2=Y1_piatti(1:PIATTI_TESTA+PIATTI_CODA);

YLK_piatticoda=YVbmoli((PIATTI_CODA+1:-1:1),LK)';
YLK_piatti=YVmoli(:,LK)';
YLK_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA+1)=YLK_piatticoda;
YLK_piatti2=YLK_piatti(1:PIATTI_TESTA+PIATTI_CODA);

YHK_piatticoda=YVbmoli((PIATTI_CODA+1:-1:1),HK)';
YHK_piatti=YVmoli(:,HK)';
YHK_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA+1)=YHK_piatticoda;
YHK_piatti2=YHK_piatti(1:PIATTI_TESTA+PIATTI_CODA);

VLK_piatticoda=Vbmoli((PIATTI_CODA:-1:1),LK)';
VLK_piatti=Vmoli(:,LK)';
VLK_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA)=VLK_piatticoda;

VHK_piatticoda=Vbmoli((PIATTI_CODA:-1:1),HK)';
VHK_piatti=Vmoli(:,HK)';
VHK_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA)=VHK_piatticoda;

V1_piatticoda=Vbmoli((PIATTI_CODA:-1:1),1)';
V1_piatti=Vmoli(:,1)';
V1_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA)=V1_piatticoda;

V2_piatticoda=Vbmoli((PIATTI_CODA:-1:1),2)';
V2_piatti=Vmoli(:,2)';
V2_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA)=V2_piatticoda;

V6_piatticoda=Vbmoli((PIATTI_CODA:-1:1),6)';
V6_piatti=Vmoli(:,6)';
V6_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA)=V6_piatticoda;

V4_piatticoda=Vbmoli((PIATTI_CODA:-1:1),4)';
V4_piatti=Vmoli(:,4)';
V4_piatti(PIATTI_TESTA+1:PIATTI_TESTA+PIATTI_CODA)=V4_piatticoda;


hold on

% plot(PIATTI,Tpiatti)
plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),YLK_piatti2)
plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),YHK_piatti2)
% plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),VHK_piatti)
% plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),V1_piatti)
% plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),V2_piatti)
% plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),V6_piatti)
% plot(PIATTI(1:PIATTI_TESTA+PIATTI_CODA),V4_piatti)



% FUNCTION PER IL CALCOLO DELLE TENSIONI DI VAPORE
function[Po]=ANTOINE(T,A,B,C) %Function che riceve la temperatura e le costanti di Antoine e ci restituisce la P° 
Po=exp(A-(B/(T+C)));
end

function[Teb]=ANTOINEINV(P,A,B,C);  %Function che riceve una pressione (quella operativa della colonna) e ci restituisce la temperatura di ebollizione dei componenti (questo i serve per avere una Tmin e Tmax)
Teb=(A*C-C*log(P)-B)/(log(P)-A);
end
function[P0]=PENG_ROBINSON(T,w1,Tc1,Pc2)
R=8.314;
%benzene
w=w1;
Tc=Tc1;% in kelvin
Pc1=Pc2;%in bar
Pc=Pc1*100000;%conv in pascal
%parametri che servono per il loop
n_loop=50;
Pmax=Pc;
Pmin=0;
% n_cicli=1000;% quante volte deve eseguire il ciclo
% Tmax=Tc;
% Tmin=393;
% delta_T=(Tmax-Tmin)/n_cicli;
% Tmax1=Tmax+delta_T;
 
% T(i)=Tmax1-delta_T;
 Tx=T;
%parametri che servono per il loop
Pmax=Pc*2.7;
Pmin=0;
    for k=1:n_loop
    Px=(Pmax+Pmin)/2;
    Tr=Tx/Tc;
    K=0.37463+1.54226*w-0.26992*w^2;
    a_1=(1+K*(1-(Tr^(1/2))))^2;
    a=0.45724*(R^2)*(Tc^2)/Pc;
    b=0.0778*R*Tc/Pc;
    a_=(1+(0.37464+1.54226*w-0.26992*w^2)*(1-Tr^(1/2)))^2;
    A=a*a_*Px/((R^2)*Tx^2);
    B=b*Px/(R*Tx);
    %forma polinomiale in Z della cubica Peng Robinson
    %a1*Z^3+a2Z^2+a3*Z+a4=0
    a1=1;
    a2=-(1-B);
    a3=A-2*B-3*B^2;
    a4=-(A*B-B^2-B^3);
    %matrice dei coefficienti della forma polinomiale
    M=[a1 a2 a3 a4];
    %soluzioni della M, (z1 z2 z3):
    Z=roots(M);
    %soluzioni della M, (v1 v2 v3), volumi in m3/mole:
    vj(:)=Z*(R*Tx)/Px;
    Zv=Z(1);
    Zl=Z(3);
    %fugacità vapore
    fv=Px*exp((Zv-1)-log(Zv-B)-(A/(2*sqrt(2)*B))*log((Zv+(1+sqrt(2))*B)/(Zv+(1-sqrt(2))*B)));
    %fugacità liquido
    fl=Px*exp((Zl-1)-log(Zl-B)-(A/(2*sqrt(2)*B))*log((Zl+(1+sqrt(2))*B)/(Zl+(1-sqrt(2))*B)));
        if fl>=fv
        Pmin=Px;
        else
        Pmax=Px;
        end
    end
P0=(Px/100000)*760;

end