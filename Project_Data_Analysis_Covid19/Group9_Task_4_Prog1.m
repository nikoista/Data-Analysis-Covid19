%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS

%% EISAGWGH TWN DEDOMENWN GIA TON COVID 19 1/1/2020 --- 13/12/2020
DATAConfirmed = importdata('Covid19Confirmed.xlsx','headerlinesIn');
DATAConfirmed = DATAConfirmed.data;

DATADeath = importdata('Covid19Deaths.xlsx','headerlinesIn');
DATADeath = DATADeath.data;

g=0;
population = DATAConfirmed(131+g,1);
if(population == 46937060)
    g = 0;
else 
    g=-1;
end


%% XWRA OMADAS : ISPANIA 131
%% XWRES EU : AUSTRIA 9 , BELGIO 14,  CROATIA 34 ,  DANIA 38, LITHUANIA 83,
%% ESTONIA 45 ,FINNLAND 48 , OLLANDIA 98, NORBIGIA 104, ELVETIA 135
nfig=0;

%% SPAIN 131

popSpa = DATAConfirmed(131+g ,1);
spac = DATAConfirmed(131+g ,2:end);
spad = DATADeath(131 +g ,2:end);

%% BELGIUM 14
popBel = DATAConfirmed(14 +g,1);
belc = DATAConfirmed(14+g  ,2:end);
beld = DATADeath(14+g ,2:end);

%% DENMARK 38
popDen = DATAConfirmed(38+g ,1);
denc = DATAConfirmed(38 +g ,2:end);
dend = DATADeath(38+g ,2:end);

%% NETHERLANDS 98
popOll = DATAConfirmed(98+g  ,1);
ollc = DATAConfirmed(98+g,2  :end);
olld = DATADeath(98+g ,2:end);

%% NORWAY 104
popNor = DATAConfirmed(104+g  ,1);
norc = DATAConfirmed(104+g ,2:end);
nord = DATADeath(104+g ,2:end);

%% SWITZERLAND 135
popElv = DATAConfirmed(135+g  ,1);
elvc = DATAConfirmed(135 +g ,2:end);
elvd = DATADeath(135 +g ,2:end);


startc= [ 57 ;  61 ;   66 ;  62  ; 62 ; 60];
endc = [ 154 ; 186 ; 186 ;   179 ; 138; 130];


startd= [ 66 ;  73  ;  75  ;  74  ; 74 ; 71 ];
endd = [ 158 ; 206 ;  194 ;  190 ; 138; 144];


bounds = [60, 158; 65 , 206; 70 , 194; 68 , 194 ; 68 , 138 ; 66 , 144];

CC = [spac ;  belc ;  denc ;  ollc ; norc ; elvc];
D =  [spad ;  beld ;  dend ;  olld ; nord ; elvd];

COUNTRY = {'SPAIN',' BELGIUM', 'DENMARK' ,'NETHERLANDS', 'NORWAY','SWITZERLAND'};


%% STATHERES 
alpha=0.05;
B1 = -20;
B2 = 20;

%% ARXIKOPOIHSH PINAKWN
N=size(CC,1);
cc = zeros(size(CC));
d = zeros(size(D));
duration = zeros(N,1);
r = zeros(N,1);
r2 = zeros(N,1);
t =  zeros(N,1);
bestT =  zeros(N,1);
max_r2 =  zeros(N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    
    duration(i) = (endd(i) - endc(i));
    
    %% DATA SMOOTHING
    [cc(i,:) , d(i,:)] =  Group9Exe1Fun1(CC(i,:),D(i,:));
    
    %% EPANALHPTIKH DIADIKASIA EURESHS XRONIKWN USTERISEWN
    c=1;
    for T=B1:1:B2
        
        firstFD = d(i,(bounds(i,1)):(bounds(i,2)));
        
        %% r(X(t),Y(t+T))
        fc = cc(i,bounds(i,1):(bounds(i,2)));
        fd = d(i,(bounds(i,1)+T):(bounds(i,2)+T));

        x = fc;
        y = fd;
        
        n = length(x);
        meanx = mean(x);
        meany = mean(y);
        sumxy = 0;
        sumx2 = 0;
        sumy2 = 0;
        for k=1:1:n
            sumxy = sumxy + x(k)*y(k);
            sumx2 = sumx2 + (x(k)^2);
            sumy2 = sumy2 + (y(k)^2);
        end
        sx = sqrt(sumx2 - n*(meanx^2));
        sy = sqrt(sumy2 - n*(meany^2));
        cor = corrcoef(x,y);
        
        r(i,c) = (sumxy - n*meanx*meany)/(sx*sy);
        r2(i,c) = r(i,c)*r(i,c)*100;
        
        t(c)=T;
        c=c+1;
    end  
    
    %% BRISKW TO MEGISTO r(Pearson) KAI THN XRONIKH USTERHSH POU ANTISTOIXH
    [max_r2(i),I] = max(r2(i,:));
    bestT(i) = t(I);
    
    %% PLOT HMERISIOUS THANATOUS ME KAI XWRIS XRONIKH USTERHSH 
    delayed_FD = d(i,(bounds(i,1)+bestT(i)):(bounds(i,2)+bestT(i)));
    
    nfig = nfig +1;
    figure(nfig)
    N = length(firstFD);
    plot(1:N,firstFD,'-.',1:N,delayed_FD)
    title(strcat(COUNTRY{i},' Time Delay given by r(Pearson) for daily Deaths Covid 19'));
    modtxt = sprintf('r = %2.5f ',max(r(i,:)));
    ax = axis;
    text(ax(1)+0.05*(ax(2)-ax(1)),ax(4)-0.15*(ax(4)-ax(3)),modtxt);
    ylabel('$Daily Deaths$','Interpreter','latex','fontsize',10);
    xlabel('$Days$','Interpreter','latex','fontsize',10);
    legend('No Time Delay' , 'Time Delay');
    
    fprintf('\n');
    fprintf(' PRESS ANY KEY TO CONTINUE \n');
    pause;
end

%% PRINTS APOTELESMATWN
fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' COUNTRIES     :   SPAIN , BELGIUM , DENMARK , NETHERLANDS , NORWAY , SWITZERLAND \n');  
fprintf(' r^2(percent)  :  [%2.4f , %2.4f  , %2.4f   , %2.4f , %2.4f  , %2.4f]\n',max_r2(1),max_r2(2),max_r2(3),max_r2(4),max_r2(5),max_r2(6));
fprintf(' TIME DELAYS   :  [%d       , %d        , %d        , %d        , %d      , %d  ]\n',bestT(1),bestT(2),bestT(3),bestT(4),bestT(5),bestT(6));
fprintf(' *******************************************************************\n');
fprintf('\n');

%% ********************** SXOLIASMOS APOTELESMATWN ************************
% STHN SUNEXIA THS ERGASIA PROSPATHISAME NA BROUME THN XRONIKH USTERUSH TWN
% HMERISIWN KROUSMATWN POU HTAN ANAMENOMENH UPOLOGIZONTAS TON SUNTELESTH 
% r PEARSON. 

%% ERWTHSH 1)
% ME TON UPOLOGISMO AUTOU TOU SUNTELESTH KATALAVAINOUME PWS ENAS 
% ARITHMOS X KROUSMATWN ENERGOPOIEI ENAN Y ARITHMO THANATWN . OPOTE OSO 
% PERISSOTERO SUSXETIZONTAI TA X,Y (MEGALHTERO r) TOSO KALUTERES PROVLEPSEIS
% THA EXOUME GIA THN METAXU TOUS USTERHSH.
% ETSI BRISKONTAS THN MEGISTH TIMH TOU SUNTELESTH r ,BRISKOUME TAUTOXRONA 
% THN XRONIKH USTERHSH POU TOU ANTISTOIXEI KAI EXOUME TA EXHS :
% USTERHSEIS     : [6       , 5     , 0    , 5        , 11       , 6  ] 
% XWRES          : [ISPANIA ,BELGIO ,DANIA , OLLANDIA , NORVIGIA , ELVETIA]
% r^2 %100:   [83.8317 , 79.0079  , 71.6224   , 86.0119 , 38.2292  , 54.8521]
% VLEPOUME PWS H ISPANIA KAI H OLLANDIA EXOUN MEGALH GRAMMIKH SUSXETHSH >
% 80 % ENW TO BELGIO EINAI EPISHS KONTA STO 80% . LIGO PIO MARKUA APO TO
% 80% EINAI H DANIA(MESEA SUSXETHSH) KAI ARKETA MAKRUA ( KAKH SUSXETISH) 
% EXEI H NORBIGIA KAI H ELVETIA.ARA MPOROUME NA POUME PWS GIA THN ISPANIA
% THN OLLANDIA KAI TO BELGIO OTI H PROSEGGISH AUTH EKTIMA SWSTA THN
% USTERHSH THS POREIAS TWN HMERHSIWN THANATWN WS PROS THN POREIA TWN
% HMERHSIWN KROUSMATWN , OMWS DEN MPOROUME NA POUME TO IDIO GIA TIS
% UPOLUPES XWRES. ARA H TROPOS AUTOS PROSEGGISHS EIXE 50% EPITUXIA.

%% ERWTHMA 2)
% STO PROIGOUMENO ZHTHMA 3 EXETASAME TO AN MPOROUME NA POUME OTI UPARXEI
% SUSXETISH METAXU THS X meras (GIA TA NEA KROUSMATA) KAI THN X+14 (GIA
% TOUS THANATOUS) OPOU KAI VGALAME ENA DIASTHMA EMPISTOSUNHS ~~ [4 , 14] .
% STO ZHTHMA OMWS AUTO THELOUME NA VROUME THN MEGISTH GRAMMIKH SUSXETISH 
% POU DHMIOURGEI H XRONIKH USTERHSH T .
% TO SPOUDAIO EINAI OTI H EKTIMHSH AUTH THS USTERHSHS ME THN EKTIMHSH TOU
% ZHTHMATOS 3 EINAI KONTA!!!
% 3) TIME DELAYS : [6  , 5 , 4 , 2  , 14 , 9  ] 
% 4) TIME DELAYS:  [6  , 5 , 0 , 5  , 11 , 6  ]


%%  PROVLIMATISMOI
% OSTOSO UPHRXAN KAI PERIPTWSEIS
% ME XWRES POU H USTERHSH EVGAZE ARNHTIKES MERES POU SUMAINEI OTI EIXAME
% XRONIKH OLISTHHSH PROS TA THETIKA ARA TO ANTITHETO THS USTERHSHS OMWS GIA
% NA APODEIXOUME THN UPOTHESH TOU ZHTHMATOS EPILEXAME MERES ME THETIKES
% XRONIKES USTERHSEIS.