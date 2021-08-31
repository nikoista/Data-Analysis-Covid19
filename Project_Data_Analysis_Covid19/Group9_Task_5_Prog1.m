%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 5
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


startc= [ 57 ;  61 ;   66 ;  58   ; 62 ; 60];
endc = [ 154 ; 186 ; 186 ;   140 ; 138; 130];

startd= [ 66 ;  73  ;  75  ;  74  ; 74 ; 71 ];
endd = [ 158 ; 206 ;  194 ;  190 ; 138; 144];


bounds = [60, 158; 65 , 206; 70 , 194; 68 , 194 ; 68 , 138 ; 66 , 144];
CC = [spac ;  belc ;  denc ;  ollc ; norc ; elvc];
D =  [spad ;  beld ;  dend ;  olld ; nord ; elvd];
COUNTRY = {'SPAIN',' BELGIUM', 'DENMARK' ,'NETHERLANDS', 'NORWAY','SWITZERLAND'};


distnames = {'Normal','LogLogistic','Logistic','Poisson',...
    'Rayleigh','tLocationScale','Exponential','ExtremeValue','Gamma'};


%% CONSTANTS
alpha=0.05;
B1 = -20;
B2 = 20;

%% MATRICES INITIALIZATION
N=size(CC,1);
cc = zeros(size(CC));
d = zeros(size(D));
duration = zeros(N,1);
r = zeros(N,B2+1);
bestT =  zeros(N,1);
beta = zeros(2,B2+1);
R2 = zeros(N,B2+1);
adjR2 = zeros(N,B2+1);
bestB =  zeros(N,2);
bestr =  zeros(N,1);
bestR2 =  zeros(N,1);
bestadjR2 =  zeros(N,1);
bestI =  zeros(N,1);
bestSumEstar =  zeros(N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    
    close all ;
    
    %% DATA SMOOTHING
    [cc(i,:) , d(i,:)] =  Group9Exe1Fun1(CC(i,:),D(i,:));
    L = length(bounds(i,1):(bounds(i,2)));
    y = d(i,bounds(i,1):(bounds(i,2)));
    duration(i) = (endd(i) - endc(i));
    
    %% ARXIKOPOIHSH PINAKWN
    b = zeros(2,B2+1);
    yhat = zeros(L,B2+1);
    estar = zeros(L,B2+1);
    x = zeros(L,B2+1);
    bestYhat = zeros(1,L);
    bestEstar = zeros(1,L);
    bestX = zeros(1,L);
    
    %% EPANALHPTIKH DIADIKASIA EURESHS XRONIKWN USTERISEWN
    c=1;
    for T=0:1:B2
        T
        
        %% X(t-T) AND Y(t)
        x(:,c) = cc(i,(bounds(i,1)-T):bounds(i,2)-T);
        n = length(x(:,c));
        
        %% tcritical Zcritical
        tcrit = tinv(1-alpha/2,n-2);
        zcrit = norminv(1-alpha/2);
        
        
        %% STATISTICAL VALUES OF x,y
        n = length(x(:,c));
        meanx = mean(x(:,c));
        meany = mean(y);
        
        sumxy = 0;
        sumx2 = 0;
        sumy2 = 0;
        for k=1:1:n
            sumxy = sumxy + x(k,c)*y(k);
            sumx2 = sumx2 + (x(k,c)^2);
            sumy2 = sumy2 + (y(k)^2);
        end
        sx = sqrt(sumx2 - n*(meanx^2));
        sy = sqrt(sumy2 - n*(meany^2));
        
        varx = (sumx2 - n*(meanx^2))/(n-1);
        vary = (sumy2 - n*(meany^2))/(n-1);
        
        stdy = std(y);
        stdx = std(x(:,c));
        
        Sxx  = varx*(n-1);
        
        %% UPOLOGISMOS APLHS GRAMMIKHS PALINDROMHSHS
        regressionModel = fitlm(x(:,c),y);
        C = table2array(regressionModel.Coefficients);
        b(:,c) = C(:,1) ;
        k1 = length(b(:,c));
        yhat(:,c) = [ones(length(x(:,c)),1) x(:,c)]*b(:,c);
        
        %% ERROR-SFALMA
        e = y - yhat(:,c)';
        
        %% TUPOPOIHMENO SFALMA PROSARMOGHS
        se2 = (1/(n-(k1+1)))*(sum((y-yhat(:,c)').^2));
        se = sqrt(se2);
        estar(:,c) = e ./ se;
        
        %% PARAMETRIKO DIASTHMA EMPISTOSUNHS GIA   bo   KAI  b1
        BCI = coefCI(regressionModel);
        
        %% PEARSON CORRELASION OF REGRESSION
        r(i,c) = b(2,c)*(sx/sy);
        
        %% R2 AND adjR2
        R2(i,c) = 1-(sum(e.^2))/(sum((y-mean(y)).^2));
        adjR2(i,c) =1-((n-1)/(n-(k1+1)))*(sum(e.^2))/(sum((y-mean(y)).^2));
        
        
        %% PLOT OF LINEAR REGRESSION FOR EVERY TIME DELAY [0,20]
        nfig = nfig+1;
        figure(nfig)
        
        subplot(2,1,1);
        scatter(x(:,c),y,5);
        hold on
        plot(x(:,c),b(1,c)+b(2,c)*x(:,c),'linewidth',1.5)
        ax = axis;
        text(ax(1)+0.7*(ax(2)-ax(1)),ax(3)+0.3*(ax(4)-ax(3)),['R^2(y)=',...
            num2str(R2(i,c),3)])
        text(ax(1)+0.7*(ax(2)-ax(1)),ax(3)+0.1*(ax(4)-ax(3)),['adjR^2(y)=',...
            num2str(adjR2(i,c),3)])
        xlabel( 'Daily Confirmed Cases')
        ylabel( 'Daily Deaths')
        title( sprintf(strcat(COUNTRY{i},' Deaths vs Confirmed Cases r=%1.6f'),r(i,c)))
        if b(2,c)<0
            modtxt = sprintf('y = %2.3f - %2.3f x',b(1,c),abs(b(2,c)));
        else
            modtxt = sprintf('y = %2.3f + %2.3f x',b(1,c),abs(b(2,c)));
        end
        text(ax(1)+0.05*(ax(2)-ax(1)),ax(4)-0.15*(ax(4)-ax(3)),modtxt)
        hold off;
        
        %% DIAGNWSTIKO DIAGRAMMA DIASPORAS TUPOPOIHMENWN SFALMATWN
        subplot(2,1,2);
        plot(y,estar(:,c),'.','linewidth',4)
        hold on
        ax = axis;
        plot([ax(1) ax(2)],zcrit*[1 1],'c--')
        plot([ax(1) ax(2)],-zcrit*[1 1],'c--')
        title(  sprintf(strcat(COUNTRY{i},' Regression standardised errors')))
        xlabel(  'Daily Deaths')
        ylabel(  'Residuals of Linear Regression')
        hold off;
        
        
        t(c)=T;
        c=c+1;
    end
    
    %% BRISKW THN KALUTERH PROSARMOGH TWN MODELWN PALINDROMHSHS MESW TOU r
    %% KALUTEROU r KAI STHN SUNEXEIA THA PLOTARW TO DIAGRAMMA DIASPORAS KAI
    %% THN EUTHEIA THS APLHS PALINDROMHSHS
    
    [~,Imax] = max(r(i,:));
    bestT(i) = t(Imax);
    bestB(i,:) =  b(:,Imax);
    bestX(i,:) = x(:,Imax);
    bestr(i) = r(i,Imax);
    bestR2(i) = R2(i,Imax);
    bestadjR2(i) = adjR2(i,Imax);
    bestYhat(i,:) = yhat(:,Imax);
    bestEstar(i,:) = estar(:,Imax);
    bestSumEstar(i) = sum(abs(bestEstar(i,:)));
    b = bestB(i,:);
    
    %% DIAGRAMMA DIASPORAS-EUTHEIAS
    nfig = nfig+1;
    figure(nfig)
    subplot(2,1,1);
    
    
    scatter( bestX(i,:),y,5);
    hold on
    plot( bestX(i,:),b(1)+b(2)* bestX(i,:),'k','linewidth',1.5)
    ax = axis;
    xlabel(  'Daily Confirmed Cases')
    ylabel(  'Daily Deaths')
    title(  sprintf(strcat(COUNTRY{i},' Best Model Describing Deaths - Confirmed Cases  | adjR^2=%1.4f'),bestadjR2(i)))
    if b(2)<0
        modtxt = sprintf('y = %2.5f - %2.3f x',b(1),abs(b(2)));
    else
        modtxt = sprintf('y = %2.5f + %2.3f x',b(1),abs(b(2)));
    end
    text(ax(1)+0.05*(ax(2)-ax(1)),ax(4)-0.15*(ax(4)-ax(3)),modtxt)
    hold off;
    
    %% DIAGNWSTIKO DIAGRAMMA DIASPORAS TUPOPOIHMENWN SFALMATWN
    subplot(2,1,2);
    scatter(y,bestEstar(i,:));
    hold on;
    plot(y,repmat(1.96,1,length(y)));
    hold on;
    plot(y,repmat(0,1,length(y)));
    hold on;
    plot(y,repmat(-1.96,1,length(y)));
    title(  sprintf(strcat(COUNTRY{i},' Regression standardised errors')))
    xlabel(  'Daily Deaths')
    ylabel( 'Residuals of Linear Regression')
    
    fprintf('\n');
    fprintf(' PRESS ANY KEY TO CONTINUE \n');
    pause;
end

%% PRINTS APOTELESMATWN
fprintf('\n');
fprintf(' ******************************************************************* \n');
fprintf(' COUNTRIES     :   SPAIN , BELGIUM , DENMARK , NETHERLANDS , NORWAY , SWITZERLAND \n');
fprintf(' TIME DELAYS   :  [%d       , %d        , %d        , %d        , %d      , %d  ] \n',bestT(1),bestT(2),bestT(3),bestT(4),bestT(5),bestT(6));
fprintf(' r^2(percent)  :  [%2.4f , %2.4f  , %2.4f   , %2.4f , %2.4f  , %2.4f] \n',(bestr(1)^2)*100,(bestr(2)^2)*100,(bestr(3)^2)*100,(bestr(4)^2)*100,(bestr(5)^2)*100,(bestr(6)^2)*100);
fprintf(' adjR^2        :  [%2.4f , %2.4f  , %2.4f   , %2.4f , %2.4f  , %2.4f] \n',bestadjR2(1),bestadjR2(2),bestadjR2(3),bestadjR2(4),bestadjR2(5),bestadjR2(6));
fprintf(' SUM(e(star))  :  [%5.2f , %5.2f  , %5.2f   , %5.2f , %5.2f  , %5.2f] \n',bestSumEstar(1),bestSumEstar(2),bestSumEstar(3),bestSumEstar(4),bestSumEstar(5),bestSumEstar(6));
fprintf(' ******************************************************************* \n');
fprintf('\n');


%% ********************** SXOLIASMOS APOTELESMATWN ************************
% FTANONTAS STO 5 ZHTHMA THS ERGASIAS PROSPATHOUME NA DIEREUNHSOUME THN
% DUNATOTHTA PROBLEPSHS TWN HMERHSIWN THANATWN APO HMERISIA KROUSMATA STO
% PRWTO KUMA MIAS XWRAS SUGKRINONTAS MODELA GRAMMIKHS PALINDROMHSH NEWN
% THANATWN SE MIA MERA APO NEA KROUSMATA ME XRONIKH USTERHSH SE ENA
% PARATHURO [ 0 , 20 ].

%% ERWTHSH 1)
% AFOU TREXOUME TO PROGRAMMA VLEPOUME PWS GIA KAPOIO T H PROSARMOGH DEN
% EINAI TO IDIO KALH GIA OLES TIS XWRES . AN LAVOUME UPOPSH TO adjR^2
% GIA TIS XWRES POU DIALEXAME KAI STO ZHTHMA 4  EXOUME TA EXHS APOTELESMATA
% XWRES  :   SPAIN  , BELGIUM , DENMARK  , NETHERLANDS , NORWAY  , SWITZERLAND
% adjR^2 :  [0.8373 , 0.8050  , 0.7116   , 0.8579      , 0.3589  , 0.5396]
% OPOTE EINAI FANERO PWS H ISPANIA TO BELGIO KAI H OLLANDIA EXOUN KALO
% SUNTELESTH PROSARMOGHS POLLAPLOU PROSDIORISMOU ENW KAI H DANIA EXEI
% IKANOPOIHTIKO APOTELESMA. OSTOSO H NORBIGIA KAI H ELVETIA EXOUN MIKRO
% SUNTELESTH PROSARMOGHS POLLAPLOU PROSDIORISMOU.

%% ERWTHSH 2)
% AN LAVOUME UPOSPH TON DIAGNWSTIKO ELEGXO  GIA KATHE USTERHSH (DIAGRAMA
% TUPOPOIHMENWN SFALMATWN) THA DOUME MESW TOU GRAFHMATOS PWS:

% H ISPANIA GIA MIKRES TIMES TOU y EXOUME SUGKEDRWSH KONTA STO 0
% ENW GIA MEGALES TIMES TOU y EXOUME APOMAKRUNSH APO TO 0 SUNEPWS H DIASPORA
% MEGALWNEI DEN EINAI STATHERH . ARA EXOUME SE MEGALO VATHMO MH
% OMOSKEDASTIKOTHTA KAI AUTO ODHGEI SE PROVLHMATA SUMPERASMATOLOGIAS
%(CI,ELEGXOS). EXOUME MIA APOMAKRH TIMH(OUTLIER) H OPOIA EPIREAZEI THN
% EKTIMHSH TOU MODELOU. TO MODELO EINAI DEN IKANOPOIHTIKO - EPARKES.

% TO BELGIO MIKRES TIMES TOU y EXOUME SUGKEDRWSH KONTA STO 0 ENW APO THN
% MESH KAI META TOU y EXOUME APOMAKRUNSH APO TO 0 SUNEPWS H DIASPORA
% MEGALWNEI KAI DEN EINAI STATHERH. DEN EXEI OMOSKEDASTIKOTHTA KAI AUTO ODHGEI
% SE PROVLHMATA SUMPERASMATOLOGIAS ALLA KAI EXOUME MERIKES APOMAKRES TIMES
% (OUTLIER) OI OPOIES EPIREAZOUN THN EKTIMHSH TOU MODELOU.
% TO MODELO DEN EINAI IKANOPOIHTIKO - EPARKES.

% H DANIA MAS PAROUSIAZEI ENA DIAGRAMA POU DEN EXEI OMOSKEDASTIKOTHTA. STHN
% ARXH MIAZEI NA EXEI STATHERH DIASPORA OMWS APO THN MESH KAI META ARXIZEI
% NA 'APLWNETAI' ME APOTELESMA H DIASPORA NA MHN EINAI STATHERH.
% ETSI ODHGOUMASTE SE PROVLHMATA SUMPERASMATOLOGIAS
% EXOUME KAI EDW MERIKES APOMAKRH APOMAKRES TIMES(OUTLIER) OI
% OPOIES EPIREAZOUN THN EKTIMHSH TOU MODELOU.
% TO MODELO DEN EINAI IKANOPOIHTIKO - EPARKES.

% H OLLANDIA GIA MIKRES TIMES TOU y EXOUME SUGKEDRWSH KONTA STO 0
% ENW GIA LIGO MEGALHTERES TIMES TOU y EXOUME APOMAKRUNSH APO TO 0 KAI
% STATHEROPOIHSH THN DIASPORAS APO EKEI KAI PERA.EXOUME IKANOPOIHTIKH
% OMOSKEDASTIKOTHTA KAI AUTO ODHGEI SE KALH SUMPERASMATOLOGIA . EPISHS
% EXOUME KAI MIA APOMAKRH TIMH(OUTLIER) H OPOIA EPIREAZEI THN EKTIMHSH TOU
% MODELOU. AKOMH PANW APO TO 95% TWN PARATHRHSEWN BRISKETAI MESA STO
% DIASTHMA [-1.96,1.96] . TO MODELO EINAI IKANOPOIHTIKO - EPARKES.

% H NORBIGIA MAS PAROUSIAZEI ENA DIAGRAMA POU DEN EXEI OMOSKEDASTIKOTHTA
% KATHWS STHN ARXH EXOUME TIMES STO ARNHTIKO KOMMATI KAI STHN SUNEXIA TIMES
% STO THETIKO KOMMATI KAI PANW. MIA MH GRAMMIKH PROSEGISH THA EFERNE KALUTERA
% APOTELESMATA .OLA AUTA ODHGOUN SE LANTHASMENH SUMPERASMATOLOGIA KATHWS
% H DIASPORA DEN EINAI STATHERH .KAI EDW EXOUME 3 OUTLIERS PRAGMA POU
% EPIREAZEI THN EKTIMHSH TOU MODELOU. AKOMH PANW APO TO 95% TWN PARATHRHSEWN
% BRISKETAI MESA STO DIASTHMA [-1.96,1.96]
% TO MODELO EINAI IKANOPOIHTIKO - EPARKES.

% H ELVETIA EXEI MH STATHERH DIASPORA H OPOIA ODHGEI SE MH OMOSKEDASTIKOTHTA
% KAI AUTO SUNEPAGETAI LANTHASMENH SUMPERASMATOLOGIA . EPISHS EXOUME KAI 3
% APOMAKRES TIMES(OUTLIERS) OI OPOIES EPIREAZOUN THN EKTIMHSH TOU
% MODELOU. AKOMH PANW APO TO 95% TWN PARATHRHSEWN BRISKETAI MESA STO
% DIASTHMA [-1.96,1.96]
% TO MODELO DEN EINAI IKANOPOIHTIKO - EPARKES.

%% ERWTHSH 3)
% SUMFWNA ME TON SUNTELESTH PROSARMOGHS POLLAPLOU PROSDIORISMOU adjR^2
% MPOROUME NA SUMPERANOUME AN UPARXEI SUSXETISH METAXU TWN NEWN KROUMATWN
% KAI THANATWN GIA MIA SUGKEKRIMENH MERA SUMFWNA ME MIA XRONIKH USTERISH T
% OSO PIO KONTA EINAI STHN MONADA O SUNTELESTHS TOSO PERISSOTERO H EKTIMHSH
% MAS PLHSIAZEI THN PRAGMATIKH TIMH . OPOTE AN TREXOUME TO PROGRAMMA
% EXOUME:  [0.8373 , 0.8050  , 0.7116   , 0.8579 , 0.3589  , 0.5396]
% EPOMENOS H ISPANIA ,BELGIO , OLLANDIA EXOUN MEGALO SUNTELESTH adjR^2 KAI
% SUNEXIZONTAS H DANIA EXEI ENAN IKANOPOIHTIKO SUNTELESTH adjR^2 . OI XWRES
% NORBIGIA KAI ELVETIA EXOUN XAMHLO SUNTELESTH POU ODHGEI SE KAKH EKTIMHSH
% MELLONTIKWN TIMWN TWN THANATWN.

%% ERWTHSH 4)
% OI BELTISTES YSTERHSEIS KATHE ZHTHMATOS EINAI H EXHS GIA KATHE XWRA:
% XWRES          :  SPAIN  , BELGIUM , DENMARK  , NETHERLANDS , NORWAY , SWITZERLAND
% 3) TIME DELAYS : [6      , 5       , 4        , 2           , 14     , 9  ]
% 4) TIME DELAYS : [6      , 5       , 0        , 5           , 11     , 6  ]
% 5) TIME DELAYS : [6      , 5       , 0        , 5           , 11     , 6  ]
% EINAI ANAMENOMENO TO ZHTHMA 4 KAI 5 NA EXOUN TIS IDIES XRONIKES
% USTERHSEIS AFOU EXETAZOUME TIS IDIES XWRES KAI O SUNTELESTHS r(Pearson)
% PAIRNEI THN MEGALUTERH TIMH STHN XRONIKH USTERISH OPOU KAI O SUNTELESTHS
% KALHS PROSARGMOGHS PAIRNEI THN MEGALUTERH TOU TIMH. OSON AFORA TO ZHTHMA 3
% EXOUME MIKRES APOKLHSEIS STIS XRONIKES USTERISEIS SE MERIKES XWRES.
