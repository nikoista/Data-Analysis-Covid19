%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA  6
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

%% FIRST WAVE BOUNDS
bounds = [60, 158; 65 , 206; 70 , 194; 68 , 194 ; 68 , 138 ; 66 , 144];

CC = [spac ;  belc ;  denc ;  ollc ; norc ; elvc];
D =  [spad ;  beld ;  dend ;  olld ; nord ; elvd];
COUNTRY = {'SPAIN',' BELGIUM', 'DENMARK' ,'NETHERLANDS', 'NORWAY','SWITZERLAND'};

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
R2wise = zeros(N,B2+1);
R2s = zeros(N,B2+1);
adjR2wise = zeros(N,B2+1);
adjR2s = zeros(N,B2+1);
bestR2 =  zeros(N,3);
bestadjR2 =  zeros(N,3);
paramStepwise =  zeros(N,1);


%% APO TO ZHTHMA 5 H KALUTERES XRONIKES USTERHSEIS
bestDelays = [6 ,5,0,5,11,6] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    i
    
    %% DATA SMOOTHING
    duration(i) = (endd(i) - endc(i));
    [cc(i,:) , d(i,:)] =  Group9Exe1Fun1(CC(i,:),D(i,:));
    
    
    
    %% 1 KUMA BOUNDS GIA THANATOUS KAI NEA KROUSMATA
    x = zeros(B2+1,length(cc(i,(bounds(i,1)):bounds(i,2))));
    y = d(i,bounds(i,1):(bounds(i,2)));
    xs = cc(i,(bounds(i,1)-bestDelays(i)):bounds(i,2)-bestDelays(i));
    
    %% EPANALHPTIKH DIADIKASIA EURESHS XRONIKWN USTERISEWN
    c=1;
    for T=0:1:B2
        x(c,:) = cc(i,(bounds(i,1)-T):bounds(i,2)-T);
        c=c+1;
    end
    xs =xs';
    x = x';
    y = y';
    n = length(y);
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UPOLOGISMOS APLHS GRAMMIKHS PALINDROMHSHS
    simpleModel = fitlm(xs,y);
    C = table2array(simpleModel.Coefficients);
    bs = C(:,1) ;
    yhats = [ones(length(xs),1) xs]*bs;
    
    %% ERROR-SFALMA
    es = y - yhats;
    
    %% BATHMOS k
    ks = length(bs);
    
    %% TUPOPOIHMENO SFALMA PROSARMOGHS
    se2s = (1/(n-(ks+1)))*(sum((y-yhats).^2));
    ses = sqrt(se2s);
    estars = es ./ ses;
    
    %% SUNTELESTES R^2 , adjR^-2
    R2s = 1-(sum(es.^2))/(sum((y-mean(y)).^2));
    adjR2s =1-((n-1)/(n-(ks+1)))*(sum(es.^2))/(sum((y-mean(y)).^2));
    
    %% DIAGNWSTIKO DIAGRAMA TUPOPOIHMENWN SFALMATWN FULL REGRESSION
    nfig = nfig+1;
    figure(nfig)
    subplot(3,1,1);
    scatter(y,estars);
    hold on;
    plot(y,repmat(1.96,1,length(y)));
    hold on;
    plot(y,repmat(0,1,length(y)));
    hold on;
    plot(y,repmat(-1.96,1,length(y)));
    title(  sprintf(strcat(COUNTRY{i},' Simple Model Regression standardised errors | adjR^2 = %2.4f'),adjR2s))
    xlabel(  'Daily Deaths')
    ylabel(  'Residuals of Simple Model')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% FULL REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REGRESSION FULL MODEL
    REG = fitlm(x,y);
    ball = table2array(REG.Coefficients);
    ball = ball(:,1);
    xall = [ones(length(x),1) x];
    yhatall = xall*ball;
    eall = y-yhatall;
    
    %% BATHMOS k
    kall = length(ball)-1;
    
    %% TUPOPOIHMENO SFALMA
    se2all = (1/(n-(kall+1)))*(sum(eall.^2));
    seall = sqrt(se2all);
    estarall = eall / seall;
    
    %% SUNTELESTES R^2 , adjR^2
    R2all = 1-(sum(eall.^2))/(sum((y-mean(y)).^2));
    adjR2all =1-((n-1)/(n-(kall+1)))*(sum(eall.^2))/(sum((y-mean(y)).^2));
    
    %% DIAGNWSTIKO DIAGRAMA TUPOPOIHMENWN SFALMATWN FULL REGRESSION
    subplot(3,1,2);
    scatter(y,estarall);
    hold on;
    plot(y,repmat(1.96,1,length(y)));
    hold on;
    plot(y,repmat(0,1,length(y)));
    hold on;
    plot(y,repmat(-1.96,1,length(y)));
    title(   sprintf(strcat(COUNTRY{i},' Full Model Regression standardised errors | adjR^2 = %2.4f'),adjR2all))
    xlabel(  'Daily Deaths')
    ylabel(  'Residuals of Full Model')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%% STEPWISE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEPWISE REGRESSION
    [b,~,~,model,stats] = stepwisefit(x,y);
    b0 = stats.intercept;
    bwise = [b0; b(model)];
    xwise = [ones(length(x),1) x(:,model)];
    yhatwise = xwise*bwise;
    ewise = y - yhatwise;
    
    %% BATHMOS k
    k1 = length(bwise);
    
    %% TUPOPOIHMENO SFALMA
    se2wise = (1/(n-(k1+1)))*(sum(ewise.^2));
    sewise = sqrt(se2wise);
    estarwise = ewise / sewise;
    
    %% SUNTELESTES R^2 , adjR^2
    R2wise = 1-(sum(ewise.^2))/(sum((y-mean(y)).^2));
    adjR2wise =1-((n-1)/(n-(k1+1)))*(sum(ewise.^2))/(sum((y-mean(y)).^2));
    
    %% DIAGNWSTIKO DIAGRAMA TUPOPOIHMENWN SFALMATWN  STEPWISE MODEL
    subplot(3,1,3);
    scatter(y,estarwise);
    hold on;
    plot(y,repmat(1.96,1,length(y)));
    hold on;
    plot(y,repmat(0,1,length(y)));
    hold on;
    plot(y,repmat(-1.96,1,length(y)));
    title(  sprintf(strcat(COUNTRY{i},' Stepwise Regression standardised errors | adjR^2 = %2.4f'),adjR2wise))
    xlabel( 'Daily Deaths')
    ylabel( 'Residuals of Stepwise')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% EVRESH KALUTEROU SUNTELESTH
    bestR2(i,:) = [R2s R2all R2wise];
    bestadjR2(i,:) = [adjR2s adjR2all adjR2wise];
    
    %% UPOLOGISMOS PARAMETRWN POU XRHSIMOPOIHTHIKAN GIA REGRESSION
    paramStepwise(i) = length(bwise);
    paramFullReg = length(ball);
    
    clc;
    fprintf('\n');
    fprintf(' PRESS ANY KEY TO CONTINUE \n');
    
    pause;
end
%% PRINTS APOTELESMATWN
fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' COUNTRIES            :   SPAIN , BELGIUM , DENMARK , NETHERLANDS , NORWAY , SWITZERLAND \n');
fprintf(' adjR^2 SIMPLE MODEL  :  [%2.4f , %2.4f  , %2.4f   , %2.4f , %2.4f  , %2.4f] \n',bestadjR2(1,1),bestadjR2(2,1),bestadjR2(3,1),bestadjR2(4,1),bestadjR2(5,1),bestadjR2(6,1));
fprintf(' adjR^2 FULL MODEL    :  [%2.4f , %2.4f  , %2.4f   , %2.4f , %2.4f  , %2.4f] \n',bestadjR2(1,2),bestadjR2(2,2),bestadjR2(3,2),bestadjR2(4,2),bestadjR2(5,2),bestadjR2(6,2));
fprintf(' adjR^2 STEPWISE      :  [%2.4f , %2.4f  , %2.4f   , %2.4f , %2.4f  , %2.4f] \n',bestadjR2(2,3),bestadjR2(2,3),bestadjR2(3,3),bestadjR2(4,3),bestadjR2(5,3),bestadjR2(6,3));
fprintf(' FULLMODEL parameters :  [%d      , %d       , %d       , %d       , %d     , %d ]\n',paramFullReg,paramFullReg,paramFullReg,paramFullReg,paramFullReg,paramFullReg);
fprintf(' STEPWISE  parameters :  [%d       , %d        , %d        , %d        , %d      , %d  ]\n',paramStepwise(1),paramStepwise(2),paramStepwise(3),paramStepwise(4),paramStepwise(5),paramStepwise(6));
fprintf(' *******************************************************************\n');
fprintf('\n');


%
%  *******************************************************************
%  COUNTRIES     :   SPAIN , BELGIUM , DENMARK , NETHERLANDS , NORWAY , SWITZERLAND
%  TIME DELAYS   :  [6       , 5        , 0        , 5        , 11      , 6  ]
%  r^2(percent)  :  [84.0611 , 80.7803  , 71.6224   , 86.0118 , 37.7241  , 55.1450]
%  adjR^2        :  [0.8390 , 0.8064  , 0.7139   , 0.8590 , 0.3682  , 0.5456]
%  SUM(e(star))  :  [69.15 , 89.60  , 86.29   , 87.19 , 48.20  , 58.95]
%  *******************************************************************

%% ********************** SXOLIASMOS APOTELESMATWN ************************
% STO ZHTHMA 6 STO OPOIO BRISKOMASTE THELOUME NA DIEREUNHSOUME AN H
% PROVLEPSH TWN HMERISIWN THANATWN STO PRWTO KUMA BELTIWNETAI AN ANTI TOU
% MONTELOU APHS GRAMMIKHS PALINDROMHSHS WS PROS TA HMERHSIA KROUSMATA THS
% IDIAS H KAPOIAS PROHGOUMENHS HMERAS , OPOTE THA XRHSIMOPOIHSOUME MODELO
% POLLAPLHS GRAMMIKHS PALINDROMHSHS HMERHSIWN KROUSMATWN SE PERISSOTERES
% APO MIA USTERHSEIS . TO MODELO POU EPILEXAME GIA TO ZHTHMA AUTO EINAI TO
% STEP WISE REGRESSION MODEL DIOTI EXEI POLU IKNOPOIHTIKA APOTELESMATA KAI
% O SUNTELESTHS PROSARMOGHS POLLAPLOU PROSDIORISMOU POU VRHKAME STO ZHTHMA
% 5 EINAI POLU KONTA SE AUTO TOU FULL MODEL POLLAPLHS GRAMMIKHS
% PALINDROMHSHS.

%% ERWTHSH 1)
% KATARXAS H PROSARMOGH POIKILH APO XWRA SE XWRA . SE MERIKES EINAI UPSHLH
% ENW SE MERIKES ALLES EINAI ARKEtA XAMHLH. THA PAROUSIASOUME TOUS
% SUNTELESTES  adjR^2 GIA KATHE MIA APO TIS 3 METHODOUS PALINDROMHSHS :
%  COUNTRIES     :    SPAIN , BELGIUM , DENMARK  , NETHERLANDS , NORWAY , SWITZERLAND
% 1) APLH GRAMM PALIND [0.8373 , 0.8050  , 0.7116   , 0.8579 , 0.3589  , 0.5396]
% 2) FULL GRAMM PALIND [0.9123 , 0.9741  , 0.8028   , 0.9111 , 0.3662  , 0.6818]
% 3) STEP WISE  PALIND [0.9736 , 0.9736  , 0.8084   , 0.9015 , 0.4547  , 0.6713]
% OPWS VLEPOUME H PROSARARMOGH BELTISTOU MODELOU GIA KATHE XWRA EINAI:
% SPAIN :       2 > 3 > 1
% BELGIUM :     2 > 3 > 1
% DENMARK :     3 > 2 > 1
% NETHERLANDS : 2 > 3 > 1
% NORWAY :      3 > 2 > 1
% SWITZERLAND : 2 > 3 > 1
% EINAI FANERO PWS H PROSARARMOGH BELTISTOU MODELOU GIA KATHE XWRA DEN
% EINAI H IDIA . EPISHS VLEPOUME OTI TO STEPWISE MODELO DINEI TA KALUTERA
% APOTELESMATA MERIKES FORES SE SUGKRISH ME TO FULL MODELO ME TIS 21 METAVLHTES
% TELOS TO APLO MODELO GRAMMIKHS PALINDROMHSHS EINAI IKANOPOIHTIKO SE MERIKES
% XWRES NA LAVOUME UPOSPH TO adjR^2.

%% ERWTHSH 2)
% EINAI PROFANES PWS OTAN XRHSIMOPOIEITAI MONTELO POLLAPLHS GRAMMIKHS
% PALINDROMHSHS BELTIWNETAI O SUNTELESTHS POLLAPLOU PROSDIORISMOU KAI
% SUNEPWS H IKANOTHTA TOU MODELOU STHN PROVLEPSH HMERHSIWN THANATWN.
% EPISHS AN PARATHRHSOUME TA DIAGRAMMATA TUPOPOIHMENWN SFALMATWN THA DOUME
% OTI H DIASPORA PLEON EINAI PIO STATHERH SE SXESH ME THN APLH GRAMMIKH
% PALINDROMHSH PRAGMA POU VOHTHAEI STHN DHMIOURGIA OMOSKEDASTIKOTHTAS TOU
% MONTELOU. AKOMA KAI ME THN KALUTERH PROSARMOGH TOU MODELOU POLLAPLHS
% PALINDROMHSHS EINAI ESTHHTO PWS OI APOMAKRISMENES TIMES (OUTLIERS) DEN
% DIORTHONWNTAI , KATHWS PARAMENOUN KAI STIS 2 PERIPTWSEIS POLLAPLHS
% GRAMMIKHS PALINDROMHSHS.

%%  PROVLIMATISMOI
% VLEPOUME OTI SE MERIKES XWRES TO STEPWISE MONTELO EXEI adjR^2 MEGALHTERO
% APO OTI TO MONTELO POLLAPLHS PALINDROMHSHS ME 21 METABLHTES . EINAI
% PERIERGO TO GEGONOS OTI ENA MONTELO ME EPILOGH MERIKWN METABLHTWN GIA
% PALINDROMHSH MPOREI NA EINAI KALUTERO APO ENA POLLAPLHS ME OLES TIS
% MERAVLHTES . TO MONO POU MPOROUME NA SKEFTOUME GIA AUTH THN PERIPTWSH
% EINAI OTI APO TON TUPO TOU adjR^2 H AUXHSH TWN BATHMWN k MIKRENEI TON
% PARANOMASTH ARA TO KLASMA (n-1)/(n-(k+1)) AUXANETAI POU SHMAINEI OTI THA
% EPREPE OLOS AUTOS O OROS ((n-1)/(n-(k+1)))*ssRes/ssTot NA AUXANETAI
% ESTW OTI THETW X ----> ((n-1)/(n-(k+1)))*ssRes/ssTot . STHN IDANIKH
% PERIPTWSH OPOU TO ssRes/ssTot EINAI IDIO STO STEPWISE KAI FULL MONTELO
% GRAMMIKHS PALINDROMHSHS TOTE :
% STEPWISE GRAMMIKH PALINDROMHSH       adjR^2(1) = 1  - A  OPOU X>A
% POLLAPLH GRAMMIKH PALINDROMHSH       adjR^2(2) = 1  - X  OPOU X>A
% OPOTE   adjR^2(1) >  adjR^2(2). MIA PROSPATHIA EXHGHSHS ME MATHIMATIKA
% TIPOTA PARAPANW .