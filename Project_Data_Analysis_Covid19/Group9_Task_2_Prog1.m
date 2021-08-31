%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 2
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
spad = DATADeath(131+g  ,2:end);

%% AUSTRIA 9
popAus = DATAConfirmed(9+g  ,1);
ausc = DATAConfirmed(9+g  ,2:end);
ausd = DATADeath(9 +g ,2:end);

%% BELGIUM 14
popBel = DATAConfirmed(14+g ,1);
belc = DATAConfirmed(14+g  ,2:end);
beld = DATADeath(14+g ,2:end);

%% CROATIA 34
popCro = DATAConfirmed(34+g  ,1);
croc = DATAConfirmed(34 +g ,2:end);
crod = DATADeath(34+g  ,2:end);

%% DENMARK 38
popDen = DATAConfirmed(38+g ,1);
denc = DATAConfirmed(38 +g ,2:end);
dend = DATADeath(38+g ,2:end);

%% LITHUANIA 83
popLit = DATAConfirmed(83+g  ,1);
litc = DATAConfirmed(83+g ,2:end);
litd = DATADeath(83+g  ,2:end);

%% ESTONIA 45
popEst = DATAConfirmed(45+g ,1);
estc = DATAConfirmed(45+g ,2:end);
estd = DATADeath(45+g  ,2:end);

%% FINNLAND 48
popFil = DATAConfirmed(48+g ,1);
filc = DATAConfirmed(48 +g ,2:end);
fild = DATADeath(48+g  ,2:end);

%% NETHERLANDS 98
popOll = DATAConfirmed(98 +g ,1);
ollc = DATAConfirmed(98+g,2  :end);
olld = DATADeath(98+g ,2:end);

%% NORWAY 104
popNor = DATAConfirmed(104+g  ,1);
norc = DATAConfirmed(104+g ,2:end);
nord = DATADeath(104+g ,2:end);

%% SWITZERLAND 135
popElv = DATAConfirmed(135+g  ,1);
elvc = DATAConfirmed(135+g  ,2:end);
elvd = DATADeath(135 +g ,2:end);


startc = [ 57  ; 63 ; 61  ; 76  ;  66 ; 74  ; 71  ; 66  ; 62  ; 62  ; 60];
endc   = [ 154 ;146 ; 186 ; 138 ; 186 ; 140 ; 123 ; 156 ; 179 ; 138 ; 130];

startd= [ 66 ; 81 ; 73  ;  89 ; 75  ; 83  ; 87  ; 83  ; 74  ; 74 ; 71 ];
endd = [ 158 ;157 ; 206 ; 147 ; 194 ; 155 ; 131 ; 156 ; 190 ; 138; 144];

bounds = [60 ,158 ; 80 ,155 ; 65 , 206; 72 ,137; 70 , 194; 78 , 153; 71 , 135; 65 , 156; 68 , 194 ; 68 , 138 ; 66 , 144];

CC = [ spac ; ausc ; belc ; croc ; denc ; litc ; estc ; filc ; ollc ; norc ; elvc];
D  =  [spad ;ausd ; beld ; crod ; dend ; litd ; estd ; fild ; olld ; nord ; elvd];
COUNTRY = {'SPAIN','AUSTRIA',' BELGIUM', 'CROATIA', 'DENMARK' ,'LITHUANIA'...
    ' ESTONIA', 'FINNLAND', 'NETHERLANDS', 'NORWAY',' SWITZERLAND'};

distnames = {'Normal','LogLogistic','Logistic','Poisson',...
    'Rayleigh','tLocationScale','Exponential','ExtremeValue','Gamma'};

N=size(CC,1);
cc = zeros(size(CC));
d = zeros(size(D));
indexc = zeros(N);
indexd = zeros(N);
MSEd = zeros(N,1);
MSEc = zeros(N,1);
RMSEd = zeros(N,1);
RMSEc = zeros(N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    
    %% DATA SMOOTHING
    [cc(i,:) , d(i,:)] =  Group9Exe1Fun1(CC(i,:),D(i,:));
    
    % MERIKA PLOTS ,STHN OUSIA EINAI TA RABDOGRAMATA GIA KATHE SENARIO
    % WSTE NA DOUME GRAFIKA PRWTA TA SHMEIA POU KATHORIZOUN THN ARXH KAI
    % TO TELOS TOUS PRWTOU KUMATOS AFOU EXOUME PROEPEXERGASTEI TA DEDOMENA
    
    
    nfig = nfig + 1;
    figure(nfig)
    plot(cc(i,:),'black')
    title(strcat(COUNTRY{i},' Covid 19 Confirmed Case from 1/1/12 to 13/12/2020'))
    ylabel('$Confirmed Cases$','Interpreter','latex','fontsize',10);
    xlabel('$Days$','Interpreter','latex','fontsize',10);
    
    nfig = nfig + 1;
    figure(nfig)
    plot(d(i,:),'red')
    title(strcat(COUNTRY{i},' Covid 19 Deaths from 1/1/12 to 13/12/2020'))
    ylabel('$Deaths$','Interpreter','latex','fontsize',10);
    xlabel('$Days$','Interpreter','latex','fontsize',10);
    
    fc = cc(i,startc(i):endc(i));
    fd = d(i,startc(i):endd(i));
    durationc = abs(startc(i) - endc(i));
    durationd = abs(startd(i) - endd(i));
    
    n = length(fc);
    fc = fc';
    fd = fd';
    
    
    xc = 1:length(fc);
    xc = xc';
    xd = 1:length(fd);
    xd = xd';
    
    LogLogisticc = fitdist(xc,'LogLogistic','Frequency',fc);
    LogLogisticd = fitdist(xd,'LogLogistic','Frequency',fd);
    
    yc = pdf(LogLogisticc,xc);
    yc = yc*sum(fc);
    errorc = ((fc - yc).^2)/length(fc);
    MSEc(i)= sum(errorc);
    ec = fc - yc;
    
    yd = pdf(LogLogisticd,xd);
    yd = yd*sum(fd);
    errord = ((fd - yd).^2)/length(fd);
    MSEd(i)= sum(errord);
    ed = fd - yd;
    
    
    RMSEc(i) = sqrt(MSEc(i));
    RMSEd(i) = sqrt(MSEd(i));
    
    %% BAR CHARTS FOR THE BEST DATA DISTRIBUTION FIT FOR CONFIRMED CASES
    nfig = nfig + 1;
    figure(nfig)
    bar(fc)
    hold on
    plot(yc,'LineWidth',2)
    title(strcat(COUNTRY{i},' LogLogisticc Distribution - Covid 19 Confirmed Cases'))
    ylabel('$Deaths$','Interpreter','latex','fontsize',10)
    xlabel('$Days$','Interpreter','latex','fontsize',10)
    
    %% BAR CHARTS FOR THE BEST DATA DISTRIBUTION FIT FOR CONFIRMED DEATHS
    nfig = nfig + 1;
    figure(nfig)
    bar(fd)
    hold on
    plot(yd,'LineWidth',2)
    title(strcat(COUNTRY{i},' LogLogisticc Distribution - Covid 19 Deaths'))
    ylabel('$Deaths$','Interpreter','latex','fontsize',10)
    xlabel('$Days$','Interpreter','latex','fontsize',10)
    
    fprintf('\n');
    fprintf(' *******************************************************************\n');
    fprintf(' TIME PERIOD 1st WAVE COVID 19 \n');
    fprintf(' COUNTRY : %s \n',COUNTRY{i});
    fprintf(' TIME PERIOD FOR NEW CASES [ %d , %d ] \n',startc(i),endc(i));
    fprintf(' DURATION : %d DAYS \n',durationc);
    fprintf(' TIME PERIOD FOR DEATHS [ %d , %d ] \n',startd(i),endd(i));
    fprintf(' DURATION : %d DAYS \n',durationd);
    fprintf(' *******************************************************************\n');
    fprintf('\n');
    fprintf(' *******************************************************************\n');
    fprintf(' BEST FIT CONTROL BY MSE FOR DAILY CASES %s \n',COUNTRY{i});
    fprintf(' LogLogisticc DISTRIBUTION MSE  CONFIRMED CASES: %.6f \n',MSEc(i));
    fprintf(' LogLogisticc DISTRIBUTION RMSE CONFIRMED CASES : %.6f \n',RMSEc(i));
    fprintf(' *******************************************************************\n');
    fprintf(' BEST FIT CONTROL BY MSE FOR DAILY DEATHS %s \n',COUNTRY{i});
    fprintf(' LogLogisticc DISTRIBUTION MSE  DEATHS   : %.6f \n',MSEd(i));
    fprintf(' LogLogisticc DISTRIBUTION RMSE DEATHS   : %.6f \n',RMSEd(i));
    fprintf(' *******************************************************************\n');
    fprintf(' \n')
    fprintf(' PRESS ANY KEY TO CONTINUE \n');
    
    pause;
end
clc
%% TAXINOMHSH TWN XWRWN ANALOGA TO RMSE SE AUXOUSA SEIRA
[~,Ic]=sort(RMSEc);
[~,Id]=sort(RMSEd);
fprintf('TAXINOMHSH TWN XWRWN OS PROS TA EPIVEVAIWMENA KROUSMATA :');
fprintf('%s ',COUNTRY{Ic});
fprintf(' \n')
fprintf('TAXINOMHSH TWN XWRWN OS PROS TOUS THANATOUS :');
fprintf('%s ',COUNTRY{Id});
fprintf(' \n')
%% ********************** SXOLIASMOS APOTELESMATWN ************************
% H LOG LOGISTIC KATANOMH POU VRHKAME OTI EFARMOZEI SE IKANOPOIHTIKO VATHMO
% STA DEDOMENA COVID 19 THS ISPANIAS EFARMOSTHKE SE AUTO TO ZHTHMA SE 10
% EUROPAIKES XWRES. YSTERA APO OLA TA PLOTS POU EIDAME MESW TIS EPANALHPTIKHS
% DIADIKASIAS PARAPANW , SUMPERANAME PWS H LOGLOGISTIC KATANOMH PROSARMOZETAI SE
% SE MEGALO VATHMO STA DEDOMENA TWN 10 EUROPAIKWN XWRWN. AN VALOUME 4 PITHANES
% TIMES    [ OXI KALA      ,    METRIA    ,        KALA     ,   POLU KALA]
% DIASTHMA [ 0-0.3         , 0.3-0.6      , 0.6 - 0.8       , 0.8-1]
% EXOUME TA EXHS APOTELESMATA GIA TA  NEA KROUSMATA  |  THANATOUS
% NEA KROUSMATA : [ 2 , 1 , 4 , 4 ] = KALH   PROSARMOGH THS LOGLOGISTIC
% THANATOI      : [ 2 , 4 , 1 , 4 ] = MESAIA PROSARMOGH THS LOGLOGISTIC
% TELOS H KALUTERH PROSARMOGH KATA AUXOUSA SEIRA GIA TA NEA KROUSMATA
% [ CROATIA , LITHUANIA , ESTHONIA  ,FILLANDIA, NORBIGIA  DANIA , AUSTRIA...
%   OLLANDIA , ELVETIA , BELGIO, ISPANIA ]
% KAI H KALUTERH PROSARMOGH KATA AUXOUSA SEIRA GIA TOUS THANATOUS
% [ ESTHONIA , LITHUANIA ,CROATIA ,DANIA, NORBIGIA ,FILLANDIA , AUSTRIA...
%   BELGIO , ELVETIA , OLLANDIA, ISPANIA ]
% POLU KALH PROSARMOGH EXEI EPISHS H GAMMA KATANOMH STIS PARAPANW XWRES
% OPOU HTAN H DEUTERH KALUTERH KATANOMH POU PROSARMOZONTAN STA DEDOMENA
% THS ISPANIAS!!!

% XWRA OMADAS : ISPANIA 131
% XWRES EU : AUSTRIA 9, BELGIO 14, CROATIA 34, DANIA 38 ,LITHOUANIA 83
% XWRES EU : ESTONIA 45, FILANDIA 48, OLLANDIA 98, NORBIGIA 104, ELVETIA 135
