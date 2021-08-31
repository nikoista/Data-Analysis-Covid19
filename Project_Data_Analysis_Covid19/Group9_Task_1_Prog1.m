%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS

%% EISAGWGH TWN DEDOMENWN GIA TON COVID 19 1/1/2020 --- 13/12/2020
%% XWRA OMADAS : ISPANIA 131
DATAConfirmed = importdata('Covid19Confirmed.xlsx','headerlinesIn');
DATAConfirmed = DATAConfirmed.data;

DATADeath = importdata('Covid19Deaths.xlsx','headerlinesIn');
DATADeath = DATADeath.data;

% EPEIDH MAS PAROUSIAZONTAN PROVLHMA ME TON TROPO ME TON OPOIO GINOTAN
% IMPORT TA DEDMENA STO PROGRAMMA ( TO INDEX TIS ISPANIAS EINAI 131 STON
% UPOLOGISTH MOU ALLA SE ALLO UPOLOGISTH MPOREI NA EINAI 130 ) DHMIOURGHSA
% ENA TROPO ME TON OPOIO ELEGXO AN O PLHTHISMOS THS XWRAS POU XOUME SAN 
% OMADA (ISPANIA ME PLHTHISMO = 46937060 EISAI SWSTOS
g=0;
population = DATAConfirmed(131+g,1);
if(population == 46937060)
    g = 0;
else 
    g=-1;
end

population = DATAConfirmed(131+g,1);
dc = DATAConfirmed(131+g,2:end);
dd = DATADeath(131+g,2:end);

nfig = 0;

% MERIKA PLOTS ,STHN OUSIA EINAI TA RABDOGRAMATA GIA KATHE SENARIO
% WSTE NA DOUME ARXIKA SE GRAFIKH APEIKONISH TA SHMEIA POU KATHORIZOUN
% THN ARXH KAI TO TELOS TOUS PRWTOU KUMATOS XWRIS NA EXOUME "KATHARISEI"
% AKOMA TO DATASET KAI XWRIS KANONIKOPOIHSH

nfig = nfig + 1;
figure(nfig)
plot(dc,'black');
title('Spain Covid 19 Confirmed Cases from 1/1/12 to 13/12/2020');
ylabel('$Confirmed Cases$','Interpreter','latex','fontsize',10);
xlabel('$Days$','Interpreter','latex','fontsize',10);

nfig = nfig + 1;
figure(nfig)
plot(dd,'red');
title('Spain Covid 19 Deaths from 1/1/12 to 13/12/2020');
ylabel('$Deaths$','Interpreter','latex','fontsize',10);
xlabel('$Days$','Interpreter','latex','fontsize',10);

old_dc = dc;
old_dd = dd;

%% START-END OF COVID 19 CONFIRMED CASES AND DEATHS
startc = 57;
endc = 157;
startd = 57;
endd = 157;

% OI KAMPULESPOU DHMIOURGOUDE GIA TA KROUSMATA KAI THANATOUS EXOUN
% EXOKEIMENES TIMES OPOTE EXOUME MIA KAMPULH SAN NA EXEI THORUVO. GIA NA
% ANTIMETOPOIHSOUME TO PROBLHMA MPOROUME NA  EFARMOSOUME ENA
% FILTRO - EURWSTO ALGORITHMO GIA NA KANOUME SMOOTHING TA DEDOMENA MAS
% XRHSIMOPOIONTAS smoothdata(...) ME variables->(dataset,MovingMeanMethod,Window)
% OSTOSO THA XRHSIMOPOIHSOUME MIA METHODO POU XWRIZEI TIS SUXNOTHTES
% SE 3 MERES ANTI GIA MIA OPWS PX STHN ISPANIA POU THN DEUTERA EXW TA
% KROUSMATA TOU SAVVATOKURIAKOU ,  ARA ANAKATATASOUME TA DEDOMENA
% KAI STIS 3 MERES ME THN METHODO POU EXOUME STO ARXEIO Group9Exe1Fun1.
% dc = smoothdata(dc,'gaussian');
% dd = smoothdata(dd,'gaussian');


% META APO THN Group9Exe1Fun1 EXOUME DIWRTHWSEI TIS ARNHTIKES TIMES
% KAI THN METAVLHTOTHTA TWN TIMWN APO 0 SE KAPOIO MEGALO ARITHMO POU
% SUMVAINEI STA DEDOMENA MAS
t=1:length(dc);
[dc , dd] = Group9Exe1Fun1(dc,dd);



%% MERIKA PLOTS GIA THN XWRA THS OMADAS 
nfig = nfig + 1;
figure(nfig)
plot(dc);
title('Smoothed Data Confirmed Cases');
ylabel('$Smoothed Data Confirmed Cases$','Interpreter','latex','fontsize',10);
xlabel('$Days$','Interpreter','latex','fontsize',10);

nfig = nfig + 1;
figure(nfig)
plot(dd);
title('Smoothed Data Deaths');
ylabel('$Smoothed Data Deaths$','Interpreter','latex','fontsize',10);
xlabel('$Days$','Interpreter','latex','fontsize',10);

nfig = nfig + 1;
figure(nfig)
plot(t,dc,'--',t,old_dc);
title('Smoothed Confirmed Cases - Noise Data Confirmed Cases')
ylabel('$Smoothed Data Confirmed Cases - Noise Data Confirmed Cases$'...
    ,'Interpreter','latex','fontsize',10)
xlabel('$Days$','Interpreter','latex','fontsize',10)
legend('Smoothed Data','Data with Noise')

nfig = nfig + 1;
figure(nfig)
plot(t,dd,'--',t,old_dd);
title('Smoothed Data Deaths - Noise Data Deaths')
ylabel('$Smoothed Data Deaths -Noise Data Deaths$'...
    ,'Interpreter','latex','fontsize',10)
xlabel('$Days$','Interpreter','latex','fontsize',10)
legend('Smoothed Data','Data with Noise')

%% BOXPLOT FOR CONFIRMED CASES AND DEATHS
nfig = nfig + 1;
figure(nfig)
boxplot(dc)
title('Box Plot Covid 19 Confirmed Cases')
ylabel('$Confirmed Cases $','Interpreter','latex','fontsize',10)
xlabel('$Counts$','Interpreter','latex','fontsize',10)

nfig = nfig + 1;
figure(nfig)
boxplot(dd)
title('Box Plot Covid 19 Deaths')
ylabel('$Deaths $','Interpreter','latex','fontsize',10)
xlabel('$Counts$','Interpreter','latex','fontsize',10)


% PARAMETRIKH KATANOMH PITHANOTHTAS ME THN KALYTERH PROSARMOGH STA
% DEDOMENA HMERISIWN KROUSMATWN TOU PRWTOU KUMATOS GIA THN XWRA = ISPANIA
% SUXNOTHTES EMFANISHS NEWN KROUSMATWN GIA KATHE MERA GIA TO PRWTO KUMA

fc = dc(startc:endc);
fd = dd(startd:endd);

%% PARAMETRIKES KATANOMES PITHANOTHTAS 8 SE ARITHMO GIA NEA KROUSMATA
xc = 1:length(fc);
xc = xc';

normalc = fitdist(xc,'Normal','Frequency',fc);
loglogistic = fitdist(xc,'Loglogistic','Frequency',fc);
logisticc = fitdist(xc,'Logistic','Frequency',fc);
poissonc = fitdist(xc,'Poisson','Frequency',fc);
rayleighc = fitdist(xc,'Rayleigh','Frequency',fc);
tlocationc = fitdist(xc,'tLocationScale','Frequency',fc);
exponentialc = fitdist(xc,'Exponential','Frequency',fc);
extremeValuec = fitdist(xc,'ExtremeValue','Frequency',fc);
gammac = fitdist(xc,'Gamma','Frequency',fc);

%% PARAMETRIKES KATANOMES PITHANOTHTAS 8 SE ARITHMO GIA THANATOUS
xd = 1:length(fd);
xd = xd';

normald = fitdist(xd,'Normal','Frequency',fd);
loglogistid = fitdist(xd,'Loglogistic','Frequency',fd);
logisticd = fitdist(xd,'Logistic','Frequency',fd);
poissond = fitdist(xd,'Poisson','Frequency',fd);
rayleighd = fitdist(xd,'Rayleigh','Frequency',fd);
tlocationd = fitdist(xd,'tLocationScale','Frequency',fd);
exponentiald = fitdist(xd,'Exponential','Frequency',fd);
extremeValued = fitdist(xd,'ExtremeValue','Frequency',fd);
gammad = fitdist(xd,'Gamma','Frequency',fd);


fc = fc';
fd = fd';
n1 = length(fc);
n2 = length(fd);

%% GNWSTES KATANOMES POU THA XRHSIMOPOITHOUN
distnamesC = {normalc,loglogistic,logisticc,poissonc,rayleighc,tlocationc,...
    exponentialc,extremeValuec,gammac};

distnamesD = {normald,loglogistid,logisticd,poissond,rayleighd,tlocationd,...
    exponentiald,extremeValued,gammad};


%% BRES TO ELAXISTO MSE GIA KATHE KATANOMH TWN NEWN KROUSMATWN
% [h,pgof,gofstats]=chi2gof(y1);
% pgof

for i=1:1:9
    
    yc(:,i) = pdf(distnamesC{i},xc);
    yc(:,i) = yc(:,i)*sum(fc);
    ec(:,i) = fc - yc(:,i);
    
    yd(:,i) = pdf(distnamesD{i},xd);
    yd(:,i) = yd(:,i)*sum(fd);
    ed(:,i) = fd - yd(:,i);
    
    MSEc(i) = sum((fc-yc(:,i)).^2);
    MSEc(i) =  MSEc(i)/n1;
    
    MSEd(i) = sum((fd-yd(:,i)).^2);
    MSEd(i) =  MSEd(i)/n2;
    
    RMSEc(i) = sqrt(MSEc(i));
    RMSEd(i) = sqrt(MSEd(i));
end

[~,ICmin] = min(MSEc);
[~,IDmin] = min(MSEd);

if ( ICmin == IDmin)
    fprintf('\n');
    fprintf(' LOWEST MSE FOR DEATHS AND CONFIRMED CASES CORESPOND TO THE SAME DISTRIBUTION \n');
end

%% PLOT TIS PDF KAI THN PRAGMATIKH GRAFIKH PARASTASH GIA NEA KROUMSATA
nfig = nfig + 1;
figure(nfig)
plot(1:length(xc),yc(:,1),'k',1:length(xc),yc(:,2),'b',1:length(xc),yc(:,3),...
    'r',1:length(xc),yc(:,4),'-y',1:length(xc),yc(:,5),'g',1:length(xc),yc(:,6),...
    'm',1:length(xc),yc(:,7),'c',1:length(xc),yc(:,8),'ob',1:length(xc),yc(:,9),...
    'ok',1:length(xc),fc,'-.')
title('SPAIN Pdfs of Every Distribution for Confirmed Cases')
ylabel('$Pdfs of Every Distribution $','Interpreter','latex','fontsize',10)
xlabel('$Counts$','Interpreter','latex','fontsize',10)
legend('Normal','LogLogistic','Logistic','Poisson',...
    'Rayleigh','tLocationScale','Exponential','ExtremeValue','Gamma','Confirmed Cases')

%% PLOT TIS PDF KAI THN PRAGMATIKH GRAFIKH PARASTASH GIA THANATOUS
nfig = nfig + 1;
figure(nfig)
plot(1:length(xd),yd(:,1),'k',1:length(xd),yd(:,2),'b',1:length(xd),yd(:,3),...
    'r',1:length(xd),yd(:,4),'-y',1:length(xd),yd(:,5),'g',1:length(xd),yd(:,6),...
    'm',1:length(xd),yd(:,7),'c',1:length(xd),yd(:,8),'ob',1:length(xd),yd(:,9),...
    'ok',1:length(xd),fd,'-.')
title('SPAIN Pdfs of Every Distribution for Deaths')
ylabel('$Pdfs of Every Distribution $','Interpreter','latex','fontsize',10)
xlabel('$Counts$','Interpreter','latex','fontsize',10)
legend('Normal','LogLogistic','Logistic','Poisson',...
    'Rayleigh','tLocationScale','Exponential','ExtremeValue','Gamma','Deaths')

%% RABDOGRAMMA GIA THN KALUTERH PROSARMOGH KATANOMHS GIA NEA KROUSMATA
nfig = nfig + 1;
figure(nfig)
bar(fc)
hold on
plot(yc(:,2),'LineWidth',2)
title('SPAIN LogLogistic Distribution - Covid 19 Confirmed Cases')
ylabel('$Confirmed Cases $','Interpreter','latex','fontsize',10)
xlabel('$Days$','Interpreter','latex','fontsize',10)


%% RABDOGRAMMA GIA THN KALUTERH PROSARMOGH KATANOMHS GIA THANATOUS
nfig = nfig + 1;
figure(nfig)
bar(fd)
hold on
plot(yd(:,2),'LineWidth',2)
title('SPAIN LogLogistic Distribution - Covid 19 Deaths')
ylabel('$Deaths$','Interpreter','latex','fontsize',10)
xlabel('$Days$','Interpreter','latex','fontsize',10)


%% PRINTS APOTELESMATWN
fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' FIND DISTRIBUTION BEST FIT WITH MSE FOR NEW CASES \n');
fprintf(' *******************************************************************\n');
fprintf(' NORMAL DISTRIBUTION           MSE : %.4f || RMSE: %.4f  \n',MSEc(1) ,RMSEc(1));
fprintf(' LOG LOGISTIC DISTRIBUTION     MSE : %.4f || RMSE: %.4f  \n',MSEc(2) ,RMSEc(2));
fprintf(' LOGISTIC DISTRIBUTION         MSE : %.4f || RMSE: %.4f  \n',MSEc(3) ,RMSEc(3));
fprintf(' POISSON DISTRIBUTION          MSE : %.3f || RMSE: %.4f  \n',MSEc(4) ,RMSEc(4));
fprintf(' RAYLEIGH DISTRIBUTION         MSE : %.4f || RMSE: %.4f  \n',MSEc(5) ,RMSEc(5));
fprintf(' t LOCATION SCALE DISTRIBUTION MSE : %.4f || RMSE: %.4f  \n',MSEc(6) ,RMSEc(6));
fprintf(' EXPONENTIAL DISTRIBUTION      MSE : %.3f || RMSE: %.4f  \n',MSEc(7) ,RMSEc(7));
fprintf(' EXTREME VALUE DISTRIBUTION    MSE : %.4f || RMSE: %.4f  \n',MSEc(8) ,RMSEc(8));
fprintf(' GAMMA DISTRIBUTION            MSE : %.4f || RMSE: %.4f  \n',MSEc(9) ,RMSEc(9));
fprintf(' WE CAN SEE THAT THE BEST DISTRIBUTION FOR DATA FIT IS LogLogistic \n');
fprintf(' AND THE SECOND BEST DESTRIBUTION IS Gamma \n');
fprintf(' *******************************************************************\n');

fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' FIND DISTRIBUTION BEST FIT WITH MSE FOR DEATHS \n');
fprintf(' *******************************************************************\n');
fprintf(' NORMAL DISTRIBUTION           MSE : %.4f || RMSE: %.4f  \n',MSEd(1) ,RMSEd(1));
fprintf(' LOG LOGISTIC DISTRIBUTION     MSE : %.4f || RMSE: %.4f  \n',MSEd(2) ,RMSEd(2));
fprintf(' LOGISTIC DISTRIBUTION         MSE : %.4f || RMSE: %.4f  \n',MSEd(3) ,RMSEd(3));
fprintf(' POISSON DISTRIBUTION          MSE : %.3f || RMSE: %.4f  \n',MSEd(4) ,RMSEd(4));
fprintf(' RAYLEIGH DISTRIBUTION         MSE : %.4f || RMSE: %.4f  \n',MSEd(5) ,RMSEd(5));
fprintf(' t LOCATION SCALE DISTRIBUTION MSE : %.4f || RMSE: %.4f  \n',MSEd(6) ,RMSEd(6));
fprintf(' EXPONENTIAL DISTRIBUTION      MSE : %.3f || RMSE: %.4f  \n',MSEd(7) ,RMSEd(7));
fprintf(' EXTREME VALUE DISTRIBUTION    MSE : %.4f || RMSE: %.4f  \n',MSEd(8) ,RMSEd(8));
fprintf(' GAMMA DISTRIBUTION            MSE : %.4f || RMSE: %.4f  \n',MSEd(9) ,RMSEd(9));
fprintf(' WE CAN SEE THAT THE BEST DISTRIBUTION FOR DATA FIT IS LogLogistic \n');
fprintf(' AND THE SECOND BEST DESTRIBUTION IS Gamma \n');
fprintf(' *******************************************************************\n');


%% ********************** SXOLIASMOS APOTELESMATWN ************************
% ARXIKA BLEPODAS TO THIKOGRAMA SUMPERENOUME PWS OI MUSTAKES EXOUN
% DIAFORETIKO MHKOS KAI H DIAMESOS VRHSKETAI KATW APO TO KEDRO THS THHKHS.
% SUNEPWS PERIMENOUME H KATALLHLH KATANOMH NA EINAI MH SUMMETRIKH KAI NA
% EXEI MAKRUES OURES . TETOIES KATANOMES EINAI px LOGISTIC , LOGLOGISTIC..
% SUMFWNA ME TA APOTELESMATA H KALUTERH-KATALLHLOTERH KATANOMH
% GIA TA DEDOMENA MAS EINAI  H LogLogistic KAI EINAI H KATALLHLOTERH KAI
% GIA TA DUO DATASET (NEA KROUSMATA , THANATOI).

% LOG LOGISTIC DISTRIBUTION     MSE : 560695.3369 || RMSE: 748.7959  
% LOG LOGISTIC DISTRIBUTION     MSE : 11469.3224 || RMSE: 107.0949  

% SUMFWNA ME THN ELAXISTOPOIHSH TOU MSE H RMSE 
% EPILEXAME THN LogLogistic H OPOIA OPWS VLEPOUME KAI STO PLOT-RABDOGRAMMA
% KATAXRHSTIKA ISTOGRAMA PROSARMOZETAI POLU KALA STA DEDOMENA. TELOS
% VLEPOUME MERIKES KATANOMES DEN PROSARMOZODE KATHOLOU KALA STA DEDOMENA.
