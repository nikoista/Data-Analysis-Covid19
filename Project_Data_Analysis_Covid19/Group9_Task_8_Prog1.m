%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 8
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


%% 1 KUMA BOUNDS
bounds = [60, 158; 65 , 206; 70 , 194; 68 , 194 ; 68 , 138 ; 66 , 144];

%% 2 KUMA BOUNDS
bounds2 = [185, 345; 200, 348; 216 , 348; 202 , 348 ; 213 , 346 ; 265 , 346];

CC = [spac ;  belc ;  denc ;  ollc ; norc ; elvc ];
D =  [spad ;  beld ;  dend ;  olld ; nord ; elvd ];

COUNTRY = {'SPAIN',' BELGIUM', 'DENMARK' ,'NETHERLANDS', 'NORWAY','SWITZERLAND'};
REGRESSION_METHODS = {' PLS',' LASSO'};

%% CONSTANTS
alpha=0.05;
B1 = -20;
B2 = 20;

%% MATRICES INITIALIZATION
N=size(CC,1);
cc = zeros(size(CC));
d = zeros(size(D));
duration = zeros(N,1);
duration2 = zeros(N,1);
R2 = zeros(N,1);
R22 =  zeros(N,B2+1);
adjR2 = zeros(N,1);
adjR22 = zeros(N,1);
RMSE =  zeros(N,1);
NMSE =  zeros(N,1);
NDEI =  zeros(N,1);
RMSE2 =  zeros(N,1);
NMSE2 =  zeros(N,1);
NDEI2 =  zeros(N,1);
paramStepwise =  zeros(N,1);
bestAdjR2 =  zeros(N,2);
bestRMSE =  zeros(N,2);
all1adjR2 =  zeros(256,1);
all2adjR2 =  zeros(256,1);
all1RMSE =  zeros(256,1);
all2RMSE =  zeros(256,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    
     
    %% DATA SMOOTHING
    [cc(i,:) , d(i,:)] =  Group9Exe1Fun1(CC(i,:),D(i,:));
    
    %% 1 KUMA BOUNDS GIA THANATOUS KAI NEA KROUSMATA
    y = d(i,bounds(i,1):(bounds(i,2)));
    y2 = d(i,bounds2(i,1):(bounds2(i,2)));
    
    x = zeros(B2+1,length(cc(i,(bounds(i,1)):bounds(i,2))));
    x2 = zeros(B2+1,length(cc(i,(bounds2(i,1)):bounds2(i,2))));
    
    %% EPANALHPTIKH DIADIKASIA EURESHS XRONIKWN USTERISEWN
    c=1;
    for T=0:1:B2
        x(c,:) = cc(i,(bounds (i,1)-T):bounds (i,2)-T);
        x2(c,:) = cc(i,(bounds2(i,1)-T):bounds2(i,2)-T);
        c=c+1;
    end
    
    x=x';
    y=y';
    x2=x2';
    y2=y2';
    
    
    
    %% ******************* EPILOGH METHODOU REGRESSION ***************** %%
    % 1 GIA PLS  
      reg_num = 1;
    
    % 2 GIA LASSO
      %reg_num = 2;
    

    %% SUNARTHSH UPOLOGISMOU REGRESSION MONTELOU 
    [Y,Yhat,Y2,Y2hat,bestAdjR2,...
    bestRMSE,bestStandarError,...
    bestStandarError2] = Group9Exe8Fun1(x,y,x2,y2,reg_num);
    
    ALL_AdjR2(i,:) = bestAdjR2;
    ALL_RMSE(i,:) = bestRMSE;
    ALL_standard_Errors{i} = bestStandarError;
    ALL_standard_Errors2{i} = bestStandarError2;
    
    
    %% DIAGRAMMATA THANATWN GIA PRAGMATIKES TIMES KAI EKTIMHSEIS GIA TO
    %% 1 KUMA KAI TO 2 KUMA
    nfig = nfig+1;
    figure(nfig)
    subplot(2,1,1);
    
    plot(Y);
    hold on
    plot(Yhat);
    title(  sprintf(strcat(COUNTRY{i},REGRESSION_METHODS{reg_num},' Regression for 1 Wave | adjR^2 = %2.4f'),bestAdjR2(1)))
    xlabel(  'Daily Deaths')
    ylabel(  'Days')
    legend('Y-real','Y-predict')
    hold off;
    
    subplot(2,1,2);
    plot(Y2);
    hold on
    plot(Y2hat);
    title(   sprintf(strcat(COUNTRY{i},REGRESSION_METHODS{reg_num},' Regression for 1 Wave | adjR^2 = %2.4f'),bestAdjR2(2)))    
    xlabel(  'Daily Deaths')
    ylabel(  'Days')
    legend('Y2-real','Y2-predict')
    hold off;
    
    %% DIAGNWSTIKO DIAGRAMA TUPOPOIHMENWN SFALMATWN  STEPWISE MODEL
    nfig = nfig+1;
    figure(nfig)
    subplot(2,1,1);
    scatter(Y,ALL_standard_Errors{i});
    hold on;
    plot(Y,repmat(1.96,1,length(Y)));
    hold on;
    plot(Y,repmat(0,1,length(Y)));
    hold on;
    plot(Y,repmat(-1.96,1,length(Y)));
    title( sprintf(strcat(COUNTRY{i},REGRESSION_METHODS{reg_num},' Regression standardised errors ')))
    xlabel( 'Daily Deaths')
    ylabel( 'Residuals of Regresion')
    
    subplot(2,1,2);
    scatter(Y2,ALL_standard_Errors2{i});
    hold on;
    plot(Y2,repmat(1.96,1,length(Y2)));
    hold on;
    plot(Y2,repmat(0,1,length(Y2)));
    hold on;
    plot(Y2,repmat(-1.96,1,length(Y2)));
    title(  sprintf(strcat(COUNTRY{i},REGRESSION_METHODS{reg_num},' Regression standardised errors ')))
    xlabel( 'Daily Deaths')
    ylabel( 'Residuals of Regresion')
    
end

%% PRINTS APOTELESMATWN
fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' COUNTRIES                 :    SPAIN , BELGIUM , DENMARK  , NETHERLANDS , NORWAY , SWITZERLAND \n');
fprintf(strcat(REGRESSION_METHODS{reg_num},' || adjR^2 || 1 WAVE :  [%2.4f , %2.4f  , %2.4f   , %2.4f      , %2.4f , %2.4f] \n'),ALL_AdjR2(1,1),ALL_AdjR2(2,1),ALL_AdjR2(3,1),ALL_AdjR2(4,1),ALL_AdjR2(5,1),ALL_AdjR2(6,1));
fprintf(strcat(REGRESSION_METHODS{reg_num},' || adjR^2 || 2 WAVE :  [%2.4f , %2.4f  , %2.4f   , %2.4f      , %2.4f , %2.4f] \n'),ALL_AdjR2(1,2),ALL_AdjR2(2,2),ALL_AdjR2(3,2),ALL_AdjR2(4,2),ALL_AdjR2(5,2),ALL_AdjR2(6,2));
fprintf(strcat(REGRESSION_METHODS{reg_num},' ||  RMSE  || 1 WAVE :  [%2.4f , %2.4f  , %2.4f   , %2.4f      , %2.4f , %2.4f] \n'),ALL_RMSE(1,1),ALL_RMSE(2,1),ALL_RMSE(3,1),ALL_RMSE(4,1),ALL_RMSE(5,1),ALL_RMSE(6,1));
fprintf(strcat(REGRESSION_METHODS{reg_num},' ||  RMSE  || 2 WAVE :  [%2.4f , %2.4f  , %2.4f   , %2.4f      , %2.4f , %2.4f] \n'),ALL_RMSE(1,2),ALL_RMSE(2,2),ALL_RMSE(3,2),ALL_RMSE(4,2),ALL_RMSE(5,2),ALL_RMSE(6,2));
fprintf(' *******************************************************************\n');
fprintf('\n'); 

% 1 KUMA!

%  COUNTRIES                 :    SPAIN , BELGIUM , DENMARK  , NETHERLANDS , NORWAY , SWITZERLAND 
%  PLS   || adjR^2 || 1 WAVE :  [0.7908 , 0.8759  , 0.6170   , 0.7083      , 0.1802 , 0.3608] 
%  LASSO || adjR^2 || 1 WAVE :  [0.9123 , 0.9741  , 0.8028   , 0.9110      , 0.3662 , 0.6818] 
%  STEP  || adjR^2 || 1 WAVE :  [0.9048 , 0.9736  , 0.8084   , 0.9015      , 0.4547 , 0.6713] 

% 2 KUMA!

%  COUNTRIES                 :    SPAIN , BELGIUM , DENMARK  , NETHERLANDS , NORWAY , SWITZERLAND 
%  PLS   || adjR^2 || 2 WAVE :  [0.2590 , 0.6546  , 0.5826   , 0.6776      , 0.1733 , 0.5100] 
%  LASSO || adjR^2 || 2 WAVE :  [0.2668 , 0.2316  , 0.6070   , 0.5882      , 0.1516 , 0.3260] 
%  STEP  || adjR^2 || 2 WAVE :  [0.3633 , 0.3072  , 0.6716   , 0.6284      , 0.2771 , 0.5117] 

%% ********************** SXOLIASMOS APOTELESMATWN ************************
% FTANONTAS STO TELOS THS ERGASIA DHLADH STO ZHTHMA 8 THA EPANALABOUME THN
% IDIA DIADIKASIA OPWS STO ZHTHMA 7 ALLA THA XRHSIMOPOIHSOUME ENA ALLO
% MONTELO MEIWSEIS DIASTASHS ( PLS, LASSO, RIDGE). KLEINONTAS THA
% SUGKRINOUME TA DUO MONTELA ME VASH TA SFALMATA PROSARMOGHS (RMSE). 
% SUNEPWS VLEPOUME TA APOTELESMATA TWN PARAPANW PINAKWN KAI SUMPERENOUME OTI
% STO SUNOLO EKPAIDEUSHS H LASSO EINAI H KALHTERH KAI AKOLOUTHEI H STEPWISE
% KAI TELOS H PLS.
% STO SUNOLO AXIOLOGHSHS H STEPWISE EINAI H KALUTERH KAI AKOLOUTHEI H PLS KAI
% TELOS H LASSO.



%%  PROVLIMATISMOI
% UPHRXAN MERIKH PROVLIMATISMOI PANW STON TROPO ME TON OPOIO THA GINEI H
% KANONIKOPOIHSH KAI THN METHODO ME THN OPOIA THA GINEI. tO SWSTO THA HTAN
% EPREPE NA KANONIKOPOIHSOUME TO TRAINING SET MAS KAI USTERA NA
% KANONIKOPOIHSOUME TO TESTING SET SUMFWNA ME THN MESH TIMH KAI TUPIKH
% APOKLHSH TOU TRAINING SET OMWS AUTO THA MAS ODHGOUSE SE PROVLHMATA WS
% PROS THN SUGKRISH TWN adjR^2 TOU TRAINING KAI TESTING SET. OSON AFORA THN
% METHODO ME THN OPOIA EGINE KANONIKOPOIHSH PROSPATHISAME NA VROUME ENAN
% TROPO ME TON OPOIO H KANONIKOPOIHSH MAS THA ODHGOUSE TON adjR^2  STHN
% MEGISTH TIMH TOU OPOU KAI EGINE.