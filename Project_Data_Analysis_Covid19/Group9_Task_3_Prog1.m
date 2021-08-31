%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 3
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



startc= [ 57 ;63 ; 61 ; 76  ;  66 ; 74  ; 71  ; 66  ; 62  ; 62 ; 60];
endc = [ 154 ;146 ; 186 ; 138 ; 186 ; 140 ; 123 ; 156 ; 179 ; 138; 130];

startd= [ 66 ; 81 ; 73  ;  89 ; 75  ; 83  ; 87  ; 83  ; 74  ; 74 ; 71 ];
endd = [ 158 ;157 ; 206 ; 147 ; 194 ; 155 ; 131 ; 156 ; 190 ; 138; 144];

bounds = [60 ,158 ; 80 ,155 ; 65 , 206; 72 ,137; 70 , 194; 78 , 153; 71 , 135; 65 , 156; 68 , 194 ; 68 , 138 ; 66 , 144];

CC = [spac ; ausc ; belc ; croc ; denc ; litc ; estc ; filc ; ollc ; norc ; elvc];
D =  [spad ;ausd ; beld ; crod ; dend ; litd ; estd ; fild ; olld ; nord ; elvd];
COUNTRY = {'SPAIN','AUSTRIA',' BELGIUM', 'CROATIA', 'DENMARK' ,'LITHUANIA'...
           ' ESTONIA', 'FINNLAND', 'NETHERLANDS', 'NORWAY','SWITZERLAND'};
       
N=size(CC,1);
cc = zeros(size(CC));
d = zeros(size(D));
indexc = zeros(N);
indexd = zeros(N);
ed = zeros(N,1);
ec = zeros(N,1);
maxdateC = zeros(N,1);
maxdateD = zeros(N,1);
cTOd = zeros(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:N
    [cc(i,:) , d(i,:)] =  Group9Exe1Fun1(CC(i,:),D(i,:));
    
    %% AN THEORISOUME PWS TO XEKINHMA KAI TELOS TOU PRWTOU KUMATOS EINAI
    %% DIAFORETIKO GIA NEA KROUSMATA KAI THANATOUS
    
    fc = cc(i,startc(i):endc(i));
    fd = d(i,startd(i):endd(i));
    durationc = abs(startc(i) - endc(i));
    durationd = abs(startd(i) - endd(i));
    
    maxdateC(i) = Group9Exe3Fun1(fc,COUNTRY{i},'Confirmed Cases',i);
    maxdateD(i) = Group9Exe3Fun1(fd,COUNTRY{i},'Deaths',i+20);
    
    maxC = startc(i)+ maxdateC(i);
    maxD = startd(i)+ maxdateD(i);
    cTOd(i) =  (maxD - maxC);
    
    fprintf('\n');
    fprintf(' PRESS ANY KEY TO CONTINUE \n');
    pause;
end

%% PARAMETRIKO - BOOTSTRAP 95% DIASTHMA EMPISTOSUNHS
close all;

M=1;
a = 0.05;
B = 1000;
X = cTOd;

%% ARITHMOS DEIGMATOS KAI MESH TIMH
n = length(cTOd);
meanCtoD = mean(cTOd);

%% Plot BOOTSTRAP
bootstrapMean = bootstrp(B,@mean,cTOd);

nfig = nfig +1;
figure(nfig)
clf;
histogram(bootstrapMean)
hold on;
plot([meanCtoD meanCtoD],ylim,'r');
title('Bootstrap-1000 for our sample');


%% Bootstrap CI
CIBoot = bootci(B,{@mean,X'},'type','percentile')';

%% Parametric CI
[h,p,CIParam,~]= ttest(X,0);

%% Histograms
nfig = nfig +1;
figure(nfig);
histogram(CIBoot(1), M);
hold on;
histogram(CIParam(1), M);
legend('Bootstrap','Parametric');
title('Confidence Interval lower bound');
hold off;

nfig = nfig +1;
figure(nfig);
histogram(CIBoot(2), M);
hold on;
histogram(CIParam(2), M);
legend('Bootstrap','Parametric');
title('Confidence Interval upper bound');


%% STATISTIKOS ELEGXOS GIA THN UPOTHESH AN H USTERHSH 14 HMERWN VRISKETAI
%% RANDOM PERMUTATION HYPOTHESYS TEST

lowerLim = (B+1)*a/2;
upperLim = B+1-lowerLim;
limits = [lowerLim upperLim];

TEST = [14,14,14,14,14,14,14,14,14,14];
CIBootstrap = zeros(M,2);
result = zeros(M,1);
for i = 1:M
    % Random permutation test
    samples = [cTOd(i,:) TEST(i,:)];
    difference = zeros(B,1);
    for j = 1:B
        samplesTemp = samples(randperm(length(samples)));
        samplesX = samplesTemp(1:n);
        samplesY = samplesTemp(n+1:end);
        mx = mean(samplesX);
        my = mean(samplesY);
        difference(j) = mx - my;
    end
    stat = mean(cTOd(i,:)) - mean(TEST(i,:));
    difference = [difference; stat];
    difference = sort(difference);
    
    rankStat = find(difference == stat);
    if( length(rankStat) > 1)
        L = length(rankStat);
        sample = randsample(L,1);
        rankStat = rankStat(sample);
    end
    
    if( rankStat < limits(1) || rankStat > limits(2) )
        result(i) = 1;
    end
    
    if rankStat > 0.5*(B+1)
        pd = 2*(1-rankStat/(B+1));
    else
        pd = 2*rankStat/(B+1);
    end
    
end


%% PRINTS APOTELESMATWN
fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' COUNTRIES             : SPAIN, AUSTRIA , BELGIUM , CROATIA , DENMARK , LITHUANIA ,ESTONIA , FINNLAND , NETHERLANDS , NORWAY , SWITZERLAND \n');
fprintf(' (Max-Max) TIME DELAYS : [%d   ,   %d    ,   %d    ,  %d      , %d       , %d         , %d      , %d        , %d           , %d     , %d         ] \n',cTOd(1),cTOd(2),cTOd(3),cTOd(4),cTOd(5),cTOd(6),cTOd(7),cTOd(8),cTOd(9),cTOd(10),cTOd(11));
fprintf(' PARAMETRIC MEAN CONFIDENCE INTERVAL FOR TIME DELAYS : [%.5f , %.5f] \n',CIParam(1),CIParam(2));
fprintf(' BOOTSTRAP  MEAN CONFIDENCE INTERVAL FOR TIME DELAYS : [%.5f , %.5f] \n',CIBoot(1),CIBoot(2));
fprintf(' *******************************************************************\n');

fprintf('\n');
fprintf(' *******************************************************************\n');
fprintf(' RANDOM PERMUTATION HYPOTHESYS TEST FOR TIME DELAY = 14  \n');
if (pd < a)
    fprintf(' p-value FOR OUR HYPOTHESYS TEST IS %.4f \n',pd);
    fprintf(' SO p-value < %.4f AND OUR HYPOTHESYS IS DENIED \n',a);
else
    fprintf(' p-value FOR OUR HYPOTHESYS TEST IS %.5f \n',pd);
    fprintf(' SO p-value > %.4f AND OUR HYPOTHESYS IS ALLOWED \n',a);
end

%% ********************** SXOLIASMOS APOTELESMATWN ************************
% SE AUTO TO ZHTHMA UPOLOGISAME TIS XRONIKES USTERHSEIS GIA TO PLHTHOS TWN
% XWRWN POU ANAFERAME KAI STA PROIGOUMENA ZHTHMATA. VRHKAME THN PIO
% KATALLHLH KATANOMH PITHANOTHTAS XEXWRISTA GIA EPIVEVAIWMENA
% KROUSMATA KAI THANATOUS .TAUTOXRONA UPOLOGISAME MESW AUTHS THN MERA
% KORIFWSEIS EPIVEVAIWMENA KROUSMATWN KAI  THANATWN ANTISOIXA. 
% TELOS UPOLOGISAME TO PARAMETRIKO KAI BOOTSTRAP DIASTHMA EMPISTOSUNHS 
% GIA TO DEIGMA TWN 11 XWRWN (10 XWRES EU + XWRA OMADAS) OPOTE EXOUME:
% PARAMETRIKO  : [3.74441 , 14.25559]
% BOOTSTRAP    : [5.27273 , 14.22727] 
% EINAI FANERO APO TA DIASTHMATA EMPISTOSUNHS PWS H UPOTHESH OTI
% H XRONIKH USTERHSH NA PAREI THN TIMH 14 APORRIPTETAI DIOTI H TIMH 14
% ANHKEI STO 0.05 DIASTHMA. EPISHS KANAME KAI ENA RANDOM PERMUTATION TEST
% GIA NA UPOGRAMISOUME TO APORRIPSH THS UPOTHESHS ELEGXOU MIAS KAI TO
% p-value < 0.05 .
