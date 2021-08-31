%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND MAX CASES/DEATHS PER DAY FUNCTION

function [maxdate] = Group9Exe3Fun1(data,COUNTRY,datatype,nfig)

f = data;
f = round(f);

%% PARAMETRIKES KATANOMES PITHANOTHTAS 8 SE ARITHMO GIA NEA KROUSMATA
x = 1:length(f);
x = x';

normalc = fitdist(x,'Normal','Frequency',f);
loglogistic = fitdist(x,'Loglogistic','Frequency',f);
logisticc = fitdist(x,'Logistic','Frequency',f);
poissonc = fitdist(x,'Poisson','Frequency',f);
rayleighc = fitdist(x,'Rayleigh','Frequency',f);
tlocationc = fitdist(x,'tLocationScale','Frequency',f);
exponentialc = fitdist(x,'Exponential','Frequency',f);
extremeValuec = fitdist(x,'ExtremeValue','Frequency',f);
gammac = fitdist(x,'Gamma','Frequency',f);

distnames = {normalc,loglogistic,logisticc,poissonc,rayleighc,tlocationc,...
    exponentialc,extremeValuec,gammac};

distnamesS = {'Normal','LogLogistic','Logistic','Poisson',...
    'Rayleigh','tLocationScale','Exponential','ExtremeValue','Gamma'};

for i=1:1:9
    
    Y(:,i) = pdf(distnames{i},x);
    Y(:,i) = Y(:,i)*sum(f);
    e(:,i) = f' - Y(:,i);
    
    n = length(f);
    
    MSE(i) = sum((f'-Y(:,i)).^2);
    MSE(i) =  MSE(i)/n;
    
    RMSE(i) = sqrt(MSE(i));
end

%% BRISKW THN MERA MEGISTOU

[~,Imin] = min(RMSE);
y = Y(:,Imin);

[~,Imax] = max(y);
maxdate = Imax;

s = strcat('DAY WITH THE MOST (',datatype,') ACCORDING TO THE DISTRIBUTION');
fprintf("\n");
fprintf(" *******************************************************************\n");
fprintf(" %s \n",s);
fprintf(" WHICH FITS BEST THE DATASET IS  = %d \n",maxdate);
fprintf(" *******************************************************************\n");


%% BARCHART FOR BEST DISTRIBUTION FIT
figure(nfig)
bar(f)
hold on
plot(y,'LineWidth',2)
title(strcat(COUNTRY,'--',distnamesS{Imin},' Distribution - Covid 19',datatype))
ylabel('$Confirmed Cases $','Interpreter','latex','fontsize',10)
xlabel('$Days$','Interpreter','latex','fontsize',10)
xline(maxdate,'red','LineWidth',3);
hold off;


end