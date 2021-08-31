%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPOLOGISMOS MONTELOU STEPWISE REGRESSION

function [Y,Yhat,Y2,Y2hat,bestAdjR2,bestRMSE,...
    bestStandarError,bestStandarError2,...
    number_of_params] = Group9Exe7Fun2(x,y,x2,y2)

n1 = length(y);
n2 = length(y2);

%% EPILEGW THN SWSTH METHODO KANONIKOPOIHSHS TRAINING SET KAI TEST SET
c=1;
for u=1:1:4
    for f=1:1:4
        for z=1:1:4
            for a=1:1:4
                
                %% NORMALIZE TA DEDOMENA EKPAIDEUSHS
                xNorm = Group9Exe7Fun1(x ,u);
                yNorm(:,c) = Group9Exe7Fun1(y ,f);
                
                %% NORMALIZE TA DEDOMENA AXIOLOGHSHS
                x2Norm = Group9Exe7Fun1(x2 ,z);
                y2Norm(:,c) = Group9Exe7Fun1(y2 ,a);
                
                
                
                
                %%%%%%%%% STEPWISE REGRESSION STO SET EKPAIDEUSHS %%%%%%%%%
                %% STEPWISE REGRESSION
                [b,~,~,model,stats] = stepwisefit(xNorm,yNorm(:,c));
                b0 = stats.intercept;
                bwise = [b0; b(model)];
                xwise = [ones(length(xNorm),1) xNorm(:,model)];
                yhatwise(:,c) = xwise*bwise;
                ewise = yNorm(:,c) -  yhatwise(:,c);
                
                %% PARAMETERS NUMBER
                k = length(bwise);
                
                %% TUPOPOIHMENO SFALMA
                se2wise = (1/(n1-(k+1)))*(sum(ewise.^2));
                sewise = sqrt(se2wise);
                estarwise = ewise ./ sewise;
                errorStarWise(:,c)  =  estarwise ;
                
                %% MSE - RMSE KAI R^2 - adjR^2
                ssRes = sum((yNorm(:,c)-yhatwise(:,c) ).^2);
                ssTot = sum((yNorm(:,c)-mean(yNorm(:,c))).^2);
                RMSE = sqrt(ssRes/length(yNorm(:,c)));
                adjR2 =1-((n1-1)/(n1-(k+1)))*ssRes/ssTot;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                
                %%%%%%%%% STEPWISE REGRESSION STO SET AXIOLOGISHS %%%%%%%%%
                %% STEPWISE REGRESSION
                x2wise = [ones(length(x2Norm),1) x2Norm(:,model)];
                yhatwise2(:,c) = x2wise*bwise;
                ewise2 = y2Norm(:,c) -  yhatwise2(:,c);

                %% TUPOPOIHMENO SFALMA
                se2wise2 = (1/(n2-(k+1)))*(sum(ewise2.^2));
                sewise2 = sqrt(se2wise2);
                estarwise2 = ewise2 ./ sewise2;
                errorStarWise2(:,c)  =  estarwise2 ;
                
                %% MSE - RMSE KAI R^2 - adjR^2
                ssRes2 = sum((y2Norm(:,c)- yhatwise2(:,c)).^2);
                ssTot2 = sum((y2Norm(:,c)-mean(y2Norm(:,c))).^2);
                RMSE2 = sqrt(ssRes2/length(y2Norm(:,c)));
                adjR22 =1-((n2-1)/(n2-(k+1)))*ssRes2/ssTot2;
                
                
                all1adjR2(c) =  adjR2;
                all2adjR2(c) =  adjR22;
                all1RMSE(c) =  RMSE;
                all2RMSE(c) =  RMSE2;
                c=c+1;
            end
        end
    end
end
[max1AdjR2,I1max]=max(all1adjR2);
[max2AdjR2,I2max]=max(all2adjR2);
Y = yNorm(:,I1max);
Yhat = yhatwise(:,I1max);
Y2 = y2Norm(:,I2max);
Y2hat = yhatwise2(:,I2max);
bestAdjR2 = [max1AdjR2 ,  max2AdjR2];
bestRMSE = [ all1RMSE(:,I1max) ,  all2RMSE(:,I2max)];
bestStandarError = errorStarWise(:,I1max);
bestStandarError2= errorStarWise2(:,I2max);


%% UPOLOGISMOS PARAMETRWN POU XRHSIMOPOIHTHIKAN GIA REGRESSION
number_of_params = length(bwise);

end