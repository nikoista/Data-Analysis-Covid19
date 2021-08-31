%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS
%%%%%%%%%% UPOLOGISMOS MONTELOU REGRESSION SUMFWNA METHN EPILOGH TOU XRHSTH

function [Y,Yhat,Y2,Y2hat,bestAdjR2,bestRMSE,...
    bestStandarError,bestStandarError2,...
    number_of_params] = Group9Exe8Fun1(x,y,x2,y2 ,REG)

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
                
                
              
                %%
                %%%%%%%%% STEPWISE REGRESSION STO SET EKPAIDEUSHS %%%%%%%%%
                
                %% EDW EPILEGETAI H METHODOS GIA REGRESSION
                if REG == 1
                    [~,~,~,~,b] = plsregress(xNorm,yNorm(:,c),1);
                    
                     %% PARAMETERS NUMBER
                     k = length(b) - 1;
                     
                elseif REG ==2
                    lambda = 0.0001; 
                    b = lasso(xNorm,yNorm(:,c),'Lambda',lambda);
                    
                    %% PARAMETERS NUMBER
                     k = length(b) ;
                     
                end
                if REG == 2
                    yhat(:,c) = xNorm*b;
                else
                    yhat(:,c) = [ones(length(xNorm),1) xNorm]*b;
                end
                
                e = yNorm(:,c) -  yhat(:,c) ;
                
               
                
                %% TUPOPOIHMENO SFALMA
                se2= (1/(n1-(k+1)))*(sum(e.^2));
                se = sqrt(se2);
                estar = e ./ se;
                errorStar(:,c)  =  estar ;
                
                %%  RMSE , adjR^2
%                 ssRes = sum((yNorm(:,c)-yhat(:,c) ).^2);
%                 ssTot = sum((yNorm(:,c)-mean(yNorm(:,c))).^2);
                RMSE =  sqrt(sum((e.^2)/length(yNorm(:,c))));
                adjR2 =1-((n1-1)/(n1-(k+1)))*(sum(e.^2))/(sum((yNorm(:,c)-mean(yNorm(:,c))).^2));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 %% SUNTELESTES R^2 , adjR^2
%                 R2wise = 1-(sum(ewise.^2))/(sum((y-mean(y)).^2));
%                 adjR2wise =1-((n-1)/(n-(k1+1)))*(sum(ewise.^2))/(sum((y-mean(y)).^2));
%                 
%                 errord = ((fd - yd).^2)/length(fd);
%                 MSEd(i)= sum(errord);
%                 ed = fd - yd;
%                 
%                 
%                 RMSEc(i) = sqrt(MSEc(i));
%                 
                
                
                %%
                %%%%%%%% STEPWISE REGRESSION STO SET AXIOLOGHSHS %%%%%%%%%%
                if REG == 2
                    y2hat(:,c) = x2Norm*b;
                else
                    y2hat(:,c) = [ones(length(x2Norm),1) x2Norm]*b;
                end
                e2 = y2Norm(:,c) -  y2hat(:,c) ;
                
                %% PARAMETERS NUMBER
                k = length(b);
                
                %% TUPOPOIHMENO SFALMA
                s2e2= (1/(n2-(k+1)))*(sum(e2.^2));
                s2e = sqrt(s2e2);
                estar2 = e2 ./ s2e;
                errorStar2(:,c)  =  estar2 ;
                
                %%  RMSE , adjR^2
                RMSE2 =  sqrt(sum((e2.^2)/length(y2Norm(:,c))));
                adjR22 =1-((n2-1)/(n2-(k+1)))*(sum(e2.^2))/(sum((y2Norm(:,c)-mean(y2Norm(:,c))).^2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                all1adjR2(c) =  adjR2;
                all2adjR2(c) =  adjR22;
                all1RMSE(c) =  RMSE;
                all2RMSE(c) =  RMSE2;
                c=c+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end
[max1AdjR2,I1max]=max(all1adjR2);
[max2AdjR2,I2max]=max(all2adjR2);
Y = yNorm(:,I1max);
Yhat = yhat(:,I1max);
Y2 = y2Norm(:,I2max);
Y2hat = y2hat(:,I2max);
bestAdjR2 = [max1AdjR2 ,  max2AdjR2];
bestRMSE = [ all1RMSE(I1max) , all2RMSE(I2max)];
bestStandarError = errorStar(:,I1max);
bestStandarError2 = errorStar2(:,I2max);

%% UPOLOGISMOS PARAMETRWN POU XRHSIMOPOIHTHIKAN GIA REGRESSION
number_of_params = length(b);

end