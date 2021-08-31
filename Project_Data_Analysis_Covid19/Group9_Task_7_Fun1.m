%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZHTHMA 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NORMALIZATION FUNCTION 
function [NormalizedData] = Group9Exe7Fun1(data,choise)

if choise == 1
    NormalizedData = normalize(data,'norm',1);
end
if choise == 2
    NormalizedData = normalize(data,'norm',2);
end
if choise == 3
    NormalizedData = normalize(data,'scale','std');
end
if choise == 4
    NormalizedData = normalize(data);
end

end