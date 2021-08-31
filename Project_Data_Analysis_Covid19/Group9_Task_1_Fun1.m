%% DATA ANALYSIS Project 2020
%% NIKOLAOS ISTATIADIS  AEM:9175
%% KYPARISSIS ODYSSEAS  AEM:8955

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ALGORITHM FOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DATA SMOOTHING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FOR WEEKEND.
function [new_confirmed ,new_dead ] = Group9Exe1Fun1(dc,dd)

NANfinder = isnan(dc);
dc(NANfinder) = 0;
NANfinder = find(isnan(dd));
dd(NANfinder) = 0;

% Smoothing Data for values smaller than zero for daily cases.
[~,Ic] = find(dc < 0);
dc(Ic) = (dc(Ic-1) + dc(Ic+1))/3.333;
dc(Ic+1) = (dc(Ic-1) + dc(Ic+1))/3.333;
dc(Ic-1) = (dc(Ic-1) + dc(Ic+1))/3.333;

% Smoothing Data for values smaller than zero for daily deaths.
[~,Id] = find(dd < 0);
dd(Id) = (dd(Id-1) + dd(Id+1))/3.333;
dd(Id+1) = (dd(Id-1) + dd(Id+1))/3.333;
dd(Id-1) = (dd(Id-1) + dd(Id+1))/3.333;

smothing = 3.3333;
new_confirmed = dc;
new_dead = dd;
[~,I] = find(new_confirmed>0);

for i = I(1):length(new_confirmed)-2
    % Check if Saturdays confirmed cases are zero
    if ( new_confirmed(i) <= 0)
        % Check if Sundays confirmed cases are zero
        if ( new_confirmed(i+1) <= 0)
            % Smoothing Mondays confirmed cases with the weekend days.
            new_confirmed(i) = round(new_confirmed(i+2)/smothing);
            new_confirmed(i+1) = round(new_confirmed(i+2)/smothing);
            new_confirmed(i+2) = round(new_confirmed(i+2)/smothing);
            
        else
            new_confirmed(i) = round(new_confirmed(i+1)/smothing);
            new_confirmed(i+1) = round(new_confirmed(i+1)/smothing);
            new_confirmed(i-1) = round(new_confirmed(i+1)/smothing);
        end
    end
    if ( new_dead(i) <= 0)
        if ( new_dead(i+1) <= 0)
            new_dead(i) = round(new_dead(i+2)/smothing);
            new_dead(i+1) = round(new_dead(i+2)/smothing);
            new_dead(i+2) = round(new_dead(i+2)/smothing);
        else
            new_dead(i) = round(new_dead(i+1)/smothing);
            new_dead(i+1) = round(new_dead(i+1)/smothing);
            new_dead(i-1) = round(new_dead(i+1)/smothing);
        end
    end
end

% Rounding values in order to be integer values.
new_confirmed = round(new_confirmed);
new_dead = round(new_dead);
end