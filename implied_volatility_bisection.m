%#ok<*GVMIS> 
format shortG

% Calculating implied volatilities for all strikes using the bisection method 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parsing input
opts = detectImportOptions('options_data.csv');
opts.SelectedVariableNames = [1, 4, 12, 13];  % Strike, quoted value/midpoint/C*, IV, option type

Table = readtable('options_data.csv',opts);

Strikes = Table.Strike;
Quoted_vals = Table.Midpoint; 
IV = Table.IV;
Option_types = Table.Type;
rows = height(Table);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global S, global T, global r, global E, global V_true, global type;
S = 4697.96;
T = 0.0329; % 12 (days) divided by 365 (days)
r = 0.0011;

volatilities = zeros(rows, 1);
computed_prices = zeros(rows, 1);
call_idx = false(rows, 1); put_idx = false(rows, 1); % indexes where option is call or put
eps = 0.00001;

for i = 1:rows
    E = Strikes(i);
    V_true = Quoted_vals(i);
    type = string(Option_types(i));
    
    if type == 'Call', call_idx(i) = 1; else, put_idx(i) = 1; end
    
    % Find a, b over which F changes sign.
    [a,b] = find_a_b();
    
    % Apply bisection:
    while b-a>=eps && a<b
        mid = (a+b)/2;
        
        if F(a)*F(mid)<0
            b = mid;
        else
            a = mid;
        end
    end

    %Solution
    sigma = (a+b)/2;
    volatilities(i) = sigma;
end

% Comparing implied volatilities to column IV:
impliedVolatilities = volatilities*100;
V = table(impliedVolatilities, IV);
disp(V)

% Plotting the results:
plot(Strikes(call_idx), volatilities(call_idx), marker='diamond', Color='blue');
hold on
plot(Strikes(put_idx), volatilities(put_idx), marker='diamond');
xline(S,'--',{'Current asset price'}, 'LabelOrientation', 'horizontal');

xlabel('Exercise price'), ylabel('Implied volatility')
title('Implied volatilities using bisection')
legend('call options','put options')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that finds [sigma_a, sigma_b]
function [a, b] = find_a_b()
    K = 0.05; a = 0;

    while F(a)*F(a+K)>0
        a=a+K;
    end
    b = a+K;
end

% The function F(sigma)
function V = F(sigma)
    global V_true, global type;

    % Compute option value:
    [C, P] = optionValue(sigma);
    
    if type == 'Call'
        V  = C - V_true; % Call
    else 
        V = P - V_true; % Put
    end
end

% Function that computes Put and Call Values:
function [C, P] = optionValue(sigma)
    global T, global S, global r, global E;

    if T > 0
        d1 = (log(S/E) + (r + 0.5*sigma^2)*(T))/(sigma*sqrt(T));
        d2 = d1 - sigma*sqrt(T);
        N1 = 0.5*(1+erf(d1/sqrt(2)));
        N2 = 0.5*(1+erf(d2/sqrt(2)));
        C = S*N1-E*exp(-r*(T))*N2;
        P = C + E*exp(-r*T) - S;
    else 
        C = max(S-E,0);
        P = max(E-S,0);
    end
end

