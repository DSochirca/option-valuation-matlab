format shortG

% Calculating implied volatilities for all strikes using the Newton's method

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

S = 4697.96;
T = 0.0329; % 12 (days) divided by 365 (days)
r = 0.0011;

volatilities = zeros(rows, 1);
call_idx = false(rows, 1); put_idx = false(rows, 1); % indexes where option is call or put

for i = 1:rows
    E = Strikes(i);
    V_true = Quoted_vals(i);
    type = string(Option_types(i));
    
    if type == 'Call', call_idx(i) = 1; else, put_idx(i) = 1; end
    
    % starting value
    sigmahat = sqrt(2*abs( (log(S/E) + r*T)/T ) );

    tol = 1e-8;
    sigma = sigmahat;
    sigmadiff = 1;
    k = 1;
    kmax = 100;
    while (sigmadiff >= tol && k < kmax)
        [C, Cvega, P, Pvega] = val(S,E,r,sigma,T);
        
        if type == 'Call' 
            increment = (C-V_true)/Cvega;
        else
            increment = (P-V_true)/Pvega;
        end
        
        sigma = sigma - increment;
        k = k+1;
        sigmadiff = abs(increment);
    end

    % The solution
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
title("Implied volatilities using Newton's method")
legend('call options','put options')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that computes Put, Call values and vega
function [C, Cvega, P, Pvega] = val(S,E,r,sigma,tau)
    if tau > 0
        d1 = (log(S/E) + (r + 0.5*sigma^2)*(tau))/(sigma*sqrt(tau));
        d2 = d1 - sigma*sqrt(tau);
        N1 = 0.5*(1+erf(d1/sqrt(2)));
        N2 = 0.5*(1+erf(d2/sqrt(2)));
        C = S*N1-E*exp(-r*(tau))*N2;
        Cvega = S*sqrt(tau)*exp(-0.5*d1^2)/sqrt(2*pi);
        P = C + E*exp(-r*tau) - S;
        Pvega = Cvega;
    else 
        C = max(S-E,0);
        Cvega = 0;
        P = max(E-S,0);
        Pvega = 0;
    end
end