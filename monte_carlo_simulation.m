randn('state',100)

% Estimating the time-zero option price using a Monte Carlo simulation.

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 52; E = 50; sigma = 0.3; r = 0.12; T = 0.25; 
Dt = 1e-3; N = T/Dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_values = 2.^(1:17); %2^1, 2^2, ... 2^17
disp("Monte Carlo errors:")

for k = 1:length(M_values)
    M = M_values(k);
    
    V = zeros(M,1);
    for i = 1:M
        Sfinal = S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*randn);
        V(i) = exp(-r*T)*max(Sfinal-E,0);
    end

    aM = mean(V);
    bM = std(V);
    conf = [aM - 1.96*bM/sqrt(M), aM + 1.96*bM/sqrt(M)];
    err = bM/sqrt(M); % Monte Carlo error
    
    % output the error
    disp("M = 2^" + string(k) + ":  " + string(err))
    
    % plot mean and confidence intervals:
    plot([k k], [conf(1) conf(2)], 'k')
    hold on
    plot(k, aM, 'rx')
    hold on
end

% Visualizing the Monte Carlo convergence:
BS_solution = 5.057;
yline(BS_solution ,'--','5.057');

ylim([0 10])
xticks(1:length(M_values))
xticklabels(M_values)
title('Monte Carlo convergence')
xlabel('Num samples'), ylabel('Approximation')

disp("Approximate option price: " + aM)