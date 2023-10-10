format long g
randn('state',100)
clf
%------------------------------------------
% Parameters:
T = 1; S0=3400; sigma = sqrt(0.0003); mu = 0.05;
%------------------------------------------

L = 10;  % no. of intervals
dt = T/L;
tvals = 0:dt:T;  % array of [0, 0.1, 0.2, ..., 0.9, 1]

%-----------------------
% Visualizing one asset path simulation for DISCRETE TIME MODEL

% here i got the s0 in front as a common factor:
S = S0*cumprod(1 + mu*dt + sigma*sqrt(dt)*randn(1,L), 2); % 2 signifies that the cum. product is done across rows
S = [S0, S];  % add S0 in front as the time 0 price

plot(tvals,S)
title('Asset path according to the discrete time model')
xlabel('t'), ylabel('S(t)')
grid on

%-----------------------
% Simulating 500 times, and getting the sample mean & variance

S = S0*cumprod(1 + mu*dt + sigma*sqrt(dt)*randn(500,L), 2); % 500 rows
S = [S0*ones(500,1), S];

S_T = transpose(S(:, end)); % get last column and transpose it to a vector

figure(2)
histogram(S_T, 20)
title('Histogram of asset values at time T=1, according to discrete time model.')

sample_mean = mean(S_T)
sample_variance = var(S_T)

%---------------------------
%  Determining a value of L so that the sample mean and sample variance are 
%  close enough to the theoretical mean and theoretical variance.


for l=9:11
    theoretical_mean = S0*exp(mu*T);
    theoretical_variance = S0^2*exp(2*mu*T)*(exp(sigma^2*T)-1);

    S = S0*cumprod(1 + mu*dt + sigma*sqrt(dt)*randn(500,l), 2); % 500 rows
    S = [S0*ones(500,1), S];
    S_T = transpose(S(:, end)); % get last column and transpose it to a vector
    
    diff_mean = mean(S_T) - theoretical_mean;
    diff_var = var(S_T) - theoretical_variance;

    disp("L=" + l)
    % print difference between sample values and theoretical values
    disp("diff mean: " + diff_mean)
    disp("diff variance: " + diff_var)
end

%Conclusion: L = 10 is already optimal

%---------------------------
% Simulating the SPX index according to CONTINUOUS TIME MODEL

exponent = (mu-0.5*(sigma^2))*dt + sigma*sqrt(dt)*randn(1, L);
S = S0*cumprod(exp(exponent), 2);
S = [S0, S]; %add S0 in the front

figure(3)
plot(tvals,S)
title('Asset path according to the continuous time model')
xlabel('t'), ylabel('S(t)')
grid on
