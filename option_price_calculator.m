format long g
S = 52;
E = 50;
r = 0.12;
sigma = 0.3;
t = 0;
T = 0.25;

if T-t == 0
   C = max(S-E,0);
   P = max(E-S,0);
else
    d1_ = d1(S, E, r, sigma, T, t);
    d2_ = d2(S, E, r, sigma, T, t);
    N_d1 = N(d1_);
    N_d2 = N(d2_);
    N_d1_minus = N(-d1_);
    N_d2_minus = N(-d2_);

    C = CallPrice(S, E, r, sigma, T, t);
    P = PutPrice(S, E, r, sigma, T, t);
end

disp("European call price: " + C);
disp("European put price: " + P);

function x = N(x) 
    x = 0.5 * (1 + erf(x/sqrt(2)));
end

function d = d1(S, E, r, sigma, T, t)
    d = (log(S/E) + (r + 0.5* sigma^2) * (T-t))  /  (sigma * sqrt(T-t));
end

function d = d2(S, E, r, sigma, T, t)
    d = (log(S/E) + (r - 0.5* sigma^2) * (T-t))  /  (sigma * sqrt(T-t));
end

function v = CallPrice(S, E, r, sigma, T, t)
    v = S * N(d1(S, E, r, sigma, T, t)) - E * exp(-r*(T-t))*N(d2(S, E, r, sigma, T, t));
end

function v = PutPrice(S, E, r, sigma, T, t)
    v = E * exp(-r*(T-t))*N(-d2(S, E, r, sigma, T, t)) - S * N(-d1(S, E, r, sigma, T, t));
end