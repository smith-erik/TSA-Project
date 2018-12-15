%% Once
close all;
clear;
clc;

% Släng 100 för att svänga in?


% Simulate a process
n = 500;
A = [1 -1.35 0.43];
sigma2 = 4;
noise = sqrt(sigma2) * randn(n + 100, 1);
y = filter(1,A,noise);
y = y(101:end);
subplot (2,1,1)
plot(y)
subplot (2,1,2)
plot( noise )

% Split and create data sets
n_est = floor(2/3 * n);
y_est = iddata( y(1: n_est) );
y_val = iddata( y(n_est + 1:end) );


NN = (1:10)';


V = arxstruc(y_est, y_val, NN);
n_order = selstruc(V, 0); % 0 for LS criterion
n_aic = selstruc(V, 'aic');


%% With a loop and histogram plot

close all;
clear;
clc;

n = 10000;
n_est = floor(2/3 * n);

A = [1 -1.35 0.43];
sigma2 = 4;

NN = (1:10)';
n_times = 100;
n_order = zeros(n_times, 1) - 1;
n_aic = zeros(n_times, 1) - 1;

for i = 1:n_times
    
    noise = sqrt(sigma2) * randn(n + 100, 1);
    y = filter(1,A,noise);
    y = y(101:end);
    
    y_est = iddata( y(1:n_est) );
    y_val = iddata( y(n_est+1:end) );
    
    V = arxstruc(y_est, y_val, NN);
    n_order(i) = selstruc(V, 0);
    n_aic(i) = selstruc(V, 'aic');
    
end

subplot(1,2,1)
histogram(n_order)
title('LS')

subplot(1,2,2)
histogram(n_aic)
title('AIC')


% Q4.1: AIC är ett  bättre / mer consistent kriterium.
% Q4.2: Order 2 sticker ivär även för LS med n=10000.


ar_model = arx(y, n_aic(end));
ar_model.NoiseVariance
ar_model.CovarianceMatrix
present(ar_model)



















