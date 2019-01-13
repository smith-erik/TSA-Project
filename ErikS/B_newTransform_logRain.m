



%% 0. Data Preprocessing
addpath('../matlab')
addpath('../data')
clear;
close all;
clc;

load('../data/proj18.mat')
load('../data/kalmanRain.mat')
kalmanRain = kalmanRain';
kalmanRain = kalmanRain(793:end);

ndvi = ElGeneina.nvdi;
ndvi_t = ElGeneina.nvdi_t;
rain = kalmanRain;
rain_t = ElGeneina.rain_t;

%% Transforms
%lambda = bcNormPlot(ndvi);
%ndvi_inv_sqr = ndvi.^(lambda);
%lambda = bcNormPlot(ndvi);
%lambda2 = bcNormPlot(rain);
%ndvi = (ndvi).^(lambda);
%rain = (rain+2).^(lambda2);

% Rescale 0 - 255 values.
% 0 -> -1 and 255 -> 1.
% Zero will be at 127.5.
ndvi_transformed = (2/255) * ndvi - 1;
fprintf('-- Transformed --\nMax: %f.\nMin: %f.\n', ...
        max(ndvi_transformed), min(ndvi_transformed))

% Log-transform the rain
rain = log(rain + 1); % Add 1 else 0 -> -Inf. Now 0 -> 0.

%% Split into model set (70%) and test set (30%)

N = length(ndvi_transformed);
test_size = 0.7;

ndvi_test = ndvi_transformed(floor(N*test_size)+1 : end);
rain_test = rain(floor(N*test_size)+1 : end);

ndvi_t_test = ndvi_t(floor(N*test_size)+1 : end);
rain_t_test = rain_t(floor(N*test_size)+1 : end);

ndvi = ndvi_transformed(1 : floor(N*test_size));
rain = rain(1 : floor(N*test_size));

ndvi_t = ndvi_t(1 : floor(N*test_size)); % OBS: We re-assign NDVI here.
rain_t = rain_t(1 : floor(N*test_size));

%% Note seasonality in ACF
% AR(1) seems reasonable.
% Also, we have a season of 36/3 = 12 months.
analyzets(ndvi, 80)

%% Remove the season

A36 = [1 zeros(1,35) -1];
ndvi_s = filter(A36, 1, ndvi);

figure(2)
plotNTdist(ndvi_s)

figure(1)
analyzets(ndvi_s, 80)




%% 1. Simple (almost) AR model, no rain as input.

data = iddata(ndvi_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );

model_init.Structure.c.Free = [1 zeros(1,35) 1];

model_ar = pem(data, model_init);
present(model_ar)

ehat = resid(ndvi_s, model_ar);

analyzets(ehat.y)

figure(2)
whitenessTest(ehat.y)




%% 2. More advanced Box-Jenkins model with rain as input.

% Remove season and select order for pre-whitening.
A36 = [1 zeros(1,35) -1];
rain_s = filter(A36, 1, rain);

% ARMA(1, 2) seems reasonable.
% Adding an extra order 36 coeficcient in C (and/or A) makes it worse.
analyzets(rain_s)

%% Pre-whiten rain with chosen polynomials

data = iddata(rain_s);
model_init = idpoly([1 0], [] , [1 zeros(1,36)] );
model_init.Structure.c.Free = [1 0 0 0 1 zeros(1,31) 1];
model_arma = pem(data,model_init);
present(model_arma)

A3 = model_arma.A;
C3 = model_arma.C;

rain_pw = filter(A3, C3, rain_s);
ndvi_pw = filter(A3, C3, ndvi_s);

analyzets(rain_pw);

figure(2)
whitenessTest(rain_pw)

%% Estimate impulse response weights

% Analyze plot and select order of B and A2.
M = 40;
stem(-M:M, crosscorr(rain_pw, ndvi_pw, M) );
title('Crosscorr function');
xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(rain_pw))*ones(1,2*M+1), '--') % Why plot for negative?
plot(-M:M, -2/sqrt(length(rain_pw))*ones(1,2*M+1), '--')
hold off

%% Estimate with selected orders and calculate x(t)

d = 2;
s = 4;
r = 1;
% idpoly(A,B,C,D,F,NoiseVariance,Ts)
% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)
A2 = [1 zeros(1,r)];
B = [zeros(1,d) 1 zeros(1,s)]; 

% Initialize and lock delay coefficients to zero.
model_init = idpoly(1, B, [], [], A2);
model_init.Structure.b.Free = [zeros(1,d) 1 ones(1,s)];

data = iddata(ndvi_pw, rain_pw);
Model_B_A2 = pem(data, model_init);
present(Model_B_A2)

%figure(2)
%resid(Model_B_A2, data)

% Analyze x to select orders for A1 and C1.
% A1 order 1, C1 order 0 seems reasonable.
x = ndvi - filter(Model_B_A2.B, Model_B_A2.F, rain);
analyzets(x)

%% Also check that CCF from u to x is white
M = 100;
stem(-M:M, crosscorr(x, rain, M) );
title('Crosscorr function');
xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(rain_pw))*ones(1,2*M+1), '--') % Why plot for negative?
plot(-M:M, -2/sqrt(length(rain_pw))*ones(1,2*M+1), '--')
hold off

%% Complete model
close all;

% Estimate the complete model with orders as selected above.
A1 = [1 0];
C1 = [1 zeros(1,36)];

% idpoly(A,B,C,D,F,NoiseVariance,Ts)
% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)

data = iddata(ndvi_s, rain_s);
model_init = idpoly(1, B, C1, A1, A2);
model_init.Structure.c.Free = [1 1 1 zeros(1,33) 1];

MboxJ = pem(data, model_init);
present(MboxJ)

% Analyze resudials
finalRes = resid(data, MboxJ);

figure(2)
whitenessTest(finalRes.y)

figure(1)
analyzets(finalRes.y)

figure(3)
plotNTdist(finalRes.y)



%% k-step prediction
close all;
clc;

k = 1;

% Predict on data with season removed.
ndvi_test_s = filter(A36, 1, ndvi_test);
rain_test_s = filter(A36, 1, rain_test);

data_test = iddata(ndvi_test_s, rain_test_s);
pred_test = predict(MboxJ, data_test, k);

data_train = iddata(ndvi_s, rain_s);
pred_train = predict(MboxJ, data, k);

% Add back season when plotting.
pred_test.y = filter(1, A36, pred_test.y);
data_test.y = filter(1, A36, data_test.y);

pred_train.y = filter(1, A36, pred_train.y);
data_train.y = filter(1, A36, data_train.y);


subplot(2,1,1)
plot(ndvi_t_test(k+1:end), data_test.y(1:end-k))
hold on
plot(ndvi_t_test(k+1:end), pred_test.y(1:end-k))
title('Test Set')
legend('Real','Estimate', 'location', 'southwest')

subplot(2,1,2)
plot(ndvi_t(k+1:end), data_train.y(1:end-k))
hold on
plot(ndvi_t(k+1:end), pred_train.y(1:end-k))
title('Training Set')
legend('Real','Estimate', 'location', 'southwest')







%% Predict and save plots for k as 1 through 10

for k = 1:10
    close all;
    clc;

    ndvi_test_s = filter(A36, 1, ndvi_test);
    rain_test_s = filter(A36, 1, rain_test);

    data_test = iddata(ndvi_test_s, rain_test_s);
    pred_test = predict(MboxJ, data_test, k);

    data_train = iddata(ndvi_s, rain_s);
    pred_train = predict(MboxJ, data, k);

    % Add back season when plotting.
    pred_test.y = filter(1, A36, pred_test.y);
    data_test.y = filter(1, A36, data_test.y);

    pred_train.y = filter(1, A36, pred_train.y);
    data_train.y = filter(1, A36, data_train.y);

    % Plot
    subplot(2,1,1)
    plot(ndvi_t_test(k+1:end), data_test.y(1:end-k))
    hold on
    plot(ndvi_t_test(k+1:end), pred_test.y(1:end-k))
    title(sprintf('Test Set with k = %d', k))
    legend('Real','Estimate', 'location', 'southwest')

    subplot(2,1,2)
    plot(ndvi_t(k+1:end), data_train.y(1:end-k))
    hold on
    plot(ndvi_t(k+1:end), pred_train.y(1:end-k))
    title(sprintf('Training Set with k = %d', k))
    legend('Real','Estimate', 'location', 'southwest')
    
    saveas(gcf,sprintf('Plots/PredPlots/k%d_log.png', k))

end





