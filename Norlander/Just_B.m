addpath('../matlab')
addpath('../data')
clear;
close all;
clc;

load('../data/proj18.mat')
load('../data/kalmanRain.mat')
kalmanRain = kalmanRain';
kalmanRain = kalmanRain(793:end);

%% B - Modelling and validation
close all
ndvi = ElGeneina.nvdi;
ndvi_t = ElGeneina.nvdi_t;
rain = kalmanRain;
rain_t = ElGeneina.rain_t;
lambda = bcNormPlot(ndvi);
ndvi_inv_sqr = ndvi.^(lambda);
%% WITH INV_square
% Rescale to [-1,1]
k = 2/(max(ndvi_inv_sqr)-min(ndvi_inv_sqr));
ndvi_transformed = ndvi_inv_sqr*k + (1 - k*max(ndvi_inv_sqr));

% Split data into 70 % modelling, 30 % test.
N = length(ndvi_transformed);
test_size = 0.7;

ndvi_test = ndvi_transformed(floor(N*test_size)+1 : end);
rain_test = rain(floor(N*test_size)+1 : end);

ndvi_t_test = ndvi_t(floor(N*test_size)+1 : end);
rain_t_test = rain_t(floor(N*test_size)+1 : end);

ndvi = ndvi_transformed(1 : floor(N*test_size));
rain = rain(1 : floor(N*test_size));

ndvi_t = ndvi_t(1 : floor(N*test_size));
rain_t = rain_t(1 : floor(N*test_size));


%% See PACF
analyzets(ndvi)
% -> AR(1)

%% Note seasonality in ACF
analyzets(ndvi, 400)
% -> Season 36/3 = 12 months, so a year.


%% Remove season
A36 = [1 zeros(1,35) -1];
ndvi_s = filter(A36,1, ndvi);
%%
analyzets(ndvi_s,400)
%%
data = iddata(ndvi_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );

model_init.Structure.c.Free = [1 zeros(1,35) 1];

model_ar = pem(data, model_init);
present(model_ar)

ehat = resid(ndvi_s, model_ar);
histogram(ehat.y)
title('Histogram over the errors')
xlabel('Deviation')
ylabel('Number of samples')
set(gca,'Fontsize',13)
%%
% If not removing outliers, doesn't pass Ljung Box Pierce!
[~,i_max] = max(ehat.y);
[~,i_min] = min(ehat.y);
ehat.y(i_max) = [];
ehat.y(i_min) = [];
%analyzets(ehat);
%whitenessTest(ehat.y)
%% Pre-whiten rain 
A36 = [1 zeros(1,35) -1];
rain_s = filter(A36, 1, rain);
data = iddata(rain_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );
model_init.Structure.c.Free = [1 1 zeros(1,34) 1];

model_arma = pem(data,model_init);
present(model_arma)

A = model_arma.A;
C = model_arma.C;

rain_pw = filter(A, C, rain_s);
ndvi_pw = filter(A, C, ndvi_s);

%analyzets(rain_pw);
%whitenessTest(rain_pw)
%% 
% Data is now white, estimate impulse response weights
% to find order of B and A2.

M = 40;
stem(-M:M, crosscorr(rain_pw, ndvi_pw, M) );
title('Crosscorr function');
xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(rain_pw))*ones(1,2*M+1), '--') % Why plot for negative?
plot(-M:M, -2/sqrt(length(rain_pw))*ones(1,2*M+1), '--')
hold off
%%
% Yields
d = 1; % 10 day delay after rain
s = 0;
r = 0;

% idpoly(A,B,C,D,F,NoiseVariance,Ts)
% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)
A2 = [1 zeros(1,r)];

B = [zeros(1,d) 1 zeros(1,s)]; 

Mi = idpoly(1, B, [], [], A2);

% Lock delay coefficients
Mi.Structure.b.Free = [ zeros(1,d) 1 ones(1,s)];

zpw = iddata(ndvi_pw, rain_pw);

Mba2 = pem(zpw, Mi);
present(Mba2)

%vhat = resid(Mba2, zpw);
%resid(zpw, Mba2)

x = ndvi_s - filter(Mba2.B, Mba2.F, rain_s);
analyzets(x) % A order 1, C order 0

%% Final Whiteness Testing
close all
A1 = [1 0];
C1 = 1;

Mi = idpoly(1, B, C, A1, A2);
z = iddata(ndvi_s, rain_s);
MboxJ = pem(z,Mi);
present(MboxJ)

resid = resid(z, MboxJ); % Nice!
whitenessTest(resid.y)
%%
plotNTdist(resid.y)
%%
plot(resid.y)
%% Prediction
close all;
clc;

z_test = iddata(ndvi_test, rain_test);
y_test = predict(MboxJ, z_test,1);

z = iddata(ndvi, rain);
y = predict(MboxJ, z, 1);


subplot(211)
plot(ndvi_t_test(1:end), z_test.y(1:end))
hold on
plot(ndvi_t_test, y_test.y)
title('TEST SET')
legend('Real','Estimate')

subplot(212)
plot(ndvi_t, z.y)
hold on
plot(ndvi_t, y.y)
title('TRAINING SET')
legend('Real','Estimate')