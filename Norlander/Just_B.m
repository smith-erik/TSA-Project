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
%lambda = bcNormPlot(ndvi);
%ndvi_inv_sqr = ndvi.^(lambda);
%% WITH INV_square
% Rescale to [-1,1]
lambda = bcNormPlot(ndvi);
lambda2 = bcNormPlot(rain);
ndvi = (ndvi).^(lambda);
rain = (rain+2).^(lambda2);

%% Split into model set and test set
k = 2/(max(ndvi)-min(ndvi));
ndvi_transformed = ndvi*k + (1 - k*max(ndvi));

%ndvi_transformed = ndvi./255;

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
%ndvi_s = ndvi; % Don't do separate, add a36 in A instead.
%%
analyzets(ndvi_s, 50)

%% 1. Simple ARMA model, no rain as input.
data = iddata(ndvi_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );

model_init.Structure.c.Free = [1 zeros(1,35) 1];
%model_init.Structure.a.Free = [1 zeros(1,35) 1];

model_arma = pem(data, model_init);
present(model_ar)
%%
ehat = resid(ndvi_s, model_arma);

%histogram(ehat.y)
%title('Histogram over the errors')
%xlabel('Deviation')
%ylabel('Number of samples')
%set(gca,'Fontsize',13)

analyzets(ehat.y)
figure(2)
whitenessTest(ehat.y)
%% Prediction with AR1-model
pred = iddata(ndvi_test);
k = 1;
z = predict(model_ar,pred,k);
figure(2); clf;
plot(ndvi_t_test, ndvi_test)
hold on
plot(ndvi_t_test(1:end-k), z.y(k+1:end))
title('Prediction')
legend('Real','Estimate')
ehat1 = resid(pred,model_ar);
sum(ehat1.y.^2)
%% 2. More advanced Box-Jenkins model with rain as input.
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
whitenessTest(rain_pw)
%% 
% Data is now white, estimate impulse response weights
% to find order of B and A2. 
figure(1)
M = 40;
stem(-M:M, crosscorr(rain_pw, ndvi_pw, M) );
title('Crosscorr function of upw to ypw');
xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(rain_pw))*ones(1,2*M+1), '--') % Why plot for negative?
plot(-M:M, -2/sqrt(length(rain_pw))*ones(1,2*M+1), '--')
set(gca,'Fontsize',11)
hold off
%%
% Yields
d = 2; % 10 day delay after rain
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
%close all
%A1 = [1 zeros(1,36)];
A1 = [1 0];
C1 = 1;

Mi = idpoly(1, B, C1, A1, A2);

%Mi.Structure.c.Free = [1 zeros(1,35) 1];

z = iddata(ndvi_s, rain_s);
MboxJ = pem(z,Mi);
present(MboxJ)

res = resid(z, MboxJ); % Nice!
whitenessTest(res.y)
%%
plotNTdist(res.y)
%%
plot(res.y)
%% Prediction
%close all;
clc;
k = 1;

A36 = [1 zeros(1,35) -1];
ndvi_test_s = filter(A36, 1, ndvi_test);

z_test = iddata(ndvi_test, rain_test);
y_test = predict(MboxJ, z_test,k);

z = iddata(ndvi, rain);
y = predict(MboxJ, z, 1);
%%
figure(2); clf;
subplot(211)
plot(ndvi_t_test, z_test.y)
hold on
plot(ndvi_t_test(1:end-k), y_test.y(k+1:end),'-')
title('Prediction')
legend('Real','Estimate')

subplot(212)
plot(ndvi_t, z.y)
hold on
plot(ndvi_t(1:end-k), y.y(k+1:end))
title('TRAINING SET')
legend('Real','Estimate')
%%
ehat = resid(z_test,MboxJ);
sum(ehat.y.^2)

figure(3)
subplot(211)
autocorr(ehat1.y)
title('ACF for the resiuduals of the ARMA(1,2) model')
set(gca,'Fontsize',13)
subplot(212)
autocorr(ehat.y)
title('ACF for the residuals of the Box-Jenkins model')
set(gca,'Fontsize',13)
%% Prediction with 4 < k < 10
%close all;
clc;
k = 4:10;
ehat_arma = zeros(length(ndvi_test),length(k));
ehat_bj = ehat_arma;
y_arma = ehat_arma;
y_bj = ehat_arma;

A36 = [1 zeros(1,35) -1];
ndvi_test_s = filter(A36, 1, ndvi_test);
z_test = iddata(ndvi_test, rain_test);
z = iddata(ndvi, rain);

one_step_arma = predict(model_arma, z_test, 1);
one_step_bj = predict(MboxJ, z_test, 1);
var_one_arma = var(one_step_arma.y - ndvi_test);
var_one_bj = var(one_step_bj.y - ndvi_test);

for i = 1:length(k)
   
   tmp_arma = predict(model_arma, z_test, k(i));
   tmp_bj = predict(MboxJ, z_test,k(i));
   y_arma(:,i) = tmp_arma.y;
   y_bj(:,i) = tmp_bj.y;
   
   ehat_arma(:,i) = tmp_arma.y - ndvi_test;
   ehat_bj(:,i) = tmp_bj.y - ndvi_test;
end

limits = [ndvi_t_test(1), ndvi_t_test(end),-2,2];

figure(1);clf;
subplot(211)
plot(ndvi_t_test, z_test.y, '.-')
hold on
plot(ndvi_t_test, one_step_arma.y)
legend('Real','Estimated','Orientation','horizontal')
title('1-step prediction with the ARMA(1,2)-model')
axis(limits)
set(gca,'Fontsize',16)

subplot(212)
plot(ndvi_t_test, z_test.y, '.-')
hold on
plot(ndvi_t_test, one_step_bj.y)
legend('Real','Estimated','Orientation','horizontal')
title('1-step prediction with the Box-Jenkins-model')
axis(limits)
set(gca,'Fontsize',16)

figure(2);clf;
subplot(211)
plot(ndvi_t_test, y_arma)
hold on
plot(ndvi_t_test, z_test.y, '.-')
legend('k=4', 'k=5', 'k=6', 'k=7', 'k=8', 'k=9', 'k=10', 'REAL','Orientation','horizontal')
title('k-step prediction with the ARMA(1,2)-model')
axis(limits)
set(gca,'Fontsize',16)

subplot(212)
plot(ndvi_t_test, y_bj)
hold on
plot(ndvi_t_test, z_test.y, '.-')
legend('k=4', 'k=5', 'k=6', 'k=7', 'k=8', 'k=9', 'k=10', 'REAL', 'Orientation','horizontal')
title('k-step prediction with the Box-Jenkins-model')
axis(limits)
set(gca,'Fontsize',16)


variance = [var_one_arma, var(ehat_arma); var_one_bj, var(ehat_bj)]

%%
figure(1)
subplot(212)
plot(ndvi_t_test, z_test.y)
hold on
plot(ndvi_t_test(1:end-k), y_test.y(k+1:end),'-')
title('Prediction using the Box-Jenkins framework')
legend('Real','Estimate')
%% Try for Kassala with the same model
ndvi = Kassala.nvdi(0.7*end:end);
ndvi_t = Kassala.nvdi_t(0.7*end:end);

k = 2/(max(ndvi)-min(ndvi));
ndvi_transformed = ndvi*k + (1 - k*max(ndvi));

z_test = iddata(ndvi_transformed, kalmanRain(0.7*end:end));
y_test = predict(MboxJ, z_test,k);

k = 1;


figure(2); clf;
plot(ndvi_t, z_test.y)
hold on
plot(ndvi_t(1:end-k), y_test.y(k+1:end))
title('Prediction')
legend('Real','Estimate')

ehat = resid(z_test,MboxJ);
sum(ehat.y.^2)