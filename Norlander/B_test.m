addpath('../matlab')
clear;
close all;
clc;

load('../data/proj18.mat')

ndvi = ElGeneina.nvdi;
ndvi_t = ElGeneina.nvdi_t;
rain = ElGeneina.rain;
rain_t = ElGeneina.rain_t;

% Throw away data from 1960 to 1982: that's 22*12*3 = 792 points.
rain = rain(793:end);

%clear Kassala ElGeneina;

ndvi_noscale = ndvi;
rain_noscale = rain;

% Rescale to [-1,1]
k = 2/(max(ndvi)-min(ndvi));
ndvi = ndvi*k + (1 - k*max(ndvi));

figure(1); clf;

subplot(211)
plot(ndvi_t,ndvi)
xlabel('Year')
ylabel('NDVI')
set(gca,'Fontsize',12)

subplot(212)
plot(rain_t,ElGeneina.rain)
xlabel('Year')
ylabel('Rain')
set(gca,'Fontsize',12)

%% Split data into 70 % modelling, 30 % test.
N = length(ndvi);

ndvi_test = ndvi(floor(N*0.7)+1 : end);
rain_test = rain(floor(N*0.7)+1 : end);

ndvi_t_test = ndvi_t(floor(N*0.7)+1 : end);
rain_t_test = rain_t(floor(N*0.7)+1 : end);

ndvi = ndvi(1 : floor(N*0.7));
rain = rain(1 : floor(N*0.7));

ndvi_t = ndvi_t(1 : floor(N*0.7));
rain_t = rain_t(1 : floor(N*0.7));

%% See PACF
analyzets(ndvi)
% -> AR(2)

%% Note seasonality in ACF
analyzets(ndvi, 400)
% -> Season 36/3 = 12 months, so a year.


%% Remove season

A36 = [1 zeros(1,35) -1];
ndvi_s = filter(A36,1, ndvi);

%analyzets(ndvi_s)

%figure(2)
%plotNTdist(ndvi_s)

data = iddata(ndvi_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );

model_init.Structure.c.Free = [1 zeros(1,35) 1];

model_ar = pem(data, model_init);
present(model_ar)

%ehat = resid(ndvi_s, model_ar);
%analyzets(ehat);

figure(1); clf;
y = predict(model_ar, ndvi_test,1);
plot(y)
hold on
plot(ndvi_test)

%% Box Jenkins 
A36 = [1 zeros(1,35) -1];
rain_s = filter(A36, 1, rain);

data = iddata(rain_s);

model_init = idpoly([1 0 0 0], [], [1 0 0]);

model_arma = pem(data,model_init);
present(model_arma)

ehat = resid(rain_s,model_arma);

whitenessTest(ehat.y)
%% Pre-whiten rain 
A = model_arma.A;
C = model_arma.C;

rain_pw = filter(A, C, rain_s);
ndvi_pw = filter(A, C, ndvi_s);

analyzets(rain_pw);

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
s = 2;
r = 1;

% idpoly(A,B,C,D,F,NoiseVariance,Ts)
% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)
A2 = [1 zeros(1,r)];

B = [zeros(1,d) 1 zeros(1,s)]; 

Mi = idpoly(1, B, [], [], A2);

% Lock delay coefficients
% Mi.Structure.b.Free = [ zeros(1,d)) 1];

zpw = iddata(ndvi_pw, rain_pw);

Mba2 = pem(zpw, Mi);
present(Mba2)

vhat = resid(Mba2, zpw);
%resid(zpw, Mba2)

x = ndvi_s - filter(Mba2.B, Mba2.F, rain_s);

analyzets(x) % A order 1, C order 0

%% 
close all
A1 = [1 0];
C1 = 1;

Mi = idpoly(1, B, C, A1, A2);
z = iddata(ndvi_s, rain_s);
MboxJ = pem(z,Mi);

present(MboxJ)

resid = resid(z, MboxJ); % Nice!
whitenessTest(resid.y)
%% Prediction
close all;
clc;


z_test = iddata(ndvi_test, rain_test);
y_test = predict(MboxJ, z_test,1);

z = iddata(ndvi_s, rain_s);
y = predict(MboxJ, z, 1);


subplot(211)
plot(ndvi_t_test, z_test.y)
hold on
plot(ndvi_t_test, y_test.y)

subplot(212)
plot(ndvi_t, z.y)
hold on
plot(ndvi_t, y.y)

