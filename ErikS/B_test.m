cd '/Users/erik/github/TSA-Project/ErikS'
addpath('../matlab')

%%
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

clear Kassala ElGeneina;

ndvi_noscale = ndvi;
rain_noscale = rain;

% Rescale to [-1,1]
k = 2/(max(ndvi)-min(ndvi));
ndvi = ndvi*k + (1 - k*max(ndvi));

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

analyzets(ndvi_s)

figure(2)
plotNTdist(ndvi_s)

%% Simple AR(1) model

data = iddata(ndvi_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );

model_init.Structure.c.Free = [1 zeros(1,35) 1];

model_ar1 = pem(data, model_init);
present(model_ar1)

ehat_ar1 = resid(ndvi_s, model_ar1);
analyzets(ehat_ar1);





%% Box Jenkins with rain as input

analyzets(rain)

%% Begin with removing seasonality

A36 = [1 zeros(1,35) -1];
rain_s = filter(A36,1, rain);

analyzets(rain_s) % ARMA(4,4), (3,2) ?

%% Pre-whiten with chosen model
close all;
clc;

data = iddata(rain_s);

model_init = idpoly([1 0 0 0], [] , [1 0 0] );
model_arma = pem(data, model_init);
present(model_arma)

ehat = resid(rain_s, model_arma);

analyzets(ehat);

figure(2);
whitenessTest(ehat.y)

% Seems white enough


%%
close all;
clc;


A3 = model_arma.A;
C3 = model_arma.C;

upw = filter(A3, C3, rain_s);
ypw = filter(A3, C3, ndvi_s);


%%

% Data is now white, estimate impulse response weights
% to find order of B and A2.

M = 40;
stem(-M:M, crosscorr(upw, ypw, M) );
title('Crosscorr function');
xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(upw))*ones(1,2*M+1), '--')
plot(-M:M, -2/sqrt(length(upw))*ones(1,2*M+1), '--')
hold off

% OK, wtf... Men säg d=1 då 10 dagar.
% Sen s = 2, r = 1, helt omotiverat.


%%
close all;
clc;

d = 1;
s = 2;
r = 0;

A2 = [1 zeros(1,r) ];

B = [zeros(1,d) 1 zeros(1,s)]; % s = 0, 1, 2?

% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)
% idpoly(A,B,C,D,F,NoiseVariance,Ts)

Mi = idpoly(1, B, [], [], A2);

zpw = iddata(ypw, upw);

Mba2 = pem(zpw, Mi);
present(Mba2)

vhat = resid(Mba2, zpw);
resid(zpw, Mba2) 

%%
x = ndvi_s - filter(Mba2.B, Mba2.F, rain_s);

analyzets(x) % A as order 1, C as order 0?

%%
close all;
clc;

A1 = [1 0];
C1 = 1;

Mi = idpoly(1, B, C1, A1, A2);

z = iddata(ndvi_s, rain_s);
MboxJ = pem(z, Mi);

present(MboxJ)

figure(2);
resids = resid(z, MboxJ);
resid(z, MboxJ)

figure(1);
analyzets(resids); % Nice!

figure(3)
whitenessTest(resids.y) % Pretty good!


%% Predict
close all;
clc;

k = 4;

predMboxJ = predict(MboxJ, z, k) % Box - Jenkins prediction

predAR1 = predict(model_ar1, data, k) % Simple AR(1) prediction


figure(4)
plot(ndvi_t, ndvi_s, 'b-');
hold on
plot(ndvi_t, predMboxJ.y, 'r-')
plot(ndvi_t, predAR1.y, 'k-')























