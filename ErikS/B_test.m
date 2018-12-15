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

ndvi = ndvi(1 : floor(N*0.7));
rain = rain(1 : floor(N*0.7));

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

%%

data = iddata(ndvi_s);

model_init = idpoly([1 0], [] , [1 zeros(1,36)] );

model_init.Structure.c.Free = [1 zeros(1,35) 1];

model_ar = pem(data, model_init);
present(model_ar)

ehat = resid(ndvi_s, model_ar);
analyzets(ehat);











%% Box Jenkins with 













