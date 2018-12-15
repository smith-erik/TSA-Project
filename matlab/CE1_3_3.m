addpath('../data')

%%

close all;
clear;
clc;

load('data.dat')
load('noise.dat')

data = iddata(data);


ar1_model = arx(data, 1);
ar2_model = arx(data, 2);
% ar3_model = arx(data, 3);

present(ar1_model)
present(ar2_model)
% present(ar3_model)


% resid(ar1_model, data);
% resid(ar2_model, data);
% resid(ar3_model, data);
rar1 = resid(ar1_model, data);
rar2 = resid(ar2_model, data);
rar3 = resid(ar3_model, data);