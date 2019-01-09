addpath('../matlab')
addpath('../data')
clear;
close all;
clc;

load('../data/proj18.mat')
%% A - Kalman filter reconstruction of data
rain = Kassala.rain_org;
rain_t = Kassala.rain_org_t;
rain_long = Kassala.rain;
rain_t_long = Kassala.rain_t;

% Simulate/set process
y = rain;
N = length(y);

% Define the state space equations
a1 = 0.7;
sig_e = 10;
A = [-a1^3 0 0; 
      a1^2 0 0; 
      a1^1 0 0];
Re = sig_e*eye(3); % Hidden state noise covariance matrix
Rw = 1; % Observation noise variance

C = [1 1 1];

% Set some initial values
Rxx_1 = 1 * eye(3);
xtt_1 = [0 0 0]';

% Store values here
xsave = zeros(3,N);

for k = 1:N
    
    % Update
    Ryy = C * Rxx_1 * C' + Rw; 
    
    Kt = (Rxx_1 * C') / Ryy; 
    
    xtt = xtt_1 + Kt*(y(k) - C * xtt_1);
    
    Rxx = (eye(3) - Kt * C) * Rxx_1;
    
    % Save
    xsave(:,k) = xtt;
    
    % Predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A * xtt;
    
end

rainsum = sum(xsave,1);
%raincum = cumsum(xsave,1);
rainflat = reshape(xsave,1,3*480);
%%
figure(1); clf;
subplot(211)
plot(rain_t,rainsum,'o')
hold on 
plot(rain_t,rain)

subplot(212)
plot(rain_t_long,rainflat,'.-');
hold on
plot(rain_t_long(2:end),rain_long(1:end-1),'.-')

figure(2); clf;
plot(rain_t,rain,'o')
hold on
plot(rain_t_long,rain_long,'.-')
title('Sum of Kalman reconstruction compared with actual rain.')

clear rain rain_t


%% B - Modelling an validation
close all
ndvi = Kassala.nvdi;
ndvi_t = Kassala.nvdi_t;
rain = rainflat;
%rain = rainlog;
rain_t = rain_t_long;
%rain = Kassala.rain;
%rain_t = Kassala.rain_t;

% Throw away data from 1960 to 1982: that's 22*12*3 = 792 points.
rain = rain(793:end);

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
%figure(2)
%plotNTdist(ndvi_s)

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
data = iddata(rain_s');

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
s = 1;
r = 0;

% idpoly(A,B,C,D,F,NoiseVariance,Ts)
% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)
A2 = [1 zeros(1,r)];

B = [zeros(1,d) 1 zeros(1,s)]; 

Mi = idpoly(1, B, [], [], A2);

% Lock delay coefficients
Mi.Structure.b.Free = [ zeros(1,d) 1 ones(1,s)];

zpw = iddata(ndvi_pw, rain_pw');

Mba2 = pem(zpw, Mi);
present(Mba2)

%vhat = resid(Mba2, zpw);
%resid(zpw, Mba2)

x = ndvi_s - filter(Mba2.B, Mba2.F, rain_s)';

analyzets(x') % A order 1, C order 0

%% 
close all
A1 = [1 0];
C1 = 1;

Mi = idpoly(1, B, C, A1, A2);
z = iddata(ndvi_s, rain_s');
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


z_test = iddata(ndvi_test, rain_test');
y_test = predict(MboxJ, z_test,1);

z = iddata(ndvi_s, rain_s');
y = predict(MboxJ, z, 1);


subplot(211)
plot(ndvi_t_test(1:end), z_test.y(1:end))
hold on
plot(ndvi_t_test, y_test.y)

subplot(212)
plot(ndvi_t, z.y)
hold on
plot(ndvi_t, y.y)

%% Prediction with season:ed test-data
close all;
clc;
A36 = [1 zeros(1,35) -1];
rain_test_s = filter(A36, 1, rain_test);
ndvi_test_s = filter(A36, 1, ndvi_test);
z_test = iddata(ndvi_test_s, rain_test_s');
y_test = predict(MboxJ, z_test,1);

z = iddata(ndvi_s, rain_s');
y = predict(MboxJ, z, 1);

subplot(211)
plot(ndvi_t_test(1:end), z_test.y(1:end))
hold on
plot(ndvi_t_test, y_test.y)
legend('Actual Value', 'Estimation')

subplot(212)
plot(ndvi_t, z.y)
hold on
plot(ndvi_t, y.y)
legend('Actual Value', 'Estimation')