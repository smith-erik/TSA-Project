cd '/Users/erik/github/TSA-Project/ErikS'
addpath('../matlab')

%%
clear;
close all;
clc;

load('../data/proj18.mat')

% Rename monthly to 'rain' since that is what we are using all the time.
rain = ElGeneina.rain_org;
rain_t = ElGeneina.rain_org_t;

intRain = ElGeneina.rain;
intRain_t = ElGeneina.rain_t;

clear Kassala ElGeneina;

figure(1)
plot(rain_t, rain)

%% A quick AR(1) fit to get a better initial guess for a1
close all;
clc;

model_ar1 = arx(rain, 1);
present(model_ar1)

model2_ar1 = arx(intRain, 1);
present(model2_ar1)

% So fine initial a guess would be inbetween -0.6 and -0.9.
% Negative at least.

% Re from B:
% var(ehat.y) is 119.3986
% std(ehat.y) is 10.9270

% What about Rw? Fine as one?


%% Kalman for an AR(1) process
clearvars -except rain rain_t intRain intRain_t
close all;
clc;

% Simulate/set process
y = rain;
N = length(y);

% Define the state space equations
a1 = -0.7;

A = [-(a1.^3) 0 0;
       a1.^2  0 0;
      -a1     0 0];

C = [1 1 1];

% Hidden state noise covariance matrix
sig_e_2 = 100;
Re = [sig_e_2 0 0;
      0 sig_e_2 0;
      0 0 sig_e_2];

%Re = [sig_e_2 50 20;
%      50 sig_e_2 30;
%      20 30 sig_e_2];


% Observation noise variance
Rw = 1;

% Set some initial values
Rxx_1 = 10 * eye(3);
xtt_1 = [0 0 0]';

% Store values here
xsave = zeros(3,N);


% Kalman filter

for k = 1:N
    
    % Update
    Ryy = C * Rxx_1 * C' + Rw; % Not t+1|t ?
    
    Kt = (Rxx_1 * C') / Ryy; % Same as above but t|t-1 ?
    
    xtt = xtt_1 + Kt*(y(k) - C * xtt_1);
    
    Rxx = (eye(3) - Kt * C) * Rxx_1;
    
    % Save
    xsave(:,k) = xtt;q
    
    % Predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A * xtt;
    
end

xsave_flat = reshape(xsave, 1, length(xsave)*length(xsave(:,1)) );


figure(2)
plot(rain_t, rain, 'k.-') % What do we do about negative?
hold on
plot(intRain_t, xsave_flat, 'b.-')

% Plot
% 1. 30-day sum of Kalman data
% 2. Interpolated 10-day rain data and Kalman 10-day estimate


rainsum = sum(xsave, 1);
figure(3)
plot(rain_t, rain, 'ko')
hold on
plot(rain_t, rainsum, 'b-x')
title('Monthly Sum of Estimations vs True Monthly Precipitation')
set(gca, 'fontsize', 14)
legend('True', 'Estimated', 'Location', 'northeast')



k = 2;
rainPred = zeros(N,1);

% From 1+k since we skip the first k due to prediction
for i = 1+k:N
    rainPred(i) = C * A^k * xsave(:,i-k);
end

rainPred = rainPred(1+k:end);
rain_cut_k = rain(1+k:end);


figure(5)
plot(rain_cut_k, 'k-')
hold on
plot(rainPred, 'b-')
legend('True', 'Estimated', 'Location', 'northeast')

% save('KalmanRain', xsave_flat)




