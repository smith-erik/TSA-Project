cd '/Users/erik/github/TSA-Project/ErikS'
addpath('../matlab')

%%
clear;
close all;
clc;

load('../data/proj18.mat')

rain = ElGeneina.rain_org;
rain_t = ElGeneina.rain_org_t;

dumb_rain = ElGeneina.rain;
dumb_rain_t = ElGeneina.rain_t;

clear Kassala ElGeneina;


plot(rain_t, rain)

%% Kalman for an AR(1) process
clearvars -except rain rain_t
close all;
clc;

% Simulate/set process
y = rain;
N = length(y);

% Define the state space equations
a1 = 0.7;
A = [-a1 0 0; 1 0 0; 0 1 0];
C = [1 1 1];
Re = [1 0 0; 0 0 0; 0 0 0]; % Hidden state noise covariance matrix
Rw = 1; % Observation noise variance

% Set some initial values
Rxx_1 = 1 * eye(3);
xtt_1 = [0 0 0]';

% Store values here
xsave = zeros(3,N);

% Kalman filter, start from k = 3 since we need old values of y.
for k = 1:N
    
    % Update
    Ryy = C * Rxx_1 * C' + Rw; % Not t+1|t ?
    
    Kt = (Rxx_1 * C') / Ryy; % Same as above but t|t-1 ?
    
    xtt = xtt_1 + Kt*(y(k) - C * xtt_1);
    
    Rxx = (eye(3) - Kt * C) * Rxx_1;
    
    % Save
    xsave(:,k) = xtt;
    
    % Predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A * xtt;
    
end

figure(3)
plot(xsave(1,:)) % What do we do about negative?


%%



rainsum = zeros(N,1);

for i = 3:N
    rainsum(i) = xsave(i) + xsave(i-1) + xsave(i-2);
end

rainsum = sum(xsave, 1);

%rainflat = [xsave() ]


figure(2)
plot(rainsum, 'b-')
hold on
plot(y, 'r-')
legend('True', 'Estimated', 'Location', 'northeast')






%%
% Do k-step prediction
k = 1;
yPred = zeros(N-k,1);

for i = 1:N-k
    yPred(i) = y(i) * -xsave(i+k);
end

figure(4)
plot(y)
hold on
plot(yPred)
legend('True', 'Estimated', 'Location', 'northeast')

















