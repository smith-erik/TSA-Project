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
A = eye(1);
Re = 0.00001; % Hidden state noise covariance matrix
Rw = 10; % Observation noise variance



% Usually C should be set here
% but in this case C is a function of time (previous y).

% Set some initial values
Rxx_1 = 1 * eye(1);
xtt_1 = 0;

% Store values here
xsave = zeros(1,N);

% Kalman filter, start from k = 3 since we need old values of y.
for k = 8:N
   
    % Since C is function of time here, set C = [?] here.
    C = -y(k-1);
    
    % Update
    Ryy = C * Rxx_1 * C' + Rw; % Not t+1|t ?
    
    Kt = (Rxx_1 * C') / Ryy; % Same as above but t|t-1 ?
    
    xtt = xtt_1 + Kt*(y(k) - C * xtt_1);
    
    Rxx = (eye(1) - Kt * C) * Rxx_1;
    
    % Save
    xsave(:,k) = xtt;
    
    % Predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A * xtt;
    
end

figure(3)
plot(-xsave(1,:)) % Why - ?


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

















