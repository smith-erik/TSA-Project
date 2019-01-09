addpath('../matlab')
addpath('../data')
clear;
close all;
clc;

load('../data/proj18.mat')

rain = ElGeneina.rain_org;
rain_t = ElGeneina.rain_org_t;
rain_long = ElGeneina.rain;
rain_t_long = ElGeneina.rain_t;



% Simulate/set process
y = rain;
N = length(y);

% Define the state space equations
a1 = 0.7;
sig_e = 10;
A = [-a1^3 0 0; a1^2 0 0; a1 0 0];
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

figure(1); clf;
subplot(211)
plot(rain_t,rainsum,'o')
hold on 
plot(rain_t,rain)

subplot(212)
plot(rain_t_long,rainflat,'.-');
hold on
plot(rain_t_long(2:end),rain_long(1:end-1),'.-')
%%
figure(2);
plot(rain_t,rain,'o')
hold on
plot(rain_t_long,rain_long,'.-')