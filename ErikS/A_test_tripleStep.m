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
sig_e_2 = 119;
Re = [sig_e_2 0 0;
      0 sig_e_2 0;
      0 0 sig_e_2];

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
    xsave(:,k) = xtt;
    
    % Predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A * xtt;
    
end


xflat = reshape(xsave, 1, length(xsave)*length(xsave(:,1)) );
kalmanRain = reshape(xsave, 1, length(xsave)*length(xsave(:,1)) );
kalmanRain(kalmanRain < 0) = 0;

% Some values are < 0, not possible.
% Difference in final result is negligeable.
xsave(xsave < 0) = 0; 

% Plot
% 1. 30-day sum of Kalman data
% 2. Interpolated 10-day rain data and Kalman 10-day estimate

figure(2)
plot(kalmanRain)

rainsum = sum(xsave, 1);
figure(3)
plot(rain, 'b-')
hold on
plot(rainsum, 'r-.')
legend('True', 'Estimated', 'Location', 'northeast')


%% Kalman reconstruction vs interpolation
close all;
clc;

xsave_flat = reshape(xsave, 1, length(xsave)*length(xsave(:,1)) );
kalmanMonthSum = sum(xsave, 1)';

intMonthSum = zeros(length(rain), 1);
k = 1;
for i = 3:3:length(intRain)
    intMonthSum(k) = intRain(i) + intRain(i-1) + intRain(i-2);
    k = k + 1;
end

% Sum the squared difference between true monthly and estimated monthly
SSR_Kalman = sum( (rain - kalmanMonthSum).^2 );
SSR_Interpolated = sum( (rain - intMonthSum).^2 );

fprintf('Squared sum of resuduals for the monthly sum from ...\n')
fprintf('Kalman: %f.\nInterpolated: %f.\n', SSR_Kalman, SSR_Interpolated)

fprintf('\nTotal rain ...\n')
fprintf('True: %f.\n', sum(rain))
fprintf('Kalman: %f.\nInterpolated: %f.\n', sum(xsave_flat), sum(intRain))



figure(4)
plot(rain_t, rain, 'ko')
hold on
plot(intRain_t, xsave_flat, 'b:')
plot(intRain_t, intRain, 'r:')
plot(rain_t, kalmanMonthSum, 'bx')
plot(rain_t, intMonthSum, 'rx')

legend('True monthly total', 'Kalman Estimated', 'Interpolated', 'Monthly sum from Kalman', ...
       'Monthly sum from interpolated', 'Location', 'northeast')
set(gca, 'fontsize', 12)
set(gcf, 'position', [100  200 1200 600])


%% Save interpolation in data folder

save('../data/kalmanRain.mat', 'kalmanRain')




