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

y = rain;
N = length(y);
a1 = -0.7;

RwCoeffs = linspace(-5, 5, 10000);
SSRs = zeros(length(RwCoeffs),1) - 1;
idx = 1;

for Rw = RwCoeffs
    

    % Define the state space equations
    A = [-(a1.^3) 0 0;
           a1.^2  0 0;
          -a1     0 0];

    C = [1 1 1];

    % Hidden state noise covariance matrix
    sig_e_2 = 100;
    Re = [sig_e_2 0 0;
          0 sig_e_2 0;
          0 0 sig_e_2];

    % Observation noise variance
    % Rw = 1;

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
    
    % Some values are < 0, not possible.
    % Difference in final result is negligeable.
    xsave(xsave < 0) = 0;
    
    % Calculate and save squared sum of residuals
    kalmanMonthSum = sum(xsave, 1)';
    SSRs(idx) = sum( (rain - kalmanMonthSum).^2 );
    idx = idx + 1;

end

figure(2)
plot(RwCoeffs, SSRs);
title('Squared sum of residuals between true monthly and monthly sum of 10-day estimates')
xlabel('R_w value')
%ylabel('Accumulated rain  [mm]')

[RwMin, RwMinIdx] = min(SSRs);

fprintf('Min: %f.\nFor R_w = %f.\n', RwMin, RwCoeffs(RwMinIdx))



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




