
addpath('../matlab')

close all;
clear;
clc;

%%
close all;
clear;
clc;

A1 = [ 1 -1.79 0.84 ];
C1 = [ 1 -0.18 -0.11 ];

A2 = [ 1 -1.79 ];
C2 = [ 1 -0.18 -0.11 ];

arma1_p = idpoly(A1, [], C1);
arma2_p = idpoly(A2, [], C2);

N = 200;
e = sqrt(2) * randn(N,1);

y1 = filter(arma1_p.c, arma1_p.a, e);
y2 = filter(arma2_p.c, arma2_p.a, e);

figure(1);
subplot(2,1,1)
plot(y1)
subplot(2,1,2)
plot(y2)

figure(2);
pzmap(arma1_p)
title('Process 1')

figure(3);
pzmap(arma2_p)
title('Process 2') % Q1: This is unstable

%%
m = 50;
sigma2 = 1.5;

rtheo = kovarians(arma1_p.c, arma1_p.a, m);
figure(4);
stem(0:m, rtheo*sigma2)
hold on;
rest = covf(y1, m+1); % Obsolete?
stem(0:m, rest, 'r')
legend('Theoretical', 'Estimated')

% Q2.1: Because cvf is biased.
% Q2.2: Rule of thumb is N/4, i.e. 50. Our bad after ~10-25.




%% 

% Plot -1? Börja på 0?

figure(5);
maxOrd = 50;

subplot(3,1,1)
stem(acf(y1, maxOrd));
title('ACF')

subplot(3,1,2)
stem(pacf(y1,maxOrd));
title('PACF')

subplot(3,1,3)
normplot(y1)


%%

data = iddata(y1);


ar_model = arx(y1, 2);
arma_model = armax(y1, [2 2]); % Ask about input order, [a b c]?

present(ar_model)
present(arma_model)

% FPE nära varanda?

%%

e_hat = filter(arma1_p.a, arma1_p.c, y1);

figure(7);
plot(e_hat)

figure(8);
maxOrd = 50;

subplot(3,1,1)
stem(0:maxOrd , acf(e_hat, maxOrd));
title('ACF')

subplot(3,1,2)
stem(0:maxOrd , pacf(e_hat,maxOrd));
title('PACF')

subplot(3,1,3)
normplot(e_hat)

% Q3: Typ samma egentligen, men borde ju vara den "rätta".


