%% Generar

clear;
close all;
clc;

A = [1 -0.5 0 0.5];
C = [1 0.9 zeros(1,10) -0.5];
A5 = [1 0 0 0 0 -1];
A_star = conv(A,A5);
e = randn(600 ,1);
y = filter(C,A_star,e);
y = y(100:end);

y_s = filter(A5,1,y);

figure(1);
plot(y)

figure(2);
plot(y_s)

%%

data = iddata(y_s);

% A(q) y(t) = [B(q)/F(q)] u(t) + [C(q)/D(q)] e(t)
% model_init = idpoly (A,B,C,D,F);

model_init = idpoly([1 0 0 0] ,[] , [1 zeros(1,12)] );

model_init.Structure.c.Free = [1 1 zeros(1,10) 1];
model_init.Structure.a.Free = [1 1 0 1];


model_armax = pem(data , model_init)

pzmap(model_armax)









