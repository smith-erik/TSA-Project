A = [1 -0.3 0 0.5 ];
C = [1 0.5 zeros(1 ,10) -0.3];
A5 = [1 zeros(1 ,4) -1];
Astar = conv(A, A5);
e = randn(600 ,1);
y = filter(C,Astar,e);
y = y(100:end);
plot(y);


y_s = filter(A5,1,y);
data = iddata(y_s);

model_init = idpoly([1 0 0 0],[],[1 zeros(1,12)]);
model_init.Structure.a.Free = [1 1 0 1];
model_init.Structure.c.Free = [1 1 zeros(1,10) 1];
model_armax = pem(data, model_init)
