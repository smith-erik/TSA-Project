addpath('/Users/eriknorlander/Documents/MATLAB/FMSN45/downloadMaterialTSA/data')
addpath('/Users/eriknorlander/Documents//MATLAB/FMSN45/downloadMaterialTSA/matlab')
EG = struct2cell(ElGeneina);
KA = struct2cell(Kassala);
%%
close all
figure(1)
subplot(311)
plot(EG{2,1},EG{1,1})

subplot(312)
plot(EG{4,1},EG{3,1})

subplot(313)
plot(EG{6,1},EG{5,1})

figure(2)
subplot(311)
plot(KA{2,1},KA{1,1})

subplot(312)
plot(KA{4,1},KA{3,1})

subplot(313)
plot(KA{6,1},KA{5,1})

% Modeling as an AR(1)-process
A = [1 1];
C = [];
AR_poly = idpoly(A, [], C);
AR_poly.a
