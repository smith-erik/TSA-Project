function res=res_id(y,s,koeff)
% res=res_id(y,s,koeff)
% r�knar ut residualen m.a.p ett givet Markov tillst�nd
arp=[1 -koeff];
res=arp*y'/s;

