function pip=filter_start(p,r)
%pip=filter_start(p,r)
%fixar initalv�rde p� filtersannolikheterna enligt den station�ra f�rdelningen
p=p';
M=[eye(r-1) zeros(r-1,1)];
p=p(1:r-1,:);
p=p-M;
p=[p ; ones(1,r)];
a=[zeros(r-1,1); 1];
pip=p\a;
pip=pip';

