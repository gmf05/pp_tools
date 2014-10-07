function table = compare_models2(modelH,ms)
% mp = parent model
% mc = child model

[M,N] =find(modelH);
[M,ind] = sort(M);
N = N(ind);
Ncol = 6;
table = cell(length(M)+1,Ncol);
% for each edge, compare pair of models
table{1,1} = 'Models';
table{1,2} = '$\Delta$df';
table{1,3} = '$\Delta$Dev';
table{1,4} = '$p_{Dev}$';
table{1,5} = '$\Delta$AIC';
table{1,6} = '$p_{AIC}$';
for m = 1:length(M)
  m1=ms{M(m)}; m2=ms{N(m)};
  D = compare_model_pair(ms{M(m)},ms{N(m)});
%   D = compare_model_pair(ms{N(m)},ms{M(m)});
  
  table{m+1,1} = [m1.name ' - ' m2.name];
  table{m+1,2} = D(1);
  table{m+1,3} = D(2);
  table{m+1,4} = D(3);
  table{m+1,5} = D(4);
  table{m+1,6} = D(5);
end

end


function D = compare_model_pair(mc,mp)
  
dDev = mp.dev - mc.dev;
dAIC = mp.AIC - mc.AIC;
P1 = size(mc.b,1); P2 = size(mp.b,1); 
F1 = dDev/(P1-P2);
F2 = dAIC/(P1-P2);
pchi1 = 1-chi2cdf(dDev,P1-P2);
pchi2 = 1-chi2cdf(dAIC,P1-P2);
pF1 = 1-fcdf(F1,P1,P1-P2);
pF2 = 1-fcdf(F2,P1,P1-P2);
D = [P1-P2,dDev,pchi1,dAIC,pchi2];

end