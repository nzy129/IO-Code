function G = Gnnd(X,n,Tn,delta, v)
%% this function is to get the nth lowest empirical CDF for lower bound
%%X is the data generated from genedata, has the 3-Dimensional size n*2*Tn,
%%the third demension is Tn

G=0;
for t = 1: Tn
    Xin=mink(X(:,2,t),n);
    if Xin+delta<=v
        G=G+1;
    end
    
end
G=G/Tn;