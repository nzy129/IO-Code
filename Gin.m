function G = Gin(X,i,Tn,v)
%% this function is to get the ith lowest empirical CDF
%%X is the data generated from genedata, has the 3-Dimensional size n*2*Tn,
%%the third demension is Tn

G=0;
for t = 1: Tn
    Xin=mink(X(:,2,t),i);
    if Xin<=v
        G=G+1;
    end
    
end
G=G/Tn;

