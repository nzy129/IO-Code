%% Created by Zheyu Ni for Take Home Exam for IO8872 
%% Replicate Haile and Tamer

n=6; Tn=25; delta=0.05; %%or 1 
v=linspace(0,200);
Times=500;
zz1=ones(Times,length(v));
llb=zeros(Times,length(v));
zzm=zz1;
%lamda=0;
lamda=0.25;



for time = 1:Times
    
    X = genedata(Tn, n,delta,lamda);
% For upper bound   
    for tt =1:length(v)
        zz1(tt)=inverse(Gin(X,1,Tn,v(tt)),1,n);%%use self written inverse function
        zz2=inverse(Gin(X,2,Tn,v(tt)),2,n);
        zz3=inverse(Gin(X,3,Tn,v(tt)),3,n);
        zz4=inverse(Gin(X,4,Tn,92),4,n);
        disp(v(tt));
        zz5=inverse(Gin(X,5,Tn,v(tt)),5,n);
        zz6=inverse(Gin(X,6,Tn,v(tt)),6,n);
        if zz2<zzm(time,tt)
            zzm(time,tt)=zz2;
        end
   % if zz3<zzm(tt)
        %zzm(tt)=zz3;
   % end
        if zz4<zzm(time,tt)
            zzm(time,tt)=zz4;
        end
    %if zz5<zzm(tt)
       % zzm(tt)=zz5;
   % end
        if zz6<zzm(time,tt)
            zzm(time,tt)=zz6;
        end
    
    end

%For lower bound    
    for ttt =1:length(v)
        ii=5;
        fun=@(x)(x.^(ii-1)).*((1-x).^(n-ii));
        Huge=@(z)(factorial(n)/(factorial(n-ii))/(factorial(ii-1))*integral(fun,0,z))-Gnnd(X,6,Tn,delta,v(ttt));
        f=fzero(Huge,[0,1]);
        llb(time,ttt)=f;
        
    end

 
    
    
end

    
    

%zz=inverse(Gin(X,i,Tn,200),i,n);
    
Upaverage=mean(zzm,1);
nfifthup=prctile(zzm,95);
Loaverage=mean(llb,1);
fifthlow=prctile(llb,5);
plot(v,Loaverage);

%plot upper bound
figure(1)

plot(v,Upaverage);
hold on
plot(v,Loaverage);
hold on
plot(v,nfifthup);
hold on
plot(v,fifthlow);
%plot lognormal cdf
plot(v,logncdf(v,4,0.5));
hold off


%other model 2 with lamda=0
Times=500;
model_1=zeros(Times,length(v));
for t = 1:Times
    X = genedata(Tn, n,delta,0);
    for mm=1:length(v)
        model_1(t,mm)=Gin(X,5,Tn,v(mm));
    end
modelm=mean(model_1,1);
model95=prctile(model_1,95);
model5=prctile(model_1,5);
end
figure(2)
subplot(3,2,1);
title('delta=0.05,lamda=0 Tn=200')
plot(v,modelm);
hold on 
plot(v,model95);
hold on
plot(v,model5);
hold on
plot(v,logncdf(v,4,0.5));
hold off

subplot(3,2,2);
title('delta=0.05,lamda=0 Tn=50')
%other model 1 with lamda 0
Times=500;
model_12=zeros(Times,length(v));
for t = 1:Times
    X = genedata(50, n,delta,0);
    for mm=1:length(v)
        model_12(t,mm)=Gin(X,5,50,v(mm));
    end
modelm2=mean(model_12,1);
model952=prctile(model_12,95);
model52=prctile(model_12,5);
end
plot(v,modelm2);
hold on 
plot(v,model952);
hold on
plot(v,model52);
hold on
plot(v,logncdf(v,4,0.5));
hold off

subplot(3,2,3);
title('delta=0.05,lamda=0.25 Tn=200')
%other model 1 with lamda =0.25
Times=500;
model_12=zeros(Times,length(v));
for t = 1:Times
    X = genedata(200, n,delta,0.25);
    for mm=1:length(v)
        model_12(t,mm)=Gin(X,5,200,v(mm));
    end
modelm2=mean(model_12,1);
model952=prctile(model_12,95);
model52=prctile(model_12,5);
end
plot(v,modelm2);
hold on 
plot(v,model952);
hold on
plot(v,model52);
hold on
plot(v,logncdf(v,4,0.5));
hold off

subplot(3,2,4);
title('delta=0.05,lamda=0.25 Tn=50')
%other model 1 with lamda =0.25
Times=500;
model_12=zeros(Times,length(v));
for t = 1:Times
    X = genedata(50, n,delta,0.25);
    for mm=1:length(v)
        model_12(t,mm)=Gin(X,5,50,v(mm));
    end
modelm2=mean(model_12,1);
model952=prctile(model_12,95);
model52=prctile(model_12,5);
end
plot(v,modelm2);
hold on 
plot(v,model952);
hold on
plot(v,model52);
hold on
plot(v,logncdf(v,4,0.5));
hold off

subplot(3,2,5);
title('delta=0.05,lamda=F(v)^4 Tn=200')
%other model 1 with F(v)^4
Times=500;
model_12=zeros(Times,length(v));
for t = 1:Times
    X = genedata2(200, n,delta);
    for mm=1:length(v)
        model_12(t,mm)=Gin(X,5,200,v(mm));
    end
modelm2=mean(model_12,1);
model952=prctile(model_12,95);
model52=prctile(model_12,5);
end
plot(v,modelm2);
hold on 
plot(v,model952);
hold on
plot(v,model52);
hold on
plot(v,logncdf(v,4,0.5));
hold off

subplot(3,2,6);
title('delta=0.05%,lamda=0.25 Tn=200')
%other model 1 with lamda =0.25
Times=500;
model_12=zeros(Times,length(v));
for t = 1:Times
    X = genedata(200, n,0.0005,0.25);
    for mm=1:length(v)
        model_12(t,mm)=Gin(X,5,200,v(mm));
    end
modelm2=mean(model_12,1);
model952=prctile(model_12,95);
model52=prctile(model_12,5);
end
plot(v,modelm2);
hold on 
plot(v,model952);
hold on
plot(v,model52);
hold on
plot(v,logncdf(v,4,0.5));
hold off
