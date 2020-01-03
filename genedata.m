function X = genedata(Tn, n,delta,lamda)
X=[ ];
for T=1:Tn
  v= lognrnd(4,0.5,[1,n]);
  v=[[1:n]',v'];  %%get the valuation for n players
  Blist=zeros(n,1);
  b=0;
  bi=1;%%assume the fist bidding player is 1
  active=1:n;

  while size(active,2)>1
      cho=active(active~=bi);
      i=cho(randi(length(cho)));
      jump=binornd(1,lamda);
      if jump==1
          if b<=v(i,2)
              bi=i;
              b=(v(i,2)-b)*rand+b;%jump to a random varaible between bid and valuation
              Blist(i,1)=b;
          
          else 
              active=active(active~=i);
          end
      end
      
          
      if jump==0
          if  b+delta<=v(i,2)
              bi=i;
              b=b+delta;
              Blist(i,1)=b;
          
          else 
              active=active(active~=i);
        
          end
          
      end
          
      
  end  
  data=[[1:n]',Blist];
    
X(:,:,T)=data;    
    
end


    


    
    

