clc
clear
param=StructCreator();
f=logspace(1,7,5000);
Z=zeros(1,length(f));
V_in=150;
for n=1:length(f)
   
    
    Z(n)=cablemodel(f(n),V_in,param);


end
plot(f,abs(Z))