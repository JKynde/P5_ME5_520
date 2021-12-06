function [DeltaZ] = Zfit(f,V_in,param,Z_ex)
Z_sim=zeros(1,length(f));
Z_dif=zeros(1,length(f));
for n=1:length(f)
    [~,~,Z_sim(n)]=Matricer2(f(n),V_in,param);
    Z_sim(n)=abs(Z_sim(n));
    Z_dif(n)=abs(Z_sim(n)-Z_ex(n));
end

DeltaZ=sum(Z_dif);

end

