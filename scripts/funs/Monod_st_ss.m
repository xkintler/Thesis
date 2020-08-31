function st = Monod_st_ss(D,param)
    mu_max = param(1); %h-1
    Ks = param(3); %g/l
    
    st = -(D.*Ks)./(D - mu_max);
end