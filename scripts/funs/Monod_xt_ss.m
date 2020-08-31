function xt = Monod_xt_ss(D,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    s = -(D.*Ks)./(D - mu_max);
    xt = -(D.*(s - sin))./(D./Yx + nu/Yp);
end