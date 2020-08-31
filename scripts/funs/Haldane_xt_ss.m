function ixt = Haldane_xt_ss(D,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    s = -(Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D);
    ixt = -(D*(s - sin))/(D/Yx + nu/Yp);
end