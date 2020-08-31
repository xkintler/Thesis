function f = hybrid_obj_sub_corr(D,us,ys,alpha,param,P)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
     
    s_us = -(Ki.*(D - mu_max + (((D.^2).*Ki - 4.*(D.^2).*Ks + Ki*mu_max^2 - 2.*D.*Ki.*mu_max)./Ki).^(1/2)))./(2.*D);
    s_hy = (D-us)*sum(P) + ys;
    s = s_us + s_hy;
    x = -(D.*(s - sin))./(D./Yx + nu/Yp);
    f = D.*(1 - alpha.*x);
end