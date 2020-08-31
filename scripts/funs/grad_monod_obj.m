function f = grad_monod_obj(D,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    s = -(D.*Ks)./(D - mu_max);
    x = -(D*(s - sin))/(D/Yx + nu/Yp);
    dsdD = (D*Ks)/(D - mu_max)^2 - Ks/(D - mu_max);
    dxdD = (((sin - s) - D*dsdD)*(D/Yx + nu/Yp) - (D*(sin - s)*(1/Yx)))/((D/Yx + nu/Yp)^2);
    f = 1 - alpha*(D*dxdD + x);
end
