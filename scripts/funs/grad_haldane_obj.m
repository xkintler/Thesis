function f = grad_haldane_obj(D,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    s = -(Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D);
    x = -(D*(s - sin))/(D/Yx + nu/Yp);
    dsdD = (Ki*((2*Ki*mu_max - 2*D*Ki + 8*D*Ks)/(2*Ki*((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)) - 1))/(2*D) + (Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D^2);
    dxdD = (((sin - s) - D*dsdD)*(D/Yx + nu/Yp) - (D*(sin - s)*(1/Yx)))/((D/Yx + nu/Yp)^2);
    f = 1 - alpha*(D*dxdD + x);
end