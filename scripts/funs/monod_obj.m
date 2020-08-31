function f = monod_obj(D,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    %Monod model
    s = -(D.*Ks)./(D - mu_max);
    x = -(D.*(s - sin))./(D./Yx + nu/Yp);
    f = D.*(1 - alpha.*x);
    
end