function [xs,ss,ps] = get_Monod_ss(D,params)
    mu_max = params(1);
    nu = params(2);
    Ks = params(3);
    Ki = params(4);
    Yx = params(5);
    Yp = params(6);
    s_in = params(7);
    
    ss = -(D.*Ks)./(D - mu_max);
    xs = -(D.*(ss - s_in))./(D./Yx + nu/Yp);
    ps = xs.*nu./D;
end