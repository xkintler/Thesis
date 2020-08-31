function f = myOdesWithInhib(y,param)

    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    s_in = param(7); %g/l
    D = param(8); %h-1

    % s = y(1)
    % x = y(2)
    % p = y(3)
    
    f = [- ((mu_max*y(1)/(Ks + y(1) + ((y(1)^2)/Ki)))*y(2)/Yx) - (nu*y(2)/Yp) + D*(s_in - y(1));...
           (mu_max*y(1)/(Ks + y(1) + ((y(1)^2)/Ki)))*y(2) - D*y(2);...
            nu*y(2) - D*y(3)];
end