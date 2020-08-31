close all
clear all
clc

% chemostat parameters
F = 1; %l/h
V = 3.33; %l
D = F/V; %h-1
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ksi = 1.2;
Ki = 350; %g/l
Yx = .4;
Yp = 1;
sin = 20; %g/l
s0 = 0; x0 = 1; p0 = 0;
param = [mu_max nu Ks Ki Yx Yp sin D Ksi];

% figure
s = linspace(0,sin,100);

figure(1)
mus_inh = (mu_max*s./(Ksi + s + ((s.^2)/Ki)));
mus = (mu_max*s./(Ks + s));

hold on
plot(s,mus,'k','LineWidth',2)
plot(s,mus_inh,'b','LineWidth',2)
xlabel('Koncetrácia substrátu s [g/L]')
ylabel('Špecifická rýchlost rastu \mu(s) [h^{-1}]')
%close(1)

%chemostat steady states
% syms s x mu_max Ks D sin Yx nu Yp Ki
% % eqs = [0 == D*(sin - (D*Ks/(mu_max - D))) - ((D/Yx) + (nu/Yp))*x];
% eqs = [0 == D*(sin - (-(Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D))) - ((D/Yx) + (nu/Yp))*x];
% solx = solve(eqs,x)
% simplify(solx)

% x_us = (D*Yp*Yx*(D*sin - mu_max*sin + D*Ks))/((D - mu_max)*(Yx*nu + D*Yp))
% x_inhb_us = (D*Yp*Yx*(sin + (Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D)))/(Yx*nu + D*Yp)

%figure
%biomass steady states with changing dilution rate D

figure(2)
d = linspace(0,0.35,10);
for i = 1:1:length(d)
    D = d(i);
    x_us(1,i) = (D*Yp*Yx*(D*sin - mu_max*sin + D*Ks))/((D - mu_max)*(Yx*nu + D*Yp));
    x_inhb_us(1,i) = (D*Yp*Yx*(sin + (Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D)))/(Yx*nu + D*Yp);
end
subplot(2,1,1)
plot(d,x_us,'bo')
subplot(2,1,2)
plot(d,x_inhb_us,'bo')
close(2)

%numerical solution
alpha = .5;

    param = [mu_max nu Ks Ki Yx Yp sin D Ksi];
    a = sdpvar(1,1);
    obj = fun(a,alpha,param);
    cons = [a >= 0; a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron');
    optimize(cons,obj,ops);
    m_opt_d = value(a);

    a = sdpvar(1,1);
    obj = inh_fun(a,1,param);
    cons = [a >= 0; a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron');
    optimize(cons,obj,ops);
    inh_opt_d = value(a);

    %disp([m_opt_d inh_opt_d])
    D = m_opt_d;
    x_us = (D*Yp*Yx*(D*sin - mu_max*sin + D*Ks))/((D - mu_max)*(Yx*nu + D*Yp));
    D = inh_opt_d;
    x_inhb_us = (D*Yp*Yx*(sin + (Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ksi + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D)))/(Yx*nu + D*Yp);
    %disp([x_us x_inhb_us])
    %opt_sol(i,:) = [m_opt_d inh_opt_d x_us x_inhb_us];



figure(3)
d = linspace(0,m_opt_d,100);
%alpha = linspace(0.01,1,5);
hold on
for j = 1:1:length(alpha)
    for i = 1:1:length(d)
        plot(d(i),fun(d(i),alpha(j),param),'.b')
        plot(d(i),inh_fun(d(i),alpha(j),param),'.k')
    end
end
hold off


function f = inh_fun(D,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(9); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    %Monod-Haldane model
    f = D*(1 - alpha*(D*Yp*Yx*(sin + (Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D)))/(Yx*nu + D*Yp));
end

function f = fun(D,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    %Monod model
    f = D*(1 - alpha*(D*Yp*Yx*(D*sin - mu_max*sin + D*Ks))/((D - mu_max)*(Yx*nu + D*Yp)));
end




