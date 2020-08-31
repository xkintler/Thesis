% uprava opt D pomocou s_us + s_dat
% -II- pomocou mus = mus(s_us) + mus(s_dat)
% mozno uprava gradientu ucelovej fcie

clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')

F = 1; %l/h
V = 3.33; %l
D = 0.33; %h-1
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 350; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l
tf = 50;
s0 = 0; x0 = 5; p0 = 0;
is0 = 0; ix0 = 5; ip0 = 0;

D = [0.2 .25 0.33];
st = []; xt = []; pt = [];
ist = []; ixt = []; ipt = [];

figure(1)
t = linspace(0,tf,50);
for i = 1:1:length(D)
    param = [mu_max nu Ks Ki Yx Yp s_in D(i)];
    iS = ode45(@(t,y) myOdesWithInhib(y,param),[0 tf],[is0;ix0;ip0]);
    ist = [ist deval(iS,t,1)];
    ixt = [ixt deval(iS,t,2)];
    ipt = [ipt deval(iS,t,3)];
    
    S = ode45(@(t,y) myOdes(y,param),[0 tf],[s0;x0;p0]);
    st = [st deval(S,t,1)];
    xt = [xt deval(S,t,2)];
    pt = [pt deval(S,t,3)];
    
    is0 = ist(end); ix0 = ixt(end); ip0 = ipt(end);
    s0 = st(end); x0 = xt(end); p0 = pt(end);
end

t = linspace(0,i*tf,i*50);    
subplot(2,1,1)
hold on
plot(t,ist,'-b')
plot(t,st,'-r')
xlabel('Cas t [h]')
ylabel({'Koncentrácia substrátu';'s(t) [g/L]'})
box on
subplot(2,1,2)
hold on
plot(t,ixt,'-b')
plot(t,xt,'-r')
xlabel('Cas t [h]')
ylabel({'Koncentrácia biomasy';'x(t) [g/L]'})
box on
hold off
close(1)

% FIR identification
ds = st - ist;
%ds = ds(51:end);
%ds = ds - ds(1);
%D = 0.33;
u = ones([50 1])*D;
u = u(:)';
ds = ds';
model_error = 0.001;
rad_p = 1;
C = zeros(rad_p);
v = 1;

while v==1
   C = napln_maticu(u',rad_p);
   s = size(C);
   n = zeros(1,s(2));
   v = isempty(linprog(n,[-C;C],[model_error-ds; model_error+ds]));
   if v == 1
       rad_p = rad_p + 1;
   else
       rad_p;   
   end
   if rad_p == 100
       v = 0;
   end
end

P = inv(C'*C)*C'*ds;
y_approx = C*P;

figure(2)
hold on
stairs(y_approx,'b')
stairs(ds,'r')
xlabel('Cas t [h]')
ylabel('\Delta s(t) [g/L]')
%close(2)

figure
hold on
stairs(st,'r')
stairs(ist,'b')
stairs(ist+y_approx','k')

figure 
hold on
mus = mu_max*st./(Ks + st);
mus_s_corr = mu_max*(ist+y_approx')./(Ks + (ist+y_approx') + ((ist+y_approx').^2./Ki));
mus_mus_corr = mu_max*ist./(Ks + ist + (ist.^2./Ki)) + mu_max*y_approx'./(Ks + y_approx' + (y_approx'.^2./Ki));
plot(mus,'b')
plot(mus_s_corr,'r')
plot(mus_mus_corr,'k')


return
% hybrid model implementation
% data-base model steady state addtion
% to substrate concentration
alpha = 0.5;

a = sdpvar(1,1);
obj = fun(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron');
optimize(cons,obj,ops);
m_opt_d = value(a);

a = sdpvar(1,1);
obj = substrate_corr(a,zeros(size(P)),1,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron');
optimize(cons,obj,ops);
inh_opt_d = value(a);

opt_sol(1,:) =  [m_opt_d inh_opt_d];

a = sdpvar(1,1);
obj = substrate_corr(a,P,1,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron');
optimize(cons,obj,ops);
inh_opt_d = value(a);

opt_sol(2,:) =  [m_opt_d inh_opt_d];

a = sdpvar(1,1);
obj = mus_corr(a,P,1,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron');
optimize(cons,obj,ops);
inh_opt_d = value(a);

opt_sol(3,:) =  [m_opt_d inh_opt_d];

disp('Comparison - normal values')
disp(opt_sol(1,:))
disp('Substrate steady state correction')
disp(opt_sol(2,:))
disp('Specific growth rate correction')
disp(opt_sol(3,:))
%-------------------------------------------------------------------
%-------------------------------------------------------------------
% syms D s sin Yx nu Yp a Ks mu_max
% s = -(D*Ks)/(D - mu_max);
% x = -(D*(s + a*D - sin))/(D/Yx + nu/Yp);
% dxdD = diff(x,D)
% %dfdD = 1 - a*(x + D*dxdD) + aD;
% %sol = solve(dfdD,D);
% sol = eval(subs(dxdD,[Yx nu Yp Ks mu_max],[.4 .5 1 1.2 .53]));
% simplify(sol)
% pretty(sol)

function f = myOdes(y,param)

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
    
    f = [- ((mu_max*y(1)/(Ks + y(1)))*y(2)/Yx) - (nu*y(2)/Yp) + D*(s_in - y(1));...
           (mu_max*y(1)/(Ks + y(1)))*y(2) - D*y(2);...
            nu*y(2) - D*y(3)];
end

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

function f = substrate_corr(D,P,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    %Monod-Haldane model
    n = size(P);
    s_us = -(Ki.*(D - mu_max + (((D.^2).*Ki - 4.*(D.^2).*Ks + Ki*mu_max^2 - 2.*D.*Ki.*mu_max)./Ki).^(1/2)))./(2.*D);
    s_hy = ones([n(2) n(1)]).*D*P;
    s = s_us + s_hy;
    x = -(D.*(s - sin))./(D./Yx + nu/Yp);
    f = D.*(1 - alpha.*x);
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

function f = mus_corr(D,P,alpha,param)
    mu_max = param(1); %h-1
    nu = param(2); %h-1 
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    Yx = param(5);
    Yp = param(6);
    sin = param(7); %g/l
    
    %Monod-Haldane model
    n = size(P);
    ds = ones([n(2) n(1)]).*D*P;
    s = -(Ki*(D*ds^2 - ds^2*mu_max + ((D^2*Ki^3*Ks^2 + 2*D^2*Ki^3*Ks*ds + D^2*Ki^3*ds^2 - 4*D^2*Ki^2*Ks^3 - 8*D^2*Ki^2*Ks^2*ds - 2*D^2*Ki^2*Ks*ds^2 + 2*D^2*Ki^2*ds^3 - 8*D^2*Ki*Ks^2*ds^2 - 8*D^2*Ki*Ks*ds^3 + D^2*Ki*ds^4 - 4*D^2*Ks*ds^4 - 2*D*Ki^3*Ks^2*mu_max - 6*D*Ki^3*Ks*ds*mu_max - 4*D*Ki^3*ds^2*mu_max + 8*D*Ki^2*Ks^2*ds*mu_max + 4*D*Ki^2*Ks*ds^2*mu_max - 6*D*Ki^2*ds^3*mu_max + 8*D*Ki*Ks*ds^3*mu_max - 2*D*Ki*ds^4*mu_max + Ki^3*Ks^2*mu_max^2 + 4*Ki^3*Ks*ds*mu_max^2 + 4*Ki^3*ds^2*mu_max^2 - 2*Ki^2*Ks*ds^2*mu_max^2 + 4*Ki^2*ds^3*mu_max^2 + Ki*ds^4*mu_max^2)/Ki)^(1/2) + D*Ki*Ks + D*Ki*ds - Ki*Ks*mu_max - 2*Ki*ds*mu_max))/(2*(D*ds^2 + D*Ki*Ks + D*Ki*ds - Ki*ds*mu_max));
    x = -(D*(s - sin))/(D/Yx + nu/Yp);
    f = D*(1 - alpha*x);
end
