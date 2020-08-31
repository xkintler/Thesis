% cost function comparison
clear all
close all
clc

addpath('C:\Users\Matej\Documents\Diplomovka\scripts\funs')

F = 1; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

params = [mu_max nu Ks Ki Yx Yp s_in];
alpha = 1;
D = linspace(0,0.42,200);

% monod data
fm = monod_obj(D,alpha,params);

a = sdpvar(1,1);
obj = monod_obj(a,alpha,params);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
m_d = value(a);
m_j = value(obj);

% haldane data
fh = haldane_obj(D,alpha,params);

a = sdpvar(1,1);
obj = haldane_obj(a,alpha,params);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
h_d = value(a);
h_j = value(obj);

figure
hold on
plot(D,fm,'-b','LineWidth',2)
plot(D,fh,'-r','LineWidth',2)
plot(m_d,m_j,'ok','LineWidth',1.5)
plot(h_d,h_j,'ok','LineWidth',1.5)
xlabel('Rýchlosť riedenia D [h^{-1}]')
ylabel('Hodnota účelovej funkcie J')
legend('Monod','Haldane','Optimálne D*','Location','Best')
xlim([D(1) D(end)])
box on

%%
clear  
close all
clc

F = 1.5; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

tf = 80;
Ts = 100;
time = linspace(0,tf,Ts);

params = [mu_max nu Ks Ki Yx Yp s_in];

[x,s,p] = generate_Monod_data(D,params,[10;0;0],tf,Ts);
[ix,is,ip] = generate_Haldane_data(D,params,[10;0;0],tf,Ts);

figure
hold on
plot(time,x,'-b','LineWidth',2)
plot(time,s,'-r','LineWidth',2)
plot(time,p,'-g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia [gL^{-1}]')
set(gca,'FontSize',12)
box on

figure
hold on
plot(time,ix,'-b','LineWidth',2)
plot(time,is,'-r','LineWidth',2)
plot(time,ip,'-g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia [gL^{-1}]')
set(gca,'FontSize',12)
box on
%% porovnanie ucelovych funkcii monod a haldane modelu
clear  
close all
clc

addpath('funs')

F = 1.5; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

tf = 80;
Ts = 100;
time = linspace(0,tf,Ts);
alpha = 0.5;

params = [mu_max nu Ks Ki Yx Yp s_in];

a = sdpvar(1,1);
obj = monod_obj(a,alpha,params);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

a = sdpvar(1,1);
obj = haldane_obj(a,alpha,params);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = value(obj);

d = linspace(0,0.42,1000);

fm = monod_obj(d,alpha,params);
fh = haldane_obj(d,alpha,params);

figure
hold on
plot(d,fm,'-b','LineWidth',2)
plot(d,fh,'-r','LineWidth',2)
plot(real_d,real_j,'ok','MarkerFaceColor','k','LineWidth',2)
plot(nom_d,nom_j,'ok','MarkerFaceColor','k','LineWidth',2)
xlabel('Rýchlosť riedenia D [h^{-1}]')
ylabel('Hodnota účelovej funkcie')
axis([0 0.42 -0.5 0.1])
box on
set(gca,'FontSize',12)

%% specificka rychlost rastu MO
clear  
close all
clc

F = 1.5; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

s = linspace(0,s_in,200);

mon = mu_max.*s./(Ks + s);
hal = mu_max.*s./(Ks + s + (s.^2)./Ki);

figure
hold on
plot(s,mon,'-k','LineWidth',2)
plot(s,hal,'--k','LineWidth',2)
xlabel('Koncentrácia substrátu [gL^{-1}]')
ylabel('Špecifická rýchlosť rastu [hod^{-1}]')
legend('Monod','Haldane')
set(gca,'FontSize',12)
box on

%%
clear
close all
clc

addpath('funs')

F = 1.5; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

tf = 50;
Ts = 10;
time = linspace(0,tf,Ts);
alpha = 0.5;

params = [mu_max nu Ks Ki Yx Yp s_in];

[x0,s0,p0] = get_Monod_ss(0.2,params);
init_c = [x0;s0;p0];

D = [0.2 0.35];

[xt,st,~] = generate_Monod_data(D,params,init_c,tf,Ts);
        
% noise addition
rng(0)
power = 0.05;
noise_bio = power*(2*rand(size(xt)) - 1);

m_xt = xt + noise_bio;

rng(1)
power = 0.025;
noise_sub = power*(2*rand(size(st)) - 1);

m_st = st + noise_sub;

tt = linspace(0,2*tf,2*Ts);

figure
hold on
stairs(tt,m_xt,'-b','LineWidth',2)
stairs(tt,m_st,'-r','LineWidth',2)
xlabel('Čas [hod^{-1}]')
ylabel('Koncentrácia [gL^{-1}]')
set(gca,'FontSize',12)
box on

