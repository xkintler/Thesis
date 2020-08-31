clear all
close all
clc

addpath('funs')

F = 1; %l/h
V = 3.33; %l
D = F/V; %h-1
mu_max = 0.8; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 350; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

params = [mu_max nu Ks Ki Yx Yp s_in];

s0 = 0; x0 = 10; p0 = 0;
init_c = [x0; s0; p0];

%step change
D = [0.2 0.3 0.4 0.5];
k = length(D);

tf = 50;
Ts = 200;
time = linspace(0,k*tf,k*Ts);

[x,s,p] = generate_Monod_data(D,params,init_c,tf,Ts);

figure
hold on
plot(time,x,'-b','LineWidth',2)
plot(time,s,'-r','LineWidth',2)
plot(time,p,'-g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia x(t), s(t), p(t) [g/L]')
legend('Biomasa','Substrát','Produkt','Location','Best')
box on

%% Comparison Monod - Haldane
clear all 
close all

D = 0.4; %h-1
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 20; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

params = [mu_max nu Ks Ki Yx Yp s_in];

s0 = 0; x0 = 10; p0 = 0;
init_c = [x0; s0; p0];

tf = 50;
Ts = 200;

k = length(D);
time = linspace(0,k*tf,k*Ts);

[x,s,p] = generate_Monod_data(D,params,init_c,tf,Ts);
[ix,is,ip] = generate_Haldane_data(D,params,init_c,tf,Ts);

st = linspace(0,s_in,200);
m_mus = mu_max.*st./(Ks + st);
h_mus = mu_max.*st./(Ks + st + ((st.^2)./Ki));

figure
hold on 
plot(st,m_mus,'-k','LineWidth',2)
plot(st,h_mus,'-b','LineWidth',2)
plot([st(1) st(end)],[D D],'--y','LineWidth',2)
xlabel('Koncentrácia substrátu s(t) [g/L]')
ylabel('Špecifická rýchlosť rastu \mu(s) [h^{-1}]')
box on

figure
hold on
plot(time,x,'-b','LineWidth',2)
plot(time,s,'-r','LineWidth',2)
plot(time,p,'-g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia x(t), s(t), p(t) [g/L]')
legend('Biomasa','Substrát','Produkt','Location','Best')
set(gca,'FontSize',15)
box on

figure
hold on
plot(time,ix,'-b','LineWidth',2)
plot(time,is,'-r','LineWidth',2)
plot(time,ip,'-g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia x(t), s(t), p(t) [g/L]')
legend('Biomasa','Substrát','Produkt','Location','Best')
set(gca,'FontSize',15)
box on

