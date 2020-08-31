clear all
close all
clc

addpath('funs')

% chemostat parameters
F = 1; %l/h
V = 3.33; %l
D = F/V; %h-1
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 22; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l
param = [mu_max nu Ks Ki Yx Yp s_in D];


s0 = 0; x0 = 10; p0 = 0;
init_c = [x0 s0 p0];

% figures
s = linspace(0,s_in,100);

figure(1)
mus_inh = (mu_max*s./(Ks + s + ((s.^2)/Ki)));
mus = (mu_max*s./(Ks + s));
opt_s_inh = sqrt(Ks*Ki);
opt_mu_inh = (mu_max*opt_s_inh./(Ks + opt_s_inh + ((opt_s_inh^2)/Ki)));

hold on
plot(s,mus,'k','LineWidth',2)
plot(s,mus_inh,'b','LineWidth',2)
plot([0 s_in],[D D],'.-.','Color',[1 102/255 0],'LineWidth',2)
plot(1.5688,D,'ks','MarkerFaceColor','r','MarkerSize',8)
plot(1.7510,D,'ks','MarkerFaceColor','r','MarkerSize',8)
plot(15.0768,D,'bo','MarkerFaceColor','b')
%legend('bez inhibície','s inhibíciou','rýchlost riedenia D','Location','Best')
xlabel('Koncetrácia substrátu s(t) [g/L]')
ylabel('Špecifická rýchlost rastu \mu(s) [h^{-1}]')
box on
%print('spec_grow_rate_comparison','-depsc')



s = [22 3 18 8 2];
x = [1.1 0.5 2 4.6 2.2];
clr = {[220/255 20/255 60/255] [1 69/255 0] [124/255 252/255 0] [0 1 1] [139/255 0 139/255]};
tf = 200; %hrs
Ts = 500;
t = linspace(0,tf,500);

figure(2)
for i = 1:length(s)
    s0 = s(i);
    x0 = x(i);
    init_c = [x0 s0 p0];
    
    [xt,st,~] = generate_Haldane_data(D,param,init_c,tf,Ts);
    
    subplot(2,1,1)
    hold on
    plot(t,st,'-','color',cell2mat(clr(i)),'LineWidth',2)
    xlabel('Čas t [h]')
    ylabel('s(t) [g/L]')
    set(gca,'FontSize',15)
    box on
    subplot(2,1,2)
    hold on
    plot(t,xt,'-','color',cell2mat(clr(i)),'LineWidth',2)
    xlabel('Čas t [h]')
    ylabel('x(t) [g/L]')
    set(gca,'FontSize',15)
    box on
end
hold off
%print('init_cond_inhb','-depsc')

figure(4)
hold on
for i = 1:length(s)
    s0 = s(i);
    x0 = x(i);
    init_c = [x0 s0 p0];
    
    [xt,st,~] = generate_Monod_data(D,param,init_c,tf,Ts);
    
    subplot(2,1,1)
    hold on
    plot(t,st,'-','color',cell2mat(clr(i)),'LineWidth',2)
    xlabel('Čas t [h]')
    ylabel('s(t) [g/L]')
    set(gca,'FontSize',15)
    axis([0 40 0 25])
    box on
    subplot(2,1,2)
    hold on
    plot(t,xt,'-','color',cell2mat(clr(i)),'LineWidth',2)
    xlabel('Čas t [h]')
    ylabel('x(t) [g/L]')
    set(gca,'FontSize',15)
    axis([0 40 0 6])
    box on
end
hold off
%print('init_cond_Monod','-depsc')


figure(3)
s = linspace(0.1,25,10);
x = linspace(0.1,5,20);
hold on
for i = 1:length(s)
    for j = 1:length(x)
        s0 = s(i);
        x0 = x(j);
        init_c = [x0 s0 p0];
    
        [xt,st,~] = generate_Haldane_data(D,param,init_c,tf,Ts);
        
        st = st./s_in;
        xt = xt./s_in;

        plot(xt,st,'-','color',[115/255 115/255 115/255]) 
    end
end

s = [22 3 18 8 2 ];
x = [1.1 0.5 2 4.6 2.2 ];
clr = {[220/255 20/255 60/255] [1 69/255 0] [124/255 252/255 0] [0 1 1] [139/255 0 139/255]};
for i = 1:1:length(x)
    s0 = s(i);
    x0 = x(i);
    init_c = [x0 s0 p0];
    
    [xt,st,~] = generate_Haldane_data(D,param,init_c,tf,Ts);
    
    st = st./s_in;
    xt = xt./s_in;

    plot(xt,st,'-','color',cell2mat(clr(i)),'LineWidth',2)
end
plot(0,1,'rs','MarkerFaceColor','r','MarkerSize',8)
plot(1.18/s_in,15.0768/s_in,'bo','MarkerFaceColor','b','MarkerSize',8)
plot(xt(end),st(end),'rs','MarkerFaceColor','r','MarkerSize',8)
ylabel('Relatívna koncentrácia substrátu')
xlabel('Relatívna koncentrácia biomasy')
axis([0 0.27 0 1.1])
axis square
set(gca,'FontSize',15)
hold off
box on
%print('phase_inhb','-depsc')

figure(5)
%D = .6;
s = linspace(0.01,25,10);
x = linspace(0.01,5,20);
hold on
for i = 1:length(s)
    for j = 1:length(x)
        s0 = s(i);
        x0 = x(j);
        init_c = [x0 s0 p0];
    
        [xt,st,~] = generate_Monod_data(D,param,init_c,tf,Ts);
        
        st = st./s_in;
        xt = xt./s_in;

        plot(xt,st,'-','color',[115/255 115/255 115/255]) 
    end
end

s = [22 3 18 8 2];
x = [1.1 0.5 2 4.6 2.2];
clr = {[220/255 20/255 60/255] [1 69/255 0] [124/255 252/255 0] [0 1 1] [139/255 0 139/255]};
for i = 1:1:length(x)
    s0 = s(i);
    x0 = x(i);
    init_c = [x0 s0 p0];
    
    [xt,st,~] = generate_Monod_data(D,param,init_c,tf,Ts);
    
    st = st./s_in;
    xt = xt./s_in;

    plot(xt,st,'-','color',cell2mat(clr(i)),'LineWidth',2)
end
plot(0,1,'rs','MarkerFaceColor','r','MarkerSize',8)
plot(xt(end),st(end),'rs','MarkerFaceColor','r','MarkerSize',8)
ylabel('Relatívna koncentrácia substrátu')
xlabel('Relatívna koncentrácia biomasy')
axis([0 0.3 0 1.1])
axis square
set(gca,'FontSize',15)
hold off
box on
%print('phase_Monod','-depsc')
