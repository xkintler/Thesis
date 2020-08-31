%% MAS w/o noise
clear all
close all
clc
addpath('funs')

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

tf = 50;
Ts = 10; % number of measured samples through the simulation time

alpha = .5;

param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];
t = linspace(0,tf,Ts);

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

kmax = 21;
D = [nom_d-0.01, nom_d];

c = [0.2 0.4 0.6 0.8 0.9] % bez sumu

clr = {'r','b','g','c','m','y','k'};
lam_matrix = zeros([length(c) kmax-1]);

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);

for k = 1:1:length(c)
    lam = 0;
    disp('------------------------------------------------------')
    for j = 2:1:kmax
        % Real data measurement
        [xt,~,~] = generate_Monod_data(D,param,[x0;s0;p0],tf,Ts);
        
        % noise addition
        rng(0)
        power = 0.03;
        noise_bio = power*(2*rand(size(xt)) - 1);
        
        m_xt = xt + noise_bio;

        % modifier adaptation scheme approach
        % cost function gradient approximation
        real_grad_app = grad_monod_obj(D(j),alpha,param);
        
        % evaluate modifier
        temp_lam = real_grad_app - grad_haldane_obj(D(j),alpha,param);
        lam = c(k)*lam + (1 - c(k))*temp_lam;
        
%         figure(3)
%         hold on
%         plot(j,real_grad_app,'o','Color',clr{k})
%         plot(j,grad_monod_obj(D(j),alpha,param),'ko')
            
        
        lam_matrix(k,j-1) = lam;
        
        % corrected nominal model optimum
        a = sdpvar(1,1);
        obj = modifier_adapt(a,alpha,param,lam);
        cons = [a >= 0.1; a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        mas_d = value(a);
        fprintf(['\nModifier adaptation scheme (Haldane) optimal D : ', num2str(mas_d), '\n'])
        
        D = [D mas_d];
        optimal_d(:,j-1) = [real_d; nom_d; mas_d];
    end
    cell_data{k} = optimal_d;
end

% Value of lambda (modifier) to push the value of nominal objective function to 
% real (plant/Monod) value of objective function
lam_limit = grad_monod_obj(real_d,alpha,param) - grad_haldane_obj(real_d,alpha,param)+ 0.07;
%% Figures
close all

figure(1)
hold on
iter = 1:1:kmax-1;
for i = 1:1:size(lam_matrix,1)
    stairs(iter,lam_matrix(i,:),clr{i},'LineWidth',2)
end
plot(iter, lam_limit*ones(size(iter)),'k:','LineWidth',1.5)
legend('c = 0.2','c = 0.4','c = 0.6','c = 0.8','c = 0.9','Location','Best')
xlabel('Iterácia')
ylabel('Hodnota modifikátora \lambda')
xlim([iter(1) iter(end)])
ylim([-2.1 0])
set(gca,'FontSize',12)
box on
hold off

figure(2)
hold on
for i = 1:1:length(c)
   stairs(iter,monod_obj(cell_data{i}(3,:),alpha,param),'Color',clr{i},'LineWidth',2)
end
plot([iter(1) iter(end)],[monod_obj(real_d,alpha,param) monod_obj(real_d,alpha,param)],'k:','LineWidth',1.5)
legend('c = 0.2','c = 0.4','c = 0.6','c = 0.8','c = 0.9','Location','Best')
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
xlim([iter(1) iter(end)])
ylim([-0.4705 -0.4625])
set(gca,'FontSize',12)
box 
hold off

%% MAS with noisy data
clear all
close all
clc

addpath('funs')

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

tf = 50;
Ts = 10; % number of measured samples through the simulation time

alpha = .5;

param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];
t = linspace(0,tf,Ts);

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

kmax = 21;
D = [0.1, nom_d];

c = 0.8; % so sumom
clr = {'r','b','g','c','y','m','k'};
lam_matrix = zeros([length(c) kmax-1]);

[x0,s0,p0] = get_Monod_ss(nom_d,param);
[ix0,is0,ip0] = get_Haldane_ss(nom_d,param);

for k = 1:1:length(c)
    lam = 0;
    grad_store = [];
    disp('------------------------------------------------------')
    for j = 2:1:kmax
        % Real data measurement
        [xt,~,~] = generate_Monod_data(D,param,[x0;s0;p0],tf,Ts);
        
        % noise addition
        rng(0)
        power = 0.03;
        noise_bio = power*(2*rand(size(xt)) - 1);
        
        m_xt = xt + noise_bio;

        % modifier adaptation scheme approach
        % cost function gradient approximation
        J0 = D(j).*(1 - alpha*m_xt(end));
        J1 = D(j-1).*(1 - alpha*m_xt(end-Ts));
        real_grad_app = (J0 - J1)/(D(j) - D(j-1));
        
        if j > 2 && abs(real_grad_app-grad_store(j-2)) > 8
            real_grad_app = 1.2*grad_store(j-2);
        end
        

        % evaluate modifier
        temp_lam = real_grad_app - grad_haldane_obj(D(j),alpha,param);
        lam = c(k)*lam + (1 - c(k))*temp_lam;
        
        figure(3)
        hold on
        plot(j,real_grad_app,'o','Color',clr{k})
        plot(j,grad_monod_obj(D(j),alpha,param),'ko')
            
        
        lam_matrix(k,j-1) = lam;
        
        % corrected nominal model optimum
        a = sdpvar(1,1);
        obj = modifier_adapt(a,alpha,param,lam);
        cons = [a >= 0.1; a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        mas_d = value(a);
        fprintf(['\nModifier adaptation scheme (Haldane) optimal D : ', num2str(mas_d), '\n'])
        
        D = [D mas_d];
        optimal_d(:,j-1) = [real_d; nom_d; mas_d];
        grad_store = [grad_store real_grad_app];
    end
end

% Value of lambda (modifier) to push the value of nominal objective function to 
% real (plant/Monod) value of objective function
lam_limit = grad_monod_obj(real_d,alpha,param) - grad_haldane_obj(real_d,alpha,param)+ 0.07;

%% Figures
close all

figure(1)
hold on
iter = 1:1:kmax-1;
plot(iter, lam_limit*ones(size(iter)),'k:')
for i = 1:1:size(lam_matrix,1)
    stairs(iter,lam_matrix(i,:),clr{i})
end
xlabel('Iteration')
ylabel('Value of \lambda')
xlim([iter(1) iter(end)])
box on
hold off

% figure(2)
% hold on
% plot(xt,':k')
% plot(m_xt,'-r')
% hold off
%% Different value of noise dispersion
clear all
close all
clc
addpath('funs')

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

tf = 150;
Ts = 10; % number of measured samples through the simulation time

alpha = .5;

param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];
t = linspace(0,tf,Ts);

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

kmax = 21;
D = [0.01, nom_d];
disper = [0 0.01 0.03];

c = 0.945; 

clr = {'r','b','g','c','m','y','k'};
lam_matrix = zeros([length(disper) kmax-1]);

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);

for k = 1:1:length(disper)
    lam = 0;
    disp('------------------------------------------------------')
    for j = 2:1:kmax
        % Real data measurement
        [xt,~,~] = generate_Monod_data(D,param,[x0;s0;p0],tf,Ts);
        
        % noise addition
        rng(0)
        power = disper(k);
        noise_bio = power*(2*rand(size(xt)) - 1);
        
        m_xt = xt + noise_bio;
        
        a_xt0 = mean(m_xt(1,end-3:end));
        a_xt1 = mean(m_xt(1,end-3-Ts:end-Ts));
        % modifier adaptation scheme approach
        % cost function gradient approximation
        J0 = D(j).*(1 - alpha*a_xt0);
        J1 = D(j-1).*(1 - alpha*a_xt1);
        real_grad_app = (J0 - J1)/(D(j) - D(j-1));
        
        % evaluate modifier
        temp_lam = real_grad_app - grad_haldane_obj(D(j),alpha,param);
        lam = c*lam + (1 - c)*temp_lam;
       
        lam_matrix(k,j-1) = lam;
        
        % corrected nominal model optimum
        a = sdpvar(1,1);
        obj = modifier_adapt(a,alpha,param,lam);
        cons = [a >= 0.1; a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        mas_d = value(a);
        fprintf(['\nModifier adaptation scheme (Haldane) optimal D : ', num2str(mas_d), '\n'])
        
        D = [D mas_d];
        optimal_d(:,j-1) = [real_d; nom_d; mas_d];
    end
    cell_data{k} = optimal_d;
end

% Value of lambda (modifier) to push the value of nominal objective function to 
% real (plant/Monod) value of objective function
lam_limit = grad_monod_obj(real_d,alpha,param) - grad_haldane_obj(real_d,alpha,param)+ 0.07;
%% Figures
close all

figure(1)
hold on
iter = 1:1:kmax-1;
for i = 1:1:size(lam_matrix,1)
    stairs(iter,lam_matrix(i,:),clr{i},'LineWidth',2)
end
plot(iter, lam_limit*ones(size(iter)),'k:','LineWidth',1.5)
legend('|e| = 0.00','|e| = 0.01','|e| = 0.03','Location','Best')
xlabel('Iterácia')
ylabel('Hodnota modifikátora \lambda')
xlim([iter(1) iter(end)])
%ylim([-2.1 0])
set(gca,'FontSize',12)
box on
hold off

figure(2)
hold on
for i = 1:1:length(disper)
   stairs(iter,monod_obj(cell_data{i}(3,:),alpha,param),'Color',clr{i},'LineWidth',2)
end
plot([iter(1) iter(end)],[monod_obj(real_d,alpha,param) monod_obj(real_d,alpha,param)],'k:','LineWidth',1.5)
legend('|e| = 0.00','|e| = 0.01','|e| = 0.03','Location','Best')
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
xlim([iter(1) iter(end)])
ylim([-0.471 -0.461])
set(gca,'FontSize',12)
box 
hold off
%% Cost Function Correction
clear all
close all
clc

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

tf = 50;
Ts = 10; % number of measured samples through the simulation time

alpha = .5;

param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];
t = linspace(0,tf,Ts);

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = value(obj);

lam = grad_monod_obj(real_d,alpha,param) - grad_haldane_obj(real_d,alpha,param);

% corrected nominal model optimum
a = sdpvar(1,1);
obj = modifier_adapt(a,alpha,param,lam);
cons = [a >= 0.1; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
mas_d = value(a);
mas_j = value(obj);

dd = 0:0.001:0.42;
for i = 1:1:length(dd)
   fm(i) = monod_obj(dd(i),alpha,param);
   fh(i) = haldane_obj(dd(i),alpha,param);
   fmas(i) = modifier_adapt(dd(i),alpha,param,lam);
end
figure(1)
hold on
plot(dd,fm,'-b','LineWidth',2)
plot(dd,fh,'-r','LineWidth',2)
plot(dd,fmas,'-g','LineWidth',2)
plot(real_d,real_j,'ok','MarkerFaceColor','k')
plot(nom_d,nom_j,'ok','MarkerFaceColor','k')
plot(mas_d,mas_j,'ok','MarkerFaceColor','k')
plot([mas_d mas_d],[-10 10],'--k','LineWidth',1.5)
xlabel('Rýchlosť riedenia [hod^{-1}]')
ylabel('Hodnota účelovej funkcie')
legend('Zariadenie (Monod)','Nominálny model (Haldane)','Upravená','Location','Best')
axis([dd(1) dd(end) -1.4 0.2])
set(gca,'FontSize',12)
box on
hold off