clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')
addpath('funs')
addpath('data')

F = 1; %l/h
V = 3.33; %l
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l
alpha = 0.5;

tf = 50;
Ts = 10; % number of measured samples through the simulation time
t = linspace(0,tf,Ts);

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = value(obj);

[x0,s0,p0] = get_Monod_ss(nom_d,param);
[ix0,is0,ip0] = get_Haldane_ss(nom_d,param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];



jmax = 3;
kmax = 100;

for k = 1:1:kmax
    D = [nom_d; nom_d];
    param = [mu_max nu Ks Ki Yx Yp s_in D(1,1) tf];
    
    for j = 1:1:jmax
        % hybrid model based on biomass concentration
        % Real data measurement
        [xt,~,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
        
        % noise addition
        rng(k)
        power = 0.05;
        noise_bio = power*(2*rand(size(xt)) - 1);
        
        m_xt = xt + noise_bio;
        
        % nominal model simulation
        [ixt,~,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
        
        % FIR model design
        dx = m_xt - ixt;
        u = ones([Ts 1])*D(1,:);
        u = u(:)';
        dx = dx';
        
        u = u - nom_d;
        dx0 = x0 - ix0;
        dx = dx - dx0;
        
        model_error = 2*power;
        
        [r_bio(j),P_bio] = gpe_fir_min_order(u,dx,model_error,1000,0,0);
        
        % hybrid model optimum
        % biomass
        a = sdpvar(1,1);
        obj = hybrid_obj_bio_corr(a,nom_d,dx0,alpha,param,P_bio);
        cons = [a >= 0.1;...
            a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        hyb_d_bio = value(a);
        
        P_us(1,j) = sum(P_bio)*(D(1,j)-nom_d) + dx0;
        
        % Pareto front data
        if j == jmax
            for r = r_bio(j):1:r_bio(j)+4
                C = get_matica_vstupov(u,r);
                P_bio_lsm = (C'*C)\C'*dx;
                [P_bio_min,P_bio_max] = gpe_fir_min_max_bound(r,u,dx,model_error,0);
                
                dx_lsm = C*P_bio_lsm;
                dx_min = C*P_bio_min;
                dx_max = C*P_bio_max;
                
                bio_acc = sqrt(sum((dx-dx_lsm).^2))/length(dx);
                bio_oest = sqrt(sum(max([(dx-dx_max).^2, (dx-dx_min).^2]')))/length(dx);
                
                % hybrid model optimum
                % biomass
                a = sdpvar(1,1);
                obj = hybrid_obj_bio_corr(a,nom_d,dx0,alpha,param,P_bio_lsm);
                cons = [a >= 0.1;...
                    a <= mu_max];
                assign(a,0.1);
                ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
                optimize(cons,obj,ops);
                hyb_d_bio = value(a);
                
                bio_pareto_data(r-r_bio(j)+1,:) = [r bio_acc bio_oest hyb_d_bio P_us(1,1) sum(P_bio_lsm)*hyb_d_bio];
            end
            bio_pareto_cell{k} = bio_pareto_data;
        end
        
        % hybrid model based on substrate
        % Real data measurement
        [~,st,~] = generate_Monod_data(D(2,:),param,init_c,tf,Ts);
        
        % noise addition
        rng(k+kmax)
        power = 0.025;
        noise_sub = power*(2*rand(size(st)) - 1);
        
        m_st = st + noise_sub;
        
        % nominal model simulation
        [~,ist,~] = generate_Haldane_data(D(2,:),param,i_init_c,tf,Ts);
        
        %FIR model design
        ds = m_st - ist;
        ds = ds';
        u = ones([Ts 1])*D(2,:);
        u = u(:)';
        
        u = u - nom_d;
        ds0 = s0 - is0;
        ds = ds - ds0;
        
        model_error = 2*power;
        
        [r_sub(j),P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);
        
        % hybrid model optimum
        % substrate
        a = sdpvar(1,1);
        obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_sub);
        cons = [a >= 0.1;...
            a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        hyb_d_sub = value(a);
        hyb_sub_j = obj;
        
        P_us(2,j) = sum(P_sub)*(hyb_d_sub-D(2,j)) + ds0;
        
        % Pareto front data
        if j == jmax
            for r = r_sub(j):1:r_sub(j)+4
                C = get_matica_vstupov(u,r);
                P_sub_lsm = (C'*C)\C'*ds;
                [P_sub_min,P_sub_max] = gpe_fir_min_max_bound(r,u,ds,model_error,0);
                
                ds_lsm = C*P_sub_lsm;
                ds_min = C*P_sub_min;
                ds_max = C*P_sub_max;
                
                sub_acc = sqrt(sum((ds-ds_lsm).^2))/length(ds);
                sub_oest = sqrt(sum(max([(ds-ds_max).^2, (ds-ds_min).^2]')))/length(ds);
                
                % hybrid model optimum
                % sub
                a = sdpvar(1,1);
                obj = hybrid_obj_sub(a,alpha,param,P_sub_lsm);
                cons = [a >= 0.1;...
                    a <= mu_max];
                assign(a,0.1);
                ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
                optimize(cons,obj,ops);
                hyb_d_sub = value(a);
                
                sub_pareto_data(r-r_sub(j)+1,:) = [r sub_acc sub_oest hyb_d_sub P_us(2,1) sum(P_sub_lsm)*hyb_d_sub];
            end
            sub_pareto_cell{k} = sub_pareto_data;
        end
        
        D(:,j+1) = [hyb_d_bio; hyb_d_sub];
        optimal_d(:,j) = [nom_d; hyb_d_bio; hyb_d_sub];
        P_us(:,j) = [sum(P_bio)*D(1,j); sum(P_sub)*D(2,j)];
    end
end

%% FIR identification (for figures)

D = [nom_d; nom_d];

for j = 1:1:3
    % hybrid model based on biomass concentration
        % Real data measurement
        [xt,~,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
        
        % noise addition
        rng(k)
        power = 0.05;
        noise_bio = power*(2*rand(size(xt)) - 1);
        
        m_xt = xt + noise_bio;
        
        % nominal model simulation
        [ixt,~,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
        
        % FIR model design
        dx = m_xt - ixt;
        u_bio = ones([Ts 1])*D(1,:);
        u_bio = u_bio(:)';
        dx = dx';
        
        u_bio = u_bio - nom_d;
        dx0 = x0 - ix0;
        dx = dx - dx0;
        
        model_error = 2*power;
        
        [r_bio(j),P_bio] = gpe_fir_min_order(u_bio,dx,model_error,1000,0,0);
        
        % hybrid model optimum
        % biomass
        a = sdpvar(1,1);
        obj = hybrid_obj_bio_corr(a,nom_d,dx0,alpha,param,P_bio);
        cons = [a >= 0.1;...
            a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        hyb_d_bio = value(a);
        
        P_us(1,j) = sum(P_bio)*(D(1,j)-nom_d) + dx0;
        
        % hybrid model based on substrate
        % Real data measurement
        [~,st,~] = generate_Monod_data(D(2,:),param,init_c,tf,Ts);
        
        % noise addition
        rng(k+kmax)
        power = 0.025;
        noise_sub = power*(2*rand(size(st)) - 1);
        
        m_st = st + noise_sub;
        
        % nominal model simulation
        [~,ist,~] = generate_Haldane_data(D(2,:),param,i_init_c,tf,Ts);
        
        %FIR model design
        ds = m_st - ist;
        ds = ds';
        u_sub = ones([Ts 1])*D(2,:);
        u_sub = u_sub(:)';
        
        u_sub = u_sub - nom_d;
        ds0 = s0 - is0;
        ds = ds - ds0;
        
        model_error = 2*power;
        
        [r_sub(j),P_sub] = gpe_fir_min_order(u_sub,ds,model_error,1000,0,0);
        
        % hybrid model optimum
        % substrate
        a = sdpvar(1,1);
        obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_sub);
        cons = [a >= 0.1;...
            a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        hyb_d_sub = value(a);
        hyb_sub_j = obj;
        
        P_us(2,j) = sum(P_sub)*(hyb_d_sub-D(2,j)) + ds0;
        
        D(:,j+1) = [hyb_d_bio; hyb_d_sub];
        optimal_d(:,j) = [nom_d; hyb_d_bio; hyb_d_sub];
        P_us(:,j) = [sum(P_bio)*D(1,j); sum(P_sub)*D(2,j)];
end

%% Data conditioning
bio_pareto_data = [];
sub_pareto_data = [];

for i = 1:1:kmax
    bio_pareto_data = [bio_pareto_data; bio_pareto_cell{i}];
    sub_pareto_data = [sub_pareto_data; sub_pareto_cell{i}];
end

avg_bio_matrix = average_if_matrix(bio_pareto_data);
avg_sub_matrix = average_if_matrix(sub_pareto_data);

rel_bio_matrix = bio_pareto_data(:,2:end)./avg_bio_matrix(3,2:end);
rel_bio_matrix = [bio_pareto_data(:,1) rel_bio_matrix];
[rel_bio_matrix,rel_bio_count] = average_if_matrix(rel_bio_matrix);

rel_sub_matrix = sub_pareto_data(:,2:end)./avg_sub_matrix(3,2:end);
rel_sub_matrix = [sub_pareto_data(:,1) rel_sub_matrix];
[rel_sub_matrix,rel_sub_count] = average_if_matrix(rel_sub_matrix);

bio_r = rel_bio_matrix(1:5,1);
bio_acc = rel_bio_matrix(1:5,2);
bio_oest = rel_bio_matrix(1:5,3);
sub_r = rel_sub_matrix(2:6,1);
sub_acc = rel_sub_matrix(2:6,2);
sub_oest = rel_sub_matrix(2:6,3);

% bio_acc = avg_bio_matrix(:,2);
% bio_oest = avg_bio_matrix(:,3);
% sub_acc = avg_sub_matrix(:,2);
% sub_oest = avg_sub_matrix(:,3);

%% Load data
% load('FIR_pareto_data.mat')
%% Figures
close all
clr = linspecer(11,'qualitative');

figure(1)
plot(bio_acc,bio_oest,'bo','MarkerFaceColor','b')
for i = 1:1:length(bio_acc)
    text(bio_acc(i)-0.002,bio_oest(i)-0.005,num2str(bio_r(i)))
end
xlabel('Relatívna presnosť odhadu')
ylabel('Maximálny relatívny rozptyl odhadu')
% axis([0.86 1.02 0.98 1.15])
grid on
box on

figure(2)
tt = linspace(0,3*tf,length(dx));
e = (2*0.05).*ones(size(dx));
C = get_matica_vstupov(u_bio,1);
P_bio = (C'*C)\C'*dx;
[P_bio_min, P_bio_max] = gpe_fir_min_max_bound(1,u_bio,dx,2*0.05,0);

fir_out = C*P_bio + dx0;
fir_min = C*P_bio_min + dx0;
fir_max = C*P_bio_max + dx0;

hold on
errorbar(tt,dx+dx0,-e,e,'o','Color','k','MarkerFaceColor','k','LineWidth',1.5)
stairs(tt,fir_min,'-','Color','r','LineWidth',2)
stairs(tt,fir_max,'-','Color','b','LineWidth',2)
stairs(tt,fir_out,'-','Color','g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Rozdiel koncentrácie biomasy [gL^{-1}]')
%legend('Namerané údaje','Výstup FIR','Location','Best')
set(gca,'FontSize',12)
box on
hold off

figure(3)
hold on
stairs(tt,m_xt,'-b','LineWidth',2)
stairs(tt,ixt,'-c','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia biomasy [gL^{-1}]')
legend('Zariadenie (Monod)','Nominálny model (Haldane)','Location','Best')
set(gca,'FontSize',12)
ylim([4.3 4.52])
box on
hold off


figure(4)
plot(sub_acc,sub_oest,'bo','MarkerFaceColor','b')
for i = 1:1:length(sub_acc)
    text(sub_acc(i)-0.002,sub_oest(i)-0.005,num2str(sub_r(i)))
end
xlabel('Relatívna presnosť odhadu')
ylabel('Maximálny relatívny rozptyl odhadu')
ylim([0.91 1.04])
xlim([0.9 1.175])
set(gca,'FontSize',12)
grid on
box on

figure(5)
e = (2*0.025).*ones(size(ds));
C = get_matica_vstupov(u_sub,2);
P_sub = (C'*C)\C'*ds;
[P_sub_min, P_sub_max] = gpe_fir_min_max_bound(2,u_sub,ds,2*0.025,0);

fir_out = C*P_sub + ds0;
fir_min = C*P_sub_min + ds0;
fir_max = C*P_sub_max + ds0;

hold on
errorbar(tt,ds+ds0,-e,e,'o','Color','k','MarkerFaceColor','k','LineWidth',1.5)
stairs(tt,fir_min,'-','Color','r','LineWidth',2)
stairs(tt,fir_max,'-','Color','b','LineWidth',2)
stairs(tt,fir_out,'-','Color','g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Rozdiel koncentrácie substrátu [gL^{-1}]')
%legend('Namerané údaje','Výstup FIR','Location','Best')
set(gca,'FontSize',12)
box on
hold off

figure(6)
hold on
stairs(tt,m_st,'-r','LineWidth',2)
stairs(tt,ist,'-m','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia substrátu [gL^{-1}]')
legend('Zariadenie (Monod)','Nominálny model (Haldane)','Location','Best')
set(gca,'FontSize',12)
box on
hold off
