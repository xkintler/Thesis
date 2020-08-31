clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')
addpath('funs')

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
nom_j = obj;

D = nom_d;
[x0,s0,p0] = get_Monod_ss(nom_d,param);
[ix0,is0,ip0] = get_Haldane_ss(nom_d,param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

kmax = 3;
%D = [nom_d 0.3772 0.3826];
for k = 1:1:kmax
        
    % hybrid model based on substrate 
    % Real data measurement
    [~,st,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
        
    % noise addition
    rng(1)
    power = 0.025;
    noise_sub = power*(2*rand(size(st)) - 1);

    m_st = st + noise_sub;
    
    % nominal model simulation
    [~,ist,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
        
    %ARX model design
    ds = m_st - ist;
    ds = ds';
    u = ones([Ts 1])*D(1,:);
    u = u(:)';
    
    u = u - nom_d;
    ds0 = s0 - is0;
    ds = ds - ds0;
    
    model_error = 2*power;
    y0 = 0;
    
    gpe_out = gpe_arx_min_order(u,ds,model_error,y0,100,20,0);
    
    na = gpe_out.na;
    nb = gpe_out.nb;
    y_arx = gpe_out.ARX_output' + ds0;
    
    % ARX model design
    % arx identification toolbox
    io_data = iddata(ds,u');
    nk = 1;
    orders = [na nb nk];
    sys = arx(io_data,orders);
    opt_A = sys.A;
    opt_B = sys.B;

    P_us(k,1) = sum(opt_B)/sum(opt_A);
    
    % hybrid model optimum 
    a = sdpvar(1,1);
    obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_us(k));
    cons = [a >= 0.1;... 
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_sub = value(a);
    hyb_sub_j = obj;
    
    D(1,k+1) = hyb_d_sub;
        
end

e = model_error.*ones(size(ds));
[A_min,B_min,A_max,B_max] = gpe_arx_min_max_bound(na,nb,u,ds,y0,model_error,50);

arx_out = generate_ARX_output(opt_A(2:end),opt_B(2:end),u,y0) + ds0;
arx_min = generate_ARX_output(A_min,B_min,u,y0) + ds0;
arx_max = generate_ARX_output(A_max,B_max,u,y0) + ds0;

%% Figures
close all

tt = linspace(0,kmax*tf,kmax*Ts);

figure
hold on
stairs(tt,m_st,'-r','LineWidth',2)
stairs(tt,ist,'-m','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia substrátu [gL^{-1}]')
legend('Zariadenie (Monod)','Nominálny model (Haldane)','Location','Best')
%axis([tt(1) tt(end) 2.85 3.55])
set(gca,'FontSize',15)
box on
hold off

figure
hold on
errorbar(tt,ds+ds0,-e,e,'LineStyle','none','Marker','o','Color','k','MarkerFaceColor','k','LineWidth',2)
stairs(tt,arx_min,'-r','LineWidth',2)
stairs(tt,arx_max,'-b','LineWidth',2)
stairs(tt,arx_out,'-g','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Rozdiel koncentrácie substrátu [gL^{-1}]')
% legend('Namerané údaje','Výstupy ARX','Location','Best')
set(gca,'FontSize',15)
box on
hold off

%%
% %% Pareto front
% clear all
% close all
% clc
% 
% addpath('e:\skola\semestralny projekt\1')
% addpath('funs')
% 
% F = 1; %l/h
% V = 3.33; %l
% mu_max = 0.53; %h-1
% nu = 0.5; %h-1
% Ks = 1.2; %g/l
% Ki = 70; %g/l
% Yx = .4;
% Yp = 1;
% s_in = 20; %g/l
% alpha = 0.5;
% s0 = 2.9291; x0 = 4.4573; p0 = 5.9277;
% is0 = 3.3115; ix0 = 4.3574; ip0 = 5.7949;
% 
% ist = []; ixt = []; ipt = [];
% st = []; xt = []; pt = [];
% 
% tf = 50;
% Ts = 10; % number of measured samples through the simulation time
% t = linspace(0,tf,Ts);
% 
% param = [mu_max nu Ks Ki Yx Yp s_in [] tf];
% init_c = [x0;s0;p0];
% i_init_c = [ix0;is0;ip0];
% 
% jmax = 5;
% kmax = 2;
% 
% % Nominal (Haldane) model optimum
% a = sdpvar(1,1);
% obj = haldane_obj(a,alpha,param);
% cons = [a >= 0; a <= mu_max];
% assign(a,0.1);
% ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
% optimize(cons,obj,ops);
% nom_d = value(a);
% nom_j = obj;
% 
% D = nom_d;
% param = [mu_max nu Ks Ki Yx Yp s_in D(1,1) tf];
% 
% for j = 1:1:jmax
%     for k = 1:1:kmax
%         
%         % hybrid model based on substrate
%         % Real data measurement
%         [~,st,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
%         
%         % noise addition
%         rng(j)
%         power = 0.025;
%         noise_sub = power*(2*rand(size(st)) - 1);
%         
%         m_st = st + noise_sub;
%         
%         % nominal model simulation
%         [~,ist,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
%         
%         %FIR model design
%         ds = m_st - ist;
%         ds = ds';
%         u = ones([Ts 1])*D(1,:);
%         u = u(:)';
%         model_error = power;
%         y0 = -0.3824;
%         
%         gpe_out = gpe_arx_min_order(u,ds,model_error,y0,100,20,0);
%         
%         na = gpe_out.na;
%         nb = gpe_out.nb;
%         
%         % Pareto front data
%         if j == jmax
%             n = 1;
%             for rb = nb:1:nb+1
%                 for ra = na:1:na+2
%                     [A_min,B_min,A_max,B_max] = gpe_arx_min_max_bound(ra,rb,u,ds,y0,model_error,2);
%                     
%                     ds_min = generate_ARX_output(A_min,B_min,u,y0);
%                     ds_max = generate_ARX_output(A_max,B_max,u,y0);
%                     
%                     % ARX LSM
%                     io_data = iddata(ds,u');
%                     nk = 1;
%                     orders = [ra rb nk];
%                     sys = arx(io_data,orders);
%                     
%                     ds_lsm = generate_ARX_output(sys.A(2:end)',sys.B(2:end)',u,y0);
%                     
%                     sub_acc = sqrt(sum((ds-ds_lsm).^2))/length(ds);
%                     sub_oest = sqrt(sum(max([(ds-ds_max).^2, (ds-ds_min).^2]')))/length(ds);
%                     
%                     P_sub_lsm = sum(sys.B)/sum(sys.A);
%                     
%                     a = sdpvar(1,1);
%                     obj = hybrid_obj_sub(a,alpha,param,P_sub_lsm);
%                     cons = [a >= 0.1;...
%                         a <= mu_max];
%                     assign(a,0.1);
%                     ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
%                     optimize(cons,obj,ops);
%                     hyb_d_sub = value(a);
%                     
%                     sub_pareto_data(n,:) = [rb ra sub_acc sub_oest hyb_d_sub P_us(1,1)*D(1,2) sum(P_sub_lsm)*hyb_d_sub];
%                     n = n + 1;
%                 end
%             end
%             sub_pareto_cell{j} = sub_pareto_data;    
%         else
%             % ARX model design
%             % arx identification toolbox
%             io_data = iddata(ds,u');
%             nk = 1;
%             orders = [na nb nk];
%             sys = arx(io_data,orders);
%             opt_A = sys.A;
%             opt_B = sys.B;
%             
%             P_us(k,1) = sum(opt_B)/sum(opt_A);
%             
%             % hybrid model optimum
%             a = sdpvar(1,1);
%             obj = hybrid_obj_sub(a,alpha,param,P_us(k));
%             cons = [a >= 0.1;...
%                 a <= mu_max];
%             assign(a,0.1);
%             ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
%             optimize(cons,obj,ops);
%             hyb_d_sub = value(a);
%             hyb_sub_j = obj;
%             
%             D(1,k+1) = hyb_d_sub;
%         end        
%     end    
% end
% 
% % Data conditioning
% sub_pareto_data = [];
% for i = 1:1:kmax
%     sub_pareto_data = [sub_pareto_data; sub_pareto_cell{i}];
% end
% 
% 
