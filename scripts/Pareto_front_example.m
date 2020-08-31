%%
clear all
close all
clc
%% import
addpath('e:\skola\semestralny projekt\1')
%% FIR DATA
% % input data - PSBR
% t = linspace(0,100,100);
% rng(0)
% y_input = round(rand(size(t)));
% %disp(y_input);
% % output data - FIR
% data_order = 8;
% P = ones(data_order,1);
% U = napln_maticu(y_input',data_order);
% y_plant = U*P;

%% TF DATA
t = linspace(0,20,100);
Gs = tf(8,[2 1]);

y_input = ones(size(t));
y_plant = step(Gs,t);

% noise
per_c = 0.1;
for k = 1:1:10
rng(k)
w_noise = max(y_plant)*per_c*(2*rand(size(y_plant))-1);
y_measured = w_noise + y_plant;
%% minimal solution
model_error = 3;
rad_p = 1;
C = zeros(rad_p);
v = 1;

while v==1
   C = napln_maticu(y_input',rad_p);
   s = size(C);
   n = zeros(1,s(2));
   v = isempty(linprog(n,[-C;C],[model_error-y_measured; model_error+y_measured]));
   if v == 1
       rad_p = rad_p + 1;
   else
       rad_p;   
   end
   if rad_p == 100
       v = 0;
   end
end
disp('Rad modelu : ');
disp(rad_p); 

P = inv(C'*C)*C'*y_measured;
[P_min,P_max] = get_coef(C,model_error,y_measured);
y_approx = C*P;
y_approx_min = C*P_min;
y_approx_max = C*P_max;
%% figures
% figure(1)
% 
% subplot(3,1,1)
% stairs(t,y_input,'k')
% axis([0 length(t) 0 2])
% grid on
% box on
% 
% subplot(3,1,2)
% hold on
% stairs(t,y_measured,'r')
% stairs(t,(y_measured-w_noise),'b')
% grid on
% box on
% 
% subplot(3,1,3)
% hold on
% stairs(t,y_approx,'b')
% stairs(t,y_approx_min,'g')
% stairs(t,y_approx_max,'y')
% grid on
% box on
% 
% close(1)
%% Pareto figure data
for r = rad_p:1:rad_p+15
    C = napln_maticu(y_input',r);
    P = inv(C'*C)*C'*y_measured;
    [P_min,P_max] = get_coef(C,model_error,y_measured);
    y_approx = C*P;
    y_approx_min = C*P_min;
    y_approx_max = C*P_max;
    
    rms1(1,r-rad_p+1) = sqrt(sum((y_measured-y_approx).^2)/length(y_measured));
    for i = 1:length(y_approx)
        rms_dis(1,i) = max([(y_approx(i)-y_approx_max(i).^2) (y_approx(i)-y_approx_min(i)).^2]);
    end
    max_val(1,r-rad_p+1) = max(rms_dis);
    rms2(1,r-rad_p+1) = sqrt(sum(rms_dis)/length(rms_dis));
    all_used_orders(1,r-rad_p+1) = r;
end
m_ord(k,:) = all_used_orders;
m_acc(k,:) = rms1;
m_oest(k,:) = rms2;
end
mean_ord = mean(m_ord);
mean_acc = mean(m_acc);
mean_oest = mean(m_oest);

disp([all_used_orders;rms1;rms2])
%disp(max_val)
%disp(sum_dis)
%disp(s_dev)
%%
figure(2)
hold on
plot(mean_acc./mean_acc(7),mean_oest./mean_oest(7),'bo','LineWidth',2,'MarkerFaceColor','b')
plot([0.4195 0.4195]./mean_acc(7),[10 120]./mean_oest(7),':k','LineWidth',1.5)
plot([0.35 1.25]./mean_acc(7),[20.1020 20.1020]./mean_oest(7),':k','LineWidth',1.5)
plot(0.4195./mean_acc(7),20.1020./mean_oest(7),'rsq','LineWidth',2,'MarkerFaceColor','r')

% for i = 1:1:length(mean_ord)
%     text(1.001*mean_acc(i),mean_oest(i),num2str(mean_ord(i)))
% end
xlabel('Vlastnosť A');
ylabel('Vlastnosť B');
axis([0.55 1.75 0.2 1.9])
box on
grid on


