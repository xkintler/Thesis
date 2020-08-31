clear all
close all
clc

nu = [20 1];
y = 5*ones(nu);
u = linspace(1,nu(1),nu(1))'; 
rng(0)
e = 2*rand(nu) - 1;

% e = ones(nu);
ve = ones(nu);
theta1 = 0:0.01:10;

figure
hold on
for i = 1:1:length(y)
    theta2_1 = (1/u(i))*(y(i)-e(i)) - (1/u(i)).*theta1;
    theta2_2 = (1/u(i))*(y(i)+e(i)) - (1/u(i)).*theta1;
    
    plot(theta1,theta2_2,'LineWidth',2,'Color',[200/255 200/255 200/255])
    plot(theta1,theta2_1,'LineWidth',2,'Color',[200/255 200/255 200/255])
    %disp(i)
    %pause
end
%data bez sumu
% axis([3.5 6.5 -0.15 0.15])
% 
% plot(3.8962,0.1050,'ro','MarkerFaceColor','r','LineWidth',1.5)
% plot([3.5 3.8962],[0.1050 0.1050],':r','LineWidth',1.5)
% plot([3.8962 3.8962],[-0.15 0.1050],':r','LineWidth',1.5)
% 
% plot(6.1079,-0.1062,'ro','MarkerFaceColor','r','LineWidth',1.5)
% plot([3.5 6.1079],[-0.1062 -0.1062],':r','LineWidth',1.5)
% plot([6.1079 6.1079],[-0.15 -0.1062],':r','LineWidth',1.5)

%v = ginput(4);
%v = [4.0030 0.0008; 3.8982 0.0997; 5.9990 -0.0004; 6.1018 -0.0994];
%f = [1 2 3 4];
%patch('Faces',f,'Vertices',v,'FaceColor',[200/255 200/255 200/255],'FaceAlpha',0.8,'EdgeColor','None')

% data so sumom
axis([4.73 5.3 -0.025 0.025])
% % v = ginput(4);
v = [4.8210    0.0108; 4.7463    0.0201; 5.1791   -0.0107; 5.2538   -0.0201];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[200/255 200/255 200/255],'FaceAlpha',0.8,'EdgeColor','None')

plot(4.7463,0.0201,'ro','MarkerFaceColor','r','LineWidth',1.5)
plot([4.73 4.7463],[0.0201 0.0201],':r','LineWidth',1.5)
plot([4.7463 4.7463],[-0.025 0.0201],':r','LineWidth',1.5)

plot(5.2538,-0.0201,'ro','MarkerFaceColor','r','LineWidth',1.5)
plot([4.73 5.2538],[-0.0201 -0.0201],':r','LineWidth',1.5)
plot([5.2538 5.2538],[-0.025 -0.0201],':r','LineWidth',1.5)

xlabel('\theta_1')
ylabel('\theta_2')
%legend('\theta_1^{(1)}','\theta_1^{(2)}')
set(gca,'FontSize',12)
box on

%print('gpe_ex_line1', '-depsc')


figure
hold on
errorbar(u,y+e,-ve,ve,'Marker','sq','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2,'Color','k','LineStyle','None')
xlabel('Vstupy {\it u}')
ylabel('Výstupy {\it y(u)}')
axis([u(1) 20 2 8])
set(gca,'FontSize',12)
box on

%print('gpe_ex_data1', '-depsc')

