%%  FIGURE 5 a
clear
clc
kz = 0;
N = 200; % wall-normal grids

flowtype = 'turbu'; 
method = 'SIOA'; % 'IOA', 'SIOA'
eddy = 'eddyon'; % 'eddyoff', 'eddyon'

NY = 50;
% kxvector = logspace(0,1,NY);
% cvector = linspace(4,14,NY);
% kxvector = logspace(0,2,NY);
% cvector = linspace(1,21,NY);
kxvector = logspace(0,2,NY);
cvector = linspace(1,20,NY);
% cvector = linspace(1,19,NY);
NormCom = zeros(NY,NY);
    
% Uc=19; %Xu
Uc=21; % Kim

% Kim
Cm = 2; 
Ck = 1*Uc*Uc; 
Cd = 0.5*Uc;
Cb = 0;
Ct = 0;

wall = 'com';
path1 = [method,'_',eddy,'_',wall,'/Norm_DNS_Kim0_kz',num2str(kz),'.mat'];
load(path1)
Norm_com = NormCom;
wall = 'rigid';
path2 = [method,'_',eddy,'_',wall,'/Norm_DNS_Kim_kz',num2str(kz),'.mat'];
load(path2)
Norm_rigid = NormRigid;

% plot
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\bar{\sigma}_c}/{\bar{\sigma}_0})$';
end
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\bar{\sigma}_c^e}/{\bar{\sigma}_0^e})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\mu_c}/{\mu_0})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\mu_c^e}/{\mu_0^e})$';
end

cormin=0;
cormax=1;
A = Norm_com./Norm_rigid;
figure;
contourf(kxvector,cvector,log10(A),[cormin:0.05:cormax],'EdgeColor','none');
load('MyColormaps.mat')
colormap(mycmap);

caxis([cormin,cormax]);
grid on
colorbar;

set(gca,'XScale','log');
xlabel('$k_{x}$','Interpreter','latex');
ylabel('$c$','Interpreter','latex');

axis([1,100,1,20]);
set(gca,'XTick',[1,10,100]);
set(gca,'YTick',[1,5,10,15,20]);

h=colorbar;
fontn =26;
title(str,'Interpreter','latex');
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
set(gcf,'unit','centimeters','position',[10 7 14 13]);
set(gca,'unit','centimeters','position',[3,3.5,8,8],'fontsize',fontn,'fontname','Times')
hold on

Ck_equ = Ck + Cb*(kxvector.^4+kz^4+2*kxvector.^2*kz^2) + Ct*(kxvector.^2+kz^2);
omegar = sqrt(Ck_equ./Cm .* (1 - 2*(Cd/2./sqrt(Ck_equ*Cm)).^2) );

x = kxvector;
y = omegar./x;
plot(kxvector,y,'k','LineWidth',2);
hold on

plot(8/3,4.78,'ok','MarkerSize',10,'LineWidth',3);

wall = 'com';
path = [method,'_',eddy,'_',wall,'/Norm_DNS_Kim0_kz',num2str(kz),'_ratio'];
saveas(gcf, [path1,'.png']);



%% FIGURE 5 b
clear
clc
kz = 0;
N = 200; % wall-normal grids

flowtype = 'turbu'; 
method = 'SIOA'; % 'IOA', 'SIOA'
eddy = 'eddyon'; % 'eddyoff', 'eddyon'

NY = 50;
kxvector = logspace(0,2,NY);
cvector = linspace(1,20,NY);
NormCom = zeros(NY,NY);
    
Uc=19; % Xia
% Uc=21; % Kim

% Xia
Cm = 1; 
% Ck = 0.2*Uc*Uc;  %case D 
Ck = 1.68*Uc*Uc; %case C B
Cd = 0.25*Uc; %case C D 
% Cd = 0.47*Uc; %case B
Cb = 1.2*10^(-4)*Uc*Uc;
Ct = 0.0033*Uc*Uc;

wall = 'com';
path1 = [method,'_',eddy,'_',wall,'/Norm_DNS_XiaC_kz',num2str(kz),'.mat'];
load(path1)
Norm_com = NormCom;
wall = 'rigid';
path2 = [method,'_',eddy,'_',wall,'/Norm_DNS_Xia_kz',num2str(kz),'.mat'];
load(path2)
Norm_rigid = NormRigid;

% plot
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\bar{\sigma}_c}/{\bar{\sigma}_0})$';
end
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\bar{\sigma}_c^e}/{\bar{\sigma}_0^e})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\mu_c}/{\mu_0})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\mu_c^e}/{\mu_0^e})$';
end

cormin=0;
cormax=1;
A = Norm_com./Norm_rigid;
figure;
contourf(kxvector,cvector,log10(A),[cormin:0.05:cormax],'EdgeColor','none');
load('MyColormaps.mat')
colormap(mycmap);

caxis([cormin,cormax]);
grid on
colorbar;

set(gca,'XScale','log');
xlabel('$k_{x}$','Interpreter','latex');
ylabel('$c$','Interpreter','latex');

axis([1,100,1,20]);
set(gca,'XTick',[1,10,100]);
set(gca,'YTick',[1,5,10,15,20]);

h=colorbar;
fontn =26;
title(str,'Interpreter','latex');
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
set(gcf,'unit','centimeters','position',[10 7 14 13]);
set(gca,'unit','centimeters','position',[3,3.5,8,8],'fontsize',fontn,'fontname','Times')
hold on

Ck_equ = Ck + Cb*(kxvector.^4+kz^4+2*kxvector.^2*kz^2) + Ct*(kxvector.^2+kz^2);
omegar = sqrt(Ck_equ./Cm .* (1 - 2*(Cd/2./sqrt(Ck_equ*Cm)).^2) );

x = kxvector;
y = omegar./x;
plot(kxvector,y,'k','LineWidth',2);
hold on

plot(4.68,5.96,'+k','MarkerSize',15,'LineWidth',3); % case C

wall = 'com';
path = [method,'_',eddy,'_',wall,'/Norm_DNS_XiaC_kz',num2str(kz),'_ratio'];
saveas(gcf, [path1,'.png']);



%% FIGURE 5 c
clear
clc
kz = 0;
N = 200; % wall-normal grids

flowtype = 'turbu'; 
method = 'SIOA'; % 'IOA', 'SIOA'
eddy = 'eddyon'; % 'eddyoff', 'eddyon'

NY = 50;
kxvector = logspace(0,2,NY);
cvector = linspace(1,20,NY);
NormCom = zeros(NY,NY);
    
Uc=19; % Xia
% Uc=21; % Kim

% Xia
Cm = 1; 
% Ck = 0.2*Uc*Uc;  %case D 
Ck = 1.68*Uc*Uc; %case C B
Cd = 0.25*Uc; %case C D 
% Cd = 0.47*Uc; %case B
% Cb = 1.2*10^(-4)*Uc*Uc;
% Ct = 0.0033*Uc*Uc;
Cb = 0;
Ct = 0;

wall = 'com';
path1 = [method,'_',eddy,'_',wall,'/Norm_DNS_Xia0_kz',num2str(kz),'.mat'];
load(path1)
Norm_com = NormCom;
wall = 'rigid';
path2 = [method,'_',eddy,'_',wall,'/Norm_DNS_Xia_kz',num2str(kz),'.mat'];
load(path2)
Norm_rigid = NormRigid;

% plot
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\bar{\sigma}_c}/{\bar{\sigma}_0})$';
end
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\bar{\sigma}_c^e}/{\bar{\sigma}_0^e})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\mu_c}/{\mu_0})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\mu_c^e}/{\mu_0^e})$';
end

cormin=0;
cormax=1;
A = Norm_com./Norm_rigid;
figure;
contourf(kxvector,cvector,log10(A),[cormin:0.05:cormax],'EdgeColor','none');
load('MyColormaps.mat')
colormap(mycmap);

caxis([cormin,cormax]);
grid on
colorbar;

set(gca,'XScale','log');
xlabel('$k_{x}$','Interpreter','latex');
ylabel('$c$','Interpreter','latex');

axis([1,100,1,20]);
set(gca,'XTick',[1,10,100]);
set(gca,'YTick',[1,5,10,15,20]);

h=colorbar;
fontn =26;
title(str,'Interpreter','latex');
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
set(gcf,'unit','centimeters','position',[10 7 14 13]);
set(gca,'unit','centimeters','position',[3,3.5,8,8],'fontsize',fontn,'fontname','Times')
hold on

Ck_equ = Ck + Cb*(kxvector.^4+kz^4+2*kxvector.^2*kz^2) + Ct*(kxvector.^2+kz^2);
omegar = sqrt(Ck_equ./Cm .* (1 - 2*(Cd/2./sqrt(Ck_equ*Cm)).^2) );

x = kxvector;
y = omegar./x;
plot(kxvector,y,'k','LineWidth',2);
hold on

plot(4.68,5.96,'+k','MarkerSize',15,'LineWidth',3); % case C

wall = 'com';
path = [method,'_',eddy,'_',wall,'/Norm_DNS_Xia0_kz',num2str(kz),'_ratio'];
saveas(gcf, [path1,'.png']);