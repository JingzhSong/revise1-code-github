clear
clc
n=10;
uvwp = 'p';
nc = 201;

nkx = 100;
nkz = 100;
flowtype = 'turbu';
method = 'IOA';
eddy = 'eddyoff';
wall = 'rigid';
Y = [];

lambdax_min = 500;
lambdaz_min = 80;

figure;
nline=3;
%% read data
N = 122;
R = 180;
% path0 = ['K:/revise version code/convection/',wall,'/',uvwp];
path0 = [wall,'/',uvwp];
path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'lambda1.7_4'];

[y,~] = chebdif(N,2);
yplus = (y+1)*R;

lambdaxVector = logspace(1.7,4,nkx); 
lambdazVector = logspace(1.7,4,nkx);
kxVector = 2.*pi./lambdaxVector .*R;
kzVector = 2.*pi./lambdazVector .*R;

kx_matrix = zeros(nkx,nkz,N);
kz_matrix = zeros(nkx,nkz,N);
for ikx = 1:nkx
    kx_matrix(ikx,:,:) = kxVector(ikx)* ones(1,nkz,N);
end
for ikz = 1:nkz
    kz_matrix(:,ikz,:) = kzVector(ikz)* ones(nkx,1,N);
end
n1=length(find(lambdaxVector > lambdax_min));
fprintf('lambdax = %d \n',lambdaxVector(n1))
n2=length(find(lambdazVector > lambdaz_min));
fprintf('lambdaz = %d \n',lambdazVector(n2))

path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');

h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix; % log integeral need kx*kz
u_vec1 = squeeze(sum(cmax((nkx-n1+1):end,(nkx-n2+1):end,N/2:end) .* h((nkx-n1+1):end,(nkx-n2+1):end,N/2:end),[1 2])) ./ squeeze(sum(h((nkx-n1+1):end,(nkx-n2+1):end,N/2:end),[1 2])); %小尺度不算

semilogx(yplus(N/2+n:end,1),u_vec1(n+1:end),'--b','linewidth',nline); hold on

%% read data
N = 122;
R = 550;

path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'lambda1.7_4'];

[y,~] = chebdif(N,2);
yplus = (y+1)*R;

lambdaxVector = logspace(1.7,4,nkx);
lambdazVector = logspace(1.7,4,nkx);
kxVector = 2.*pi./lambdaxVector .*R;
kzVector = 2.*pi./lambdazVector .*R;

kx_matrix = zeros(nkx,nkz,N);
kz_matrix = zeros(nkx,nkz,N);
for ikx = 1:nkx
    kx_matrix(ikx,:,:) = kxVector(ikx)* ones(1,nkz,N);
end
for ikz = 1:nkz
    kz_matrix(:,ikz,:) = kzVector(ikz)* ones(nkx,1,N);
end
n1=length(find(lambdaxVector > lambdax_min));
fprintf('lambdax = %d \n',lambdaxVector(n1))
n2=length(find(lambdazVector > lambdaz_min));
fprintf('lambdaz = %d \n',lambdazVector(n2))

path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');
h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix; 
u_vec2 = squeeze(sum(cmax((nkx-n1+1):end,(nkx-n2+1):end,N/2:end) .* h((nkx-n1+1):end,(nkx-n2+1):end,N/2:end),[1 2])) ./ squeeze(sum(h((nkx-n1+1):end,(nkx-n2+1):end,N/2:end),[1 2])); %小尺度不算

semilogx(yplus(N/2+n:end,1),u_vec2(n+1:end),'--','Color',[0.196 0.804 0.196],'linewidth',nline); hold on

%% read data
N = 200;
R = 2000;

path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'lambda1.7_4'];

[y,~] = chebdif(N,2);
yplus = (y+1)*R;

lambdaxVector = logspace(1.7,4,nkx); 
lambdazVector = logspace(1.7,4,nkx);
kxVector = 2.*pi./lambdaxVector .*R;
kzVector = 2.*pi./lambdazVector .*R;

kx_matrix = zeros(nkx,nkz,N);
kz_matrix = zeros(nkx,nkz,N);
for ikx = 1:nkx
    kx_matrix(ikx,:,:) = kxVector(ikx)* ones(1,nkz,N);
end
for ikz = 1:nkz
    kz_matrix(:,ikz,:) = kzVector(ikz)* ones(nkx,1,N);
end
n1=length(find(lambdaxVector > lambdax_min));
fprintf('lambdax = %d \n',lambdaxVector(n1))
n2=length(find(lambdazVector > lambdaz_min));
fprintf('lambdaz = %d \n',lambdazVector(n2))

path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');
h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix;
u_vec3 = squeeze(sum(cmax((nkx-n1+1):end,(nkx-n2+1):end,N/2:end) .* h((nkx-n1+1):end,(nkx-n2+1):end,N/2:end),[1 2])) ./ squeeze(sum(h((nkx-n1+1):end,(nkx-n2+1):end,N/2:end),[1 2])); %小尺度不算

semilogx(yplus(N/2+n:end,1),u_vec3(n+1:end),'--r','linewidth',nline); hold on

%% u mean
[y,~] = chebdif(N,2);
k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
DUDy = @(y) R.*(-y)./NuT(y);
U0 = zeros(N,1);
for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end


%% plot
nsize=14;
grid on

load('DNS1993p.mat','DNS1993p')
load('DNS1993u.mat','DNS1993u')
load('DNS1993v.mat','DNS1993v')
load('DNS1993w.mat','DNS1993w')
if strcmp(uvwp,'u')
semilogx(DNS1993u(1:2:end,1),DNS1993u(1:2:end,2),'bo','linewidth',nline,'Markersize',nsize);hold on
end
if strcmp(uvwp,'v')
semilogx(DNS1993v(1:2:end,1),DNS1993v(1:2:end,2),'bo','linewidth',nline,'Markersize',nsize);hold on
end
if strcmp(uvwp,'w')
semilogx(DNS1993w(1:2:end,1),DNS1993w(1:2:end,2),'bo','linewidth',nline,'Markersize',nsize);hold on
end
if strcmp(uvwp,'p')
semilogx(DNS1993p(1:2:end,1),DNS1993p(1:2:end,2),'bo','linewidth',nline,'Markersize',nsize);hold on
end

% DNS
load('DNSp.mat','DNSp')
load('DNSu.mat','DNSu')
load('DNSv.mat','DNSv')
load('DNSw.mat','DNSw')
if strcmp(uvwp,'u')
semilogx(DNSu(1:2:end,1),DNSu(1:2:end,2),'r^','linewidth',nline,'Markersize',nsize);hold on
end
if strcmp(uvwp,'v')
semilogx(DNSv(1:2:end,1),DNSv(1:2:end,2),'r^','linewidth',nline,'Markersize',nsize);hold on
end
if strcmp(uvwp,'w')
semilogx(DNSw(1:2:end,1),DNSw(1:2:end,2),'r^','linewidth',nline,'Markersize',nsize);hold on
end
if strcmp(uvwp,'p')
semilogx(DNSp(1:2:end,1),DNSp(1:2:end,2),'r^','linewidth',nline,'Markersize',nsize);hold on
end

semilogx(yplus(N/2:end,1),U0(N/2:end,1),':k','linewidth',nline-1);hold on

fontn = 28;
axis([10^(0) 10^3.5 0 25])

xlabel('$y^+$','Interpreter','latex');

set(gcf,'unit','centimeters','position',[10 7 26 17]);
set(gca,'unit','centimeters','position',[3,3,22,13.5],'fontsize',fontn,'fontname','Times')
if strcmp(uvwp,'u')
%     ylabel('$\bar{c}^+_u,U^+$','Interpreter','latex');
ylabel('$\bar{c}_u,U$','Interpreter','latex');
end
if strcmp(uvwp,'v')
%     ylabel('$\bar{c}^+_v,U^+$','Interpreter','latex');
ylabel('$\bar{c}_v,U$','Interpreter','latex');
end
if strcmp(uvwp,'w')
% ylabel('$\bar{c}^+_w,U^+$','Interpreter','latex');
ylabel('$\bar{c}_w,U$','Interpreter','latex');
end
if strcmp(uvwp,'p')
%     ylabel('$\bar{c}^+_p,U^+$','Interpreter','latex');
ylabel('$\bar{c}_p,U$','Interpreter','latex');
end
% legend('$Re_\tau=180$','$Re_\tau=550$','$Re_\tau=2000$','DNS - $Re_\tau=180$','DNS - $Re_\tau=2000$','$U^+$','Interpreter','latex','location','SouthEast')
legend('$Re_\tau=180$','$Re_\tau=550$','$Re_\tau=2000$','DNS - $Re_\tau=180$','DNS - $Re_\tau=2000$','$U$','Interpreter','latex','location','SouthEast')


path = ['cov_Re',path1,num2str(lambdax_min)];

print('-dpng','-r300', [path,'.png'])
