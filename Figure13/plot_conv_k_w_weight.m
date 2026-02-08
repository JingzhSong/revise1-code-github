clear
clc

% N = 400;
% R = 6700;
% N = 300;
% R = 3300;
N = 200;
R = 2000;
nc = 201;

nkx = 100;
nkz = 100;
method = 'IOA';
eddy = 'eddyoff';
wall = 'rigid';
Y = [];
uvwp = 'p';

path0 = [wall,'/',uvwp];

if strcmp(wall, 'rigid')
    %     path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];
    path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'lambda1.7_4'];
else
    Cm = 0.46; Ck = 181; Cd=0.091;
    %     path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'kplus-70_final'];
    path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'lambda1.7_4'];
end

path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');

lambdaxVector = logspace(1.7,4,nkx);
lambdazVector = logspace(1.7,4,nkx);
% kxVector = linspace(10,160,nkx);
% kzVector = linspace(10,160,nkz);

[y,~] = chebdif(N,2);
yplus = (y+1)*R;

if strcmp(wall, 'rigid')
    I = eye(N);
    k = 0.426;
    alpha = 25.4;
    NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
    nuT0 = NuT(y);
    DUDy = @(y) R.*(-y)./NuT(y);
    U1 = DUDy(y);
    U0 = zeros(N,1);
    for j=1:N
        U0(j) = integral(DUDy,-1,y(j));
    end
end
if strcmp(wall, 'com')
    if R==3300
        load('COMUmean3300.mat','Umean');
        U0 = Umean;
    end
    if R==6700
        load('COMUmean6700.mat','Umean');
        U0 = Umean;
    end
    if R==8900
        load('COMUmean8900.mat','Umean');
        U0 = Umean;
    end
end

%%
figure;
set(gcf,'unit','centimeters','position',[10 7 40 13]);

cormin = 0;
cormax = 2;

fontn = 24;

subplot(1,3,1)
% yi = length(find(yplus > 5));
yi = length(find(yplus > 0))+1;
fprintf('y=%d \n',yplus(yi));
c = cmax(:,:,yi);

% Divided by the model speed
if strcmp(wall, 'rigid')
    pcolor(lambdaxVector,lambdazVector,c'./12); hold on
else
    if R == 3300
        uc=26;
    else
        uc= 28;
    end
    pcolor(lambdaxVector,lambdazVector,c'./0.53./uc); hold on
end

set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+ $','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');

if R == 8900
    title('$ y^+ \approx 4 $','Interpreter','latex');
else
    title('$ y^+ \approx 7 $','Interpreter','latex');
end

title('$ y^+ \approx 0 $','Interpreter','latex');

shading interp
load('MyColormaps.mat')
colormap(mymap);
% colormap(redblue1(cormin,cormax));
% caxis([cormin,cormax]);


set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
% axis([50 5*10^5 50 5*10^5])
axis([60 10^4 60 10^4])
caxis([cormin,cormax]);

if strcmp(wall, 'rigid')
    % Draw large scale boundaries
    % hold on
    % xline(2*R,'--k','linewidth',2)
    % hold on
    % yline(0.4*R,'--k','linewidth',2)
    hold on
    plot([60,10^4],[60,10^4],'-.k','linewidth',2)
    
    hold on
    plot([60,250],[250,250],'--k','linewidth',2)
    hold on
    plot([250,250],[60,250],'--k','linewidth',2)
    
    hold on
    plot([R,10^4],[R,R],'--k','linewidth',2)
    hold on
    plot([R,R],[R,10^4],'--k','linewidth',2)
end
if strcmp(wall, 'com')
    plot(1000,100,'ko','MarkerSize',10,'LineWidth',3)
    plot(1000,300,'ks','MarkerSize',12,'LineWidth',3)
    plot(2200,1000,'k^','MarkerSize',10,'LineWidth',3)
    % plot(2100,1000,'kd','MarkerSize',10,'LineWidth',3)
end

set(gca,'unit','centimeters','position',[3.2,2.7,9,9],'fontsize',fontn,'fontname','Times')


%%
subplot(1,3,2)
% yi = length(find(yplus > 15));
yi = length(find(yplus > 24));
fprintf('y=%d \n',yplus(yi));
c = cmax(:,:,yi);

% Divided by the local speed

if strcmp(wall, 'rigid')
    pcolor(lambdaxVector,lambdazVector,c'./U0(yi)); hold on
else
    if R == 3300
        uc=26;
    else
        uc = 28;
    end
    pcolor(lambdaxVector,lambdazVector,c'./0.53./uc); hold on
end

set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');

if strcmp(wall, 'rigid')
    title('$ y^+ \approx 25 $','Interpreter','latex');
end
if strcmp(wall, 'com')
    title('$ y^+ \approx 26 $','Interpreter','latex');
end
shading interp
load('MyColormaps.mat')
colormap(mymap);
% colormap(redblue1(cormin,cormax));
% caxis([cormin,cormax]);

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
% axis([50 5*10^5 50 5*10^5])
axis([60 10^4 60 10^4])
caxis([cormin,cormax]);

if strcmp(wall, 'rigid')
    % Draw large scale boundaries
    % hold on
    % xline(2*R,'--k','linewidth',2)
    % hold on
    % yline(0.4*R,'--k','linewidth',2)
    hold on
    plot([60,10^4],[60,10^4],'-.k','linewidth',2)
    
    hold on
    plot([60,250],[250,250],'--k','linewidth',2)
    hold on
    plot([250,250],[60,250],'--k','linewidth',2)
    
    hold on
    plot([R,10^4],[R,R],'--k','linewidth',2)
    hold on
    plot([R,R],[R,10^4],'--k','linewidth',2)
end
set(gca,'unit','centimeters','position',[15,2.7,9,9],'fontsize',fontn,'fontname','Times')
set(gca,'fontsize',fontn,'fontname','Times')

%%
subplot(1,3,3)
% yi = length(find(yplus > 105));
yi = length(find(yplus > 95));
% yi = length(find(yplus > 160));
fprintf('y=%d \n',yplus(yi));
c = cmax(:,:,yi);

% Divided by the local speed
pcolor(lambdaxVector,lambdazVector,c'./U0(yi)); hold on

set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');

% if R == 3300
%     title('$ y^+ \approx 113 $','Interpreter','latex');
% else
%     title('$ y^+ \approx 110 $','Interpreter','latex');
% end

if strcmp(wall, 'rigid')
    title('$ y^+ \approx 99 $','Interpreter','latex');
end
if strcmp(wall, 'com')
    title('$ y^+ \approx 96 $','Interpreter','latex');
end
shading interp
load('MyColormaps.mat')
colormap(mymap);
% colormap(redblue1(cormin,cormax));
% caxis([cormin,cormax]);


h=colorbar;
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;

if strcmp(wall, 'com')
    set(get(h, 'Label'),'String', '$c_{p} / c_{p,c}^M$','Interpreter', 'latex','FontSize', fontn);
    h.Ticks=[1,5,10,15,20];
else
    set(get(h, 'Label'),'String', '$c_{p} / c_{p,0}^M$','Interpreter', 'latex','FontSize', fontn);
    h.Ticks=[1,2,3,4,5];
end
h.Ticks=0:0.2:2;


set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
% axis([50 5*10^5 50 5*10^5])
axis([60 10^4 60 10^4])
caxis([cormin,cormax]);

if strcmp(wall, 'rigid')
    % Draw large scale boundaries
    % hold on
    % xline(2*R,'--k','linewidth',2)
    % hold on
    % yline(0.4*R,'--k','linewidth',2)
    
    hold on
    plot([60,10^4],[60,10^4],'-.k','linewidth',2)
    
    hold on
    plot([60,370],[370,370],'--k','linewidth',2)
    hold on
    plot([370,370],[60,370],'--k','linewidth',2)
    
    hold on
    plot([R,10^4],[R,R],'--k','linewidth',2)
    hold on
    plot([R,R],[R,10^4],'--k','linewidth',2)
end

set(gca,'unit','centimeters','position',[26.7,2.7,9,9],'fontsize',fontn,'fontname','Times')

print('-dpng','-r300', ['cmaxk',path1,'_rev.png'])





