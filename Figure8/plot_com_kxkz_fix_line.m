%% 画比值
clear
clc
N = 200;
R = 2000;
Nc = 200;
kxsel = 12;

method = 'SIOA';
eddy = 'eddyon';
wall = 'com';
selc = 'SIN'; %用什么做优化
effect = 'good';

kx_vec = [12,1,6.2832];
kz_vec = [120,10,0.1257];

for ifig = 1:3
    kx = kx_vec(ifig);
    kz = kz_vec(ifig);
    path = [method,'_',eddy,'_','com','/singular_fixkxkz_line_N=',num2str(N),'_R=',num2str(R),'_kxsel=',num2str(kxsel),'_Nc=',num2str(Nc),selc,'_',effect,'_kx',num2str(kx),'_kz=',num2str(kz)];
    load([path,'.mat'],'singular_com');
    path = [method,'_',eddy,'_','rigid','/singular_fixkxkz_line_N=',num2str(N),'_R=',num2str(R),'_kxsel=',num2str(kxsel),'_Nc=',num2str(Nc),selc,'_',effect,'_kx',num2str(kx),'_kz=',num2str(kz)];
    load([path,'.mat'],'singular_rigid');
    cvector = linspace(0,24,Nc);
    
    figure;
    
    h1=plot(cvector(2:end),singular_com(2:end),'-','LineWidth',2);
    hold on
    h2=plot(cvector(2:end),singular_rigid(2:end),'--','LineWidth',3);
    
    % 画共振频率
    [Cm,Cd,Ck,Y] = readcom_final(method,eddy,'com',selc,effect,kxsel,1);
    omegar = sqrt(Ck/Cm*(1-2*(Cd/2/sqrt(Ck*Cm))^2)); %resonance
    hold on
    h3 = xline(omegar./kx,'k','LineWidth',0.8);
    omegar./kx
    if kx == 6.2832
    legend([h1,h2],{'Compliant','Rigid'}, 'Location', 'northwest')
    else
    legend([h1,h2],{'Compliant','Rigid'})
    end
    
    set(gca,'yscale','log')
    
    if kx ==12
        axis([0,24,1.5,6])
        yticks([2:6])
        xticks(linspace(0,24,7))
    else
        % axis([0,24,1,2000])
        % yticks([1,10,100,1000])
        axis([0,24,1,100])
        yticks([1,10,100])
        yticklabels({'10^0', '10^1', '10^2'});
        xticks(linspace(0,24,7))
    end
    
    grid on
    
    % xlabel('$k_x$','Interpreter','latex');
    
    xlabel('$c$','Interpreter','latex');
    ylabel('$\mu_0^e, \, \mu_c^e$','Interpreter','latex');
    
    fontn = 20;
    set(gcf,'unit','centimeters','position',[10 7 14 13.5]);
    set(gca,'unit','centimeters','position',[3.5,3,9.5,9],'fontsize',fontn,'fontname','Times')
    hold on
    grid on
    
    path = [method,'_',eddy,'_','com','/singular_fixkxkz_line_N=',num2str(N),'_R=',num2str(R),'_kxsel=',num2str(kxsel),'_Nc=',num2str(Nc),selc,'_',effect,'_kx',num2str(kx),'_kz=',num2str(kz)];
    saveas(gcf, [path,'.png']);
    
end

% %% This section is to compute the data
% clear
% clc

% clear
% clc
% N = 200;
% R = 2000;
%
% % kz = 0.1257; kx = 6.2832;
% % kz = 10; kx = 1;
% kz = 120; kx = 12;
%
% Nc = 200;
%
% flowtype = 'turbu';
% method = 'SIOA';
% eddy = 'eddyon';
% wall = 'com';
% selc = 'SIN'; 
% effect = 'good';
% kxsel = 1;
%
% % lambdaxVector = logspace(0.5,6.5,Nkx);
% % kxvector = logspace(-0.5,2.5,Nkx);
% cvector = linspace(0,24,Nc);
%
% singular_com = zeros(1,Nc);
% for cc = 2:Nc 
%     c = cvector(cc)
%     omega = c*kx;
%     [Cm,Cd,Ck,Y] = readcom_final(method,eddy,wall,selc,effect,kxsel,omega);
%     singular_com(1,cc) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
%
% end
%
% path = [method,'_',eddy,'_',wall,'/singular_fixkxkz_line_N=',num2str(N),'_R=',num2str(R),'_kxsel=',num2str(kxsel),'_Nc=',num2str(Nc),selc,'_',effect,'_kx',num2str(kx),'_kz=',num2str(kz)];
% save([path,'.mat'],'singular_com');
%
%
% wall = 'rigid';
% Y = 0;
% singular_rigid = zeros(1,Nc);
% for cc = 2:Nc 
%
%     c = cvector(cc)
%     omega = c*kx;
%     singular_rigid(1,cc) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
%
% end
% path = [method,'_',eddy,'_',wall,'/singular_fixkxkz_line_N=',num2str(N),'_R=',num2str(R),'_kxsel=',num2str(kxsel),'_Nc=',num2str(Nc),selc,'_',effect,'_kx',num2str(kx),'_kz=',num2str(kz)];
% save([path,'.mat'],'singular_rigid');
