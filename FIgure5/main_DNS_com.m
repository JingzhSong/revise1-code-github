clear
clc

kz = 0;
N = 200;
% R = 137;  %Renolds number, Kim
R = 117;  %Renolds number, Xia
k = 0.426; alpha = 25.4;
i = sqrt(-1);

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyon';
wall = 'com';

% Uc=21; %Kim
Uc=19;  %Xia

% Xia
Cm = 1; 
% Ck = 0.2*Uc*Uc;  %case D 
Ck = 1.68*Uc*Uc; %case C B
Cd = 0.25*Uc; %case C D 
% Cd = 0.47*Uc; %case B
Cb = 1.2*10^(-4)*Uc*Uc;
Ct = 0.0033*Uc*Uc;

% Kim
% Cm = 1; 
% Ck = 1*Uc*Uc; 
% Cd = 0.5*Uc;
% Cb = 0;
% Ct = 0;
% % Cb = 1.2*10^(-4)*Uc*Uc;
% % Ct = 0.0033*Uc*Uc;

NY = 50;
kxvector = logspace(0,2,NY);
cvector = linspace(1,20,NY); 
NormCom = zeros(NY,NY);
Omega_r = zeros(NY,NY);

%% 计算不同波数
% p = parpool(50);
% parfor cc = 1:NY
for cc = 1:NY
    Normtmp = zeros(1,NY);
    omega_r = zeros(1,NY);
    c = cvector(cc);
    cc
    for kxx = 1:NY
        kx = kxvector(kxx);
        kxx
        omega = c*kx;
        Ck_equ = Ck + Cb*(kx^4+kz^4+2*kx^2*kz^2) + Ct*(kx^2+kz^2);
        ReY = -(omega^2*Cd)/((Ck_equ-omega^2*Cm)^2+(omega*Cd)^2);
        ImY = omega*(Ck_equ-omega^2*Cm)/((Ck_equ-omega^2*Cm)^2+(omega*Cd)^2);
        Y = ReY + i*ImY;
        
        Normtmp(kxx) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
        
        zeta = Cd/2/sqrt(Ck_equ*Cm);
        omega_r(kxx) = sqrt(Ck_equ/Cm * (1 - 2*(Cd/2/sqrt(Ck_equ*Cm))^2) );
    end
    NormCom(cc,:) = Normtmp;
    Omega_r(cc,:) = omega_r;
end
% delete(p);
path = [method,'_',eddy,'_',wall,'/Norm_DNS_XiaC_kz',num2str(kz)];
% save([path,'.mat'],'NormCom');


% contour(Omega_r)
% omega_r(imag(omega_r)~=0) = -1;
% loglog(kxvector, omega_r,'LineWidth',2)
% grid on
% set(gcf,'unit','centimeters','position',[10 7 16 13]);
% set(gca,'unit','centimeters','position',[3,3,10,9],'fontsize',20,'fontname','Times')
% % axis([1,50,0,700]);
% axis([1,100,10,10000]);
% set(gca,'XTick',[1,10,100]);
% xlabel('$k_{x}$','Interpreter','latex');
% ylabel('$\omega_r$','Interpreter','latex');
% % path = [method,'_',eddy,'_',wall,'/Omegar_DNS_Kim_kz0'];
% path = [method,'_',eddy,'_',wall,'/Omegar_DNS_Kim_kz',num2str(kz)];
% saveas(gcf, [path,'.png']);








