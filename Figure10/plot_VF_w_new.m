clc
clear
% kx = 12;kz = 120;c = 10;  % near-wall cycle
kx = 1;kz = 10;c = 16;   % VLSMs

if kx == 12
    N = 300; % wall-normal grids
else
    N = 200;
end

R = 2000;  %Renolds number

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyon';
wall_vec = {'rigid','com'}; % 'com', 'rigid'
selc = 'SIN'; % 'SIN', 'RS' - objective function
effect = 'good';
% scaleeta = 300; % Magnification of wall displacement

i = sqrt(-1);
omega = c*kx;
T = 2*pi/omega;

[y,DM] = chebdif(N,2);
D1 = DM(1:N,1:N,1);
D2 = DM(1:N,1:N,2);
I = eye(N);
grd = [i*kx.*I; D1; i*kz.*I];
lplc = D2 - kx.^2.*I - kz.^2.*I;

k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
nuT0 = NuT(y);
nuT1 = D1*nuT0;
DUDy = @(y) R.*(-y)./NuT(y);
U1 = DUDy(y);
U0 = zeros(N,1);

for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end

for ifig = 1:2
    wall = wall_vec{ifig};
    
    if strcmp(method,'IOA')
        if kx == 12
            if strcmp(wall,'rigid')
                if strcmp(eddy,'eddyon')
                    phasemove =5;
                    fig1max=10; % amplitude limit
                    fig2max=8; % RS limit
                    t0 = 0.1*T; % fix time for snap
                end
                if strcmp(eddy,'eddyoff')
                    phasemove =0;fig1max=10; fig2max=4;t0 = 0.8*T;
                end
            end
            if strcmp(wall,'com')
                if strcmp(selc,'SIN')
                    phasemove = 1.9;fig1max=10;fig2max=4;t0 = 0;
                else
                    phasemove =1.5;fig1max=0.06;fig2max=4;t0 = 0.73*T;
                end
            end
        end
        if kx == 1
            if strcmp(wall,'rigid')
                phasemove =4;fig1max=0.04;fig2max=0.001;t0 = 0.4*T;
            end
            if strcmp(wall,'com')
                if strcmp(selc,'SIN')
                    phasemove =1;fig1max=10;fig2max=15;t0 = 0.85*T;
                else
                    phasemove =1;fig1max=0.06;fig2max=10;t0 = 0.73*T;
                end
            end
        end
    end
    
    if strcmp(method,'SIOA')
        if kx == 12
            if strcmp(wall,'rigid')
                if strcmp(eddy,'eddyoff')
                    phasemove =5;fig1max=0.04;fig2max=10;t0 = 0;
                end
                if strcmp(eddy,'eddyon')
                    phasemove =6;fig1max=0.04;fig2max=4;t0 = 0.2*T;
                end
            end
            if strcmp(wall,'com')
                if strcmp(eddy,'eddyoff')
                    if strcmp(selc,'SIN')
                        phasemove = 1.9;fig1max=0.04;fig2max=4;t0 = 0;
                    else
                        phasemove =5.5;fig1max=0.06;fig2max=4;t0 = 0.2*T;
                    end
                end
                if strcmp(eddy,'eddyon')
                    if strcmp(selc,'SIN')
                        phasemove = 2.3;fig1max=0.04;fig2max=4;t0 = 0.94*T;
                    else
                        phasemove =4.5;fig1max=0.06;fig2max=4;t0 = 0.24*T;
                    end
                end
            end
        end
        if kx == 1
            if strcmp(wall,'rigid')
                phasemove =1.5;fig1max=0.06;fig2max=6;t0 = 0*T;
            end
            if strcmp(wall,'com')
                if strcmp(selc,'SIN')
                    phasemove =2.5;fig1max=0.06;fig2max=6;t0 = 0.98*T;
                else
                    phasemove =1;fig1max=0.1;fig2max=0.5;t0 = 0.78*T;
                end
            end
        end
    end
    
    [Cm,Cd,Ck,Y] = readcom_final(method,eddy,wall,selc,effect,kx,omega);
    
    [y,DM] = chebdif(N,2);
    yplus = (y+1).*R;
    
    fontn = 24;
    if kx ==1
        ymax = 400;
    else
        ymax=80;
    end
    
    [u,v,w,p,f1,f2,f3,RSy,RS] = UF_up_w_new(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
    
    % UVW & RS
    wid = 36.5;
    hei = 11.5;
    figure;
    set(gcf,'unit','centimeters','position',[2 10 wid hei])
    
    subplot(1,4,1);
    plot(abs(u),yplus,'k','LineWidth',2);hold on
    plot(abs(v)*5,yplus,'--r','LineWidth',2);hold on
    plot(abs(w)*5,yplus,':b','LineWidth',2);hold on
    % plot(abs(p),yplus,'-m','LineWidth',1);hold on
    grid on
    
    axis([0,fig1max,0,ymax]);
    % legend('$\hat{u}^+$','$\hat{v}^+$','$\hat{w}^+$','$\hat{p}$','Interpreter','latex')
    legend('$\hat{u}$','$\hat{v}$','$\hat{w}$','Interpreter','latex')
    % xlabel('$|\widehat{\mathbf u}|$','Interpreter','latex');
    xlabel('$|\widehat{u}|,5|\widehat{v}|,5|\widehat{w}|$','Interpreter','latex');
    ylabel('$y^+$','Interpreter','latex');
    % xticks(0:0.02:0.06)
    
    hold on
    % yline(14,'--','color',[0.5, 0.5, 0.5]);
    
    set(gca,'unit','centimeters','position',[3,3,5,7],'fontsize',fontn,'fontname','Times');
    
    subplot(1,4,2);
    au = angle(u);
    av = angle(v);
    aw = angle(w);
    ap = angle(p);
    for j = 1:N
        au(j) = au(j) + pi + phasemove;
        if au(j) > 2 * pi
            au(j) = au(j) - 2*pi;
        end
    end
    for j = 1:N
        av(j) = av(j) + pi + phasemove;
        if av(j) > 2 * pi
            av(j) = av(j) - 2*pi;
        end
    end
    for j = 1:N
        aw(j) = aw(j) + pi + phasemove;
        if aw(j) > 2 * pi
            aw(j) = aw(j) - 2*pi;
        end
    end
    for j = 1:N
        ap(j) = ap(j) + pi + phasemove;
        if ap(j) > 2 * pi
            ap(j) = ap(j) - 2*pi;
        end
    end
    plot(au./2./pi,yplus,'k','LineWidth',2);hold on
    plot(av./2./pi,yplus,'--r','LineWidth',2);hold on
    plot(aw./2./pi,yplus,':b','LineWidth',2);hold on
    % plot(ap./2./pi,yplus,'-m','LineWidth',1);hold on
    grid on
    
    axis([0,1,0,ymax]);
    % legend('$\hat{u}$','$\hat{v}$','$\hat{w}$','$\hat{p}$','location','NorthWest','Interpreter','latex')
    xlabel('$(\angle{\widehat{\textbf{u}}})/2\pi$','Interpreter','latex');
    ylabel('$y^+$','Interpreter','latex');
    set(gca,'unit','centimeters','position',[11,3,5,7],'fontsize',fontn,'fontname','Times');
    
    subplot(1,4,3);
    if strcmp(method,'SIOA')
        plot(-1*real(RSy)*100000,yplus,'k','LineWidth',2);hold on
    else
        plot(-1*real(RSy),yplus,'k','LineWidth',2);hold on
    end
    
    grid on
    axis([-fig2max,fig2max,0,ymax]);
    xlabel('$-{\rm Re}(\hat{u}^\ast \hat{v})$','Interpreter','latex');
    ylabel('$y^+$','Interpreter','latex');
    set(gca,'unit','centimeters','position',[19,3,5,7],'fontsize',fontn,'fontname','Times');
    grid on
    
    %
    subplot(1,4,4);
    lambda_x = 2*pi/kx;
    lambda_z = 2*pi/kz;
    x = linspace(0,lambda_x,20);
    z = linspace(0,lambda_z,20);
    
    z0 = 0;
    x0 = 0;
    u_yz = real(u*exp(1i*(kx.*x0-omega.*t0)))*cos(kz.*z);
    v_yz = real(v*exp(1i*(kx.*x0-omega.*t0)))*cos(kz.*z);
    w_yz = -imag(w*exp(1i*(kx.*x0-omega.*t0)))*sin(kz.*z);
    
    pcolor(z.*R,(y+1).*R,v_yz)
    shading interp
    axis([0,lambda_z*R,-inf,ymax]);
    colormap(redblue);
    
    h=colorbar;
    h.Title.String = '$v$';
    h.Title.Interpreter = 'latex';
    h.Title.FontSize = fontn;
    
    set(gca,'Layer','top');
    if strcmp(wall,'com')
        set(gca,'unit','centimeters','position',[27,2.6,6,7.45],'fontsize',fontn,'fontname','Times');
    else
        set(gca,'unit','centimeters','position',[27,3,6,7],'fontsize',fontn,'fontname','Times');
    end
    
    xlabel('${z^+}$','Interpreter','latex');
    ylabel('${y^+}$','Interpreter','latex');
    grid on
    
    hold on
    quiver(z.*R,(y+1).*R,w_yz,v_yz,0.6,'Color','k','LineWidth',0.7);
    
    % % plot wall displacement
    % scaleeta = 10000;
    % if strcmp(wall,'com')
    %
    %     eta = 1i/omega*v(end,1);
    %
    %     eta_yz = real(eta*exp(1i*(kx.*x0-omega.*t0)))*cos(kz.*z);
    %     % eta_yz = real((eta*exp(1i*(kx.*x0+kz*z-omega.*t0))));  % wall displacement
    %     hold on
    %     plot(z.*R,eta_yz.*R*scaleeta,'k','LineWidth',2)
    % end
    
    %
    path = [method,'_',eddy,'_',wall,'/uvw_RS_kx=',num2str(kx),'_kz=',num2str(kz),'_c=',num2str(c),'_',selc,'_',effect,'_rev_H0'];
    print(gcf, [path,'.png'], '-dpng', '-r300');
end
