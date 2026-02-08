% fix kx,kz，psd varies with c at different y
clear
clc

N = 300;
R = 3300;
nc = 201;

nkx = 120;
nkz = 120;
method = 'IOA';
eddy = 'eddyoff';
wall = 'com';
Cm = 0.46; Ck = 181; Cd=0.091;
Y = [];
uvwp = 'p';
lambdax_vec = [2100,1000,1000];
lambdaz_vec = [1000,100,300];
% lambdax_vec = [2100];
% lambdaz_vec = [1000];

for ifig = 1:size(lambdax_vec,2)
    lambdax = lambdax_vec(ifig);
    lambdaz = lambdaz_vec(ifig);
    
    kx = 2.*pi./lambdax .*R;
    kz = 2.*pi./lambdaz .*R;
    
    [y,~] = chebdif(N,2);
    yplus = (y+1)*R;
 
    % Small scales are not considered
    
    cVector = linspace(0,30,nc);
    
    if strcmp(wall, 'com')
        if R ==3300
            uc = 26;
        end
        if R == 6700
            uc = 28;
        end
        if R ==8900
            uc = 28;
        end
        Ck1 = Ck*uc;
        Cd1 = Cd*uc;
    end
    
    %% Computation process
    wall = 'com';
    flowtype = 'com3300';
    psd_vector_com = zeros(N,nc);
    for ic=1:nc
        c = cVector(ic);
        omega = c*kx;
        if strcmp(wall, 'com')
            [ReY,ImY] = C2Y(Cm,Ck1,Cd1,omega);
            Y = ReY + sqrt(-1)*ImY;
        else
            Y=0;
        end
        psd_vector_com(:,ic) = PSD_w_weight(kx,kz,c,N,R,flowtype,method,eddy,wall,Y,uvwp);
    
    end
    path0 = [wall,'/',uvwp,'/psd_c_y_N=',num2str(N),'_R=',num2str(R),'_lambdax=',num2str(lambdax),'_lambdaz=',num2str(lambdaz)];
    save([path0,'.mat'],'psd_vector_com');    
    
    wall = 'rigid';
    flowtype = 'turbu';
    psd_vector_rigid = zeros(N,nc);
    for ic=1:nc
        c = cVector(ic);
        omega = c*kx;
        if strcmp(wall, 'com')
            [ReY,ImY] = C2Y(Cm,Ck1,Cd1,omega);
            Y = ReY + sqrt(-1)*ImY;
        else
            Y=0;
        end
        psd_vector_rigid(:,ic) = PSD_w_weight(kx,kz,c,N,R,flowtype,method,eddy,wall,Y,uvwp);
    end
    
    path0 = [wall,'/',uvwp,'/psd_c_y_N=',num2str(N),'_R=',num2str(R),'_lambdax=',num2str(lambdax),'_lambdaz=',num2str(lambdaz)];
    save([path0,'.mat'],'psd_vector_rigid');
    
    % plot figure
    kx = 2*pi*R/lambdax;
    kz = 2*pi*R/lambdaz;
    wall = 'com';
    path0 = [wall,'/',uvwp,'/psd_c_y_N=',num2str(N),'_R=',num2str(R),'_lambdax=',num2str(lambdax),'_lambdaz=',num2str(lambdaz)];
    load([path0,'.mat'],'psd_vector_com');
    
    wall = 'rigid';
    path0 = [wall,'/',uvwp,'/psd_c_y_N=',num2str(N),'_R=',num2str(R),'_lambdax=',num2str(lambdax),'_lambdaz=',num2str(lambdaz)];
    load([path0,'.mat'],'psd_vector_rigid');
    
    % y_location = [0,27,98];
    y_location = [0];
    for i_y = 1:size(y_location)
        y_loc = y_location(i_y);
        yi = length(find(yplus > y_loc))+1;
        fprintf('y=%d \n',yplus(yi));
        
        figure;
        psd_com_tmp = psd_vector_com(yi,2:end);
        h1=plot(cVector(2:end),psd_com_tmp,'-','LineWidth',2);
        psd_rigid_tmp = psd_vector_rigid(yi,2:end);
        hold on
        h2=plot(cVector(2:end),psd_rigid_tmp,'--','LineWidth',3);
        
        omegar = sqrt(Ck1/Cm*(1-2*(Cd1/2/sqrt(Ck1*Cm))^2)); %resonance
        hold on
        h3 = xline(omegar./kx,'k','LineWidth',1);
        fprintf('omegar=%d \n',omegar./kx);
        
        % convective velocity
        [psdmaxTmp,maxci] = max(psd_rigid_tmp);
        cmax_rigid = cVector(maxci);
        [psdmaxTmp,maxci] = max(psd_com_tmp);
        cmax_com = cVector(maxci);
        
        % h4 = xline(cmax_com,'-','Color','b','LineWidth',1.5);
        % hold on
        % h4 = xline(cmax_rigid,'--','Color','r','LineWidth',1.5);
        legend([h1,h2],{'Compliant','Rigid'})
        set(gca,'yscale','log')
        
        if i_y ==1
            if lambdaz==100
                axis([0,30,0.5*10^(-7),2*10^(-7)])
                yticks(linspace(0,2*10^(-7),5))
            end
            if lambdaz==300
                axis([0,30,1*10^(-7),1*10^(-4)])
                yticks([0,10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2)])
            end
            if lambdaz==1000
                axis([0,30,1*10^(-7),2*10^(-3)])
                yticks([0,10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2)])
            end
        end
        
        xticks(linspace(0,30,7))
        grid on
        
        xlabel('$c$','Interpreter','latex');
        ylabel('$\Phi_p$','Interpreter','latex');
        
        fontn = 24;
        set(gcf,'unit','centimeters','position',[10 7 14 13.5]);
        set(gca,'unit','centimeters','position',[3,3,10.5,9],'fontsize',fontn,'fontname','Times')
        hold on
        grid on
        
        wall = 'com';
        path0 = [wall,'/',uvwp,'/psd_c_y_N=',num2str(N),'_R=',num2str(R),'_yplus=',num2str(yplus(yi)),'_lambdax=',num2str(lambdax),'_lambdaz=',num2str(lambdaz)];
        print('-dpng','-r300', [path0,'_rev_n400.png'])
    end
end