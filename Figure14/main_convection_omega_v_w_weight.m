% fix the weight w bug
clear
clc

N = 300;
R = 3300;
nc = 201;
nkx = 100;
nkz = 100;
flowtype = 'com3300';
method = 'IOA';
eddy = 'eddyoff';
wall = 'com';

uvwp = 'v';


if strcmp(wall, 'com')
%     Cm = 0.57; Ck = 45.86; Cd=0.77;
%    Cm = 0.46; Ck = 73.03; Cd=0.116;
    Cm = 0.46; Ck = 181; Cd=0.091;
%     [y,~] = chebdif(N,2);
%     yplus = (y+1)*R;
%     I = eye(N);
%     k=0.426;
%     alpha=25.4;
%     NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
%     DUDy = @(y) R.*(-y)./NuT(y);
%     U0 = zeros(N,1);
%     for j=1:N
%         U0(j) = integral(DUDy,-1,y(j));
%     end
    
if R ==3300
    uc = 26;
    load('COMUmean3300.mat','Umean');
    uc = Umean(N/2);
end
if R == 6700
    uc = 28;
    load('COMUmean6700.mat','Umean');
    uc = Umean(N/2);
end
if R ==8900
    uc = 28;
    load('COMUmean8900.mat','Umean');
    uc = Umean(N/2);
end

    Ck1 = Ck*uc^2;
    Cd1 = Cd*uc;

end

omegaVector = linspace(0.1,3000,nc);

psd = zeros(nkx,nkz,nc,N);
psdmax = zeros(N,nkz,nkx);
cmax = zeros(N,nkz,nkx);

% parpool(60);
% parfor x = 1:nkx
for x = 1:nkx
    kxVector = linspace(0,200,nkx);
    kzVector = linspace(0,200,nkz);
    psd_z = zeros(nkz,nc,N);
    psdmax_z = zeros(N,nkz);
    cmax_z = zeros(N,nkz);
    
    kx = kxVector(x);
    for z = 1:nkz
        
        kz = kzVector(z);
        psdTmp = zeros(nc,N);
        
        for cc = 1:nc
            fprintf('x=%d,z=%d,c=%d\n',x,z,cc)
            omega = omegaVector(cc);

            c = omega / kx;
            if strcmp(wall, 'com')
            [ReY,ImY] = C2Y(Cm,Ck1,Cd1,omega); % Y是omega的函数
            Y = ReY + sqrt(-1)*ImY;
            else
                Y=0;
            end
            psdTmp(cc,:) = PSD_w_weight(kx,kz,c,N,R,flowtype,method,eddy,wall,Y,uvwp);
        end
        
        psd_z(z,:,:) = psdTmp;
        [psdmaxTmp,maxci] = max(psdTmp);
        cmaxTmp = omegaVector(maxci)/kx;
        psdmax_z(:,z) = psdmaxTmp';
        cmax_z(:,z) = cmaxTmp';
    end
    psd(x,:,:,:) = psd_z;
    psdmax(:,:,x) = psdmax_z;
    cmax(:,:,x) = cmax_z;
end
% delete(gcp('nocreate'));

psdmax = permute(psdmax,[3,2,1]);
cmax = permute(cmax,[3,2,1]);
if strcmp(wall, 'com')
    path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'_w_weight'];
    save([path,'.mat'],'psd','-v7.3');
    path = [wall,'/',uvwp,'/psdmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'_w_weight'];
    save([path,'.mat'],'psdmax');
    path = [wall,'/',uvwp,'/cmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'_w_weight'];
    save([path,'.mat'],'cmax');
else
    path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_w_weight'];
    save([path,'.mat'],'psd','-v7.3');
    path = [wall,'/',uvwp,'/psdmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_w_weight'];
    save([path,'.mat'],'psdmax');
    path = [wall,'/',uvwp,'/cmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_w_weight'];
    save([path,'.mat'],'cmax');
end

