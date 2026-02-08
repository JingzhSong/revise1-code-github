clear
clc

kz = 0;
N = 200;  % wall-normal grids
R = 117;  %Renolds number, Xia
% R = 137; %Renolds number, Kim
k = 0.426; alpha = 25.4;
i = sqrt(-1);

flowtype = 'turbu';
method = 'SIOA'; % 'IOA', 'SIOA'
eddy = 'eddyon'; % 'eddyoff', 'eddyon'
wall = 'rigid';

NY = 50;
kxvector = logspace(0,2,NY);
cvector = linspace(1,20,NY);
NormRigid = zeros(NY,NY);


%% 计算不同波数
% p = parpool(50);
% parfor cc = 1:NY
for cc = 1:NY
    Normtmp = zeros(1,NY);
    c = cvector(cc);
    cc
    for kxx = 1:NY
        kx = kxvector(kxx);
        kxx
        Y = 0;
        
        Normtmp(kxx) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
    end
    NormRigid(cc,:) = Normtmp;
end
% delete(p);
path = [method,'_',eddy,'_',wall,'/Norm_DNS_Xia_kz',num2str(kz)];
% save([path,'.mat'],'NormRigid');
