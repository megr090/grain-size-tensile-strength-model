function [SRmat,Tmat,Vmat,SMBmat,Xfix,Yfix,R] = readData(icestream)
% Velocity Contours
[Esing,R] = geotiffread('strainrates_in_seconds.tif');
E=double(Esing);

[Vsing,R_V] = geotiffread('surface_velocity_2014_2015_lowerres_fixedextent.tif');
V=double(Vsing);

[Tsing,R] = geotiffread('bedmachinev01_thickness_fixedres.tif');
T = double(Tsing);

[SMB,R] = geotiffread('smb_racmo_23_ant27_mean_fixedextent.tif');

if strcmp(icestream,'BM') == 1
    xtot = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(1));
    ytot = linspace(R.YWorldLimits(1),R.YWorldLimits(2),R.RasterSize(2));
    [X,Y] = meshgrid(xtot,ytot);
    
    X = rot90(X');
    Y = rot90(Y');
    Xfix = X(3500:5000,2000:4000);
    Xfix = Xfix(350:1200,200:1200);
    Xfix = Xfix(1:500,150:700);
    
    Yfix = Y(3500:5000,2000:4000);
    Yfix = Yfix(350:1200,200:1200);
    Yfix = Yfix(1:500,150:700);
    
    BM_sr = E(3500:5000,2000:4000);
    BM_sr = BM_sr(350:1200,200:1200);
    BM_sr = BM_sr(1:500,150:700);
    
    BM_t = T(3500:5000,2000:4000);
    BM_t = BM_t(350:1200,200:1200);
    BM_t = BM_t(1:500,150:700);
    
    BM_v = V(3500:5000,2000:4000);
    BM_v = BM_v(350:1200,200:1200);
    BM_v = BM_v(1:500,150:700);
    
    BM_smb = SMB(3500:5000,2000:4000);
    BM_smb = BM_smb(350:1200,200:1200);
    BM_smb = BM_smb(1:500,150:700);

    SRmat = BM_sr;
    Vmat = BM_v;
    Tmat = BM_t;
    SMBmat = BM_smb;

elseif strcmp(icestream,'Byrd') == 1
    xtot = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(1));
    ytot = linspace(R.YWorldLimits(1),R.YWorldLimits(2),R.RasterSize(2));
    [X,Y] = meshgrid(xtot,ytot);
    
    Xfix = X(4100:4400,3600:3900);
    Xfix = Xfix(100:300,50:200);
    
    Yfix = Y(4100:4400,3600:3900);
    Yfix = Yfix(100:300,50:200);
    
    Byrd_sr = E(4100:4400,3600:3900);
    Byrd_sr = Byrd_sr(100:300,50:200);
    
    Byrd_t = T(4100:4400,3600:3900);
    Byrd_t = Byrd_t(100:300,50:200);
    
    Byrd_v = V(4100:4400,3600:3900);
    Byrd_v = Byrd_v(100:300,50:200);
    
    Byrd_smb = SMB(4100:4400,3600:3900);
    Byrd_smb = Byrd_smb(100:300,50:200);
    
    SRmat = Byrd_sr;
    Vmat = Byrd_v;
    Tmat = Byrd_t;
    SMBmat = Byrd_smb;

elseif strcmp(icestream,'Recovery') == 1
    xtot = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(1));
    ytot = linspace(R.YWorldLimits(1),R.YWorldLimits(2),R.RasterSize(2));
    [X,Y] = meshgrid(xtot,ytot);
    
    Xfix = X(2200:2700,2700:3200);
    Yfix = Y(2200:2700,2700:3200);
    R_sr = E(2200:2700,2700:3200);
    R_t = T(2200:2700,2700:3200);
    R_v = V(2200:2700,2700:3200);
    R_smb = SMB(2200:2700,2700:3200);
    
    SRmat = R_sr;
    Vmat = R_v;
    Tmat = R_t;
    SMBmat = R_smb;

elseif strcmp(icestream,'Amery') == 1
    xtot = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(1));
    ytot = linspace(R.YWorldLimits(1),R.YWorldLimits(2),R.RasterSize(2));
    [X,Y] = meshgrid(xtot,ytot);
    
    Xfix = X(2000:3000,4500:6000);
    Xfix = Xfix(400:800,200:1200);
    Xfix = Xfix(:,200:800);
    
    Yfix = Y(2000:3000,4500:6000);
    Yfix = Yfix(400:800,200:1200);
    Yfix = Yfix(:,200:800);

    Am_sr = E(2000:3000,4500:6000);
    Am_sr = Am_sr(400:800,200:1200);
    Am_sr = Am_sr(:,200:800);
    
    Am_t = T(2000:3000,4500:6000);
    Am_t = Am_t(400:800,200:1200);
    Am_t = Am_t(:,200:800);
    
    Am_v = V(2000:3000,4500:6000);
    Am_v = Am_v(400:800,200:1200);
    Am_v = Am_v(:,200:800);

    Am_smb = SMB(2000:3000,4500:6000);
    Am_smb = Am_smb(400:800,200:1200);
    Am_smb = Am_smb(:,200:800);

    SRmat = Am_sr;
    Vmat = Am_v;
    Tmat = Am_t;
    SMBmat = Am_smb;

elseif strcmp(icestream,'PIG') == 1
    xtot = linspace(R_V.XWorldLimits(1),R_V.XWorldLimits(2),R_V.RasterSize(1));
    ytot = linspace(R_V.YWorldLimits(1),R_V.YWorldLimits(2),R_V.RasterSize(2));
    [X,Y] = meshgrid(xtot,ytot);
    
    X = rot90(X');
    Y = rot90(Y');
    Xfix = X(3000:4000,1000:2000);
    Xfix = Xfix(300:800,600:900);
    Xfix = Xfix(1:400,:);
    
    Yfix = Y(3000:4000,1000:2000);
    Yfix = Yfix(300:800,600:900);
    Yfix = Yfix(1:400,:);
    
    PIG_sr = E(3000:4000,1000:2000);
    PIG_sr = PIG_sr(300:800,600:900);
    PIG_sr = PIG_sr(1:400,:);
    
    PIG_t = T(3000:4000,1000:2000);
    PIG_t = PIG_t(300:800,600:900);
    PIG_t = PIG_t(1:400,:);
    
    PIG_v = V(3000:4000,1000:2000);
    PIG_v = PIG_v(300:800,600:900);
    PIG_v = PIG_v(1:400,:);

    PIG_smb = SMB(3000:4000,1000:2000);
    PIG_smb = PIG_smb(300:800,600:900);
    PIG_smb = PIG_smb(1:400,:);

    SRmat = PIG_sr;
    Vmat = PIG_v;
    Tmat = PIG_t;
    SMBmat = PIG_smb;
    
else
    disp('Incorrect Ice Stream Choice')
    SRmat = 0;
    Tmat = 0;
    Vmat = 0;
    SMBmat = 0;
    Xfix = 0;
    Yfix = 0;
    return
end

end

