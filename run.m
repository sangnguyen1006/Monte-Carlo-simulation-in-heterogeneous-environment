%***********MONTE**CARLO**********/
% 2/24/2022
%
% run.m	
% Usage:  run.m 
% which loads myname_H.mci, and saves myname_F.bin, mynameRd_yx.bin
% mynameTt_yx.bin and myname_folder
% 
% created 2010, 2012 by
% Steven L. JACQUES
% Ting LI
% Oregon Health & Science University
%
% USAGE   run.m 
%         myname is the user's choice. 
%         The program reads two files prepared by user:
%                 myname_H.mci    = header input file for mcxyz
%                 myname_T.bin    = tissue structure file
%         The output will be written to 2 files:
%                 myname_OP.m     = optical properties  (mua, mus, g for each tissue type)
%                 myname_F.bin    = fluence rate output F[i] [W/cm^2 per W delivered]
% 
% The MATLAB program maketissue.m can create the two input files (myname_H.mci, myname_T.bin).
% 
% The MATLAB program lookmcxyz.m can read the output files and display
%         1. Fluence rate F [W/cm^2 per W delivered]
%         2. Deposition rate A [W/cm^3 per W delivered].
% 
% Log:
% Written by Ting based on Steve's mcsub.c., 2010.
%     Use Ting's FindVoxelFace().
% Use Steve's FindVoxelFace(), Dec. 30, 2010.
% Reorganized by Steve. May 8, 2012:
%     Reads input files, outputs binary files.
% **********/
 

%% INPUT FILES 
% IMPORT myname_H.mci 
filename = input('Please enter myname: ');  % myname
fileID   = fopen(strcat(filename,'_H.mci'),'r');
    % run parameters
    tline = fgetl(fileID);
    Nphotons = sscanf(tline,'%d'); % number of photons used for simulation
    
    tline = fgetl(fileID);
    Nx = sscanf(tline,'%d');       % # of bins
    tline = fgetl(fileID);
    Ny = sscanf(tline,'%d');       % # of bins
    tline = fgetl(fileID);
    Nz = sscanf(tline,'%d');       % # of bins
    tline = fgetl(fileID);
    Nr = sscanf(tline,'%d');       % # of bins
    
    tline = fgetl(fileID);
    dx = sscanf(tline,'%f');       % size of bins [cm]
    tline = fgetl(fileID);
    dy = sscanf(tline,'%f');       % size of bins [cm]
    tline = fgetl(fileID);
    dz = sscanf(tline,'%f');       % size of bins [cm]
    tline = fgetl(fileID);
    dr = sscanf(tline,'%f');       % size of bins [cm]
    
    % launch parameters
    tline = fgetl(fileID);
    mcflag = sscanf(tline,'%d');       % mcflag, 0 = uniform, 1 = Gaussian, 2 = iso-pt
    tline = fgetl(fileID);
    launchflag = sscanf(tline,'%d');   % launchflag, 0 = ignore, 1 = manually set
    tline = fgetl(fileID);
    boundaryflag = sscanf(tline,'%d'); % 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only
    
    tline = fgetl(fileID);
    xs = sscanf(tline,'%f');       % initial launch point
    tline = fgetl(fileID);
    ys = sscanf(tline,'%f');       % initial launch point
    tline = fgetl(fileID);
    zs = sscanf(tline,'%f');       % initial launch point
    
    tline = fgetl(fileID);
    xfocus = sscanf(tline,'%f');   % xfocus
    tline = fgetl(fileID);
    yfocus = sscanf(tline,'%f');   % yfocus
    tline = fgetl(fileID);
    zfocus = sscanf(tline,'%f');   % zfocus
    
    tline = fgetl(fileID);
    zref = sscanf(tline,'%f');     % plane position recorded reflectance in z
    tline = fgetl(fileID);
    ztra = sscanf(tline,'%f');     % transmittance recorded plane position in z
    
    tline = fgetl(fileID);
    ux0 = sscanf(tline,'%f');      % ux trajectory
    tline = fgetl(fileID);
    uy0 = sscanf(tline,'%f');      % uy trajectory
    tline = fgetl(fileID);
    uz0 = sscanf(tline,'%f');      % uz trajectory
    
    tline = fgetl(fileID);
    radius = sscanf(tline,'%f');   % radius
    tline = fgetl(fileID);
    waist = sscanf(tline,'%f');    % waist
    
    % tissue optical properties
    tline = fgetl(fileID);
    Nt = sscanf(tline,'%d');       % # of tissue types in tissue list
    muav = zeros(1,Nt);
    musv = zeros(1,Nt);
    gv = zeros(1,Nt);
    for i=1:Nt
        tline = fgetl(fileID);
        muav(i) = sscanf(tline,'%f'); % absorption coeff [cm^-1]
        tline = fgetl(fileID);
        musv(i) = sscanf(tline,'%f'); % scattering coeff [cm^-1]
        tline = fgetl(fileID);
        gv(i) = sscanf(tline,'%f');   % anisotropy of scatter [dimensionless]
    end
    
fclose(fileID);

fprintf('name = %s\n ',filename);
fprintf('number of photons = %d\n',Nphotons);
fprintf('Nx = %d, dx = %f [cm]\n',Nx,dx);
fprintf('Ny = %d, dy = %f [cm]\n',Ny,dy);
fprintf('Nz = %d, dz = %f [cm]\n',Nz,dz);

fprintf('xs = %f [cm]\n',xs);
fprintf('ys = %f [cm]\n',ys);
fprintf('zs = %f [cm]\n',zs);
fprintf('mcflag = %d\n',mcflag);
if mcflag == 0
    fprintf('launching uniform flat-field beam\n'); 
end
if mcflag == 1
    fprintf('launching Gaissian beam\n'); 
end 
if mcflag == 2
    fprintf('launching isotropic point source\n'); 
end
if mcflag == 3
    fprintf('launching square source\n');
end
fprintf('xfocus = %f [cm]\n',xfocus);
fprintf('yfocus = %f [cm]\n',yfocus);
fprintf('zfocus = %0.2e [cm]\n',zfocus);
fprintf('zref = %f [cm]\n',zref);
fprintf('ztra = %f [cm]\n',ztra);
if launchflag == 1 
    fprintf('Launchflag ON, so launch the following:\n');
	fprintf('ux0 = %f [cm]\n',ux0);
	fprintf('uy0 = %f [cm]\n',uy0);
	fprintf('uz0 = %f [cm]\n',uz0);
else 
	fprintf('Launchflag OFF, so program calculates launch angles.\n');
	fprintf('radius = %f [cm]\n',radius);
	fprintf('waist  = %f [cm]\n',waist);
end
if boundaryflag == 0
    fprintf('boundaryflag = 0, so no boundaries.\n');
else
    if boundaryflag == 1
	    fprintf('boundaryflag = 1, so escape at all boundaries.\n');
    else
        if boundaryflag == 2
		    fprintf('boundaryflag = 2, so escape at surface only.\n');    
        else
            fprintf('improper boundaryflag. quit.\n');
        end
    end
end
fprintf('# of tissues available, Nt = %d\n',Nt);
for i=1:Nt
    fprintf('muav[%ld] = %0.4f [cm^-1]\n',i,muav(i));
    fprintf('musv[%ld] = %0.4f [cm^-1]\n',i,musv(i));
    fprintf('  gv[%ld] = %0.4f [--]\n\n',i,gv(i));
end

% SAVE optical properties, for later use by MATLAB
fileID = fopen(strcat(filename,'_props.m'),'w');
for i=1:Nt
    fprintf(fileID,'muav(%ld) = %0.4f;\n',i,muav(i));
    fprintf(fileID,'musv(%ld) = %0.4f;\n',i,musv(i));
    fprintf(fileID,'  gv(%ld) = %0.4f;\n\n',i,gv(i));
end
fclose(fileID);

% IMPORT BINARY TISSUE FILE
% read binary file
fileID = fopen(strcat(filename,'_T.bin'),'rb');
v = fread(fileID,[1 inf],'uint8'); % tissue structure
fclose(fileID);

% Show tissue on screen, along central z-axis, by listing tissue type #'s
iy = Ny/2;
ix = Nx/2;
fprintf('central axial profile of tissue types:\n');
for iz=1:Nz
    i = round((iz-1)*Ny*Nx + (ix-1)*Ny + iy);
    fprintf('%d',v(i));
end
fprintf('\n\n');

%% SAVE INPUT PARAMETERS
% "parameters.mat" and "btissue.mat"
% run parameters
% launch parameters
% tissue optical properties
save parameters.mat Nx Ny Nz Nr dx dy dz dr mcflag launchflag boundaryflag...
    xs ys zs xfocus yfocus zfocus zref ztra ux0 uy0 uz0 radius waist muav musv gv
% tissue structure
save bintissue.mat v

%% PERFORM PARALLEL COMPUTATIONS
% Perform parallel computations on multicore computers
% Using Parallel Computing Toolbox
% RUN function F = mcxyz(Nphotons_worker)
now = datetime('now');
fprintf('%s\n', now);
start_time = tic;

NN = Nx*Ny*Nz;
F          = zeros(1,NN);            % ensure F_total[] starts empty
Ftemp      = zeros(1,NN);     
Q_total    = zeros(1,NN);            % ensure Q_total[] starts empty
Rd_total   = zeros(Ny,Nx);           % diffuse reflection weight
Tt_total   = zeros(Ny,Nx);           % transmitted weight
Rd_r_total = zeros(1,Nr);            % diffuse reflectance weight with radius
Tt_r_total = zeros(1,Nr);            % transmitted weight with radius

M = maxNumCompThreads;               % maximum number of workers
Nphotons_worker = round(Nphotons/M); % number of simulated photons per worker
ticBytes(gcp);
parfor (i = 1:M,M)
    [Q_temp, Rd_temp, Tt_temp, Rd_r_temp, Tt_r_temp] = mcxyz(Nphotons_worker);
    Q_total    = Q_total    + Q_temp;     % tong trong luong photon hap thu moi voxel [-]
    Rd_total   = Rd_total   + Rd_temp;    % tong trong luong photon phan xa khuech tan yx [-]
    Tt_total   = Tt_total   + Tt_temp;    % tong trong luong photon truyen qua yx [-]
    Rd_r_total = Rd_r_total + Rd_r_temp;  % tong trong luong photon phan xa khuech tan theo phuong r [-]
    Tt_r_total = Tt_r_total + Tt_r_temp;  % tong trong luong photon truyen qua theo phuong r [-]
end
tocBytes(gcp)

fprintf('------------------------------------------------------\n');
finish_time = toc(start_time);
time_min = finish_time/60;
fprintf('Elapsed Time for %0.3e photons = %5.4f min\n',Nphotons,time_min);
fprintf('%0.2e photons per minute\n', Nphotons/time_min);

%% RAW DATA PROCESSING
% Ma tran [1,Ndr],[1,Nda] can chuyen thanh [0,Ndr-1],[0,Nda-1] bang cach tru 1, 
% dung tinh dien tich hinh khuyen va goc thoat alpha cua photon
ir = 0:1:Nr-1;
delta_area = 2*pi*(ir+0.5)*(dr^2);          % dien tich hinh khuyen (delta a) [cm2]
temp = dx*dy*dz*Nphotons;

Rd = sum(sum(Rd_total))/Nphotons;           % do phan xa [-]
Tt = sum(sum(Tt_total))/Nphotons;           % do truyen qua [-]

Rd_yx = Rd_total/(dy*dx*Nphotons);          % xac suat photon phan xa khuech tan tren mot don vi dien tich xy [cm-2]
Tt_yx = Tt_total/(dy*dx*Nphotons);          % xac suat photon truyen qua tren mot don vi dien tich xy [cm-2]

Rd_r = Rd_r_total./(delta_area*Nphotons);   % xac suat photon phan xa khuech tan tren mot don vi dien tich theo phuong r [cm-2]
Tt_r = Tt_r_total./(delta_area*Nphotons);   % xac suat photon truyen qua tren mot don vi dien tich theo phuong r [cm-2]

% Convert data to relative fluence rate [cm^-2]
% Normalize deposition (A) to yield fluence rate (F)
A_z = zeros(1,Nz);
F_z = zeros(1,Nz);
for i=1:NN
    F(i)  = Q_total(i)/(temp*muav(v(i)));   % thong luong [cm^2]
    Ftemp(i) = Q_total(i)/muav(v(i));       % Ayxz[-]/mua[cm-1]
end 
A_yxz = reshape(Q_total,Ny,Nx,Nz);          % A(y,x,z) tong trong luong photon hap thu moi voxel [-]
Ftemp = reshape(Ftemp,Ny,Nx,Nz);            % Ftemp(y,x,z) [cm]
for i = 1:Nz
    A_z(i) = sum(sum(A_yxz(:,:,i)));        % tong trong luong photon hap thu trong moi phan tu luoi theo do sau z [-] 
    F_z(i) = sum(sum(Ftemp(:,:,i)));        % tong trong luong photon hap thu trong moi phan tu luoi theo do sau z [-] 
end
A = sum(A_z)/Nphotons;                      % xac suat hap thu trong mo [-]
A_z = A_z/(dz*Nphotons);                    % xac suat photon hap thu tren moi don vi chieu dai theo huong z [cm-1]
F_z = F_z/(dz*Nphotons);                    % do truyen qua theo do sau z [-]

%% SAVE RESULTS
% Save the binary file
fprintf('saving %s\n',strcat(filename,'_F.bin'));
fileID = fopen(strcat(filename,'_F.bin'), 'wb');      % 3D voxel output 
fwrite(fileID,F,'float');
fclose(fileID);

fprintf('saving %s\n',strcat(filename,'Rd_yx.bin'));
fileID = fopen(strcat(filename,'Rd_yx.bin'), 'wb');   % 2D output 
fwrite(fileID,reshape(Rd_yx,1,Ny*Nx),'float');
fclose(fileID);

fprintf('saving %s\n',strcat(filename,'Tt_yx.bin'));
fileID = fopen(strcat(filename,'Tt_yx.bin'), 'wb');   % 2D output 
fwrite(fileID,reshape(Tt_yx,1,Ny*Nx),'float');
fclose(fileID);

fprintf('saving %s\n',strcat(filename,'_folder'));
currentFolder = pwd;              % duong dan cua thu muc hien tai              
r = dr/2:dr:dr*(Nr-0.5);          % ban kinh r [cm]
z = dz/2:dz:dz*(Nz-0.5);          % chieu sau z [cm]
[~, ~, ~] = mkdir (filename);     % tao thu muc filename
cd(filename);                     % con tro den thu muc filename 

% ghi ket qua Rd, Tt, A
fileID = fopen(strcat(filename,'_RTA','.txt'),'w');
fprintf(fileID,'%0.4f\t %26s\n',Rd,'## Diffuse Reflectivity');
fprintf(fileID,'%0.4f\t %19s\n',Tt,'## Transmittance');
fprintf(fileID,'%0.4f\t %16s\n',A,'## Absorbance');
fclose(fileID);

% ghi ket qua Rr[1/cm2]
fileID = fopen(strcat(filename,'.Rr'),'w');
fprintf(fileID,'%s\t %13s\t\n','r[cm]','Rr[1/cm2]');
Data = [r;Rd_r];
fprintf(fileID,'%0.4E\t %0.4E\t\n',Data);
fclose(fileID);
% ghi ket qua Tr[1/cm2]
fileID = fopen(strcat(filename,'.Tr'),'w');
fprintf(fileID,'%s\t %13s\t\n','r[cm]','Tr[1/cm2]');
Data = [r;Tt_r];
fprintf(fileID,'%0.4E\t %0.4E\t\n',Data);
fclose(fileID);
% ghi ket qua Az[1/cm]
fileID = fopen(strcat(filename,'.Az'),'w');
fprintf(fileID,'%s\t %12s\t\n','z[cm]','Az[1/cm]');
Data = [z;A_z];
fprintf(fileID,'%0.4E\t %0.4E\t\n',Data);
fclose(fileID);
% ghi ket qua Fz[-]
fileID = fopen(strcat(filename,'.Fz'),'w');
fprintf(fileID,'%s\t %9s\t\n','z[cm]','Fz[-]');
Data = [z;F_z];
fprintf(fileID,'%0.4E\t %0.4E\t\n',Data);
fclose(fileID);

cd(currentFolder);                % thoat khoi thu muc filename

% save reflectance 
% NOT READY: 
% strcpy(filename,myname);
% strcat(filename,"_Ryx.bin");
% printf("saving %s\n",filename);
% fid = fopen(filename, "wb");   % 2D voxel output 
% int Nyx = Ny*Nx;
% fwrite(R, sizeof(float), Nyx, fid);
% fclose(fid);
% printf("%s is done.\n",myname);

fprintf('------------------------------------------------------\n');
now = datetime('now');
fprintf('%s\n', now);

%delete parameters.mat bintissue.mat    % delete parameters