% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz.c 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   Makes figures:
%       myname_tissue.jpg   = tissue structure (shows tissue types)
%       myname_Fzx.jpg      = fluence rate vs z,x
%       myname_Fzy.jpg      = fluence rate vs z,y
%   Uses:
%       myname_H.mci    = input file from maketissue.m
%       myname_T.bin    = tissue input file from maketissue.m
%       myname_F.bin    = fluence rate output from Monte Carlo
%       reportH_mci.m   = lists input parameters in myname_H.mci
%       makecmap.m      = makes colormap for tissue types
%       makec2f.m       = makes colormap for fluence rate
%
%   This example sets myname = 'skinvessel'.
%
% 7/feb/2017, add boundaryflag (see A(10)).
% 1/june/2017 , no major changes, just clean up display outputs.
% Steven L Jacques
home; clear
format compact
commandwindow

SAVEPICSON = 1;
if SAVEPICSON
    sz = 10; fz = 7; fz2 = 5; % to use savepic.m
else
    sz = 12; fz = 9; fz2 = 7; % for screen display
end

%%%% USER CHOICES <---------- you must specify -----
myname = 'AppleTissue_1'; nm = 650; 
%%%%


disp(sprintf('------ mcxyz %s -------',myname))

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
Nphotons = A(1);
Nx = A(2);
Ny = A(3);
Nz = A(4);
Nr = A(5);
dx = A(6);
dy = A(7);
dz = A(8);
dr = A(9);
mcflag = A(10);
launchflag = A(11);
boundaryflag = A(12);
xs = A(13);
ys = A(14);
zs = A(15);
xfocus = A(16);
yfocus = A(17);
zfocus = A(18);
zref = A(19);
ztra = A(20);
ux0 = A(21);
uy0 = A(22);
uz0 = A(23);
radius = A(24);
waist = A(25);
Nt = A(26);
j = 26;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
end

reportHmci(myname)

%% Load Fluence rate F(y,x,z) 
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

%% Load tissue structure in voxels, T(y,x,z) 
filename = sprintf('%s_T.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)

%% Load Diffuse Reflection Rd_yx [cm^-2]
filename = sprintf('%sRd_yx.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx, 'float');
    fclose(fid);
toc
Rd_yx = reshape(Data,Ny,Nx); % in plane yx perpendicular to the light source

%% Load Transmission Tt_yx [cm^-2]
filename = sprintf('%sTt_yx.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx, 'float');
    fclose(fid);
toc
Tt_yx = reshape(Data,Ny,Nx); % in plane yx perpendicular to the light source

clear Data

%%
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

%% Look at structure, Tzx
Tzx = reshape(T(Ny/2,:,:),Nx,Nz)';
tissue = makeTissueList(nm);
Nt = length(tissue);
figure(1);clf
imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue','fontweight','normal','fontsize',fz2)
for i=1:Nt
    yy = zmin + (Nt-i)/(Nt-1)*zdiff;
    text(xmax*1.4,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',fz2)
end

%%% draw launch
N = 20; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'b-')
        end

%     case 1 % Gaussian
%         for i=0:N
%             plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'b-')
%         end

    case 1 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'b-')
        end
        
    case 2 % rectangle
        zz = max(z);
        for i=1:N
            xx = -radius + 2*radius*i/20;
            plot([xx xx],[zs zz],'b-')
        end
        
    case 3 % infinitely narrow beam
        plot([0 0],[zs max(z)],'b-')  
    
    case 4 % Gaussian
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'b-')
        end
end

%%% draw the diffuse and transmitted reflection acquisition plane
plot([min(x) max(x)],[zref zref],'r-','LineWidth',0.75)
plot([min(x) max(x)],[ztra ztra],'r-','LineWidth',0.75)

axis equal image

if SAVEPICSON
    name = sprintf('%s_tissueZX.jpg',myname);
    savepic(1,[4 3],name)
end

%% Look at structure, Tzy
Tzy = reshape(T(:,Nx/2,:),Ny,Nz)';
tissue = makeTissueList(nm);
Nt = length(tissue);
figure(2);clf
imagesc(y(uy),z(uz),Tzy(uz,uy),[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('y [cm]')
ylabel('z [cm]')
title('Tissue','fontweight','normal','fontsize',fz2)
for i=1:Nt
    yy = zmin + (Nt-i)/(Nt-1)*zdiff;
    text(xmax*1.4,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',fz2)
end

%%% draw launch
N = 20; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'b-')
        end

%     case 1 % Gaussian
%         for i=0:N
%             plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'b-')
%         end

    case 1 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'b-')
        end
        
    case 2 % rectangle
        zz = max(z);
        for i=1:N
            xx = -radius + 2*radius*i/20;
            plot([xx xx],[zs zz],'b-')
        end
        
    case 3 % infinitely narrow beam
        plot([0 0],[zs max(z)],'b-')  
    
    case 4 % Gaussian
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'b-')
        end
end

%%% draw the diffuse and transmitted reflection acquisition plane
plot([min(x) max(x)],[zref zref],'r-','LineWidth',0.75)
plot([min(x) max(x)],[ztra ztra],'r-','LineWidth',0.75)

axis equal image

if SAVEPICSON
    name = sprintf('%s_tissueZY.jpg',myname);
    savepic(2,[4 3],name)
end

%% Look at Diffuse Reflection Rd_yx [cm^-2]
figure(3);clf
imagesc(y,x,log10(Rd_yx),[-1 2])
hold on
text(max(x)*1.2,min(y)-0.1*max(y),'log_{10}( Rd_y_x )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('y [cm]')
title('Diffuse Reflection [1/cm^2] ','fontweight','normal','fontsize',fz)
colormap(gray)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(y)-0.15*max(y),sprintf('The number of photons = %0.1f',Nphotons),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%sRd_yx.jpg',myname);
    savepic(3,[4 3],name)
end

%% Look at Transmission Tt_yx [cm^-2]
figure(4);clf
imagesc(y,x,log10(Tt_yx),[-1 2])
hold on
text(max(x)*1.2,min(y)-0.1*max(y),'log_{10}( Tt_y_x )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('y [cm]')
title('Transmission [1/cm^2] ','fontweight','normal','fontsize',fz)
colormap(gray)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(y)-0.15*max(y),sprintf('The number of photons = %0.1f',Nphotons),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%sTt_yx.jpg',myname);
    savepic(4,[4 3],name)
end

%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source

figure(5);clf
imagesc(x,z,log10(Fzx),[-1 2.8])
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(gray)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('The number of photons = %0.1f',Nphotons),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzx.jpg',myname);
    savepic(5,[4 3],name)
end

%% look Fzy
Fzy = reshape(F(:,Nx/2,:),Ny,Nz)';

iy = round((dy*Ny/2 + 0.15)/dy);
iz = round(zs/dz);
zzs  = zs;
%Fdet = mean(reshape(Fzy(iz+[-1:1],iy+[0 1]),6,1));

figure(6);clf
imagesc(y,z,log10(Fzy),[-1 2.8])
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('y [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(gray)
axis equal image
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('The number of photons = %0.1f',Nphotons),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzy.jpg',myname);
    savepic(6,[4 3],name)
end

%% look Azx
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
mua = muav(reshape(T(Ny/2,:,:),Nx,Nz)');
Azx = Fzx.*mua;

figure(7);clf
imagesc(x,z,log10(Azx),[-1 2.8])
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( A )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Deposition A [W/cm^3/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(gray)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('The number of photons = %0.1f',Nphotons),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Azx.jpg',myname);
    savepic(7,[4 3],name)
end

drawnow

disp('done')


