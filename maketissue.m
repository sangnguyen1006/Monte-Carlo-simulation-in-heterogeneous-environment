% maketissue_18apr17.m
% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.m.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.m, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.m.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.m can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%       

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname      = 'AppleTissue_1';% name for files: myname_T.bin, myname_H.mci  
Nphotons    = 10000000; % number of photons used for simulation
nm          = 650;   	% desired wavelength of simulation (Wavelength from 600 nm to 1000 nm)
Nx          = 400;    	% # of bins in x dimension of cube 
Ny          = 400;    	% # of bins in y dimension of cube 
Nz          = 400;    	% # of bins in z dimension of cube 
Nr          = 400;    	% # of bins in r dimension of cube 
dx          = 0.00125;  % size of x bin, eg. [cm] or [mm]
dy          = 0.00125;  % size of y bin, eg. [cm] or [mm]
dz          = 0.00125;  % size of z bin, eg. [cm] or [mm]
dr          = 0.00125;  % size of r bin, eg. [cm] or [mm]

% Set Monte Carlo launch flags
mcflag      = 4;     	% launch: 0 = uniform beam, 1 = isotropic pt. 
                        % 2 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
                        % 3 = infinitely narrow beam (x = 0, y = 0)
                        % 4 = gaussian beam
launchflag  = 0;        % 0 = let mcxyz.m calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source
xs          = 0;      	% x of source [cm] or [mm]
ys          = 0;        % y of source [cm] or [mm]
zs          = 0.01;  	% z of source [cm] or [mm]
% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% Depth of the upper and lower of the tissue
% to get diffuse reflectance and transmission spectra
zref        = 0;        % Depth of plane where diffuse reflection is recorded [cm] or [mm],
                        % measured from the air/tissue interface
ztra        = 0.15;     % Depth of plane where transmission is recorded [cm] or [mm],
                        % measured from the air/tissue interface  
% only used if mcflag == 0 or 2 or 4 (not 1 = isotropic pt.)
radius      = 0.0500;   % 1/e radius of beam at tissue surface
waist       = 0.0500;   % 1/e radius of beam at focus
%radius      = 0.1500;   % 1/e radius of beam at tissue surface
%waist       = 0.1500;   % 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.7;      % trajectory projected onto x axis
uy0         = 0.4;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

% Create tissue properties
tissue = makeTissueList(nm); % also --> global tissue(1:Nt).s
Nt = length(tissue);
for i=1:Nt
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
end

% Specify Monte Carlo parameters    
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = double(zeros(Ny,Nx,Nz)); 
T = T + 1;                       % fill background with air

radius_apple = 3.89;             % cm
center_ellipsoid = 0.1;          % the depth of bruised tissue z [cm]
%center_ellipsoid = 0.03;          % the depth of bruised tissue z [cm]
a_ellipsoid = 0.1;               % [cm] a corresponds to x (-a, a)
b_ellipsoid = 0.03;              % [cm] b corresponds to y (-b, b)
c_ellipsoid = 0.03;              % [cm] c corresponds to z (-c, c)
zsurf = 0.01;                    % position of air/apple tissue [cm] 

xc    = (Nx*dx)/2;               % cm - circle center
yc    = (Ny*dy)/2;               % cm
zc    = radius_apple + zsurf;    % cm

xce    = (Nx*dx)/2;              % cm - center of ellipsoid
yce    = (Ny*dy)/2;              % cm
zce    = center_ellipsoid + zsurf; % cm

for iz = 1:Nz
    for iy = 1:Ny
        for ix = 1:Nx
            x_di = ix*dx;     % cm
            y_di = iy*dy;     % cm
            z_di = iz*dz;     % cm
            radius_temp = sqrt((xc - x_di)^2+(yc - y_di)^2+(zc - z_di)^2);
            temp1 = ((xce - x_di)/a_ellipsoid)^2;
            temp2 = ((yce - y_di)/b_ellipsoid)^2;
            temp3 = ((zce - z_di)/c_ellipsoid)^2;
            temp = temp1 + temp2 + temp3;
            if (radius_temp <= radius_apple)
                if (iy <= Ny) && (ix <= Nx) && (iz <= Nz)
                    T(iy,ix,iz) = 2; % Normal Tissue
                end
            end
            if temp <= 1
                if (iy <= Ny) && (ix <= Nx) && (iz <= Nz)
                    T(iy,ix,iz) = 3; % Bruised Tissue
                end
            end
        end
    end
end

%%
if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.0f\n',Nphotons);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%d\n'   ,Nr);
        fprintf(fid,'%0f\n',dx);
        fprintf(fid,'%0f\n',dy);
        fprintf(fid,'%0f\n',dz);
        fprintf(fid,'%0f\n',dr);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,launchflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0f\n',xs);
        fprintf(fid,'%0f\n',ys);
        fprintf(fid,'%0f\n',zs);
        fprintf(fid,'%0f\n',xfocus);
        fprintf(fid,'%0f\n',yfocus);
        fprintf(fid,'%0f\n',zfocus);
        fprintf(fid,'%0f\n',zref + zsurf);
        fprintf(fid,'%0f\n',ztra + zsurf);
        fprintf(fid,'%0f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0f\n',uy0);
        fprintf(fid,'%0f\n',uz0);
        fprintf(fid,'%0f\n',radius);
        fprintf(fid,'%0f\n',waist);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%

for i=1:Nt
    yy = (Nt-i)/(Nt-1)*Nz*dz;
    text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
end

text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

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
plot([min(x) max(x)],([zref zref] + zsurf),'r-','LineWidth',1)
plot([min(x) max(x)],([ztra ztra] + zsurf),'r-','LineWidth',1)

disp('done')

