function tissue = makeTissueList(nm)
% function tissueProps = makeTissueList(nm)
%   Wavelength from 600 nm to 1000 nm
%
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
%       global tissuenames(i).s
%
%   Uses 
%       OpticalParameters.mat

%% Load spectral library
load OpticalParameters.mat

% wavelength        401 1 double
% mua_NormalReDeli  401 1 double
% mus_NormalReDeli  401 1 double
% mua_BruisedReDeli 401 1 double
% mus_BruisedReDeli 401 1 double

MU(:,1) = interp1(wavelength, mua_NormalReDeli, nm);
MU(:,2) = interp1(wavelength, mus_NormalReDeli, nm);

MU(:,3) = interp1(wavelength, mua_BruisedReDeli, nm);
MU(:,4) = interp1(wavelength, mus_BruisedReDeli, nm);

LOADED = 1;

%% Create tissueList
j=1;
tissue(j).name  = 'air';
tissue(j).mua   = 0.0001; % Negligible absorption yet still tracks, 
tissue(j).mus   = 1.0;    % but take steps in air
tissue(j).g     = 1.0;    % and don't scatter.

j=2;
tissue(j).name  = 'Normal Tissue';
tissue(j).mua   = MU(1);
tissue(j).mus   = MU(2);    
tissue(j).g     = 0.9547;   

j=3;
tissue(j).name  = 'Bruised Tissue';
tissue(j).mua   = MU(3);
tissue(j).mus   = MU(4);   
tissue(j).g     = 0.9547;

disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(tissue)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,...
        tissue(i).mus*(1-tissue(i).g)))
end
disp(' ')

