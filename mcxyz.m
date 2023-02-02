function [Q, Rd, Tt, Rd_r, Tt_r] = mcxyz(Nphotons_worker)
% MAJOR CYCLE
% Get parameters from parameters.mat and bintissue.mat, 
% perform simulation and return results:
%       (1) Q: total weight of absorbed photons voxel [-]
%       (2) Rd: diffuse reflection weight [-]
%       (3) Tt: transmitted weight [-]
%       (4) Rd_r: diffuse reflectance weight with radius [-]
%       (5) Tt_r: transmitted weight with radius [-]

% DEFINE
% ALIVE     = 1   		     % if photon not yet terminated
% DEAD      = 0    		     % if photon is to be terminated 
ls          = 1.0E-7;        % Moving photon a little bit off the voxel face 
THRESHOLD   = 0.01;		     % used in roulette 
CHANCE      = 0.1;  		 % used in roulette 
ONE_MINUS_COSZERO = 1.0E-12; % If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad
                  
% LOAD INPUT PARAMETERS IN "parameters.mat" AND "btissue.mat"
% run parameters
% launch parameters
% tissue optical properties
load parameters.mat Nx Ny Nz Nr dx dy dz dr mcflag launchflag boundaryflag...        
    xs ys zs xfocus yfocus zfocus zref ztra ux0 uy0 uz0 radius waist muav musv gv
% tissue structure
load bintissue.mat v

% INITIALIZATIONS
% relative fluence rate [W/cm^2/W.delivered]
% NOT READY: R  = (float *)malloc(NN*sizeof(float)); escaping flux [W/cm^2/W.delivered] 
NN   = Nx*Ny*Nz;
Q    = zeros(1,NN);    % ensure Q[] starts empty [-]
Rd   = zeros(Ny,Nx);   % diffuse reflection weight [-]
Tt   = zeros(Ny,Nx);   % transmitted weight [-]
Rd_r = zeros(1,Nr);    % diffuse reflectance weight with radius [-]
Tt_r = zeros(1,Nr);    % transmitted weight with radius [-]

% RUN
% Launch N photons, initializing each one before progation.
i_photon = 0;
while i_photon < Nphotons_worker
    % LAUNCH 
	% Initialize photon position and trajectory
    inside = 1;                 % photon inside tissue = 1, diffuse reflection or tissue transmission = 0
    i_photon = i_photon + 1;	% increment photon count 
	W = 1.0;                    % set photon weight to one 
	photon_status = 1;          % Launch an ALIVE photon 
    
    % SET SOURCE 
    % Launch collimated beam at x,y center
    % Initial position	
	% trajectory
    if launchflag == 1   % manually set launch
		x	= xs; 
		y	= ys;
		z	= zs;
		ux	= ux0;
		uy	= uy0;
		uz	= uz0;
    else                 % use mcflag
        if mcflag == 0   % uniform beam (chum tia dong nhat)
			% set launch point and width of beam
			r		= radius*sqrt(rand);  % radius of beam at launch point
			phi		= rand*2.0*pi;
			x		= xs + r*cos(phi);
			y		= ys + r*sin(phi);
			z		= zs;
			% set trajectory toward focus
			r		= waist*sqrt(rand);   % radius of beam at focus
			phi		= rand*2.0*pi;
			xfocus	= r*cos(phi);
			yfocus	= r*sin(phi);
			temp	= sqrt((x - xfocus)*(x - xfocus) + (y - yfocus)*(y - yfocus) + zfocus*zfocus);
			ux		= -(x - xfocus)/temp;
			uy		= -(y - yfocus)/temp;
			uz		= sqrt(1 - ux*ux - uy*uy);
        else
            if mcflag == 1                % isotropic pt source (nguon dang huong)
                costheta = 1.0 - 2.0*rand;
				sintheta = sqrt(1.0 - costheta*costheta);
				psi = 2.0*pi*rand;
				cospsi = cos(psi);
                if (psi < pi)
					sinpsi = sqrt(1.0 - cospsi*cospsi); 
				else
					sinpsi = -sqrt(1.0 - cospsi*cospsi);
                end
                x = xs;
				y = ys;
				z = zs;
				ux = sintheta*cospsi;
				uy = sintheta*sinpsi;
				uz = costheta;
            else
                if mcflag == 2 % rectangular source collimated (nguon chuan truc hinh chu nhat)
				    x = radius*(rand*2-1); % use radius to specify x-halfwidth of rectangle
				    y = radius*(rand*2-1); % use radius to specify y-halfwidth of rectangle
				    z = zs;
				    ux = 0.0;
				    uy = 0.0;
				    uz = 1.0; % collimated beam
                else
                    if mcflag == 3 % infinitely narrow beam
                        x  = 0.0;
				        y  = 0.0;
				        z  = zs;
				        ux = 0.0;
				        uy = 0.0;
				        uz = 1.0; % collimated beam
                    else
                        if  mcflag == 4 % gaussian beam
                            r     = radius*sqrt(-log(1-rand));
                            alpha = 2*pi*rand;
                            x     = cos(alpha)*r;
				            y     = sin(alpha)*r;
				            z     = zs;
				            ux    = 0.0;
				            uy    = 0.0;
				            uz    = 1.0; % collimated beam
                        end
                    end
                end
            end
        end
    end  % end  use mcflag
    
    % Get tissue voxel properties of launchpoint.
	% If photon beyond outer edge of defined voxels, 
	% the tissue equals properties of outermost voxels
	% Therefore, set outermost voxels to infinite background value
    ix = round(Nx/2 + x/dx);
	iy = round(Ny/2 + y/dy);
	iz = round(z/dz);        
	if (ix > Nx)  ix = Nx; end
	if (iy > Ny)  iy = Ny; end
	if (iz > Nz)  iz = Nz; end
	if (ix < 1)    ix = 1; end 
	if (iy < 1)    iy = 1; end
	if (iz < 1)    iz = 1; end
    % Get the tissue type of located voxel 
    i		= round((iz-1)*Ny*Nx + (ix-1)*Ny + iy);
	type	= v(i);
	mua 	= muav(type);
	mus 	= musv(type);
	g 		= gv(type);
    
    bflag = 1; % initialize as 1 = inside volume, but later check as photon propagates
    
    % HOP_DROP_SPIN_CHECK
    % Propagate one photon until it dies as determined by ROULETTE.
    
    while photon_status == 1  % Launch an ALIVE photon 
        % if ALIVE, continue propagating 
		% If photon DEAD, then launch new photon
        
        % HOP
		% Take step to new position
		% s = dimensionless stepsize
		% ux, uy, uz are cosines of current photon trajectory
		sleft	= -log(rand);		% dimensionless step 
        
        while sleft>0
            s     = sleft/mus;		% Step size [cm]
			tempx = x + s*ux;		% Update positions. [cm] 
			tempy = y + s*uy;	
			tempz = z + s*uz;
            
            sv = SameVoxel(x,y,z, tempx, tempy, tempz, dx,dy,dz);
            
            if sv == 1 % photon in same voxel   
			    x = tempx;			% Update positions
				y = tempy;
				z = tempz;
                
				% DROP
				% Drop photon weight (W) into local bin
                absorb = W*(1 - exp(-mua*s));	% photon weight absorbed at this step 
                W = W - absorb;					% decrement WEIGHT by amount absorbed 
				% If photon within volume of heterogeneity, deposit energy in F[] 
				% Normalize F[] later, when save output. 
                if bflag == 1 
                    Q(i) = Q(i) + absorb;	    % only save data if blag==1, i.e., photon inside simulation cube
                end
                % Update sleft 
				sleft = 0;		                % dimensionless step remaining 
            else % photon has crossed voxel boundary 
				% step to voxel face + "littlest step" so just inside new voxel
				s = ls + FindVoxelFace2(x,y,z, tempx,tempy,tempz, dx,dy,dz, ux,uy,uz);
					 
				% DROP
				% Drop photon weight (W) into local bin

				absorb = W*(1-exp(-mua*s));     % photon weight absorbed at this step 
				W = W - absorb;                 % decrement WEIGHT by amount absorbed 
				% If photon within volume of heterogeneity, deposit energy in F[] 
				% Normalize F[] later, when save output
                if bflag == 1 
                    Q(i) =  Q(i) + absorb;	
                end
				% Update sleft 
				sleft = sleft - s*mus;          % dimensionless step remaining 
				if sleft <= ls 
                    sleft = 0;
                end
				
				% Update positions
				x = x + s*ux;
				y = y + s*uy;
				z = z + s*uz;
	
				% pointers to voxel containing optical properties
                ix = round(Nx/2 + x/dx);
                iy = round(Ny/2 + y/dy);
                iz = round(z/dz);
                
                % Record diffuse reflectance weight with radius r
                if z < zref
                    if inside == 1
                        r = sqrt(x^2+y^2); % radius r [cm]
                        ir = round(r/dr+0.5); 
                        if ir <= Nr 
                            Rd_r(ir) = Rd_r(ir) + W; 
                        end
                        if (iy >= 1) && (iy <= Ny) && (ix >= 1) && (ix <= Nx)
                            Rd(iy,ix) = Rd(iy,ix) + W;
                        end
                        inside = 0;
                    end
                end
                % Record the transmitted weight in terms of radius r
                if z > ztra
                    if inside == 1
                        r = sqrt(x^2+y^2); % radius r [cm]
                        ir = round(r/dr+0.5); 
                        if ir <= Nr 
                            Tt_r(ir) = Tt_r(ir) + W; 
                        end
                        if (iy >= 1) && (iy <= Ny) && (ix >= 1) && (ix <= Nx)
                            Tt(iy,ix) = Tt(iy,ix) + W;
                        end
                        inside = 0;
                    end
                end
				
                bflag = 1;  % Boundary flag. Initialize as 1 = inside volume, then check
                if (boundaryflag ==0 ) % Infinite medium
			        % Check if photon has wandered outside volume
                    % If so, set tissue type to boundary value, but let photon wander
                    % Set blag to zero, so DROP does not deposit energy
					if (iz > Nz) iz = Nz; bflag = 0; end
					if (ix > Nx) ix = Nx; bflag = 0; end
					if (iy > Ny) iy = Ny; bflag = 0; end
					if (iz < 1)   iz = 1; bflag = 0; end
					if (ix < 1)   ix = 1; bflag = 0; end
					if (iy < 1)   iy = 1; bflag = 0; end
                else
                    if (boundaryflag == 1) % Escape at boundaries
						if (iz > Nz) iz = Nz; photon_status = 0; sleft = 0; end
						if (ix > Nx) ix = Nx; photon_status = 0; sleft = 0; end
						if (iy > Ny) iy = Ny; photon_status = 0; sleft = 0; end
						if (iz < 1)   iz = 1; photon_status = 0; sleft = 0; end
						if (ix < 1)   ix = 1; photon_status = 0; sleft = 0; end
						if (iy < 1)   iy = 1; photon_status = 0; sleft = 0; end
                    else
                        if (boundaryflag == 2) % Escape at top surface, no x,y bottom z boundaries
						    if (iz > Nz) iz = Nz; bflag = 0; end
						    if (ix > Nx) ix = Nx; bflag = 0; end
						    if (iy > Ny) iy = Ny; bflag = 0; end
						    if (iz < 1)   iz = 1; photon_status = 0; sleft=0; end
						    if (ix < 1)   ix = 1; bflag = 0; end
						    if (iy < 1)   iy = 1; bflag = 0; end
                        end
                    end
                end
                
                % update pointer to tissue type
				i    = round((iz-1)*Ny*Nx + (ix-1)*Ny + iy);
				type = v(i);
                mua  = muav(type);
                mus  = musv(type);
                g    = gv(type);
            end
        end
        
        % SPIN 
		% Scatter photon into new trajectory defined by theta and psi
		% Theta is specified by cos(theta), which is determined 
		% based on the Henyey-Greenstein scattering function
		% Convert theta and psi into cosines ux, uy, uz
        
		% Sample for costheta
		if g == 0.0
		    costheta = 2.0*rand - 1.0;
		else 
		    temp = (1.0 - g*g)/(1.0 - g + 2*g*rand);
		    costheta = (1.0 + g*g - temp*temp)/(2.0*g);
        end
		sintheta = sqrt(1.0 - costheta*costheta); % sqrt() is faster than sin()
        
        % Sample psi
		psi = 2.0*pi*rand;
		cospsi = cos(psi);
		if (psi < pi)
		    sinpsi = sqrt(1.0 - cospsi*cospsi);   % sqrt() is faster than sin()
		else
			sinpsi = -sqrt(1.0 - cospsi*cospsi);
        end
        
        % New trajectory
		if (1 - abs(uz) <= ONE_MINUS_COSZERO)      % close to perpendicular
		    uxx = sintheta * cospsi;
			uyy = sintheta * sinpsi;
			uzz = costheta * SIGN(uz);             % SIGN() is faster than division
        else  % usually use this option 
			temp = sqrt(1.0 - uz * uz);
			uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
			uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
			uzz = -sintheta * cospsi * temp + uz * costheta;
        end
        
        % Update trajectory 
		ux = uxx;
		uy = uyy;
		uz = uzz;
        
        % CHECK ROULETTE 
		% If photon weight below THRESHOLD, then terminate photon using Roulette technique
		% Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
		% and 1-CHANCE probability of terminating
        if W < THRESHOLD 
		    if rand <= CHANCE
			    W = W / CHANCE;
            else
                photon_status = 0;
            end
        end
			
    end % end STEP_CHECK_HOP_SPIN
                               
end % end RUN 

end