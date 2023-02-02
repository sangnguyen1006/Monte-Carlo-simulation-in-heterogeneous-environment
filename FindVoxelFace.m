function s=FindVoxelFace(x1,y1,z1,x2,y2,z2,dx,dy,dz,ux,uy,uz)
% the boundary is the face of some voxel
% find the the photon's hitting position on the nearest face of the voxel and update the step size
x_1 = x1/dx;
y_1 = y1/dy;
z_1 = z1/dz;
x_2 = x2/dx;
y_2 = y2/dy;
z_2 = z2/dz;
fx_1 = floor(x_1);
fy_1 = floor(y_1);
fz_1 = floor(z_1);
fx_2 = floor(x_2);
fy_2 = floor(y_2);
fz_2 = floor(z_2);
x = 0; y = 0; z = 0; x0 = 0; y0 = 0; z0 = 0; s = 0;

if (fx_1 ~= fx_2) && (fy_1 ~= fy_2) && (fz_1 ~= fz_2) % #1
    fx_2=fx_1+SIGN(fx_2-fx_1); % added
    fy_2=fy_1+SIGN(fy_2-fy_1); % added
    fz_2=fz_1+SIGN(fz_2-fz_1); % added
        
    x = (max(fx_1,fx_2)-x_1)/ux;
    y = (max(fy_1,fy_2)-y_1)/uy;
    z = (max(fz_1,fz_2)-z_1)/uz;
    if (x == min([x y z])) 
        x0 = max(fx_1,fx_2);
        y0 = (x0-x_1)/ux*uy+y_1;
        z0 = (x0-x_1)/ux*uz+z_1;
    else
        if (y == min([x y z]))
            y0 = max(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        else
            z0 = max(fz_1,fz_2);
            y0 = (z0-z_1)/uz*uy+y_1;
            x0 = (z0-z_1)/uz*ux+x_1;
        end
    end
else
    if (fx_1 ~= fx_2) && (fy_1 ~= fy_2) % #2
        fx_2=fx_1+SIGN(fx_2-fx_1); % added
        fy_2=fy_1+SIGN(fy_2-fy_1); % added
        x = (max(fx_1,fx_2)-x_1)/ux;
        y = (max(fy_1,fy_2)-y_1)/uy;
        if (x == min(x,y)) 
            x0 = max(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        else 
            y0 = max(fy_1, fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        end
    else
        if (fy_1 ~= fy_2) &&(fz_1 ~= fz_2) % #3
            fy_2 = fy_1+SIGN(fy_2-fy_1); % added
            fz_2 = fz_1+SIGN(fz_2-fz_1); % added
            y = (max(fy_1,fy_2)-y_1)/uy;
            z = (max(fz_1,fz_2)-z_1)/uz;
            if (y == min(y,z)) 
                y0 = max(fy_1,fy_2);
                x0 = (y0-y_1)/uy*ux+x_1;
                z0 = (y0-y_1)/uy*uz+z_1;
            else 
                z0 = max(fz_1, fz_2);
                x0 = (z0-z_1)/uz*ux+x_1;
                y0 = (z0-z_1)/uz*uy+y_1;
            end
        else
            if (fx_1 ~= fx_2) && (fz_1 ~= fz_2) % #4
                fx_2=fx_1+SIGN(fx_2-fx_1); % added
                fz_2=fz_1+SIGN(fz_2-fz_1); % added
                x = (max(fx_1,fx_2)-x_1)/ux;
                z = (max(fz_1,fz_2)-z_1)/uz;
                if (x == min(x,z)) 
                    x0 = max(fx_1,fx_2);
                    y0 = (x0-x_1)/ux*uy+y_1;
                    z0 = (x0-x_1)/ux*uz+z_1;
                else 
                    z0 = max2(fz_1, fz_2);
                    x0 = (z0-z_1)/uz*ux+x_1;
                    y0 = (z0-z_1)/uz*uy+y_1;
                end
            else
                if (fx_1 ~= fx_2) % #5
                    fx_2=fx_1+SIGN(fx_2-fx_1); % added
                    x0 = max(fx_1,fx_2);
                    y0 = (x0-x_1)/ux*uy+y_1;
                    z0 = (x0-x_1)/ux*uz+z_1;
                else
                    if (fy_1 ~= fy_2) % #6
                        fy_2=fy_1+SIGN(fy_2-fy_1); % added
                        y0 = max(fy_1, fy_2);
                        x0 = (y0-y_1)/uy*ux+x_1;
                        z0 = (y0-y_1)/uy*uz+z_1;
                    else % #7 
                        z0 = max(fz_1, fz_2);
                        fz_2=fz_1+SIGN(fz_2-fz_1); % added
                        x0 = (z0-z_1)/uz*ux+x_1;
                        y0 = (z0-z_1)/uz*uy+y_1;
                    end
                end
            end
        end
    end
end
    % s = ( (x0-fx_1)*dx + (y0-fy_1)*dy + (z0-fz_1)*dz )/3;
    % s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz) );
    % s = sqrt(SQR(x0-x_1)*SQR(dx) + SQR(y0-y_1)*SQR(dy) + SQR(z0-z_1)*SQR(dz));
    s = sqrt(((x0-x_1)*dx)^2 + ((y0-y_1)*dy)^2 + ((z0-z_1)*dz)^2);
end