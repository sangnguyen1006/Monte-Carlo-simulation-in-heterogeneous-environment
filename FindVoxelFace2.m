function s=FindVoxelFace2(x1,y1,z1,x2,y2,z2,dx,dy,dz,ux,uy,uz)
% my version of FindVoxelFace for no scattering
% s = ls + FindVoxelFace2(x,y,z, tempx, tempy, tempz, dx, dy, dz, ux, uy, uz)

ix1 = floor(x1/dx);
iy1 = floor(y1/dy);
iz1 = floor(z1/dz);

if (ux>=0)
    ix2 = ix1+1;
else
    ix2 = ix1;
end

if (uy>=0)
    iy2 = iy1+1;
else
    iy2 = iy1;
end
    
if (uz>=0)
    iz2 = iz1+1;
else
    iz2 = iz1;
end    
    
xs = abs((ix2*dx - x1)/ux);
ys = abs((iy2*dy - y1)/uy);
zs = abs((iz2*dz - z1)/uz);   

s = min([xs ys zs]);
end