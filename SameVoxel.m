function sv=SameVoxel(x1,y1,z1,x2,y2,z2,dx,dy,dz)
% Determine if the two position are located in the same voxel
% Returns 1 if same voxel, 0 if not same voxel

xmin = min(floor(x1/dx),floor(x2/dx))*dx;
ymin = min(floor(y1/dy),floor(y2/dy))*dy;
zmin = min(floor(z1/dz),floor(z2/dz))*dz;
xmax = xmin+dx;
ymax = ymin+dy;
zmax = zmin+dz;
sv=0;
if (x1<=xmax)&&(x2<=xmax)&&(y1<=ymax)&&(y2<=ymax)&&(z1<zmax)&&(z2<=zmax)
    sv=1;
end
end