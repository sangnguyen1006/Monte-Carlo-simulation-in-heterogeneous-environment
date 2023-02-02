function reportHmci(myname)
% function reportHmci(myname)
%   Lists the values of the input file  myname_H.mci.
%   Updated Feb 8, 2017. slj, adding boundaryflag B(10) (see s(10).s)

home
fid = fopen(sprintf('%s_H.mci',myname),'r');
B = fscanf(fid,'%f');
fclose(fid);

s(1).s = 'Nphotons';
s(2).s = 'Nx';
s(3).s = 'Ny';
s(4).s = 'Nz';
s(5).s = 'Nr';
s(6).s = 'dx';
s(7).s = 'dy';
s(8).s = 'dz';
s(9).s = 'dr';
s(10).s = 'mcflag';
s(11).s = 'launch';
s(12).s = 'boundary';
s(13).s = 'xs';
s(14).s = 'ys';
s(15).s = 'zs';
s(16).s = 'xfocus';
s(17).s = 'yfocus';
s(18).s = 'zfocus';
s(19).s = 'zref';
s(20).s = 'ztra';
s(21).s = 'ux0';
s(22).s = 'uy0';
s(23).s = 'uz0';
s(24).s = 'radius';
s(25).s = 'waist';
s(26).s = 'Nt';
for i=1:26
    disp(sprintf('%d\t%10s = %0f',i,s(i).s,B(i)))
end

for j=1:B(26)
    i=i+1;
    disp(sprintf('---'))
    disp(sprintf('%d\tmua = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tmus = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tg   = %0.4f',i,B(i)))
end
