function r=RFresnel(n1,n2,ca1) 	  
% RFresnel(n1,n2,ca1) 
% n1: incident refractive index
% n2: transmit refractive index
% ca1=cos(radian): cosine of the incident (angle a1, 0<a1<90 degrees)
% ca2_Ptr: pointer to the cosine of the transmission angle a2, a2>0

% FRESNEL REFLECTANCE
% Computes reflectance as photon passes from medium 1 to 
% medium 2 with refractive indices n1,n2. Incident
% angle a1 is specified by cosine value ca1 = cos(a1).
% Program returns value of transmitted angle a1 as
% value in *ca2_Ptr = cos(a2).

if n1 == n2                            % matched boundary
    ca2_Ptr = ca1;
    r = 0.0;
else
    if ca1>(1.0 - 1.0e-12)             % normal incidence
        ca2_Ptr = ca1;
        r = (n2-n1)/(n2+n1);
        r = r*r;
    else
        if ca1< 1.0e-6                 % very slanted
            ca2_Ptr = 0.0;
            r = 1.0;
        else  			               % general
            sa1 = sqrt(1-ca1*ca1);     % sine of incident and transmission angles
            sa2 = n1*sa1/n2;           % sine of incident and transmission angles
            if sa2>=1.0                % double check for total internal reflection
                ca2_Ptr = 0.0;
                r = 1.0;
            else 
                ca2 = sqrt(1-sa2*sa2); % cosine of transmission angle
                ca2_Ptr = ca2;
                % cosines of sum ap or diff am of the two
                % angles: ap = a1 + a2, am = a1 - a2
                cap = ca1*ca2 - sa1*sa2; 
                cam = ca1*ca2 + sa1*sa2; 
                % sines
                sap = sa1*ca2 + ca1*sa2; 
                sam = sa1*ca2 - ca1*sa2; 
                r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); % rearranged for speed
            end
        end
    end
end
end


