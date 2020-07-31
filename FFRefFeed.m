function [EFx, EFy, EFz] = FFRefFeed(freq, er, jf, a, r, th, phi)
    %Calculating required parameters from given
    %Speed of light
    c = 3e8; 
    
    %Wavelength (in m)
    lam = c/freq;
    k0 = 2*pi/lam;
    
    %Impedance (in Ohm)
    eps_0 = 8.854187817e-12;
    mu_0 = 1.2566370614e-6;
    zeta = round(sqrt(mu_0/(eps_0*er)));
    
    %Defining spectral numbers
    kx = k0.*sin(th).*cos(phi);
    ky = k0.*sin(th).*sin(phi);
    kz = -1j*sqrt(-(k0^2-kx.^2-ky.^2));
    %kz = k0.*cos(th);
    %kz = sqrt(k0^2-kx.^2-ky.^2);
    
    %Calculate JFT
    J = JFT(k0, a, th, jf);

    %Calling SGF
    FFSGF = createSGF(k0, kx, ky, zeta, th);
    
    %Finding Electric field; Assuming obs pi in far field
    EFx = 2j.*kz.*(squeeze(FFSGF(1,1,:,:)).*squeeze(J(1,:,:)) ...
        + squeeze(FFSGF(1,2,:,:)).*squeeze(J(2,:,:)) ...
        + squeeze(FFSGF(1,3,:,:)).*squeeze(J(3,:,:))).*(exp(-1j*k0*r)./(2*pi*r));
    
    EFy = 2j.*kz.*(squeeze(FFSGF(2,1,:,:)).*squeeze(J(1,:,:)) ...
        + squeeze(FFSGF(2,2,:,:)).*squeeze(J(2,:,:)) ...
        + squeeze(FFSGF(2,3,:,:)).*squeeze(J(3,:,:))).*(exp(-1j*k0*r)./(2*pi*r));
    
    EFz = 2j.*kz.*(squeeze(FFSGF(3,1,:,:)).*squeeze(J(1,:,:)) ...
        + squeeze(FFSGF(3,2,:,:)).*squeeze(J(2,:,:)) ...
        + squeeze(FFSGF(3,3,:,:)).*squeeze(J(3,:,:))).*(exp(-1j*k0*r)./(2*pi*r));
    
    %FFR = EFx.*sin(th).*cos(phi) + EFy.*sin(th).*sin(phi) - EFz.*cos(th);
    %FFTheta = EFx.*cos(th).*cos(phi) + EFy.*cos(th).*sin(phi) - EFz.*sin(th);
    %FFTheta = EFy.*cos(th) - EFz.*sin(th); 
    %FFPhi = -EFx.*sin(phi) + EFy.*cos(phi);
end