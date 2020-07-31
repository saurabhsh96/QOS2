%Spillover efficiency

function [taperS, etaS, Area] = Spillover(freq, er, jf, af, r, th_dash, phi_dash, ...
    th, phi, rho, fRef, aRef)
    %% EtaS calculations
    %Char impedance of free space
    zeta = 377;
    k0 = 2*pi*freq/(3e8);
    
    %Field at complete angle pi/2
    [EFxF, EFyF, EFzF] = FFRefFeed(freq, er, jf, af, r, th, phi);
    EmagF = sqrt(EFxF.^2 + EFyF.^2 + EFzF.^2);
    
    %Field at theta0
    [EFx0, EFy0, EFz0] = FFRefFeed(freq, er, jf, af, r, th_dash, phi_dash);
    Emag0 = sqrt(EFx0.^2 + EFy0.^2 + EFz0.^2);
    
    %Constant term to negate the phase effect
    %Why has this effect taking place? Phase negation gives wrong answer
    Constant_term = (r)./(exp(-1j * k0 * r));
    
%     Constant_term = 1;
    % Converting theta0 components in spherical
    EF0r = Constant_term*(EFx0.*sin(th_dash).*cos(phi_dash)) + (EFy0.*sin(th_dash).*sin(phi_dash)) + (EFz0.*cos(th_dash));
    EF0th = Constant_term*(EFx0.*cos(th_dash).*cos(phi_dash)) + (EFy0.*cos(th_dash).*sin(phi_dash)) - (EFz0.*sin(th_dash));
    EF0phi = Constant_term*(-EFx0.*sin(phi_dash)) + (EFy0.*cos(phi_dash));

    %Calculate field intensity
    IntensityF = r^2.*abs(EmagF).^2./(2*zeta);
    Intensity0 = r^2.*abs(Emag0).^2./(2*zeta);
    
    %Calculate prad or integral involved in etaS
    %Infinitsimal change in th and phi i.e it is dth and dphi
    dth_dash = th_dash(1, 2) - th_dash(1, 1);
    dph_dash = phi_dash(2, 1) - phi_dash(1, 1);
    Prad0 = (sum(Intensity0.*sin(th_dash), 'all'))*dth_dash*dph_dash;
    
    dth = th(1, 2) - th(1, 1);
    dph = phi(2, 1) - phi(1, 1);
    PradF = (sum(IntensityF.*sin(th), 'all'))*dth*dph;
    
    etaS = Prad0/PradF;
    
    %% Taper efficiency calculations
    dRho = rho(1, 2) - rho(1, 1);
    %Finding the total area of the ap 
    %Or this area can also be just pi.r^2?
    Area = sum(rho, 'all')*dRho*dph_dash;
    %Area = pi*aRef^2;
    %Finding Aeff
    %Phase and Amplitude term of the field at the aperture
    %Cpa = exp(-2j*k0*fRef).*(cos(th_dash/2)).^2./(fRef);
    Cpa = (cos(th_dash/2)).^2./(fRef);
    %Cpa = 1;
    ErRho = -EF0th.*Cpa;
    ErPhi = -EF0phi.*Cpa;
    
    Ex = ErRho.*cos(phi_dash) - ErPhi.*sin(phi_dash);
    Ey = ErRho.*sin(phi_dash) + ErPhi.*cos(phi_dash);
    
    IntX = sum(Ex.*rho, 'all')*dRho*dph_dash;
    IntY = sum(Ey.*rho, 'all')*dRho*dph_dash;
    
    IntAbs = abs(IntX).^2+abs(IntY).^2;
    
    Eabs = sqrt(abs(Ex).^2 + abs(Ey).^2);
    %Eabs = sqrt(abs(ErRho).^2 + abs(ErPhi).^2);
    
    denom = sum(((Eabs).^2).*rho, 'all')*dRho*dph_dash;
    %numer = (sum(Eabs.*rho, 'all')*dRho*dph_dash).^2;
    
    taperS = IntAbs/(denom*Area); 
end
