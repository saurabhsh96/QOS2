%Finding the Current Distribution of the reflector
%clear;
%% Defining inputs of feed
function EFRMag = Q2(zeta, freq, Df, jf, fRef, DRef, th_obs, theta_obs, phi_obs)
    %Zeta
    % zeta = 377;

    %Freq of operation
    % freq = 60e9;

    %Speed of light
    c = 3e8;

    %Wavelength
    lam = c/freq;

    %Wavenumber
    k0 = 2*pi/lam;

    %Diameter of the feed
    % Df = 3*lam;

    %Radius of feed
    af = Df/2;

    %Defining the current distribution of feed
    % jf = [0 1 0]'; %Only across Y why?

    %Defining the meshgrid
    drad = pi/180;
    [th, phi] = meshgrid(eps:drad:pi/2-drad, eps:drad:2*pi);

    %Jf i.e. Magnitude of current
    Jf = 1;

    %% Defining inputs of the reflector

    %Focal length of the reflector 
    % fRef = 100e-2; 

    %Diameter of the reflector f/D = 0.5 in this case
    % DRef = fRef/5;

    %Radius of the reflector
    aRef = DRef/2;

    %Parameterization %Drad is same as before
    %r' th'
    dRho = 0.0001;
    rho_dash1 = eps:dRho:aRef;
    theta0 = 2*atan(DRef/(4*fRef));
    [rho_dash, phi_dash] = meshgrid(rho_dash1, eps:drad:2*pi);
    th_dash = 2*atan(rho_dash./(2*fRef));
    r_dash = fRef.*(1+(tan(th_dash/2)).^2);

    %rho' z

    %z = -fRef + rho_dash.^2./(4*fRef);

    %% Calculating the Far field of the feed
    %Obs_pt for feed is r'
    r_temp = 10000*lam; %Finding field at this point temperorily to negate the phase effect later
    Constant_term = (r_temp)./(exp(-1j * k0 * r_temp)); %To negate initial phase effect; 
    %We can multiply be 2*pi as well

    %Should this j be multiplied by 2? Coz, Schelkunoff's formulatin
    [EFfx, EFfy, EFfz] = FFRefFeed(freq, 1, jf, af, r_temp, th_dash, phi_dash);
    Emagf = sqrt(abs(EFfx).^2 + abs(EFfy).^2 + abs(EFfz).^2); %Magnitude of E
    Emaxf = max(max(Emagf));

    %Getting Spherical Coordinates
    EFfr = (EFfx.*sin(th_dash).*cos(phi_dash)) + (EFfy.*sin(th_dash).*sin(phi_dash)) + (EFfz.*cos(th_dash));
    EFfth = (EFfx.*cos(th_dash).*cos(phi_dash)) + (EFfy.*cos(th_dash).*sin(phi_dash)) - (EFfz.*sin(th_dash));
    EFfphi = (-EFfx.*sin(phi_dash)) + (EFfy.*cos(phi_dash));

    %Plotting Emag and Theta
%     plot(th_dash(1,:)*180/pi, mag2db(Emagf(1,:)/Emaxf), 'LineWidth', 2);
%     title('Magnitude of E-field Phi = 0');
%     xlabel('Theta(in deg)');
%     ylabel('Normalized E-field [dBV]');

    %Plotting in Polar coordinates
%     figure();
%     polarplot(th_dash(1,:), Emagf(1,:)/Emaxf, 'LineWidth', 2);
%     title('Normalized E-field at Phi = 0');

    %% Equivalent currents
    Js_rho = EFfth.*exp(-2j*k0*fRef).*(cos(th_dash/2)).^2./(fRef*zeta).*Constant_term;
    Js_phi = EFfphi.*exp(-2j*k0*fRef).*(cos(th_dash/2)).^2./(fRef*zeta).*Constant_term;

    %Js = sqrt(abs(JsX).^2 + abs(JsY).^2);

    Jx = (Js_rho.*cos(phi_dash)) - (Js_phi.*sin(phi_dash));
    Jy = (Js_rho.*sin(phi_dash)) + (Js_phi.*cos(phi_dash));
    Jz = zeros(size(rho_dash));
    %[JsX, JsY] = pol2cart(Js_phi, Js_rho);

    Axx = (rho_dash).*cos(phi_dash);
    Ayy = rho_dash.*sin(phi_dash);

%     figure();
%     surface(Axx, Ayy, abs(Jx)./max(max(abs(Jx))), 'linestyle','none');
%     title('Normalized Current Distribution at Equivalent Aperture, Jx');
%     xlabel('x [m]');
%     ylabel('y [m]');
% 
%     figure();
%     surface(Axx, Ayy, abs(Jy)./max(max(abs(Jy))), 'linestyle','none');
%     title('Normalized Current Distribution at Equivalent Aperture, Jy');
%     xlabel('x [m]');
%     ylabel('y [m]');

    %% Far field of the reflector
    r_obs = 1000*lam;

%     th_obs = linspace( -5*(lam/DRef), 5*(lam/DRef), 181) ;
%     [theta_obs, phi_obs] = meshgrid(th_obs, eps:pi/4:2*pi);

    [JFTx, JFTy, JFTz] = ApertureJFT(Jx, Jy, Jz, ...
        dRho, drad, k0, rho_dash, phi_dash,theta_obs, phi_obs);

    %JFTs = sqrt(abs(JFTx).^2 + abs(JFTy).^2 + abs(JFTz).^2);
    JFTa = zeros([3, size(JFTx)]);
    JFTa(1,:,:)=JFTx;
    JFTa(2,:,:)=JFTy;
    JFTa(3,:,:)=JFTz;

    %plot(th_obs(1,:), JFTs(1,:));
    % 
    %Increaments in X and Y
    % dx = -dRho*drad.*sin(phi_dash).*cos(phi_dash);
    % dy = dRho*drad.*sin(phi_dash).*cos(phi_dash);

    % th_f = eps:8*drad*dRho:8*drad-8*drad*dRho;
    % th_f = eps:5*drad/1000:5*drad;
    % [th_ff, phi_ff] = meshgrid(th_f, eps:drad:2*pi);
    % sz = size(th_ff(1,:));
    % szph = size(phi_ff(:,1));
    % JFTax = zeros(szph(1), sz(2));
    % JFTay = zeros(szph(1), sz(2));
    % %JFTF = zeros(szph(1), sz(2));
    % phi_f = 0;
    % for indth = 1:sz(2)
    %     kx = k0.*sin(th_f(1,indth)).*cos(phi_ff);
    %     ky = k0.*sin(th_f(1,indth)).*sin(phi_ff);
    %     kz = k0.*cos(th_f(1,indth));
    %     
    % %     phase = exp(1j.*kx.*rho_dash.*cos(phi_ff)).*exp(1j.*ky.*rho_dash.*sin(phi_ff));
    % %     JFTax(:,indth) = (sum(Jx.*phase.*rho_dash, 'all'))*dx*dy;
    % %     JFTay(:,indth) = (sum(Jy.*phase.*rho_dash, 'all'))*dx*dy;
    %     phasex = exp(1j.*kx.*rho_dash.*cos(phi_ff));
    %     phasey = exp(1j.*ky.*rho_dash.*sin(phi_ff));
    %     JFTax(:,indth) = sum(sum(Jx.*phasex.*phasey.*rho_dash))*dRho*drad;
    %     JFTay(:,indth) = sum(sum(Jy.*phasey.*phasex.*rho_dash))*dRho*drad;
    %     %JFTa(:,ind) = sum(JFTax.*JFTay, 'all')*dx*dy;
    % end
    % 
    % % JFTF = JFTax.*JFTay;
    % JFTA = zeros(3, szph(1), sz(2));
    % JFTA(1,:,:) = JFTax;
    % JFTA(2,:,:) = JFTay;
    % JFTA(3,:,:) = 0;
    % 
    % JFTaMag = sqrt(abs(JFTax).^2 + abs(JFTay).^2);
    % JFTaMax = max(max(JFTaMag));
    % % 
    % % JFTaMagF = (abs(JFTF));
    % % JFTaMaxF = max(max(JFTaMagF));
    % 
    % % 
    %figure();
    %plot(th_f*180/pi, mag2db(JFTaMag(1,:)/JFTaMax));

    %Field by Parabolic Aperture
    [EFRx, EFRy, EFRz] = FFRef(freq, 1, r_obs, th_obs, phi_obs, JFTa);
    EFRMag = sqrt(abs(EFRx).^2 + abs(EFRy).^2 + abs(EFRz).^2); 
    EFRMagMax = max(max(EFRMag));

    %Field by unform aperture
    [EFxu, EFyu, EFzu] = FFRefFeed(freq, 1, jf, aRef, r_obs, theta_obs, phi_obs);
    EFRMagu = sqrt(abs(EFxu).^2 + abs(EFyu).^2 + abs(EFzu).^2); 
    EFRMagMaxu = max(max(EFRMagu));

%     figure();
%     plot(th_obs(1,:)*180/pi, mag2db(EFRMag(1,:)/EFRMagMax)); 
%     hold on;
%     plot(th_obs(1,:)*180/pi, mag2db(EFRMagu(1,:)/EFRMagMaxu), '--');
%     title('Normalized E-field at plane phi = 0');
%     xlabel('Normalized E-field [dBV]');
%     ylabel('Observation angle, Theta (in deg)');
%     ylim([-35, 0]);
%     legend('Parabolic Aperture', 'Uniform Aperture');
% 
%     figure();
%     plot(th_obs(1,:)*180/pi, mag2db(EFRMag(3,:)/EFRMagMax)); 
%     hold on;
%     plot(th_obs(1,:)*180/pi, mag2db(EFRMagu(3,:)/EFRMagMaxu), '--');
%     title('Normalized E-field at plane phi = 90');
%     xlabel('Normalized E-field [dBV]');
%     ylabel('Observation angle, Theta (in deg)');
%     ylim([-35, 0]);
%     legend('Parabolic Aperture', 'Uniform Aperture');

    % plot([-th_ff(1,1000:-1:1)*180/pi, th_ff(1,1:1:1000)*180/pi], ...
    %     [mag2db(EFRMag(1,1000:-1:1)/EFRMagMax), mag2db(EFRMag(1,1:1000)/EFRMagMax)]);

    %% Calculating the Efficiencies

    %Spillover efficiency

    %Range of radius %fRef is constant as above i.e. 1m, keeping r_temp const
    %as above
    %Why NaN values?

    Arange = 0.2*fRef/2:0.01:2*fRef/2;
    etaS = zeros(size(Arange));
    taperS = zeros(size(Arange));
    Area = zeros(size(Arange));
    iter = 1;
    for aRefR = Arange   
        rho_dash1R = eps:dRho:aRefR;
        theta0R = 2*atan(aRefR*2/(4*fRef));
        [rho_dashR, phi_dashR] = meshgrid(rho_dash1R, eps:drad:2*pi);
        th_dashR = 2*atan(rho_dashR./(2*fRef));
        r_dashR = fRef.*(1+(tan(th_dashR/2)).^2);
        [taperS(1, iter), etaS(1, iter), Area(1, iter)] = Spillover(freq, 1, jf, af, r_temp, ...
            th_dashR, phi_dashR, th, phi, rho_dashR, fRef, aRefR);
        iter = iter+1;
    end

    etaA = taperS.*etaS;

%     figure();
%     plot(Arange*2, taperS(1,:),'LineWidth', 1.5);
%     hold on;
%     plot(Arange*2, etaS(1,:), 'LineWidth', 1.5);
%     plot(Arange*2, etaA(1,:), 'LineWidth', 1.5);
%     ylim([0, 1]);
%     ylabel('Efficiencies');
%     xlabel('Diameter [m]');
%     legend('Taper', 'Spillover', 'Aperture');
%     title('Efficiencies (Taper, Spillover and Aperture)');

    %% Calculating Directivity

    %Maximum directiviy
    Dirm = (4*pi/lam^2).*Area;

    %Directivity
    Dir = Dirm.*taperS;

    %Gain
    G = Dirm.*etaA;

%     figure();
%     plot(Arange*2, pow2db(Dirm(1,:)),'LineWidth', 1.5);
%     hold on;
%     plot(Arange*2, pow2db(Dir(1,:)), 'LineWidth', 1.5);
%     plot(Arange*2, pow2db(G(1,:)), 'LineWidth', 1.5);
%     %ylim([0, 1]);
%     ylabel('Dmax, Dir and Gain (dBi)');
%     xlabel('Diameter [m]');
%     legend('Dmax', 'Dir', 'Gain');
%     title('Directivity and Gain of the reflector');
end 