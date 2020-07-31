%Finding the far field of reflector feed

%% Defining inputs
%Freq of operation
freq = 60e9;

%Speed of light
c = 3e8;

%Wavelength
lam = c/freq;

%Focal length
f = 100e-2;

%Diameter of the feed
D = 3*lam;

%Radius of feed
a = D/2;

%Defining the current distribution
j = [0 1 0]'; %Only across Y why?

%Defining the meshgrid
drad = pi/180;
[th, phi] = meshgrid(10e-7:drad:pi+10e-7-drad, 0:drad:2*pi-drad);

%Jf i.e. Magnitude of current
Jf = 1;

%Obs_pt
r = 10000*lam;

%% Calculating the Far field 
%Should this j be multiplied by 2? Coz, Schelkunoff's formulatin
[EFx, EFy, EFz] = FFRefFeed(freq, 1, j, a, r, th, phi);
Emag = sqrt(abs(EFx).^2 + abs(EFy).^2 + abs(EFz).^2); %Magnitude of E
Emax = max(max(Emag));

%Plotting Emag and Theta
figure;
plot([-th(1,90:-1:1).*180/pi, th(1,1:90).*180/pi], [mag2db(Emag(1,90:-1:1)/Emax), mag2db(Emag(1,1:90)/Emax)], 'LineWidth', 2);
title('Magnitude of E-field at plane Phi = 0');
xlabel('Theta(in deg)');
ylabel('Normalized E-field [dBV]');
ylim([-40, 0]);
grid on;

figure;
plot([-th(1,90:-1:1).*180/pi, th(1,1:90).*180/pi], [mag2db(Emag(90,90:-1:1)/Emax), mag2db(Emag(90,1:90)/Emax)], 'LineWidth', 2);
title('Magnitude of E-field at plane Phi = 90');
xlabel('Theta(in deg)');
ylabel('Normalized E-field [dBV]');
ylim([-40, 0]);
grid on;

%Plotting in Polar coordinates
figure;
polarplot([th(1,1:90), -th(1,90:-1:1)], [Emag(1,1:90)/Emax, Emag(1,90:-1:1)/Emax], 'LineWidth', 1.5);
title('Normalized E-field at Phi = 0');