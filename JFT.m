%Calculating the Fourier Transform of the Current
%Why the assumption that the current is only across Y?
%Why bessel function of first kind and first order?
function J = JFT(k0, a, th, j)
    Z = k0.*a.*sin(th); %Bessel of this is needed
    Const = pi*a^2./(Z); %Constant part of bessel Airy Function %2 pi?
    Bessel = besselj(1, Z); 
    J(1,:,:) = Const.*Bessel*j(1); %Airy Function
    J(2,:,:) = Const.*Bessel*j(2);
    J(3,:,:) = Const.*Bessel*j(3);
end
    