%JFT of the aperture current %Thanks to Tworit!
function [J_ft_x, J_ft_y, J_ft_z] = ApertureJFT(J_cx, J_cy, J_cz, drho, dphi, k0, rho, ph, theta_obs, phi_obs)

    J_ft_x = zeros(size(theta_obs, 1), size(theta_obs, 2));
    J_ft_y = zeros(size(theta_obs, 1), size(theta_obs, 2));
    J_ft_z = zeros(size(theta_obs, 1), size(theta_obs, 2));
    
     for i = 1:size(theta_obs, 1)
         for j = 1:size(theta_obs, 2) %for reflector

            J_i_x = J_cx .* exp(1j * k0 .* rho .* sin(theta_obs(i, j)) .* cos(phi_obs(i, j) - ph)) .* rho ;
            J_i_y = J_cy .* exp(1j * k0 .* rho .* sin(theta_obs(i, j)) .* cos(phi_obs(i, j) - ph)) .* rho ;
            J_i_z = J_cz .* exp(1j * k0 .* rho .* sin(theta_obs(i, j)) .* cos(phi_obs(i, j) - ph)) .* rho;

            J_ft_x(i , j) = sum(sum(J_i_x))*dphi*drho;
            J_ft_y(i , j) = sum(sum(J_i_y))*dphi*drho;
            J_ft_z(i , j) = sum(sum(J_i_z))*dphi*drho;
         end
     end
end