function [ tildew_j ] = ICconvolve( C_0, tildeu, tildesigma, sigma_0, u_0, P, T, dt, nt, x, y, nx, ny)
%This function solves approximately a wave equation of the form (12) by 
% convolving the initial condition with the free space Green's function for 
% the wave equation with constant sound speed C_0. 
%   Note this assumes Gamma = 1. 
%   In the convolution, the space integration is approximated by one-point 
%   Gaussian quadrature on the rectangular domain. There is no time
%   integration, but the (distributional) derivative in the Green's
%   function is approximated with second-order finite difference.
%   
%
%   Inputs:
%       C_0 = scalar, background sound speed used in computing pj_0
%       tildeu = nx by ny matrix of perturbed u on Cartesian grid.
%       tildesigma = nx by ny matrix of perturbed absorption on Cartesian
%           grid.
%       sigma_0 = nx by ny matrix of background absorption on Cartesian
%           grid.
%       u_0 = nx by ny matrix of background u on Cartesian grid.
%       P,T = triangular mesh data used for space integration and
%           interpolation to triangle midpoints
%       dt = scalar, length of time discretization
%       nt = scalar, number of time points spaced by dt
%       x,y = 1 x (2/[dx or dy] +1) vector, contains x and y coordinates
%           which form the Cartesian grid.
%       nx, ny = scalars, number of x and y nodes respectively.
%
%   Outputs:
%       tildew_j = 1 x nt cell, contains a matrix nx by ny at each time
%           representing the pressure on Cartesian grid. 
%
%   Created by Sarah Vallélian 10/17/13
%       Last updated 10/19/13
%%%%%%


% Assuming a regular triangular FEM mesh of equally sized elements whose
% vertices coincide with Cartesian grid points
nelem = size(T, 2);
sizeelem = 4 / nelem;
sizeP = size(P, 2);


% Identify midpoints of each triangular element
mid = zeros(2, nelem);
for j = 1:nelem
   mid(:,j) = (P(:,T(1,j)) + P(:,T(2,j)) + P(:,T(3,j))) ./ 3;
end

% % same without loop:
% mid = ( P(:,T(1,:)) + P(:,T(2,:)) + P(:,T(3,:)) ) / 3;
% % the contents of T serve as indices of P.

% % simpler example:
% % AAA = repmat(1:5, [6 1]);
% % indexvec = [1 3 5 2 4 5 4 3 2 1 1];
% % AAA(:, indexvec)




% Reshape tildeu, tildesigma, sigma_0, u_0 to be on mesh P rather than 
% Cartesian grid
tildeu = reshape(tildeu, [sizeP, 1]);
tildesigma = reshape(tildesigma, [sizeP, 1]);
sigma_0 = reshape(sigma_0, [sizeP, 1]);
u_0 = reshape(u_0, [sizeP, 1]);

% Interpolate tildeu, tildesigma, sigma_0, u_0 on mesh to triangle
% midpoints
tildeu = pdeintrp(P, T, tildeu);
tildesigma = pdeintrp(P, T, tildesigma);
sigma_0 = pdeintrp(P, T, sigma_0);
u_0 = pdeintrp(P, T, u_0);

% % if you don't want to type the same command 4 times, you could use alternatively:
% temp = cellfun( @(x) pdeintrp(P, T, reshape(x, [sizeP, 1])), ...
%     {tildeu, tildesigma, sigma_0, u_0} );
% tildeu = temp{1}; tildesigma = temp{2}; sigma_0 = temp{3}; u_0 = temp{4};
% clear temp



% Loop over times and fill in values of solution w_j at Cartesian grid
% points
tildew_j = zeros(nx,ny,nt); % better to prealloate arrays than structures or cells
k_vec = (1:nt-2)' * dt;
tic
for i = 1:nx
    for j = 1:ny
        %compute the tildew_j at time tk=(k-1)*dt at grid point xij = (x(i),y(j))
        xij = [x(i); y(j)];
        
        %approximate convolution
        % ts = tic;
        tempmat = C_0 / (2 * pi) / dt * ...
            ( Greenswave_new(mid - repmat(xij, [1 nelem]), k_vec + dt/2, 1/sqrt(C_0)) - ...
            Greenswave_new(mid - repmat(xij, [1 nelem]), k_vec - dt/2, 1/sqrt(C_0)) ) .* ...
            repmat(tildesigma .* u_0 + sigma_0 .* tildeu, [nt-2 1]); % returns nelem-by-(nt-2) matrix to be summed over nelem
        % c entered difference for \partial_t G: second order
        
        % te = toc(ts);
        %disp(['tildew_j at time k at point xij: convolution done in ' num2str(te) ' seconds']);
        tildew_j(i,j,2:nt-1) = sizeelem * sum(tempmat, 2); %space integration, note elements are all the same size
    end
    
    if mod(i,10) == 0
        toc
        disp([' >> time point: j = ', num2str(j), ' << ']);
        pause(0.5)
    end
end


% if k==1
%     tempsum = tempsum + C_0 / (2 * pi) * ...
%         ( -Greenswave(mid(:,p)-xij, tk+2*dt, 1/sqrt(C_0)) + ...
%             4 * Greenswave(mid(:,p)-xij, tk+dt, 1/sqrt(C_0)) - ...
%             3 * Greenswave(mid(:,p)-xij, tk, 1/sqrt(C_0)) ) / ...
%         (2*dt) * (tildesigma(p) * u_0(p) + sigma_0(p) * tildeu(p));
%     % forward difference for \partial_t G: second order
%     
% elseif k == nt
%     tempsum = tempsum + C_0 / (2 * pi) * ...
%         ( -Greenswave(mid(:,p)-xij, tk-2*dt, 1/sqrt(C_0)) + ...
%             4 * Greenswave(mid(:,p)-xij, tk-dt, 1/sqrt(C_0)) - ...
%             3 * Greenswave(mid(:,p)-xij, tk, 1/sqrt(C_0)) ) / ...
%         (2*dt) * (tildesigma(p) * u_0(p) + sigma_0(p) * tildeu(p));
%     % backward difference for \partial_t G: second order
% else
    
end


