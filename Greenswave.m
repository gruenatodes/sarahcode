function [ value ] = Greenswave( z, t, c)
%Greenswave is a function which computes the value of the fundamental
%solution to the wave equation (with constant sound speed c) at the point
%(z,t) in R^2 x R+
%   See paper for derivation of the Green's function.
%
%   Inputs:
%       z = 2xN column vector, a collection of points in R^2
%       t = Mx1 column vector, times in R+
%       c = scalar, constant sound speed
%
%   Outputs:
%       value = MxN matrix of values; each row indicates a certain time,
%           each column indicates a certain point in R^2
%
%   Note that constant coefficient in front of Green's function is missing;
%   different coefficients must be included depending on if the delta is in
%   the source or the initial condition.

ts=tic;

N=size(z,2);
M=size(t,1);
value=zeros(M,N);

% for i=1:N
%     for j=1:M
%         value(j,i)=heaviside(t(j)-norm(z(:,i))/c)/sqrt(t(j)^2 - (norm(z(:,i))^2)/(c^2));
%     end
% end


z_norms = sqrt(sum(z .^ 2, 1)); % is a row vector containing 2-norms of columns of z
z_norms_mat = repmat(z_norms, [M 1]);
t_mat = repmat(t, [1 N]);

heavis_arg = (t_mat - z_norms_mat) ./ c; % arg of heaviside
num = zeros(M,N);
num(heavis_arg > 0) = 1;
num(heavis_arg == 0) = 0.5;

denom = (t_mat .^ 2 - z_norms_mat .^ 2) ./ c ^ 2;

value = num ./ denom;

% if any(isinf(value)) || any(isnan(value))
%     error
% end

te=toc(ts);
% disp(['Greenswave finished in ', num2str(te), ' seconds']);


end

