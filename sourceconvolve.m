function [ w_j ] = sourceconvolve( C_0, pj_0, tildeC, P, T, dt, nt, x, y, nx, ny)
%This function solves approximately a wave equation of the form (7) by 
% convolving the source term with the free space Green's function for 
% the wave equation with constant sound speed C_0. 
%   In the convolution, the space integration is approximated by one-point 
%   Gaussian quadrature on the rectangular domain, while the time 
%   integration is approximated by trapezoidal rule on the time interval 
%   [0, (nt-1)*dt].
%   
%
%   Inputs:
%       C_0 = scalar, background sound speed used in computing pj_0
%       pj_0 = 1 x nt cell, contains a matrix nx by ny at each time
%           representing the pressure on Cartesian grid. This is computed
%           using background values C_0 as well as sigma_0, uj_0.
%       tildeC = nx by ny matrix of perturbed sound speed on Cartesian
%           grid.
%       P,T = triangular mesh data used for space integration and
%           interpolation to triangle midpoints
%       dt = scalar, length of time discretization
%       nt = scalar, number of time points spaced by dt
%       x,y = 1 x (2/[dx or dy] +1) vector, contains x and y coordinates
%           which form the Cartesian grid.
%       nx, ny = scalars, number of x and y nodes respectively.
%
%   Outputs:
%       w_j = 1 x nt cell, contains a matrix nx by ny at each time
%           representing the pressure on Cartesian grid. 
%
%   Created by Sarah Vallélian 10/14/13
%       Last updated 10/19/13
%%%%%%


% Assuming a regular triangular FEM mesh of equally sized elements whose
% vertices coincide with Cartesian grid points
nelem=size(T,2);
sizeelem=4/nelem;
sizeP=size(P,2);


% Identify midpoints of each triangular element
mid=zeros(2,nelem);
for j=1:nelem
   mid(:,j)=(P(:,T(1,j))+P(:,T(2,j))+P(:,T(3,j)))./3;
end


% Construct an approximation to \partial_t^2 pj_0 and reshape to be on mesh
% P rather than Cartesian grid; also reshape tildeC
d_ttpj_0=cell(1, nt);

d_ttpj_0{1}=(pj_0{3}-2.*pj_0{2}+pj_0{1})/(dt^2);
d_ttpj_0{1}=reshape(d_ttpj_0{1}, sizeP, 1);

d_ttpj_0{nt}=(pj_0{nt}-2.*pj_0{nt-1}+pj_0{nt-2})/(dt^2);
d_ttpj_0{nt}=reshape(d_ttpj_0{nt}, sizeP, 1);

tildeC=reshape(tildeC, sizeP, 1);

%centered differences
for j=2:nt-1
    d_ttpj_0{j}=(pj_0{j-1}-2.*pj_0{j}+pj_0{j+1})/(dt^2);
    d_ttpj_0{j}=reshape(d_ttpj_0{j}, sizeP, 1);
end


% Interpolate \partial_t^2 pj_0 and tildeC on mesh to triangle midpoints
for j=1:nt
   d_ttpj_0{j}=pdeintrp(P,T,d_ttpj_0{j});
end
tildeC=pdeintrp(P,T,tildeC);


% Loop over times and fill in values of solution w_j at Cartesian grid
% points
w_j=cell(1,nt);
for k=1:nt
   w_j{k}=zeros(nx,ny);
   for i=1:nx
       for j=1:ny
           %compute the w_j at time tk=(k-1)*dt at grid point xij=(x(i),y(j))
           tempsum=0;
           xij=[x(i);y(j)];
           tk=(k-1)*dt;
           
           %approximate convolution
           for p=1:nelem
               sum1=0;
               for q=1:nt
                  if q==1 || q==nt
                      sum1=sum1+1/(2*pi)*Greenswave(mid(:,p)-xij, (q-1)*dt-tk, 1/sqrt(C_0))*tildeC(p)*d_ttpj_0{q}(p);
                  else
                      sum1=sum1+2/(2*pi)*Greenswave(mid(:,p)-xij, (q-1)*dt-tk, 1/sqrt(C_0))*tildeC(p)*d_ttpj_0{q}(p);
                  end
               end
               tempsum=tempsum+sum1*dt/2; %time integration 
           end
           
           w_j{k}(i,j)=sizeelem*tempsum; %space integration, note elements are all the same size
       end
   end
   stepcount=mod(k,100);
   if stepcount==0
      disp([' >> time point: ', num2str(k), ' << ']); 
   end 
end




end

