clear;
restoredefaultpath();
addpath("SPIDER/SPIDER_functions/");

%% 
load("data1e8_.mat");
x1 = 100;
x2 = 355;
z1 = 1;
z2 = 128;
t1 = 1;
t2 = 50;

U = data(x1:x2,z1:z2,t1:t2,1);
W = data(x1:x2,z1:z2,t1:t2,2);
P = data(x1:x2,z1:z2,t1:t2,3);
T = data(x1:x2,z1:z2,t1:t2,4);

x = x(x1:x2);
z = z(z1:z2);
t = t(t1:t2);

dxU = diff_dim(U,x,1);
dzU = diff_dim(U,z,2);
dtU = diff_dim(U,t,3);
dxW = diff_dim(W,x,1);
dzW = diff_dim(W,z,2);
dtW = diff_dim(W,t,3);
dxP = diff_dim(P,x,1);
dzP = diff_dim(P,z,2);
dtP = diff_dim(P,t,3);
dxT = diff_dim(T,x,1);
dzT = diff_dim(T,z,2);
dtT = diff_dim(T,t,3);

number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 512; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 3;   %how many dimensions does our data have?
envelope_power          = 8;   %weight is (1-x^2)^power
size_vec                = [32,32,32,32]; %how many gridpoints should we use per integration?
buffer                  = 0; %Do not use points this close to boundary

%define shorthand notation
nl = number_of_library_terms;
nw = number_of_windows;
dof= degrees_of_freedom;

%Make important objects for integration
pol      = envelope_pol( envelope_power, dimension );
G        = zeros( dof*nw, nl );
labels   = cell(nl, 1);
scales   = zeros(1,nl);

%we also need pol to be odd in time. multiply it by t
pol0 = pol; %save for normalization

size_of_data = size(U, 1:dimension);

seed = 1;
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

grid = {x,z,t};

%Nondimensionalize and compute mean and variation of data
%TODO;

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;

labels{a} =  "dTdt";
G(:,a)    = SPIDER_integrate( T, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d2Tdt2";
G(:,a)    = SPIDER_integrate( T, [3,3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T";
G(:,a)    = SPIDER_integrate( T, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T^2";
G(:,a)    = SPIDER_integrate( T.*T, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d(T^2)dt";
G(:,a)    = SPIDER_integrate( T.*T, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "U^2";
G(:,a)    = SPIDER_integrate( (U.*U + W.*W), [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d(U^2)dt";
G(:,a)    = SPIDER_integrate( (U.*U + W.*W), [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "U^2.P";
G(:,a)    = SPIDER_integrate( (U.*U + W.*W).*P, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "Grad u";
G(:,a)    = SPIDER_integrate( dxU + dzW, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "p (Grad u)";
G(:,a)    = SPIDER_integrate( P.*(dxU + dzW), [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d(UT)dx";
G(:,a)    = SPIDER_integrate( U.*T, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d(WT)dz";
G(:,a)    = SPIDER_integrate( W.*T, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T_xx";
G(:,a)    = SPIDER_integrate( T, [1,1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T_zz";
G(:,a)    = SPIDER_integrate( T, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*U + 1, [], grid, corners, size_vec, pol );      
%norm_vec = repmat( norm_vec, [dof*numel(u0s),1] );

G = G./norm_vec;
G = G./scales;

function [dF] = diff_dim(F, y, dim)
% 3-point forward/center/backwards differencing to approximate first derivative

    dF = zeros(size(F));

    % boundary points
    y0 = y(1); y1 = y(2); y2 = y(3);
    c0 = (2*y0-y1-y2)/((y0-y1)*(y0-y2));
    c1 = (2*y0-y0-y2)/((y1-y0)*(y1-y2));
    c2 = (2*y0-y0-y1)/((y2-y0)*(y2-y1));
    if dim == 1
        dF(1, :, :) = c0*F(1, :, :)+c1*F(2, :, :)+c2*F(3, :, :);
    elseif dim == 2
        dF(:, 1, :) = c0*F(:, 1, :)+c1*F(:, 2, :)+c2*F(:, 3, :);
    else
        dF(:, :, 1) = c0*F(:, :, 1)+c1*F(:, :, 2)+c2*F(:, :, 3);
    end

    for i=2:length(y)-1
        y0 = y(i-1); y1 = y(i); y2 = y(i+1);
        c0 = (2*y1-y1-y2)/((y0-y1)*(y0-y2));
        c1 = (2*y1-y0-y2)/((y1-y0)*(y1-y2));
        c2 = (2*y1-y0-y1)/((y2-y0)*(y2-y1));
        if dim == 1
            dF(i, :, :) = c0*F(i-1, :, :)+c1*F(i, :, :)+c2*F(i+1, :, :);
        elseif dim == 2
            dF(:, i, :) = c0*F(:, i-1, :)+c1*F(:, i, :)+c2*F(:, i+1, :);
        else
            dF(:, :, i) = c0*F(:, :, i-1)+c1*F(:, :, i)+c2*F(:, :, i+1);
        end
    end

    y0 = y(end-2); y1 = y(end-1); y2 = y(end);
    c0 = (2*y2-y1-y2)/((y0-y1)*(y0-y2));
    c1 = (2*y2-y0-y2)/((y1-y0)*(y1-y2));
    c2 = (2*y2-y0-y1)/((y2-y0)*(y2-y1));
    if dim == 1
        dF(end, :, :) = c0*F(end-2, :, :)+c1*F(end-1, :, :)+c2*F(end, :, :);
    elseif dim == 2
        dF(:, end, :) = c0*F(:, end-2, :)+c1*F(:, end-1, :)+c2*F(:, end, :);
    else
        dF(:, :, end) = c0*F(:, :, end-2)+c1*F(:, :, end-1)+c2*F(:, :, end);
    end

end
