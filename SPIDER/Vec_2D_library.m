clear;
restoredefaultpath();
addpath("SPIDER/SPIDER_functions/");

%%
load("data1e8_.mat");

% Truncate the domain
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

clear data;

%%
number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 512; %number of domains we integrate over - 3D - START SMALL?
degrees_of_freedom      = 2;   %vectors have one degree of freedom DOES THIS NEED TO BE ALTERED?
dimension               = 3;   %how many dimensions does our data have? TED FOR 3D
envelope_power          = 6;   %weight is (1-x^2)^power
size_vec                = [64,64,32]; % [64,64,64];  %how many gridpoints should we use per integration?
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

%% 
%%%Coefficients

%Nondimensionalize and compute mean and variation of data
%TODO;

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;

labels{a} = "U";
%G(:,a)    =[SPIDER_integrate( U, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( W, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( U,W, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "U_t";
G(:,a)    = SPIDER_integrate_vector( U,W, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\nabla_i P";
G(:,a)    =[SPIDER_integrate( P, [1], grid, corners, size_vec, pol );
            SPIDER_integrate( P, [2], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "\nabla_i T";
G(:,a)    =[SPIDER_integrate( T, [1], grid, corners, size_vec, pol );
            SPIDER_integrate( T, [2], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "P.U";
%G(:,a)    =[SPIDER_integrate( P.*U, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( P.*W, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( P.*U, P.*W, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "T.U";
%G(:,a)    =[SPIDER_integrate( T.*U, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( T.*W, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( T.*U, T.*W, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "(U_j.\nabla_j) U_i"; % = \nabla_j (U_j U_i) via conservative form 
%G(:,a)    =[SPIDER_integrate( U.*dxU + W.dzU, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( U.*dxW + W.dzW, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( U.*U,  U.*W, [1], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate_vector( W.*U,  W.*W, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "\nabla^2 u_i";
G(:,a)    = SPIDER_integrate_vector( U, W, [1,1], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate_vector( U, W, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "U_tt";
G(:,a)    = SPIDER_integrate_vector( U,W, [3,3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "U^2.U";
U_sq = U.^2 + W.^2;
%G(:,a)    =[SPIDER_integrate( (U.^2 + W.^2).*U, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( (U.^2 + W.^2).*W, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( U_sq.*U, U_sq.*W, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "P^2.U";
%G(:,a)    =[SPIDER_integrate( (P.^2).*U, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( (P.^2).*W, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( (P.^2).*U, (P.^2).*W, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "T^2.U";
%G(:,a)    =[SPIDER_integrate( (T.^2).*U, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( (T.^2).*W, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( (T.^2).*U, (T.^2).*W, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "(\nabla_i P)_t";
%G(:,a)    =[SPIDER_integrate( dxP, [3], grid, corners, size_vec, pol );
%            SPIDER_integrate( dzP, [3], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( dxP, dzP, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "(\nabla_i T)_t";
%G(:,a)    =[SPIDER_integrate( dxT, [3], grid, corners, size_vec, pol );
%            SPIDER_integrate( dzT, [3], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( dxT, dzT, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "P \nabla_i P";
%G(:,a)    =[SPIDER_integrate( P.*dxP, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( P.*dzP, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( P.*dxP, P.*dzP, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "T \nabla_i T";
%G(:,a)    =[SPIDER_integrate( T.*dxT, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( T.*dzT, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( T.*dxT, T.*dzT, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "U.(\nabla.U)";
%G(:,a)    =[SPIDER_integrate( U.*dxU + U.*dzW, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( W.*dxU + W.*dzW, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( U.*dxU + U.*dzW, W.*dxU + W.*dzW, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

%%% ADD U.(\nabla U) %%%

labels{a} = "\nabla(\nabla.U)";
G(:,a)    =[SPIDER_integrate( dxU + dzW, [1], grid, corners, size_vec, pol );
            SPIDER_integrate( dxU + dzW, [2], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} =  "P U_t";
%G(:,a)    =[SPIDER_integrate( P.*dtU, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( P.*dtW, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( P.*dtU, P.*dtW, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T U_t";
%G(:,a)    =[SPIDER_integrate( T.*dtU, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( T.*dtW, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( T.*dtU, T.*dtW, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "U P_t";
%G(:,a)    =[SPIDER_integrate( U.*dtP, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( W.*dtP, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( U.*dtP, W.*dtP, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "U T_t";
%G(:,a)    =[SPIDER_integrate( U.*dtT, [], grid, corners, size_vec, pol );
%            SPIDER_integrate( W.*dtT, [], grid, corners, size_vec, pol );];
G(:,a)    = SPIDER_integrate_vector( U.*dtT, W.*dtT, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "P e_z";
G(:,a)    =[0*SPIDER_integrate( P, [], grid, corners, size_vec, pol );
            SPIDER_integrate( P, [], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "T e_z";
G(:,a)    =[0*SPIDER_integrate( T, [], grid, corners, size_vec, pol );
            SPIDER_integrate( T, [], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

%%%% Advection of pressure %%%%
%labels{a} = "(U_j \nabla_j) P"; % = \nabla_j (u_j P) via conservative form
%G(:,a)    = SPIDER_integrate_vector( U.*P, W.*P, [1], grid, corners, size_vec, pol ) + ...
%            SPIDER_integrate_vector( U.*P, W.*P, [2], grid, corners, size_vec, pol );
%scales(a) = 1;
%a = a+1;

%%%% Advection of temperature %%%%
%labels{a} = "(U_j \nabla_j) T"; % = \nabla_j (u_j T); via conservative form
%G(:,a)    = SPIDER_integrate_vector( U.*T, W.*T, [1], grid, corners, size_vec, pol ) + ...
%            SPIDER_integrate_vector( U.*T, W.*T, [2], grid, corners, size_vec, pol );
%scales(a) = 1;
%a = a+1;

labels{a} = "\nabla_i U^2";
G(:,a)    =[SPIDER_integrate( U_sq, [1], grid, corners, size_vec, pol );
            SPIDER_integrate( U_sq, [2], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "\nabla_i P^2";
G(:,a)    =[SPIDER_integrate( P.*P, [1], grid, corners, size_vec, pol );
            SPIDER_integrate( P.*P, [2], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "\nabla_i T^2";
G(:,a)    =[SPIDER_integrate( T.*T, [1], grid, corners, size_vec, pol );
            SPIDER_integrate( T.*T, [2], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;
%%%

%% normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*U + 1, [], grid, corners, size_vec, pol );      
norm_vec = repmat( norm_vec, [dof,1] );

G = G./norm_vec;
G = G./scales;

%% 
function c=SPIDER_integrate_vector( U,W, deriv, grid, corners, size_vec, pol )
   a = SPIDER_integrate( U, deriv, grid, corners, size_vec, pol );
   b = SPIDER_integrate( W, deriv, grid, corners, size_vec, pol );
   c = [a;b];
end

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
