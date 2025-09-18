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

number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 256; %number of domains we integrate over - 3D - START SMALL?
degrees_of_freedom      = 1;   %vectors have one degree of freedom DOES THIS NEED TO BE ALTERED?
dimension               = 3;   %how many dimensions does our data have? TED FOR 3D
envelope_power          = 8;   %weight is (1-x^2)^power

%%% WHY DOES THIS HAVE 4 DIMENSIONS FOR 2D? LEAVE WITH FOUR DIMENSIONS FOR 3D AS PER OTHER EXAMPLES
size_vec                = [32,32,32]; %how many gridpoints should we use per integration?

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

labels{a} =  "U_x";
G(:,a)    = SPIDER_integrate( U, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "W_z";
G(:,a)    = SPIDER_integrate( W, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "P_x";
G(:,a)    = SPIDER_integrate( P, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "P_z";
G(:,a)    = SPIDER_integrate( P, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "U_xx";
G(:,a)    = SPIDER_integrate( U, [1,1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "W_zz";
G(:,a)    = SPIDER_integrate( W, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T";
G(:,a)    = SPIDER_integrate( T, [], grid, corners, size_vec, pol ); % 2d and 3d
scales(a) = 1;
a = a+1;

labels{a} =  "T_x";
G(:,a)    = SPIDER_integrate( T, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T_z";
G(:,a)    = SPIDER_integrate( T, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*U + 1, [], grid, corners, size_vec, pol );      
%norm_vec = repmat( norm_vec, [dof*numel(u0s),1] );

G = G./norm_vec;
G = G./scales;

