function output = nonlocalMeansDenoise( I, patch_rad, search_rad, patch_sigma, search_sigma, selfsim )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  I            : grayscale 2D or 3D image to be filtered
%  patch_rad    : radius of similarity window, increase to improve
%  denoising performance at the cost of increased calculation time
%  search_rad   : radius of search window, increase to denoise more at the
%  cost of increased calculation time and, possibly, some blur
%  patch_sigma  : sigma for the center-periphery Gausssian weights of a patch
%  (defaults to 1)
%  search_sigma : sigma for the patch-patch similarity Gausssian weights
%  (defaults to 10, increase to increase denoising at the cost of increased blur)
%  selfsim      : Self-patch input to the similarity mean
%  (defaults to 0)
%
%  Note:
%    if selfsim = 0, then w(i,i) = max_{j neq i} w(i,j), for all i
%
%  Author: Yury Petrov
%  Date: May, 2021
%
%  Extension of Christian Desrosiers's simple_nlm() (Date: 07-07-2015) to 3D images
% 
%  For details see:
%     A. Buades, B. Coll and J.M. Morel, "A non-local algorithm for image denoising"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error( 'Need a 2D or a 3D image as the first input argument!' );
else
    nd = ndims( I );
    if nd < 2 || nd > 3
        error( 'Need a 2D or a 3D image as the first input argument!' );
    end
end

if nargin < 6
    selfsim = 0;
end
if nargin < 5
    search_sigma = 10;
end
if nargin < 4
    patch_sigma = 1;
end
if nargin < 3
    if nd == 2
        search_rad = 3;
    else
        search_rad = 1;
    end
end
if nargin < 2
    if nd == 2
        patch_rad = 2;
    else
        patch_rad = 1;
    end
end

[ m, n, r ] = size( I );
s = m * n * r;
pixels = I(:);

psize = 2 * patch_rad + 1; % similarity patch size
nsize = 2 * search_rad + 1; % neighborhood search size

% Compute patches
if nd == 2
    padInput = padarray( I, [ patch_rad patch_rad ], 'symmetric' );
    filter = fspecial( 'gaussian', psize, patch_sigma );
    % image patches with intensity falling-off as a Gaussian from a patch's center
    patches = repmat( sqrt( filter(:) )', [ s 1 ] ) .* im2col( padInput, [ psize psize ], 'sliding' )'; % npatches x patchlength
else
    padInput = padarray( I, [ patch_rad patch_rad patch_rad ], 'symmetric' );
    filter = fspecial3( 'gaussian', psize, patch_sigma );
    patches = repmat( sqrt( filter(:) )', [ s 1 ] ) .* im2col3( padInput, [ psize psize psize ], [ 1 1 1 ] )';
end

% Compute list of edges (pixel pairs within the same search window)
indexes = reshape( 1 : s, m, n, r );
if nd == 2
    padIndexes = padarray( indexes, [ search_rad search_rad ] );
    neighbors = im2col( padIndexes, [ nsize nsize ], 'sliding' ); % indices within each block in columns
    target = repmat( 1:s, [ nsize^2 1 ] ); % target voxel indices in columns
else
    padIndexes = padarray( indexes, [ search_rad search_rad search_rad ] ); % zeros go to non-existing neighbors
    neighbors = im2col3( padIndexes, [ nsize nsize nsize ], [ 1 1 1 ] ); 
    target = repmat( 1:s, [ nsize^3 1 ] );
end
edges = [ target(:) neighbors(:) ];
edges( target(:) >= neighbors(:), : ) = [];  % all possible connections from a block start voxel to voxels within a block
clear target neighbors padIndexes
% Compute weight matrix (using weighted Euclidean distance)
pixPerPatch = size( patches, 2 );
if s * pixPerPatch <= 3000000
    diff = patches( edges( :, 1 ), : ) - patches( edges( :, 2 ), : ); % nedges x patchlength, difference vector between two image patches
    V = exp( -sum( diff .* diff, 2 ) / search_sigma^2 ); % similarity scalar between two patches in an edge
    clear diff
else % loop over each pixel in the patch to fit into memory
    disp(1)
    diff2 = zeros( size( edges, 1 ), 1 );
    for p = 1 : pixPerPatch % loop over patch pixels
        diff = patches( edges( :, 1 ), p ) - patches( edges( :, 2 ), p ); % nedges x 1, difference vector between corresponding pixels in image patches
        diff2 = diff2 + diff .* diff;
    end
    V = exp( -diff2 / search_sigma^2 ); % similarity scalar between two patches in an edge
    clear diff2 diff
end
disp(size(edges))
disp(size(V))
W = sparse( edges( :, 1 ), edges( :, 2 ),double(V), s, s ); % store similarities as a sparse matrix

clear V edges patches
 

% Make matrix symetric and set diagonal elements
if selfsim > 0
    W = W + W' + selfsim * speye( s );
else
    maxv = max( W, [], 2 );
    W = W + W' + spdiags( maxv, 0, s, s );
end

% Normalize weights
W = spdiags( 1 ./ sum( W, 2 ), 0, s, s ) * W;

% Compute denoised image
% W=full(W);%从稀疏矩阵还原为普通矩阵
disp(class(pixels))
disp(class(W))
output = W * pixels; % average image values over all pixels with similar patch neighborhoods
output = reshape( output, m, n, r );
clear W
end


% function output=im2col3(I,psize,type)
% first=im2col(I(:,:,1),[psize psize],'sliding');
% output=zeros(size(first,1),size(first,2),size(I,3));
% output(:,:,1)=first;
% for i=2:size(I,3)
%     output(:,:,i)=im2col(I(:,:,i),[psize psize],'sliding');
% end
% end