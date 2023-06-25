function out = im2col3( I, blocksize, stepsize )
%
% 3D analog of im2col( I, blocksize ) with specified step, e.g., stepsize =
% [ 1 1 1 ] corresponds to im2col( I, blocksize, 'sliding' ), while
% stepsize = blocksize corresponds to im2col( I, blocksize, 'distinct' )
%
% Author: Yury Petrov
% Date: May, 2021
%

% block size
nrows = blocksize( 1 );
ncols = blocksize( 2 );
nslcs = blocksize( 3 );

% step size
d_row = stepsize( 1 );
d_col = stepsize( 2 );
d_slc = stepsize( 3 );

% image size
[ m, n, r ] = size( I );

% indices for each block's starting voxel
ds = ( 0 : d_slc : r - nslcs ) * m * n;
start_ind = reshape( bsxfun( @plus, ...
    ( bsxfun( @plus, ( 1 : d_row : m - nrows + 1 )', ( 0 : d_col : n - ncols ) * m ) ), ... % indices of the block starting voxel for the first slice
    reshape( ds, [ 1 1 length( ds ) ] ) ), [], 1 ); % all in one column: nblocks x 1

% expand with the block's remaining row voxels
% block row indices added to start of each block then flipped into depth dimension: nrows x 1 x nblocks
lin_row = permute( bsxfun( @plus, start_ind, ( 0 : nrows - 1 ) )', [ 1 3 2 ] ); 

% expand with the block's remaining column voxels
% block column indices added to start of each block: nrows x ncols x nblocks, then reshaped to (nrows x ncols) x nblocks
lidx_2D = reshape( bsxfun( @plus, lin_row, ( 0 : ncols - 1 ) * m ), nrows * ncols, [] ); 

% expand with the block's remaining slice voxels
% blocks flipped into depth dimension, then block slice indices added to each block: (nrows x ncols) x nslcs x nblocks
lidx_3D = bsxfun( @plus, permute( lidx_2D, [ 1 3 2 ] ), m * n * ( 0 : nslcs - 1 ) );

 % reshape into (nrows x ncols x nslcs ) x nblocks
lidx_2D_final = reshape( lidx_3D, [], size( lidx_2D, 2 ) );

% Get the corresponding blocks from the input image
out = I( lidx_2D_final );
