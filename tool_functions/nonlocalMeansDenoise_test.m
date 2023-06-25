function denoised=nonlocalMeansDenoise_test( dim,image )
% synthetic example by Jose Vicente Manjon-Herrera extended to a 3D image

if dim == 2
    % Create synthetic 2D image
    image = repelem( [ 100 200; 50 100 ], 50, 50 );

    % Add some noise
    sigma = 10;
    noisy = image + sigma * randn( size( image ) );

    % Denoising parameters
    search_rad = 3;
    patch_rad = 2;
    patch_sigma = 1;
    search_sigma = 10;
    selfsim = 0;
    
    tic
    denoised = nonlocalMeansDenoise( noisy, patch_rad, search_rad, patch_sigma, search_sigma, selfsim );
    toc

    figure,
    subplot( 2, 2, 1 ), imagesc( image ), title( 'original' );
    subplot( 2, 2, 2 ), imagesc( noisy ), title( 'noisy' );
    subplot( 2, 2, 3 ), imagesc( denoised ), title( 'filtered' );
    subplot( 2, 2, 4 ), imagesc( noisy - denoised), title( 'residuals' );
    colormap( gray )
    
elseif dim == 3
    % Create synthetic 3D image
%     image = cat( 3, [ 100 200; 50 100 ], [ 100 50; 200 100 ] );
%     image = repelem( image, 96, 96, 96 );
% 
%     % Add some noise
%     sigma = 10;
%     noisy = image + sigma * randn( size( image ) );
    
    
    
    % Denoising parameters
    search_rad = 1;
    patch_rad = 1;
    patch_sigma = 1;
    search_sigma = 10;
    selfsim = 0;

    tic
    denoised = nonlocalMeansDenoise( image, patch_rad, search_rad, patch_sigma, search_sigma, selfsim );
    toc

%     figure,
%     subplot( 2, 2, 1 ), imagesc( image( :, :, 25 ) ), title( 'original' );
%     subplot( 2, 2, 2 ), imagesc( noisy( :, :, 25 ) ), title( 'noisy' );
%     subplot( 2, 2, 3 ), imagesc( denoised( :, :, 25 ) ), title( 'filtered' );
%     subplot( 2, 2, 4 ), imagesc( noisy( :, :, 25 ) - denoised( :, :, 25 ) ), title( 'residuals' );
%     colormap( gray )

else
    error( 'Input must be either 2 or 3' );
end

mse = norm( image(:) - denoised(:), 'fro' ) / numel( image );
fprintf( 'Frobenius norm of the residuals: %.3f\n', mse );

