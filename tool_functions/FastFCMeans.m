function [C, U, LUT, H] = FastFCMeans(im, c, q, opt)
    % Partition N-dimensional grayscale image into c classes using a memory
    % efficient implementation of the fuzzy c-means (FCM) clustering algorithm.
    % Computational efficiency is achieved by using the histogram of image
    % intensities during clustering instead of the raw image data.
    %
    % INPUT:
    %   - im  : N-dimensional grayscale image in integer format.
    %   - c   : positive integer greater than 1 specifying the number of
    %           clusters. c=2 is the default setting. Alternatively, c can be
    %           specified as a k-by-1 array of initial cluster (aka prototype)
    %           centroids.
    %   - q   : fuzzy weighting exponent. q must be a real number greater than
    %           1.1. q=2 is the default setting. Increasing q leads to an
    %           increased amount of fuzzification, while reducing q leads to
    %           crispier class memberships. Note that while in principle
    %           setting q==1 is equivalent to using a classical c-means
    %           algorithm, this setting cannot be used in practice because it
    %           produces an infinite exponent in the membership update formula.
    %   - opt : optional logical argument used to indicate how to initialize
    %           cluster centroids. If opt=true {default} then centroids are
    %           initialized are by sampling the intensity range at uniform
    %           intervals. If opt=false then the initial centroids are set
    %           using the c-means algorithm.
    %
    % OUTPUT  :
    %   - C   : 1-by-k array of cluster centroids.
    %   - U   : L-by-k array of fuzzy class memberships, where k is the number
    %           of classes and L is the intensity range of the input image,
    %           such that L=numel(min(im(:)):max(im(:))).
    %   - LUT : L-by-1 array that specifies the intensity-class relations,
    %           where L is the dynamic intensity range of the input image.
    %           Specifically, LUT(1) corresponds to the (class) label assigned to
    %           min(im(:)) and LUT(L) corresponds to the label assigned
    %           to max(im(:)). LUT is used as input to 'apply_LUT' function to
    %           create a label image.
    %   - H   : image histogram. If I=min(im(:)):max(im(:)) are the intensities
    %           present in the input image, then H(i) is the number of
    %           pixels/voxels with intensity I(i).

    % Default input arguments
    if nargin < 2 || isempty(c), c = 2; end
    if nargin < 3 || isempty(q), q = 2; end
    if nargin < 4 || isempty(opt), opt = true; end

    % Basic error checking
    if nargin < 1 || isempty(im)
        error('Insufficient number of input arguments')
    end

    msg = 'Revise variable used to specify class centroids. See function documentation for more info.';

    if ~isnumeric(c) || ~isvector(c)
        error(msg)
    end

    if numel(c) == 1 && (~isnumeric(c) || round(c) ~= c || c < 2)
        error(msg)
    end

    if ~isnumeric(q) || numel(q) ~= 1 || q < 1.1
        error('3rd input argument (q) must be a real number > 1.1')
    end

    if ~islogical(opt) || numel(opt) > 1
        error('4th input argument (opt) must a Boolean')
    end

    % Check image format
    if isempty(strfind(class(im), 'int'))
        error('Input image must be specified in integer format (e.g. uint8, int16)')
    end

    % Intensity range
    Imin = double(min(im(:)));
    Imax = double(max(im(:)));
    I = (Imin:Imax)';
    % Initialize cluster centroids
    if numel(c) > 1 % user-defined centroids
        C = c;
        opt = true;
    else % automatic initialization

        if opt
            dI = (Imax - Imin) / c;
            C = Imin + dI / 2:dI:Imax;
        else
            [C, ~, H] = FastCMeans(im, c);
        end

    end

    % Compute intensity histogram
    if opt
        H = hist(double(im(:)), I);
        H = H(:);
    end

    clear im

    % Update fuzzy memberships and cluster centroids
    a = 1 / (q - 1);
    dC = Inf;
    while dC > 1E-3

        C0 = C;

        % Distance to the centroids
        D = abs(bsxfun(@minus, I, C));
        D = D .^ (2 / (q - 1)) + eps;

        % Compute fuzzy memberships
        U = bsxfun(@times, D, sum(1 ./ D, 2));
        U = 1 ./ (U + eps);

        % Update the centroids
        UH = bsxfun(@times, U .^ q, H);
        C = sum(bsxfun(@times, UH, I), 1) ./ sum(UH, 1);
        C = sort(C, 'ascend'); % enforce natural order

        % Change in centroids
        dC = max(abs(C - C0));

    end

    % improve part 
    eta = [1 1 1 1];
    eta = ((ones(length(I), 1) * eta) .^ a);
    D_2 = power(abs(bsxfun(@minus, I, C)), 1);
    
    % Update U based on eta and C
    den = eta + D_2 .^ a;
    U = eta ./ den;

    % Defuzzify
    [~, LUT] = max(U, [], 2);
