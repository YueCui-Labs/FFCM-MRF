function [C, LUT, H] = FastCMeans(im, c)
    % Partition N-dimensional grayscale image into c classes using a memory
    % efficient implementation of the c-means (aka k-means) clustering
    % algorithm. Computational efficiency is achieved by using the histogram
    % of image intensities during clustering instead of the raw image data.
    %
    % INPUT:
    %   - im  : N-dimensional grayscale image in integer format.
    %   - c   : positive integer greater than 1 specifying the number of
    %           clusters. c=2 is the default setting. Alternatively, c can be
    %           specified as a k-by-1 array of initial cluster (aka prototype)
    %           centroids.
    %
    % OUTPUT  :
    %   - C   : 1-by-k array of cluster centroids.
    %   - LUT : L-by-1 array that specifies the intensity-class relations,
    %           where L is the dynamic intensity range of the input image.
    %           Specifically, LUT(1) corresponds to the (class) label assigned to
    %           min(im(:)) and LUT(L) corresponds to the label assigned
    %           to max(im(:)). LUT is used as input to 'apply_LUT' function to
    %           create a label image.
    %   - H   : image histogram. If I=min(im(:)):max(im(:)) are the intensities
    %           present in the input image, then H(i) is the number of
    %           pixels/voxels with intensity I(i).
    %
    % AUTHOR    : Anton Semechko (a.semechko@gmail.com)
    %

    % Default input arguments
    if nargin < 2 || isempty(c), c = 2; end

    % Basic error checking
    if nargin < 1 || isempty(im)
        error('Insufficient number of input arguments')
    end

    msg = 'Revise variable used to specify class centroids. See function documentaion for more info.';

    if ~isnumeric(c) || ~isvector(c)
        error(msg)
    end

    if numel(c) == 1 && (~isnumeric(c) || round(c) ~= c || c < 2)
        error(msg)
    end

    % Check image format
    if isempty(strfind(class(im), 'int'))
        error('Input image must be specified in integer format (e.g. uint8, int16)')
    end

    % Intensity range
    Imin = double(min(im(:)));
    Imax = double(max(im(:)));
    I = (Imin:Imax)';

    % Compute intensity histogram
    H = hist(double(im(:)), I);
    H = H(:);

    % [~, Imax] = max(H);
    % c = [Imax / 4 Imax 4 * Imax 5 * Imax];

    % Initialize cluster centroids
    if numel(c) > 1 % uder-defined
        C = c;
        c = numel(c);
    else % automatic
        dI = (Imax - Imin) / c;
        C = Imin + dI / 2:dI:Imax;
    end

    % Update cluster centroids
    IH = I .* H;
    dC = Inf;

    while dC > 1E-3

        C0 = C;

        % Distance to centroids
        D = abs(bsxfun(@minus, I, C));
        D = D .^ 2;
        % Classify by proximity
        [~, LUT] = min(D, [], 2);

        for j = 1:c
            C(j) = sum(IH(LUT == j)) / sum(H(LUT == j));
        end

        C = sort(C, 'ascend'); % enforce natural order

        % Change in centroids
        dC = max(abs(C - C0));

    end
