function [MRF_seg] = MRF(FCM_seg, IMG, vesselness, U, beta, W, IL)
    % Vessel segmentation refinement based on Markov Random Field
    % for FCM clustering result.

    % INPUT:
    %   - FCM_seg  : binary segmentation of FCM clustering algorithm.
    %   - IMG      : TOF-MRA image after histogram specification.
    %   - vesselness   : the result of multi-scalar vessel enhancement without skull
    %   - U        : L-by-k array of fuzzy class memberships, where k is the number
    %                of classes and L is the intensity range of the input image.
    %   - beta     : doubleton potential β，As β increases, regions become more homogenous.
    %   - W        : ratio of each class derived from FCM clustering result.
    %   - IL       : iterations of MRF.
    %
    % OUTPUT:
    %   - MRF_seg  : the final segmentation result.

    Iout = vesselness;
    Mask = double(IMG > 0);
    sortIout = sort(Iout(:), 'desc');
    vessel_ratio = sum(FCM_seg(:)) / sum(IMG(:) > 0);

    % stage 1, MRF with 26 neighborhood to remove most false positive, mainly strip artifacts
%     thresh = sortIout(round(numel(IMG(IMG > 0)) * vessel_ratio)); % thresh = T
%     disp('MRF stage 1 ...')
%     disp(['thresh (T) = ', num2str(thresh)]);
% 
%     tmp_IMG = double(Iout > (0.5 * thresh));
%     tmp_IMG = tmp_IMG .* Mask;
% 
%     flag_IMG = (FCM_seg + tmp_IMG) > 0; % S = S1+S2
%     flag_IMG = double(flag_IMG);
% 
%     Ni = find(flag_IMG ~= 0); % S = S1+S2
%     clear flag_IMG tmp_IMG
%     Dx_init = FCM_seg;
%     NB = 26;
% 
%     for t = 1:IL
%         Dx_ratio = length(find(Dx_init == 1)) / numel(IMG(IMG > 0));
%         T = sortIout(round(numel(IMG(IMG > 0)) * Dx_ratio));
%         Dx = ICM_estimation(IMG, U, Dx_init, Ni, T, Iout, beta, NB, W);
%         disp(['Dx ratio is ' num2str(Dx_ratio) '; T is ' num2str(T) '; ICM iteration ' num2str(t) ' times...']);
%         Dx_init = Dx;
%     end

    % stage 2, MRF with 6 neighborhood to remove false negative, mainly thin vessels
%     vessel_ratio = sum(Dx_init(:)) / sum(IMG(:) > 0);%使用26邻域时解开注释
    thresh = sortIout(round(numel(IMG(IMG > 0)) * vessel_ratio)); % thresh = T
    disp('MRF stage ...')
%     disp(['thresh (T) = ', num2str(thresh)]);

    tmp_IMG2 = double(Iout > (0.5 * thresh));
    tmp_IMG2 = tmp_IMG2 .* Mask;

    flag_IMG2 = (FCM_seg + tmp_IMG2) > 0; % S = S1+S2：仅使用6邻域时解开注释
%     flag_IMG2 = (Dx_init + tmp_IMG2) > 0; % S = S3+S4
    flag_IMG2 = double(flag_IMG2);

    Ni2 = find(flag_IMG2 ~= 0); % S = S3+S4
    NB = 6;
    clear flag_IMG2 tmp_IMG2
    Dx_init = FCM_seg;% 只用6邻域时解开注释

    for t = 1:IL
        Dx_ratio = length(find(Dx_init == 1)) / numel(IMG(IMG > 0));
        T = sortIout(round(numel(IMG(IMG > 0)) * Dx_ratio));
%         disp(['Dx ratio is ' num2str(Dx_ratio) '; T is ' num2str(T) '; ICM iteration ' num2str(t) ' times...']);
        Dx = ICM_estimation(IMG, U, Dx_init, Ni2, T, Iout, beta, NB, W);
        Dx_init = Dx;
    end

    MRF_seg = Dx;
    disp('----------- All FINISHED -------------------------');
end


function [Dnew] = ICM_estimation(IMG, U, D, Ni, T, Iout, beta, NB, W)

    sizeD = size(D);
    Dnew = zeros(sizeD);
    [a, b, c] = ind2sub(sizeD, Ni);
    s = [a b c];
    FN = (NB == 6) * 1 + (NB == 18) * 2 + (NB == 26) * 1;
    f = {@clique @clique_2};
    Hessian_map = Iout;
    Imin = min(int32((IMG(IMG > 0))));

    for n = 1:length(Ni)
        [pB, pV] = f{FN}(T, s(n, :), D, beta, Hessian_map, NB); % get energy function
        idx = int32(IMG(s(n, 1), s(n, 2), s(n, 3))) - Imin + 1;
        u = U(idx, :);
        u_vr = (u(3) * W(3) + u(4) * W(4)) / (W(3) + W(4));
        u_br = (u(1) * W(1) + u(2) * W(2)) / (W(1) + W(2));
        post_V = pV * u_vr;
        post_B = pB * u_br;
        Dnew(Ni(n)) = 1 * (post_V > post_B);
    end

end


function [pB, pV, NsigmaV, NsigmaB] = clique(T, s, D, beta, H, NB)

    [i_max, j_max, n_max] = size(D);
    i = s(1); j = s(2); n = s(3);
    ip = (i + 1 <= i_max) * (i + 1) + (i + 1 > i_max); im = (i - 1 >= 1) * (i - 1) + (i - 1 < 1) * i_max; %----行加和减
    jp = (j + 1 <= j_max) * (j + 1) + (j + 1 > j_max); jm = (j - 1 >= 1) * (j - 1) + (j - 1 < 1) * j_max; %----列加和减
    np = (n + 1 <= n_max) * (n + 1) + (n + 1 > n_max); nm = (n - 1 >= 1) * (n - 1) + (n - 1 < 1) * n_max; %----层加和减
    % define the matrix structure of the spatial neighborhood

    D_nb26 = [D(im, jm, nm) D(i, jm, nm) D(ip, jm, nm) D(im, j, nm) D(i, j, nm) D(ip, j, nm) D(im, jp, nm) D(i, jp, nm) D(ip, jp, nm) ... %3x3x3立方�? -1层标记�?
              D(im, jm, n) D(i, jm, n) D(ip, jm, n) D(im, j, n) D(i, j, n) D(ip, j, n) D(im, jp, n) D(i, jp, n) D(ip, jp, n) ... %3x3x3立方�?  0层标记�?
                  D(im, jm, np) D(i, jm, np) D(ip, jm, np) D(im, j, np) D(i, j, np) D(ip, j, np) D(im, jp, np) D(i, jp, np) D(ip, jp, np)]; %3x3x3立方�? +1层标记�?
    Hessian_nb26 = [H(im, jm, nm) H(i, jm, nm) H(ip, jm, nm) H(im, j, nm) H(i, j, nm) H(ip, j, nm) H(im, jp, nm) H(i, jp, nm) H(ip, jp, nm) ... %3x3x3立方�? -1层标记�?
                        H(im, jm, n) H(i, jm, n) H(ip, jm, n) H(im, j, n) H(i, j, n) H(ip, j, n) H(im, jp, n) H(i, jp, n) H(ip, jp, n) ... %3x3x3立方�?  0层标记�?
                        H(im, jm, np) H(i, jm, np) H(ip, jm, np) H(im, j, np) H(i, j, np) H(ip, j, np) H(im, jp, np) H(i, jp, np) H(ip, jp, np)]; %3x3x3立方�? +1层标记�?

    if NB == 6
        A = [5 11 13 15 17 23];
    else
        A = [1:13 15:27];
    end

    NsigmaV = sum(D_nb26(A) == 1);
    NsigmaB = sum(D_nb26(A) ~= 1);

    Uv6_1 = NB - NsigmaV; %
    Ub6_1 = NB - NsigmaB; %
    Uv6_2 = NB - sum(Hessian_nb26(A));
    Ub6_2 = (1 / (T) - 1) * sum(Hessian_nb26(A));

    alpha1 = 0.5;
    alpha2 = 0.5;
    Uv6 = alpha1 * Uv6_1 + alpha2 * Uv6_2;
    Ub6 = alpha1 * Ub6_1 + alpha2 * Ub6_2;

    Uv = beta * Uv6;
    pV = exp(-Uv);
    Ub = beta * Ub6;
    pB = exp(-Ub);
end
