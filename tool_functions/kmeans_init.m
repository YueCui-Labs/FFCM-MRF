function [C_init] = kmeans_init(IMG,K,dataset,sub)
    disp('FCM kmeans init...')
    img = IMG(IMG>0);
    [N,X] = hist(img(:),0:max(img(:)));
    [~,Imax] = max(N);
    
    % K-Means core function
    tic;
    disp(['Imax is:',num2str(Imax)]);

    if K>3
        
        [idx,ctrls] = kmeans(img(:),K,'Start',[Imax/4;Imax;2*Imax;4*Imax]);
        [idx,ctrls] = kmeans(img(:),K);
    else
        [idx,ctrls] = kmeans(img(:),K,'Start',[Imax/4;Imax;2*Imax]);
    end
    [Idx,Ctrls] = kmeans_reorder(idx,ctrls);
    t=toc;
    disp(['Center after kmeans are: ',num2str(Ctrls')]);
    disp(['FCM K-Means using ' num2str(t) ' seconds']);
    
    C_init = Ctrls';
end

function [Idx,Ctrs] = kmeans_reorder(idx,ctrs)
K = length(ctrs);
Ctrs_index = [ctrs,(1:length(ctrs))'];
sort_Ctrs = sortrows(Ctrs_index,1);
K_index = cell(1,K);
for k = 1:K 
    K_index{k} = find(idx==k);
end
Idx = zeros(size(idx));
Ctrs = zeros(size(ctrs));
for i = 1:K
    Ctrs(i) = sort_Ctrs(i,1);
    Idx(K_index{sort_Ctrs(i,2)}) = i;
end
end