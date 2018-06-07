%-------------------------------------------------------------------------%
% Input:
% Image: input RBG image              --> MxNx3 matrix;
% foreground: user defined foreground --> 2xn matrix,
%   where the first row is X coordinate, the seconde row is Y coordinate;
% background: user defined background --> 2xm matrix,
%   where the first row is X coordinate, the seconde row is Y coordinate;
% K: number of component in GMM the number of component is 5 as infered in paper "GrabCut"               
% Output:
% MattedImage: Processd image         --> MxN matrix;
%-------------------------------------------------------------------------%
function MattedImage = BP_matte(Image, foreground, background,K)
%% initially define Uc & Un based on user-defined foreground & background
s = size(Image);
s = s([1,2]);
y = 1:s(1);
x = 1:s(2);
[X,Y] = meshgrid(y,x);
X = X.';    Y = Y.';

level = 25;                                 % number of alpha level
alphaK = linspace(0,1,level);                  % discrete alpha into 25 levels

alpha = 0.5*ones(s);                        % initialize alpha value, which is in the region [0,1];
uncert = ones(s);
for i1 = 1:length(foreground)
    uncert(foreground(1,i1),foreground(2,i1)) = 0;
end

for i1 = 1:length(background)
    uncert(background(1,i1),background(2,i1)) = 0;% initialize uncertainty;
end

Un = uncert;                                % initialize group Un, which include all unknown pixels
U = sum(uncert(:));
U_new = 0;
Uc_tilde = zeros(s);                        % initialize group Uc~

estimated = zeros(s);

for i1 = 1:length(foreground)
    estimated(foreground(1,i1),foreground(2,i1)) = 1;                       % measure the pixel have been estimted or not
    alpha(foreground(1,i1),foreground(2,i1)) = 1;                           % set value of user-defined foreground to be 1
end

for i1 = 1:length(background)
    estimated(background(1,i1),background(2,i1)) = 1;                       % measure the pixel have been estimted or not
    alpha(background(1,i1),background(2,i1)) = 0;                           % set value of user-defined foreground to be 1
end

counter = 1;

%% build GMM based on user-defined foreground & background and assign pixels to single Gaussain
options = statset('MaxIter',500);
gmm_fore = fitgmdist(foreground.', K,'Options',options);        % build Gaussian mix model for foreground
gmm_back = fitgmdist(background.', K,'Options',options);        % build Gaussian mix model for background

k_fore = ceil(K*rand(1,length(foreground)));   % assign each pixel in foreground to one single Gaussian
k_back = ceil(K*rand(1,length(background)));   % assign each pixel in background to one single Gaussian

r1 = 15;                                    % radius of region to be added into Uc_tilde
r2 = 20;                                    % radius of region considered when calculate weight

%%
% the loop judge condition needs to be further determined.
while max(Un(:))~=0 & U>U_new % while Un is not null
    U = sum(uncert(:));                              % update the total uncertainty.
    %% if Un is not null, transfer pixelswithin 15 pixels of Uc from Un to Uc_tilde
    if max(Un(:))~=0
        X_c = X(uncert==0);                 % obtain pixels in Uc
    	Y_c = Y(uncert==0);
        for i1 = 1:length(X_c)
           Uc_tilde((X-X_c(i1)).^2+(Y-Y_c(i1)).^2<=r1^2 & Un==1) = 1;
           Un((X-X_c(i1)).^2+(Y-Y_c(i1)).^2<=r1^2 & Un==1) = 0;
        end
    end
    
    %% build MRF based on pixels in Uc
    MRF = Uc_tilde;                         % treat each pixel in Uc_tilde as a node in the MRF
    X_tilde = X(Uc_tilde==1);
    Y_tilde = Y(Uc_tilde==1);
    for i1 = 1:length(Y_tilde)                % adjacent pixels in Uc are also in MRF
        if X_tilde(i1)>1
            if  uncert(X_tilde(i1)-1,Y_tilde(i1))==0
                MRF(X_tilde(i1)-1,Y_tilde(i1)) = 1;
            end
        end
        
        if X_tilde(i1)<s(1)
            if uncert(X_tilde(i1)+1,Y_tilde(i1))==0
                MRF(X_tilde(i1)+1,Y_tilde(i1)) = 1;
            end
        end
        
        if Y_tilde(i1)>1
            if uncert(X_tilde(i1),Y_tilde(i1)-1)==0
                MRF(X_tilde(i1),Y_tilde(i1)-1) = 1;
            end
        end
        
        if Y_tilde(i1)<s(2)
            if uncert(X_tilde(i1),Y_tilde(i1)+1)==0
                MRF(X_tilde(i1),Y_tilde(i1)+1) = 1;
            end
        end
    end
    
    %% For each node in MRF, compute its data cost
    X_mrf = X(MRF==1);                      % obtain pixels in MRF
    Y_mrf = Y(MRF==1);
    
    N = 12;                                 % number of samples
    
    Vd = ones([s,level]);                     % initialize data cost
    for i1 = 1:length(X_mrf)        
        %% sample foreground & background and calculate weight
        foreSample = zeros(2,N);            % foreground sample set coordinates
        backSample = zeros(2,N);            % background sample set coordinates
        
        
        % foreground
        foreSampleX_temp = X((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 ...
            & alpha>alpha(X_mrf(i1),Y_mrf(i1)) & estimated==1);
        foreSampleY_temp = Y((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 ...
            & alpha>alpha(X_mrf(i1),Y_mrf(i1)) & estimated==1);
        
        % check whether there are N samples
        if length(foreSampleX_temp)>=N
            % select the N samples with largest weight
            w_temp = weight([X_mrf(i1),Y_mrf(i1)], foreSampleX_temp, foreSampleY_temp, uncert, r2);
            [wF,Fidx] = sort(w_temp,'descend');     % sort
            wF = wF(1:N);       % select N samples
            foreSample(1,:) = foreSampleX_temp(Fidx(1:N));
            foreSample(2,:) = foreSampleY_temp(Fidx(1:N));
        else
            % select N samples based on GMM above
            mu = gmm_fore.mu;
            sigma = gmm_fore.Sigma;
            for i2 = 1:N
                component = k_fore(ceil(length(k_fore)*rand())); % random select a single Gaussian
                muK = mu(component,:);
                sigmaK = sigma(:,:,component); % obtain parameter fo Gaussian
                r = mvnrnd(muK,sigmaK,1);
                [~,idx] = min(sum((foreground-r.').^2));    % find the most closest foreground as sample
                foreSample(:,i2) = foreground(:,idx);
            end
            wF = weight([X_mrf(i1),Y_mrf(i1)],foreSample(1,:).',foreSample(2,:).',uncert, r2);
        end
        
        % background
        backSampleX_temp = X((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 ...
            & alpha<alpha(X_mrf(i1),Y_mrf(i1)) & estimated==1);
        backSampleY_temp = Y((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 ...
            & alpha<alpha(X_mrf(i1),Y_mrf(i1)) & estimated==1);
        
        % check whether there are N samples
        if length(backSampleX_temp)>=N
            % select the N samples with largest weight
            w_temp = weight([X_mrf(i1),Y_mrf(i1)], backSampleX_temp, backSampleY_temp, uncert, r2);
            [wB,Bidx] = sort(w_temp,'descend');     % sort
            wB = wB(1:N);       % select N samples
            backSample(1,:) = backSampleX_temp(Bidx(1:N));
            backSample(2,:) = backSampleY_temp(Bidx(1:N));
        else
            % select N samples based on GMM above
            mu = gmm_back.mu;
            sigma = gmm_back.Sigma;
            for i2 = 1:N
                component = k_back(ceil(length(k_back)*rand())); % random select a single Gaussian
                muK = mu(component,:);
                sigmaK = sigma(:,:,component); % obtain parameter fo Gaussian
                r = mvnrnd(muK,sigmaK,1);
                [~,idx] = min(sum((background-r.').^2)); % find the most closest background as sample
                backSample(:,i2) = background(:,idx);
            end
            wB = weight([X_mrf(i1),Y_mrf(i1)],backSample(1,:).',backSample(2,:).',uncert, r2);
        end
        
       %% Compute likelihood
        Lk = zeros(1,level);
        
        idx = 1:N;
        [IdxF,IdxB] = meshgrid(idx, idx);
        
        Cp = Image(X_mrf(i1),Y_mrf(i1),:);   % obtain node's RGB value
        Cp = repmat(Cp,N,N);
        
        Fp = zeros(1,length(foreSample),3); % 
        for i2 = 1:length(foreSample)
            Fp(1,i2,:) = Image(foreSample(1,i2),foreSample(2,i2)); % obtain foreground samples' RGB value
        end
        Fp_rep = repmat(Fp,N,1);
        
        Bp = zeros(length(backSample),1,3);
        for i2 = 1:length(backSample)
            Bp(i2,1,:) = Image(backSample(1,i2),backSample(2,i2)); % obtain background samples' RGB value
        end
        Bp_rep = repmat(Bp,1,N);
        
        dF = squeeze(sqrt(sum(Fp.^2,3)));                  % variance of foreground sample set
        sigmaF = mean((dF-mean(dF)).^2);
        dB = squeeze(sqrt(sum(Bp.^2,3)));                  % variance of background sample set
        sigmaB = mean((dB-mean(dB)).^2);
        
        for i2 = 1:level
            sigma = alphaK(i2)*sigmaF+(1-alphaK(i2))*sigmaB;
            dc = sum((Cp-(alphaK(i2)*Fp_rep+(1-alphaK(i2))*Bp_rep)).^2,3);
            
            Lkp = wF(IdxF).*wB(IdxB).*exp(-dc/(2*sigma^2));
            Lk(i2) = 1/N^2*sum(Lkp(:));                     % calculate likelihood
        end
        
        %% Compute data cost
        Lks = sum(Lk);
        if Lks ~= 0
            Vd(X_mrf(i1),Y_mrf(i1),:) = 1-Lk/sum(Lk);                          % obtain data cost
        else
            Vd(X_mrf(i1),Y_mrf(i1),:) = 1/level;
        end
        
    end
    Vs = 1-exp(-(repmat(alphaK,level,1)-repmat(alphaK.',1,level)).^2/0.2^2);    % sigmaS = 0.2
    
    %% apply Belief Propagation algorithm
    alpha_new = BP(MRF, Vd, Vs,alphaK,level);
    alpha(MRF==1)=alpha_new(MRF==1);
    
    for i1 = 1:length(foreground)
        alpha(foreground(1,i1),foreground(2,i1)) = 1;                           % set value of user-defined foreground to be 1
    end

    for i1 = 1:length(background)
        alpha(background(1,i1),background(2,i1)) = 0;                           % set value of user-defined foreground to be 1
    end
    
    %% update uncertianty, foreground & background
    uncert((alpha==1 |alpha==0) & Uc_tilde==1)=0;                      % assigning new foreground & background uncertainty to 0;
    Uc_tilde(alpha==1 |alpha==0) = 0;
    
    foreground_new = zeros(2,length(X(alpha==1)));
    background_new = zeros(2,length(X(alpha==0)));
    
    foreground_new(1,:) = X(alpha==1).';
    foreground_new(2,:) = Y(alpha==1).';
    background_new(1,:) = X(alpha==0).';
    background_new(2,:) = Y(alpha==0).';
    
    F_opt = Image.*alpha;
    B_opt = Image-F_opt;
    
    foreValue = zeros(3,length(foreground_new));
    for i1 = 1:length(foreground_new)
        foreValue(:,i1) = squeeze(Image(foreground_new(1,i1),foreground_new(2,i1),:));  % obtain RGB value of foreground
    end
    
    backValue = zeros(3,length(background_new));
    for i1 = 1:length(background_new)
        backValue(:,i1) = squeeze(Image(background_new(1,i1),background_new(2,i1),:));  % obtain RGB value of backround
    end
    
    for i1 = 1:s(1)
       for i2 = 1:s(2)
           if(Uc_tilde(i1,i2)==1)
              
              [~,minF] = min(sum((squeeze(F_opt(i1,i2,:))-foreValue).^2),[],2);                       % obtain the index of the smallest fitting error sample
              [~,minB] = min(sum((squeeze(B_opt(i1,i2,:))-backValue).^2),[],2);                       % obtain the index of the smallest fitting error sample
              
              wF_star = weight([i1,i2],foreground_new(1,minF),foreground_new(2,minF),uncert, r2);
              wB_star = weight([i1,i2],background_new(1,minB),background_new(2,minB),uncert, r2);
              uncert(i1,i2)= 1-sqrt(wF_star.*wB_star);
           end
       end
    end
    
%     Uc_tilde =zeros(s); Uc_tilde(uncert~=0 & uncert~=1)=1;              % re-define Uc_tilde
    Un = zeros(s);  Un(uncert==1)=1;
    U_new = sum(uncert(:));                                              % compute total uncertainty
    
    estimated(MRF==1)=1;                                                % record pixels which have been estimated
%     fprintf('the %d iteration',counter);



    counter = counter+1;
end

MattedImage = alpha*255;
end