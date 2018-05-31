%-------------------------------------------------------------------------%
% Input:
% Image: input RBG image              --> MxNx3 matrix;
% foreground: user defined foreground --> 2xn matrix,
%   where the first row is X coordinate, the seconde row is Y coordinate;
% background: user defined background --> 2xm matrix,
%   where the first row is X coordinate, the seconde row is Y coordinate;
%                     
% Output:
% MattedImage: Processd image         --> MxN matrix;
%-------------------------------------------------------------------------%
function MattedImage = BP_matte(Image, foreground, background)
%% initially define Uc & Un based on user-defined foreground & background
s = size(Image);
y = 1:s(1);
x = 1:s(2);
[X,Y] = meshgrid(x,y);
X = X.';    Y = Y.';

alphaK = linspace(0,1,25);                  % discrete alpha into 25 levels

alpha = 0.5*ones(s);                        % initialize alpha value, which is in the region [0,1];
uncert = ones(s);
for i1 = 1:length(foreground)
    uncert(foreground(1,i1),foreground(2,i1)) = 0;
end

for i1 = 1:length(background)
    uncert(background(1,i1),background(2,1i)) = 0;% initialize uncertainty;
end

Uc = 1-uncert;                              % initialize group Uc, which include all known pixels
Un = uncert;                                % initialize group Un, which include all unknown pixels
U = sum(uncert(:));
U_new = 0;
Uc_tilde = zeros(s);                        % initialize group Uc~

alpha(Uc==1) = 1;                           % set value of user-defined foreground to be 1
alpha(Un==1) = 0;                           % set value of user-defined background to be 0

%% build GMM based on user-defined foreground & background and assign pixels to single Gaussain
K = 5;                                      % the number of component is 5 as infered in paper "GrabCut"
gmm_fore = fitgmdist(foreground, K);        % build Gaussian mix model for foreground
gmm_back = fitgmdist(background, K);        % build Gaussian mix model for background

k_fore = round(K*rand(1,length(foreground)));   % assign each pixel in foreground to one single Gaussian
k_back = round(K*rand(1,length(background)));   % assign each pixel in background to one single Gaussian

r1 = 15;                                    % radius of region to be added into Uc_tilde
r2 = 20;                                    % radius of region considered when calculate weight

%%
% the loop judge condition needs to be further determined.
while max(Un)~=0 && U>U_new % while Un is not null
    U = sum(uncert(:));                              % update the total uncertainty.
    %% if Un is not null, transfer pixelswithin 15 pixels of Uc from Un to Uc_tilde
    if max(Un)~=0
        X_tilde = X(Uc==1);                 % obtain pixels in Uc
    	Y_tilde = Y(Uc==1);
        for i1 = 1:length(X_tilde)
           Uc_tilde((X-X_tilde(i1)).^2+(Y-Y_tilde(i1)).^2<=r1^2 & Un==1) = 1;
           Un((X-X_tilde(i1)).^2+(X-X_tilde(i1)).^2<=r1^2 & Un==1) = 0;
        end
    end
    
    %% build MRF based on pixels in Uc
    MRF = Uc_tilde;                         % treat each pixel in Uc_tilde as a node in the MRF
    X_tilde = X(Uc_tilde==1);
    Y_tilde = Y(Uc_tilde==1);
    for i1 = 1:length(Y_tilde)                % adjacent pixels in Uc are also in MRF
        if Uc(X_tilde(i1)-1,Y_tilde(i1))==1
           MRF(X_tilde(i1)-1,Y_tilde(i1)) = 1;
        end
        
        if Uc(X_tilde(i1)+1,Y_tilde(i1))==1
           MRF(X_tilde(i1)+1,Y_tilde(i1)) = 1;
        end
        
        if Uc(X_tilde(i1),Y_tilde(i1)-1)==1
           MRF(X_tilde(i1),Y_tilde(i1)-1) = 1;
        end
        
        if Uc(X_tilde(i1),Y_tilde(i1)+1)==1
           MRF(X_tilde(i1),Y_tilde(i1)+1) = 1;
        end
    end
    
    %% For each node in MRF, compute its data cost
    X_mrf = X(MRF==1);                      % obtain pixels in MRF
    Y_mrf = Y(MRF==1);
    
    N = 12;                                 % number of samples
    
    Vd = zeros([s,25]);                     % initialize data cost
    for i1 = 1:length(X_mrf)        
        %% sample foreground & background and calculate weight
        foreSample = zeros(2,N);            % foreground sample set coordinates
        backSample = zeros(2,N);            % background sample set coordinates
        
        
        % foreground
        foreSampleX_temp = X((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 & alpha==1);
        foreSampleY_temp = Y((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 & alpha==1);
        
        % check whether there are N samples
        if length(foreSampleX_temp)>=N
            % select the N samples with largest weight
            w_temp = weight([X_mrf(i1),Y_mrf(i1)], foreSampleX_temp, foreSampleY_temp, uncert);
            [wF,Fidx] = sort(w_temp,'descend');     % sort
            wF = wF(1:N);       % select N samples
            foreSample(1,:) = foreSampleX_temp(Fidx(1:N));
            foreSample(2,:) = foreSampleY_temp(Fidx(1:N));
        else
            % select N samples based on GMM above
            mu = gmm_fore.mu;
            sigma = gmm_fore.Sigma;
            for i2 = 1:N
                component = k_fore; % random select a single Gaussian
                muK = mu(component,:);
                sigmaK = sigma(component,:); % obtain parameter fo Gaussian
                r = mnvrnd(muK,sigmaK,1);
                [~,idx] = min(sum((foreground-r.').^2));
                foreSample(:,i2) = foreground(:,idx);
            end
            wF = weight([X_mrf(i1),Y_mrf(i1)],foreSample(1,:),foreSample(2,:),uncert);
        end
        
        % foreground
        backSampleX_temp = X((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 & alpha==1);
        backSampleY_temp = Y((X-X_mrf(i1)).^2+(Y-Y_mrf(i1)).^2<=r2^2 & alpha==1);
        
        % check whether there are N samples
        if length(backSampleX_temp)>=N
            % select the N samples with largest weight
            w_temp = weight([X_mrf(i1),Y_mrf(i1)], backSampleX_temp, backSampleY_temp, uncert);
            [wB,Bidx] = sort(w_temp,'descend');     % sort
            wB = wB(1:N);       % select N samples
            backSample(1,:) = backSampleX_temp(Bidx(1:N));
            backSample(2,:) = backSampleY_temp(Bidx(1:N));
        else
            % select N samples based on GMM above
            mu = gmm_back.mu;
            sigma = gmm_back.Sigma;
            for i2 = 1:N
                component = k_back(i2); % random select a single Gaussian
                muK = mu(component,:);
                sigmaK = sigma(component,:); % obtain parameter fo Gaussian
                r = mnvrnd(muK,sigmaK,1);
                [~,idx] = min(sum((background-r.').^2));
                backSample(:,i2) = background(:,idx);
            end
            wB = weight([X_mrf(i1),Y_mrf(i1)],backSample(1,:),backSample(2,:),uncert);
        end
        
       %% Compute likelihood
        Lk = zeros(1,25);
        
        idx = 1:N;
        [IdxF,IdxB] = meshgrid(idx, idx);
        
        Cp = Image(X_mrf(i1),Y_mrf(i1),:);   % obtain node's RGB value
        Cp = repmat(Cp,N,N);
        
        Fp = Image(foreSample(1,:),foreSample(2,:)); % obtain foreground samples' RGB value
        Fp_rep = repmat(Fp,N,1);
        
        Bp = Image(backSample(1,:),backSample(2,:)); % obtain background samples' RGB value
        Bp_rep = repmat(Bp.',1,N);
        
        dF = squeeze(sum(Fp.^2,3));                  % variance of foreground sample set
        sigmaF = cov(dF);
        dB = squeeze(sum(Bp.^2,3));                  % variance of background sample set
        sigmaB = cov(dB);
        
        for i2 = alphaK
            sigma = i2*sigmaF+(1-i2)*sigmaB;
            dc = sum((Cp-(i2*Fp_rep+(1-i2)*Bp_rep)).^2,3);
            
            Lk = wF(IdxF).*wB(IdxB).*exp(-dc/(2*sigma^2));
            Lk = 1/N^2*sum(Lk(:));                  % calculate likelihood
        end
        
        %% Compute data cost
        Vd(X_mrf(i1),Y_mrf(i1),:) = 1-Lk/sum(Lk);                          % obtain data cost
        
        
        
    end
    Vs = 1-exp(-(repmat(alphaK,N,1)-repmat(alphaK.',1,N)).^2/0.2^2);    % sigmaS = 0.2
    
    %% apply Belief Propagation algorithm
    alphaK = BP(MRF, Vd, Vs,alpha);
    
    %% update uncertianty, foreground & background
    uncert((alpha==1 |alpha==0)& Uc_tilde==1)=0;                      % assigning new foreground & background uncertainty to 0;
    Uc_tilde(alpha==1 |alpha==0) = 0;
    
    foreground_new(1,:) = X(alpha==1).';
    foreground_new(2,:) = Y(alpha==1).';
    background_new(1,:) = X(alpha==0).';
    background_new(2,:) = Y(alpha==0).';
    
    for i1 = 1:s(1)
       for i2 = 1:s(2)
           if(Uc_tilde==1)
              F_opt = Image(i1,i2)*alpha(i1,i2);
              B_opt = Image(i1,i2)*(1-alpha(i1,i2));
              
              foreValue = squeeze(Image(foreground_new(1,:),foreground_new(2,:)));  % obtain RGB value of foreground
              [~,minF] = min(sum((F_opt-foreValue).^3,3));                       % obtain the index of the smallest fitting error sample
              
              backValue = squeeze(Image(background_new(1,:),background_new(2,:)));  % obtain RGB value of backround
              [~,minB] = min(sum((B_opt-backValue).^3,3));                       % obtain the index of the smallest fitting error sample
              
              wF_star = weight([i1,i2],foreground_new(1,minF),foreground_new(2,minF),uncert);
              wB_star = weight([i1,i2],foreground_new(1,minB),foreground_new(2,minB),uncert);
              uncert(i1,i2)= 1-sqrt(wF_star, wB_star);
           end
       end
    end
    
    Uc_tilde =zeros(s); Uc_tilde(uncert~=0 & uncert~=1)=1;              % re-define Uc_tilde
    Un = zeros(s);  Un(uncert==1)=1;
    U_new = sum(uncert(:));                                              % compute total uncertainty
end

MattedImage = alpha*255;
end