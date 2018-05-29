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
uncert(foreground(1,:),foreground(2,:)) = 0;
uncert(background(1,:),background(2,:)) = 0;% initialize uncertainty;

Uc = 1-uncert;                              % initialize group Uc, which include all known pixels
Un = uncert;                                % initialize group Un, which include all unknown pixels
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
while max(Un)~=0                            % while Un is not null
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
                component = round(K*rand()); % random select a single Gaussian
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
                component = round(K*rand()); % random select a single Gaussian
                muK = mu(component,:);
                sigmaK = sigma(component,:); % obtain parameter fo Gaussian
                r = mnvrnd(muK,sigmaK,1);
                [~,idx] = min(sum((background-r.').^2));
                backSample(:,i2) = background(:,idx);
            end
            wB = weight([X_mrf(i1),Y_mrf(i1)],backSample(1,:),backSample(2,:),uncert);
        end
        
        %% Compute data cost
        idx = 1:N;
        [IdxF,IdxB] = meshgrid(idx, idx);
        
        Cp = squeeze(Image(X_mrf(i1),Y_mrf(i1),:));   % obtain node's RGB value
        Cp = repmat(Cp,1,N);
        
        FP = squeeze(Image(foreSample(1,:),foreSample(2,:))); % obtain foreground samples' RGB value
        FP = FP.';
        
        Bp = squeeze(Image(backSample(1,:),backSample(2,:))); % obtain background samples' RGB value
        Bp = Bp.';
        
        for i2 = alphaK
            
            dc = Cp-(i2*Fp+(1-i2)*Bp)
            Lk = wF(IdxF).*wB(IdxB);
            Lk = 1/N^2*sum(Lk(:))*exp;                  % calculate likelihood
        end
        
        
    end
    
end


end