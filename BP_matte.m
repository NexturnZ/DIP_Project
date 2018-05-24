%-------------------------------------------------------------------------%
% Input:
% Image:      input image             --> M-by-N matrix;
% foreground: user defined foreground --> 2-by-n matrix,
%   where the first row is X coordinate, the seconde row is Y coordinate;
% background: user defined background --> 2-by-m matrix,
%   where the first row is X coordinate, the seconde row is Y coordinate;
%                     
% Output:
% MattedImage: Processd image         --> M-by-N matrix;
%-------------------------------------------------------------------------%
function MattedImage = BP_matte(Image, foreground, background)
%% initially define Uc & Un based on user-defined foreground & background
s = size(Image);
y = 1:s(1);
x = 1:s(2);
[X,Y] = meshgrid(x,y);

alpha = 0.5*ones(s);                        % initialize alpha value, which is in the region [0,1];
uncert = ones(s);
uncert(foreground(1,:),foreground(2,:)) = 0;
uncert(background(1,:),background(2,:)) = 0;% initialize uncertainty;

Uc = 1-uncert;                              % initialize group Uc, which include all known pixels
Un = uncert;                                % initialize group Un, which include all unknown pixels
Uc_tilde = zeros(s);                        % initialize group Uc~

alpha(Uc==1) = 1;                           % set value of user-defined background to be 0
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
           Uc_tilde((X-X_tilde(i1)).^2+(Y-Y_tilde(i1)).^2<r1^2 && Un==1) = 1;
           Un((X-X_tilde(i1)).^2+(X-X_tilde(i1)).^2<r1^2 && Un==1) = 0;
        end
    end
    
    %% build MRF based on pixels in Uc
    MRF = Uc_tilde;                         % treat each pixel in Uc_tilde as a node in the MRF
    X_mrf = X(Uc_tilde==1);
    Y_mrf = Y(Uc_tilde==1);
    for i1 = 1:length(Y_mrf)                % adjacent pixels in Uc are also in MRF
        if Uc(X_mrf(i1)-1,Y_mrf(i1))==1
           MRF(X_mrf(i1)-1,Y_mrf(i1)) = 1;
        end
        
        if Uc(X_mrf(i1)+1,Y_mrf(i1))==1
           MRF(X_mrf(i1)+1,Y_mrf(i1)) = 1;
        end
        
        if Uc(X_mrf(i1),Y_mrf(i1)-1)==1
           MRF(X_mrf(i1),Y_mrf(i1)-1) = 1;
        end
        
        if Uc(X_mrf(i1),Y_mrf(i1)+1)==1
           MRF(X_mrf(i1),Y_mrf(i1)+1) = 1;
        end
    end
    
    %% For each node in MRF, compute its data cost
    
end


end