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
s = size(Image);
y = 1:s(1);
x = 1:s(2);
[X,Y] = meshgrid(x,y);
Uc = zeros(s);                              % initialize group Uc, which include all known pixels
Un = ones(s);                               % initialize group Un, which include all unknown pixels
alpha = 0.5*ones(s);                        % initialize alpha value, which is in the region [0,1];

Uc(foreground(1,:),foreground(2,:)) = 1;
Uc(background(1,:),background(2,:)) = 1;    % put user-defined foreground & background into Uc;
Un(Uc==1) = 0;                              % put other pixels into Un;
alpha(Uc==1) = 1;                           % set value of user-defined background to be 0
alpha(Un==1) = 0;                           % set value of user-defined background to be 0



end