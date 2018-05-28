%-------------------------------------------------------------------------%
% Input:
% node: position of node        --> 1x2 vector
% sampleX: x value of samples   --> 1xN vector
% sampleY: y value of samples   --> 1xN vector
% uncert: uncertainty of whole image
%                     
% Output:
% w: weight of each samples     --> 1xN vector
%-------------------------------------------------------------------------%
function w = weight(node, sampleX, sampleY,uncert)

uFSample = uncert(sampleX, sampleY);    % uncertainty of samples
sFsample = (sampleX-node(1)).^2+(sampleY-node(2)).^2; % spatial distance between samples and node
w = (1-uFSample)*exp(-sFsample/(r2/2).^2);

end