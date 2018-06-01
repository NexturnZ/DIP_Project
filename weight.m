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
function w = weight(node, sampleX, sampleY,uncert, r2)

uFSample = zeros(1,length(sampleX));
for i1 = 1:length(sampleX)
    uFSample(i1) = uncert(sampleX(i1), sampleY(i1));    % uncertainty of samples
end

sFsample = (sampleX-node(1)).^2+(sampleY-node(2)).^2; % spatial distance between samples and node
w = (1-uFSample).*exp(-sFsample.'/(r2/2).^2);

end