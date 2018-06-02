function alphaK_new = BP(MRF, Vd, Vs,alphaK, level)
s = size(MRF);          % obtain the size of image
m = zeros([s,4,level]);     % initialize message
m_new = ones([s,4,level]);

x = 1:s(1); y = 1:s(2);
[X,Y] = meshgrid(x,y);
X = X.';    Y = Y.';

X_mrf = X(MRF==0);
Y_mrf = Y(MRF==0);
iteration = 20*length(X_mrf);         % iteration times
counter = 1;
% while sum(m_new(:) == m(:))~=numel(m) & counter < iteration
while counter < iteration
    m = m_new;
    idx = ceil(length(X_mrf)*rand(1));
    node = [X_mrf(idx),Y_mrf(idx)];
    m_new = mesCompute(m, Vd, Vs, node);
    counter = counter+1;
    
end

fprintf('%d iteration for BP algorithm\n',counter);

% compute belief vector, b should be a 3D matrix
b = Vd.*squeeze(sum(m_new,3));  
% compute alpha level, alphaK_level should be a 2D matrix
[~,alphaK_level] = max(b,[],3);
% compute new alpha matrix
alphaK_new = alphaK(alphaK_level);
end

