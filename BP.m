function alphaK_new = BP(MRF, Vd, Vs,alpha, level)
iteration = 30;         % iteration times
s = size(MRF);          % obtain the size of image
m = zeros([s,4,level]);     % initialize message
m_new = ones([s,4,level]);
counter = 1;

x = 1:s(2); y = 1:s(1);
[X,Y] = meshgrid(x,y);
X = X.';    Y = Y.';

while sum(m_new(:) == m(:))~=numel(m) & counter < iteration
    m = m_new;
    
    for i1 = 1:level
        Vs_rep = ones(1,1,level);
        Vs_rep(1,1,:)= Vs(i1,:);
        Vs_rep = repmat(Vs_rep,s(1),s(2),1);

        m_up = squeeze(prod(m(:,:,[2,3,4],:),3));    % MxNx25 matrix
        m_up = max(Vd.*m_up.*Vs_rep,[],4);
        m_new(1:s(1)-1,:,1,i1) = m_up(2:s(1),:,i1);

        m_right = squeeze(prod(m(:,:,[1,3,4],:),3));
        m_right = max(Vd.*m_right.*Vs_rep,[],4);
        m_new(:,2:s(2),2,i1) = m_right(:,1:s(2)-1,i1);

        m_down = squeeze(prod(m(:,:,[1,2,4],:),3));
        m_down = max(Vd.*m_down.*Vs_rep,[],4);
        m_new(2:s(1),:,3,i1) = m_down(1:s(1)-1,:,i1);

        m_left = squeeze(prod(m(:,:,[1,2,3],:),3));
        m_left = max(Vd.*m_left.*Vs_rep,[],4);
        m_new(:,1:s(2)-1,4,i1) = m_left(:,2:s(2),i1);

    end
    
    X_mrf = X(MRF==0);
    Y_mrf = Y(MRF==0);
    for i1 = 1:length(X_mrf)
        m_new(X_mrf(i1),Y_mrf(i1),:,:)= ones(1,1,4,level);
    end
    
    
    fprintf('The %d iteration\n',counter);
    counter = counter+1;
    
end

% compute belief vector, b should be a 3D matrix
b = Vd.*squeeze(sum(m_new,3));  
% compute alpha level, alphaK_level should be a 2D matrix
[~,alphaK_level] = max(b,[],3);
% compute new alpha matrix
alphaK_new = alpha(alphaK_level);
end

