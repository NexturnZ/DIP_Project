function alphaK_new = BP(MRF, Vd, Vs,alphaK, level)
iteration = 30;         % iteration times
s = size(MRF);          % obtain the size of image
m = zeros([s,4,level]);     % initialize message
m_new = ones([s,4,level]);
counter = 1;

x = 1:s(1); y = 1:s(2);
[X,Y] = meshgrid(x,y);
X = X.';    Y = Y.';
X_mrf = X(MRF==0);          % pixels that are not in MRF
Y_mrf = Y(MRF==0);
X_n = X(MRF==1);            % pixels that are in MRF
Y_n = Y(MRF==1);

while sum(m_new(:) == m(:))~=numel(m) & counter < iteration
    m = m_new;
    
    for i1 = 1:level
        Vs_rep = ones(1,1,level);
        Vs_rep(1,1,:)= Vs(i1,:);
        Vs_rep = repmat(Vs_rep,s(1),s(2),1);

        m_up = squeeze(prod(m(:,:,[2,3,4],:),3));    % MxNx25 matrix
        m_up = max(Vd.*m_up.*Vs_rep,[],3);
        m_new(1:s(1)-1,:,1,i1) = m_up(2:s(1),:);

        m_right = squeeze(prod(m(:,:,[1,3,4],:),3));
        m_right = max(Vd.*m_right.*Vs_rep,[],3);
        m_new(:,2:s(2),2,i1) = m_right(:,1:s(2)-1);

        m_down = squeeze(prod(m(:,:,[1,2,4],:),3));
        m_down = max(Vd.*m_down.*Vs_rep,[],3);
        m_new(2:s(1),:,3,i1) = m_down(1:s(1)-1,:);

        m_left = squeeze(prod(m(:,:,[1,2,3],:),3));
        m_left = max(Vd.*m_left.*Vs_rep,[],3);
        m_new(:,1:s(2)-1,4,i1) = m_left(:,2:s(2));
    end
    
    if i1==1
    cm = ones(length(x1),4,25);
        for i2 = 1:length(x1) 
            cm(i2,:,:) = m_new(x1(i2),y1(i2),:,:); 
        end
    c = mean(cm(:));
    end
    
    
    for i1 = 1:length(X_mrf)
        m_new(X_mrf(i1),Y_mrf(i1),:,:)= ones(1,1,4,level);
    end
    
    m_new = m_new/c;
    
    fprintf('The %d iteration\n',counter);
    counter = counter+1;
    
end

% compute belief vector, b should be a 3D matrix
b = Vd.*squeeze(sum(m_new,3));  
% compute alpha level, alphaK_level should be a 2D matrix
[~,alphaK_level] = max(b,[],3);
% compute new alpha matrix
alphaK_new = alphaK(alphaK_level);
end

