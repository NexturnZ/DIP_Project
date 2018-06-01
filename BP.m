function alphaK_new = BP(MRF, Vd, Vs,alpha, level)
iteration = 20;         % iteration times
s = size(MRF);          % obtain the size of image
m = ones([s,4,25]);     % initialize message


for i1 = 1:iteration
    m_new = m;
    
    for i2 = 1:level
        Vs_rep = ones(1,1,level);
        Vs_rep(1,1,:)= Vs(i2,:);
        Vs_rep = repmat(Vs_rep,s(1),s(2),1);

        m_up = squeeze(prod(m(:,:,[2,3,4],:),3));    % MxNx25 matrix
        m_up = max(Vd.*m_up.*Vs_rep,[],4);
        m_new(1:s(1)-1,:,1,i2) = m_up(2:s(1),:,i2);

        m_right = prod(m(:,:,[1,3,4],:),3);
        m_right = max(Vd.*m_right.*Vs_rep,[],4);
        m_new(:,1:s(2)-1,2,i2) = m_right(:,2:s(2),i2);

        m_down = prod(m(:,:,[1,2,4],:),3);
        m_down = max(Vd.*m_down.*Vs_rep,[],4);
        m_new(2:s(1),:,3,i2) = m_down(1:s(1)-1,:,i2);

        m_left = prod(m(:,:,[1,2,3],:),3);
        m_left = max(Vd.*m_left.*Vs_rep,[],4);
        m_new(:,1:s(2)-1,4,i2) = m_left(:,2:s(2),i2);

    end
    
%     for i1 = 2:s(1)-1
%        for i2 = 2:s(2)-1
%            
%           % compute &pass message for 4 neighbor respectively
%           if MRF(i1,i2)==1 && MRF(i1-1,i2)==1           
%               m = mesCompute(m, Vd, Vs, [i1,i2], 1);
%           end
%           
%           if MRF(i1,i2)==1 && MRF(i1,i2+1)==1
%               m = mesCompute(m, Vd, Vs, [i1,i2], 2);
%           end
%           
%           if MRF(i1,i2)==1 && MRF(i1+1,i2)==1
%               m = mesCompute(m, Vd, Vs, [i1,i2], 3);
%           end
%           
%           if MRF(i1,i2)==1 && MRF(i1,i2-1)==1
%               m = mesCompute(m, Vd, Vs, [i1,i2], 4);
%           end
%           
%        end
%     end
end

% compute belief vector, b should be a 3D matrix
b = Vd.*squeeze(sum(m,3));  
% compute alpha level, alphaK_level should be a 2D matrix
[~,alphaK_level] = max(b,[],3);
% compute new alpha matrix
alphaK_new = alpha(alphaK_level);
end

