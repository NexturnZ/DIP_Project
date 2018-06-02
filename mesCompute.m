%-------------------------------------------------------------------------%
% Input:
% m: message matrix                 --> MxNx25 matrix
% Vd
% Vs
% node: index of node               --> 1x2 vector
% dir: the position of another node --> scale
%       1:up, 2:right, 3:down, 4:left
%                     
% Output:
% m_new: new message matrix         --> MxNx25 matrix;
%-------------------------------------------------------------------------%
function m_new = mesCompute(m, Vd, Vs, node)
m_new = m;
s = size(Vd(:,:,1));

% up
if node(1) ~= 1
    m_up = squeeze(prod(m(node(1),node(2)+1,[2,3,4],:),3));
    for i1 = 1:25   % for everys
        m_new(node(1)-1,node(2),3,i1) = max(squeeze(Vd(node(1),node(2))).*Vs(:,i1).*m_up);
    end
end

% right
if node(2) ~= s(2)
    m_right = squeeze(prod(m(node(1)-1,node(2),[1,3,4],:),3));          % message from up node
    for i1 = 1:25   % for everys
        m_new(node(1),node(2)+1,4,i1) = max(squeeze(Vd(node(1),node(2))).*Vs(:,i1).*m_right);
    end
end

% down
if node(1) ~= s(1)
    m_down = squeeze(prod(m(node(1)-1,node(2),[1,2,4],:),3));          % message from up node
    for i1 = 1:25   % for everys
        m_new(node(1)+1,node(2),1,i1) = max(squeeze(Vd(node(1),node(2))).*Vs(:,i1).*m_down);
    end
end

% left
if node(2) ~= 1
    m_left = squeeze(prod(m(node(1)-1,node(2),[1,2,3],:),3));          % message from up node      
    for i1 = 1:25   % for everys
        m_new(node(1),node(2)-1,2,i1) = max(squeeze(Vd(node(1),node(2))).*Vs(:,i1).*m_left);
    end
end

end
