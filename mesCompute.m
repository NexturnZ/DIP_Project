%-------------------------------------------------------------------------%
% Input:
% m: message matrix                 --> MxNx25 matrix
% Vd
% Vs
% node: index of node               --> 1x2 vector
% dir: the position of another node --> scale
%       0:up, 1:right, 2:down, 3:left
%                     
% Output:
% m_new: new message matrix         --> MxNx25 matrix;
%-------------------------------------------------------------------------%
function m_new = mesCompute(m, Vd, Vs, node, dir)
m_new = m;

switch dir
    case 0
        m_pre = squeeze(m(node(1),node(2)+1));
        m_pre = m_pre.*squeeze(m(node(1)+1,node(2)));
        m_pre = m_pre.*squeeze(m(node(1),node(2)-1));
    case 1
        m_pre = squeeze(m(node(1)-1,node(2)));
        m_pre = m_pre.*squeeze(m(node(1)+1,node(2)));
        m_pre = m_pre.*squeeze(m(node(1),node(2)-1));
    case 2
        m_pre = squeeze(m(node(1)-1,node(2)));
        m_pre = m_pre.*squeeze(m(node(1),node(2)+1));
        m_pre = m_pre.*squeeze(m(node(1),node(2)-1));
    case 3
        m_pre = squeeze(m(node(1)-1,node(2)));
        m_pre = m_pre.*squeeze(m(node(1),node(2)+1));
        m_pre = m_pre.*squeeze(m(node(1)+1,node(2)));
end
for i1 = 1:25   % for everys
    m_new(node(1),node(2),i1) = max(squeeze(Vd(node(1),node(2))).*Vs(:,i1).*m_pre);
end
end
