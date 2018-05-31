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
function m_new = mesCompute(m, Vd, Vs, node, dir)
m_new = m;

switch dir
    case 1          % up
        m_pre = squeeze(m(node(1),node(2)+1,4,:));          % message from right node
        m_pre = m_pre.*squeeze(m(node(1)+1,node(2),1,:));   % message from down node
        m_pre = m_pre.*squeeze(m(node(1),node(2)-1,2,:));   % message from left node
    case 2          % right
        m_pre = squeeze(m(node(1)-1,node(2),3,:));          % message from up node
        m_pre = m_pre.*squeeze(m(node(1)+1,node(2),1,:));   % message from down node
        m_pre = m_pre.*squeeze(m(node(1),node(2)-1,2,:));   % message from left node
    case 3          % down
        m_pre = squeeze(m(node(1)-1,node(2),3,:));          % message from up node
        m_pre = m_pre.*squeeze(m(node(1),node(2)+1,4,:));   % message from right node
        m_pre = m_pre.*squeeze(m(node(1),node(2)-1,2,:));   % message from left node
    case 4          % left
        m_pre = squeeze(m(node(1)-1,node(2),3,:));          % message from up node           
        m_pre = m_pre.*squeeze(m(node(1),node(2)+1,4,:));   % message from right node
        m_pre = m_pre.*squeeze(m(node(1)+1,node(2),1,:));   % message from down node
end
for i1 = 1:25   % for everys
    m_new(node(1),node(2),dir,i1) = max(squeeze(Vd(node(1),node(2))).*Vs(:,i1).*m_pre);
end
end
