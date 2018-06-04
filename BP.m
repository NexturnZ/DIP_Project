function alphaK_new = BP(MRF, Vd, Vs,alphaK, level)
iteration = 40;         % iteration times
s = size(MRF);          % obtain the size of image

% initialize message 
% 3rd dimension: 1: message to up node; 
%                2: message to right node;
%                3: message to down node;
%                4: message to left node.
m = zeros([s,4,level]);     
m_new = ones([s,4,level]);
counter = 1;

x = 1:s(1); y = 1:s(2);
[X,Y] = meshgrid(x,y);
X = X.';    Y = Y.';
X_mrf = X(MRF==1);          % pixels that are in MRF
Y_mrf = Y(MRF==1);
X_n = X(MRF==0);            % pixels that are not in MRF
Y_n = Y(MRF==0);

b = zeros([s,level])

while sum(m_new(:) == m(:))~=numel(m) & counter < iteration
    m = m_new;
    
    for i1 = 1:level
        Vs_rep = ones(1,1,level);
        Vs_corner(1,1,:)= Vs(i1,:);
        Vs_edge_hor = repmat(Vs_rep,1,s(2)-2,1);
        Vs_edge_ver = repmat(Vs_rep,s(1)-2,1,1);
        Vs_cen = repmat(Vs_rep,s(1),s(2),1);

        %% message for nodes at four corners
        % left-up corner
        m_new(1,1,2,i1) = max(Vd(1,1,:).*Vs_corner.* squeeze(m(2,1,1,:)),[],3);
        m_new(1,1,3,i1) = max(Vd(1,1,:).*Vs_corner.* squeeze(m(1,2,1,:)),[],3);
        
        % right-up corner
        m_new(1,s(2),4,i1) = max(Vd(1,s(2),:).*Vs_corner.* squeeze(m(2,s(2),1,:)),[],3);
        m_new(1,s(2),3,i1) = max(Vd(1,s(2),:).*Vs_corner.* squeeze(m(1,s(2)-1,1,:)),[],3);
        
        % left-down corner
        m_new(s(1),1,1,i1) = max(Vd(s(1),1,:).*Vs_corner.* squeeze(m(s(1),2,1,:)),[],3);
        m_new(s(1),1,2,i1) = max(Vd(s(1),1,:).*Vs_corner.* squeeze(m(s(1)-1,1,1,:)),[],3);
        
        % right-down corner
        m_new(s(1),s(2),1,i1) = max(Vd(s(1),s(2),:).*Vs_corner.* squeeze(m(s(1),s(2)-1,1,:)),[],3);
        m_new(s(1),s(2),4,i1) = max(Vd(s(1),s(2),:).*Vs_corner.* squeeze(m(s(1)-1,s(2),1,:)),[],3);
        
        
        
        %% message for nodes at four edges of image
        % up edge
        m_right = squeeze(m(1,1:s(2)-2,2,:).*m(2,2:s(2)-1,1,:));
        m_new(1,2:s(2)-1,2,i1) = max(m_right.*Vd(1,2:s(2)-1,:).*Vs_edge_hor,[],3);
        
        m_down = squeeze(m(1,1:s(2)-2,2,:).*m(1,3:s(2),4,:));
        m_new(1,2:s(2)-1,3,i1) = max(m_down.*Vd(1,2:s(2)-1,:).*Vs_edge_hor,[],3);
        
        m_left = squeeze(m(2,2:s(2)-1,1,:).*m(1,3:s(2),4,:));
        m_new(1,2:s(2)-1,4,i1) = max(m_left.*Vd(1,2:s(2)-1,:).*Vs_edge_hor,[],3);
        
        
        % right edge
        m_up = squeeze(m(3:s(1),s(2),1,:).*m(2:s(1)-1,s(2)-1,2,:));
        m_new(2:s(1)-1,s(2),1,i1) = max(m_up.*Vd(2:s(1)-1,s(2),:).*Vs_edge_ver,[],3);
        
        m_down = squeeze(m(1:s(1)-2,s(2),3,:).*m(2:s(1)-1,s(2)-1,2,:));
        m_new(2:s(1)-1,s(2),3,i1) = max(m_down.*Vd(2:s(1)-1,s(2),:).*Vs_edge_ver,[],3);
        
        m_left = squeeze(m(1:s(1)-2,s(2),3,:).*m(3:s(1),s(2),2,:));
        m_new(2:s(1)-1,s(2),4,i1) = max(m_left.*Vd(2:s(1)-1,s(2),:).*Vs_edge_ver,[],3);
        
        
        % down edge
        m_up = squeeze(m(s(1),1:s(2)-2,2,:).*m(s(1),3:s(2),4,:));
        m_new(s(1),2:s(2)-1,1,i1) = max(m_up.*Vd(s(1),2:s(2)-1,:).*Vs_edge_hor,[],3);
        
        m_right = squeeze(m(s(1)-1,2:s(2)-1,3,:).*m(s(1),1:s(2)-2,2,:));
        m_new(s(1),2:s(2)-1,2,i1) = max(m_right.*Vd(s(1),2:s(2)-1,:).*Vs_edge_hor,[],3);
        
        m_left = squeeze(m(s(1)-1,2:s(2)-1,3,:).*m(s(1),3:s(2),4,:));
        m_new(s(1),2:s(2)-1,4,i1) = max(m_left.*Vd(s(1),2:s(2)-1,:).*Vs_edge_hor,[],3);
        
        
        % left edge
        m_up = squeeze(m(2:s(1)-1,2,4,:).*m(3:s(1),1,1,:));
        m_new(2:s(1)-1,1,1,i1) = max(m_up.*Vd(2:s(1)-1,1,:).*Vs_edge_ver,[],3);
        
        m_right = squeeze(m(1:s(1)-2,1,3,:).*m(3:s(1),1,1,:));
        m_new(2:s(1)-1,1,2,i1) = max(m_right.*Vd(2:s(1)-1,1,:).*Vs_edge_ver,[],3);
        
        m_down = squeeze(m(1:s(1)-2,1,3,:).*m(2:s(1)-1,2,1,:));
        m_new(2:s(1)-1,1,3,i1) = max(m_down.*Vd(2:s(1)-1,1,:).*Vs_edge_ver,[],3);
        
        %% message for nodes at center
        m_up = squeeze(m(2:s(1)-1,3:s(2),4,:).*m(3:s(1),2:s(2)-1,1,:).*m(2:s(1)-1,1:s(2)-1,2,:));
        m_new(2:s(1)-1,2:s(2)-1,1,:) = max(m_up.*Vd(2:s(1)-1,2:s(2)-1,:).*Vs_cen,[],3);
        
        m_right = squeeze(m(1:s(1)-2,2:s(2)-1,3,:).*m(3:s(1),2:s(2)-1,1,:).*m(2:s(1)-1,1:s(2)-1,2,:));
        m_new(2:s(1)-1,2:s(2)-1,2,:) = max(m_right.*Vd(2:s(1)-1,2:s(2)-1,:).*Vs_cen,[],3);
        
        m_down = squeeze(m(1:s(1)-2,2:s(2)-1,3,:).*m(2:s(1)-1,3:s(2),4,:).*m(2:s(1)-1,1:s(2)-1,2,:));
        m_new(2:s(1)-1,2:s(2)-1,2,:) = max(m_down.*Vd(2:s(1)-1,2:s(2)-1,:).*Vs_cen,[],3);
        
        m_left = squeeze(m(1:s(1)-2,2:s(2)-1,3,:).*m(2:s(1)-1,3:s(2),4,:).*m(3:s(1),2:s(2)-1,1,:));
        m_new(2:s(1)-1,2:s(2)-1,2,:) = max(m_left.*Vd(2:s(1)-1,2:s(2)-1,:).*Vs_cen,[],3);
        
    end
    
    cm = ones(length(X_mrf),4,25);
    for i1 = 1:length(X_mrf) 
        cm(i1,:,:) = m_new(X_mrf(i1),Y_mrf(i1),:,:); 
    end
    c = max(cm(:));         % obtain normalization factor
    m_new = m_new/c;        % normalize message
    
    for i1 = 1:length(X_n)
        m_new(X_n(i1),Y_n(i1),:,:)= ones(1,1,4,level);
    end
    
    fprintf('The %d iteration\n',counter);
    counter = counter+1;
    
end

% compute belief vector, b should be a 3D matrix
% corner nodes
b(1,1,:) = Vd(1,1,:).*squeeze(m(1,2,4,:)+m(2,1,1,:));
b(1,s(1),:) = Vd(1,s(1),:).*squeeze(m(2,s(1),1,:)+m(1,s(1)-1,2,:));
b(s(1),s(2),:) = Vd(s(1),s(2),:).*squeeze(m(s(1)-1,s(2),3,:)+m(s(1),s(2)-1,2,:));
b(s(1),1,:) = Vd(s(1),1,:).*squeeze(m(s(1)-1,1,3,:)+m(s(1),2,4,:));
% edge nodes
b(1,2:s(2)-1,:) = Vd(1,2:s(2)-1,:).*squeeze(m(1,3:s(2),4,:)+m(2,2:s(1)-1,1,:)+m(1,1:s(2)-2,2,:));
b(2:s(1)-1,s(2),:) = Vd(2:s(1)-1,s(2),:).*squeeze(m(1:s(1)-2,s(2),3,:)+m(3:s(1),s(2),1,:)+m(2:s(1)-1,s(2)-1,2,:));
b(s(1),2:s(2)-1,:) = Vd(s(1),2:s(2)-1,:).*squeeze(m(s(1)-1,2:s(2)-1,3,:)+m(s(1),3:s(2),4,:)+m(s(1),1:s(2)-2,2,:));
b(2:s(1)-1,1,:) = Vd(2:s(1)-1,1,:).*squeeze(m(1:s(1)-2,1,3,:)+m(2:s(2)-1,2,4,:)+m(3:s(1),1,1,:));
% center nodes
b(2:s(1)-1,2:s(2)-1,:) = Vd(2:s(1)-1,2:s(2)-1,:).*...
    squeeze(m(1:s(1)-2,2:s(2)-1,3,:)+m(2:s(1)-2,3:s(2),4,:)+m(3:s(1),2:s(2)-1,1,:)+m(2:s(1)-1,1:s(2)-2,2,:));


% compute alpha level, alphaK_level should be a 2D matrix
[~,alphaK_level] = max(b,[],3);
% compute new alpha matrix
alphaK_new = alphaK(alphaK_level);
end

