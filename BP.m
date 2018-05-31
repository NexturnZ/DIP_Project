function alphaK_new = BP(MRF, Vd, Vs,alpha)
iteration = 20;         % iteration times
counter = 0;            % counter for iteration
s = size(MRF);          % obtain the size of image
m = ones([s,4,25]);     % initialize message

while(counter<=iteration)
    for i1 = 2:s(1)-1
       for i2 = 2:s(2)-1
           
          % compute &pass message for 4 neighbor respectively
          if MRF(i1,i2)==1 && MRF(i1-1,i2)==1           
              m = mesCompute(m, Vd, Vs, node, 1);
          end
          
          if MRF(i1,i2)==1 && MRF(i1,i2+1)==1
              m = mesCompute(m, Vd, Vs, node, 2);
          end
          
          if MRF(i1,i2)==1 && MRF(i1+1,i2)==1
              m = mesCompute(m, Vd, Vs, node, 3);
          end
          
          if MRF(i1,i2)==1 && MRF(i1,i2-1)==1
              m = mesCompute(m, Vd, Vs, node, 4);
          end
          
       end
    end
end

% compute belief vector, b should be a 3D matrix
b = Vd.*squeeze(sum(m,3));  
% compute alpha level, alphaK_level should be a 2D matrix
[~,alphaK_level] = max(b,[],3);
% compute new alpha matrix
alphaK_new = alpha(alphaK_level);
end

