function BP(MRF, Vd, Vs)
iteration = 20;         % iteration times
counter = 0;            % counter for iteration
s = size(MRF);          % obtain the size of image
m = ones([s,25]);     % initialize message

while(counter<=iteration)
    for i1 = 2:s(1)-1
       for i2 = 2:s(2)-1
           
          if MRF(i1,i2)==1 && MRF(i1-1,i2)==1
              m = mesCompute(m, Vd, Vs, node, 0);
          end
          
          if MRF(i1,i2)==1 && MRF(i1,i2+1)==1
              m = mesCompute(m, Vd, Vs, node, 1);
          end
          
          if MRF(i1,i2)==1 && MRF(i1+1,i2)==1
              m = mesCompute(m, Vd, Vs, node, 2);
          end
          
          if MRF(i1,i2)==1 && MRF(i1,i2-1)==1
              m = mesCompute(m, Vd, Vs, node, 3);
          end
          
       end
    end
end