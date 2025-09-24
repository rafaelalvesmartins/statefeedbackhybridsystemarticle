function [vecNorm] = normalizeSizeVec(vec,normLen)
    
    lenCurVec = length(vec);
    vecNorm = vec;
    
    if(lenCurVec ~= normLen)       
       
        if(normLen>lenCurVec)
        diff = normLen - lenCurVec;

        for j=1:diff
            vecNorm(end+1,:) = vecNorm(end,:);
        end
        
       end
       elseif(normLen<lenCurVec)
           
           diff = lenCurVec - normLen;

           for j=1:diff
               vecNorm(end-(j-1),:)=[];   
           end    
    end
end


               
          