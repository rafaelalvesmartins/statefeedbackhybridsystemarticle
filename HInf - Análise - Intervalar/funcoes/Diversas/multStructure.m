function [res] = multStructure(sys, num)
 
 if(isstruct(sys)) % Is struct
    names = fieldnames(sys);
    
    for i=1:length(names)
        res.(names{i}) = sys.(names{i})*num;
    end
    
 else % Not struct
     res = sys*num;
 end
 
end