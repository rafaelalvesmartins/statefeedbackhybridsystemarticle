function [res] = sumStructure(sys1, sys2)
 
 if(isstruct(sys1)) % Is struct
    names = fieldnames(sys1);
    
    for i=1:length(names)
        res.(names{i}) = sys1.(names{i})+sys2.(names{i});
    end
    
 else % Not struct
     res = sys1+sys2;
 end
 
end