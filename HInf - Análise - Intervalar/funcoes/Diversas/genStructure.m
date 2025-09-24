function [res] = genStructure(sys)
 
 if(isstruct(sys{1})) % Is struct
    names = fieldnames(sys{1});
    
    for i=1:length(names)
        res.(names{i}) = zeros(size(sys{1}.(names{i})));
    end
    
 else % Not struct
     res = zeros(size(sys{1}));
 end
 
end