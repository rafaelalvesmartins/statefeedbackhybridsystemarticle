function [outPut] =  discMainPart(opt,A,B2,B1,C,D2,D1,tol,precision,combPoly,folderExample,APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,h)


%% ------------------Mid Matrix----------------------------------------
    if(contains(opt,'Mid'))

        %        Int Approach
        midA = infsup(mid(A),mid(A));
        midB2 = infsup(mid(B2),mid(B2));
        midB1 = infsup(mid(B1),mid(B1));
        midC = infsup(mid(C),mid(C));
        midD2 = infsup(mid(D2),mid(D2));
        midD1 = infsup(mid(D1),mid(D1));    
        
        clear paramInput;
        paramInput.tol = tol;
        synMidInt = synHInfKDiscIntWG(midA,midB2,midB1,midC,midD2,midD1,paramInput);
      
         
        clear paramInput;
        paramInput.disc = 1;
        hinfNormInt = calcHinfnorm(mid(A),mid(B2),mid(B1),mid(C),mid(D2),mid(D1),synMidInt.K,paramInput);
        if( trunckValue(synMidInt.norm, precision)   ~= trunckValue(hinfNormInt, precision) )
            error('ERROR - different h inf J Int :/');
        end  
        
        %       Poly Approach        
        midAPoly{1} = mid(A);
        midB2Poly{1} = mid(B2);
        midB1Poly{1} = mid(B1);
        midCPoly{1} = mid(C);
        midD2Poly{1} = mid(D2);
        midD1Poly{1} = mid(D1);
               
        clear paramInput;
        paramInput.tol = tol;
        synMidPoly = synHInfKDiscPoly(midAPoly,midB2Poly,midB1Poly, midCPoly,midD2Poly,midD1Poly,paramInput);
           
        clear paramInput;
        paramInput.disc = 1;
        hInfNormPoly = calcHinfnorm(mid(A),mid(B2),mid(B1),mid(C),mid(D2),mid(D1),synMidPoly.K,paramInput);
        if( trunckValue(synMidPoly.norm, precision)   ~= trunckValue(hInfNormPoly, precision) )
            error('ERROR - different h inf J Poly :/');
        end  
        
        % %        Check poly and int approaches
        if(trunckValue(hInfNormPoly, precision)~= trunckValue(hinfNormInt, precision))
            error('ERROR - different h inf J between Poly and Int :/');
        end
        
        
%         Riccati Approach
        clear paramInput;
        paramInput.disc = 1;
        riccatiHinf =  calcKRiccatiHinf(mid(A),mid(B2),mid(B1),mid(C),mid(D2),mid(D1),synMidInt.norm,paramInput);
        if sum(eig(riccatiHinf.X)<0)
             error('ERROR - No solution for Riccati Equation :/');
        end 
        
        clear paramInput;
        paramInput.disc = 1;
        hInfNormRiccat = calcHinfnorm(mid(A),mid(B2),mid(B1),mid(C),mid(D2),mid(D1),riccatiHinf.K,paramInput);
        if( trunckValue(synMidPoly.norm, precision)   ~= trunckValue(hInfNormRiccat, precision) )
            error('ERROR - different h inf J Riicati :/');
        end  
     
        
        
 % Simulate K from the mid system for the whole system    
 %         Int
        clear paramInput;
        paramInput.disc = 1;
        paramInput.h = h;
        paramInput.imageNameFile = [folderExample '/' 'midInt'] ;
        midIntSimul = simulatesSys(combPoly,synMidInt.K,paramInput); 
        
%         Poly
        clear paramInput;
        paramInput.disc = 1;
        paramInput.h = h;
        paramInput.imageNameFile = [folderExample '/' 'midPoly'];
        midPolySimul = simulatesSys(combPoly,synMidPoly.K,paramInput); 
                
%         Riccati
        clear paramInput;
        paramInput.disc = 1;
        paramInput.h = h;
        paramInput.imageNameFile = [folderExample '/' 'midRiccati'];
        midRiccatiSimul = simulatesSys(combPoly,riccatiHinf.K,paramInput);         
              
               
        % Save parameters to outPut
        poly.syn = synMidPoly;
        poly.simul = midPolySimul;
        
        int.syn = synMidInt;
        int.simul = midIntSimul;
        
        riccati.syn = riccatiHinf;
        riccati.simul = midRiccatiSimul;
         
        midOutPut.poly = poly;
        midOutPut.int = int;
        midOutPut.riccati = riccati;
       
        outPut.midOutPut = midOutPut;
        save([folderExample '/dataOutPutBackup'],'outPut');
    end
    %% ------------------Guaranteed Cost Int------------------------------
    if(contains(opt,'Int'))    
        
        clear paramInput;
        paramInput.tol = tol;
        synInt = synHInfKDiscIntWG(A,B2,B1,C,D2,D1,paramInput);
                 
        if(synInt.feas)
             intK = synInt.K;
             clear paramInput;
             paramInput.disc = 1;
             paramInput.h = h;
             paramInput.imageNameFile = [folderExample '/' 'Int'];
             intSimul = simulatesSys(combPoly,intK,paramInput); 
        else
            intSimul = 0;
        end
        
        % Save parameters to outPut
        outPut.intOutput.syn = synInt;
        outPut.intOutput.simul = intSimul;
        save([folderExample '/dataOutPutBackup'],'outPut');
    end
    
    %% ------------------Guaranteed Cost Poly------------------------------
    if(contains(opt,'Poly'))
        
        clear paramInput;
        paramInput.tol = tol;
        synPoly = synHInfKDiscPoly(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,paramInput);
                  
        if(synPoly.feas)
             polyK = synPoly.K;
             clear paramInput;
             paramInput.disc = 1;
             paramInput.h = h;
             paramInput.imageNameFile = [folderExample '/' 'poly'];
             polySimul = simulatesSys(combPoly,polyK,paramInput);       
        else
            polySimul = 0;
        end
        
        % Save parameters to outPut
        outPut.polyOutput.syn = synPoly;
        outPut.polyOutput.simul = polySimul;
        save([folderExample '/dataOutPutBackup'],'outPut');
    end
    
end