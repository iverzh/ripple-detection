%% Take array X and value val and return the percentile of val in X

function pctle = findPercentile(X,val,mode)

    if size(X,1) == 1
        [A,~] = sort(X(:));
        nLess = sum(A < val);
        nX = length(X(:));
        pctle = nLess/nX;
    else
        pctle = nan(1,size(X,2));
        
        for ii = 1:size(X,2)
    
            [A,~] = sort(X(:,ii));
            switch mode
                
                case 'pcn' % compute percentile
                    n = sum(A < val(ii));
                    nX = length(X(:,ii));
                    pctle(ii) = n/nX;
                case 'pval' % compute pval
                    n = sum(A >= val(ii));
                    nX = length(X(:,ii));
                    pctle(ii) = n/nX;
                case 'pval2side'
                    n1 = sum(A >= val(ii));
                    n2 = sum(A <= val(ii));
                    
                    nX = length(X(:,ii));
                    pctleAll = [n1/nX,1 - n2/nX];
                    pctle(ii) = min(pctleAll);
            end

            
        end
        
    end





return 