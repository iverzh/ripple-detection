
function ripple = recenterRipples(ripple,rippleWindow, win)

yes = 0;
no = 0;
% re-compute ripple centers based on largest positive peak in ripple band
for rip = 1:numel(ripple.ind)
    
    rippleWindowStart = rippleWindow;
    rippleWindowEnd = rippleWindow;
    
    if rip > 1 
        checkNeighbor = ripple.ind(rip) - ripple.ind(rip-1);
        
        if checkNeighbor < rippleWindow
            rippleWindowStart = checkNeighbor - 2;
        end
    end
    
    if rip < numel(ripple.ind)
        checkNeighbor = ripple.ind(rip + 1) - ripple.ind(rip);
        
        if checkNeighbor < rippleWindow
            rippleWindowEnd = checkNeighbor - 2;
        end
    end

    rip_epoch_start = ripple.ind(rip) - rippleWindowStart ;
    rip_epoch_end = ripple.ind(rip) + rippleWindowEnd;
    rip_epoch = ripple.rippleband(rip_epoch_start:rip_epoch_end);

    % find positive peaks in rippleband 
    [pks, locs] = findpeaks(rip_epoch);

    % find index of max peak
    [~,max_ind] = max(pks);

    % only update ripple center if max peak exists
    if ~isempty(max_ind)

        % if multiple max peaks use first
        max_loc = locs(max_ind(1));
        new_loc = rip_epoch_start+max_loc-1;
        ripple.ind(rip) = new_loc;

        
    end

    
    
    if max(ripple.rippleband(ripple.ind(rip)-rippleWindowStart:ripple.ind(rip)+rippleWindowEnd)) ~= ripple.rippleband(ripple.ind(rip)) 

        RBzFullWindow = ripple.RBzscore(ripple.ind(rip)-win:ripple.ind(rip)+win);
        center = round(size(RBzFullWindow,2)/2);
        score = Inf;
        winStart = 0;
        while score > 0.75 && winStart < 200
            score = abs(RBzFullWindow(center-winStart));
            winStart = winStart + 1;
        end

        score = Inf;
        winEnd = 0;
        while score > 0.75 && winEnd < 200 
            score = abs(RBzFullWindow(center+winEnd));
            winEnd = winEnd + 1;
        end
      
        if rip > 1
            checkNeighbor = ripple.ind(rip) - ripple.ind(rip-1);

            if checkNeighbor < winStart
                winStart = checkNeighbor - 5;
            end
        end

        if rip < numel(ripple.ind)
            checkNeighbor = ripple.ind(rip + 1) - ripple.ind(rip);

            if checkNeighbor < winEnd
                winEnd = checkNeighbor - 5;
            end
        end

        if winEnd > 2 && winStart > 2
           

            % recenter again based on new ripple duration
            rip_epoch_newCenter = ripple.rippleband(ripple.ind(rip)-winStart:ripple.ind(rip)+winEnd); 
            rip_epoch_newCenter_start = ripple.ind(rip)-winStart;
            [pks, locs] = findpeaks(rip_epoch_newCenter);

            % find index of max peak
            [~,max_ind] = max(pks);

            % only update ripple center if max peak exists
            if ~isempty(max_ind)

                % if multiple max peaks use first
                max_loc = locs(max_ind(1));
                new_loc = rip_epoch_newCenter_start+max_loc-1;
                ripple.ind(rip) = new_loc;
            end

            if max(ripple.rippleband(ripple.ind(rip)-winStart:ripple.ind(rip)+winEnd)) == ripple.rippleband(ripple.ind(rip)) 
                yes = yes + 1;
            else 
                no = no + 1; 
            end

        end
    end
end

% fprintf([num2str(yes), ' '])
% fprintf([num2str(no), ' '])








