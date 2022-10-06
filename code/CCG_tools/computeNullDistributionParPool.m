function nullDistributionMat = computeNullDistributionParPool(eventMask, ripple_inds, nIter, NREM_mask_ind, win, edges, fs)
%%
% just for testing purposes
% load('/home/iverzh/ripple/CortRipple/nullDistributionInputs.mat')
% eventMask=event_mask;
% nIter = 200;

% tic

% !!!! Run this at the beginning of your session !!!!
% parpool(16);
% if you run into issues with parpool you can delete sessions using:
% delete(gcp('nocreate'));

% Note: this function assumes 1000 Hz data

% mfilename('fullpath')

% instead make these inputs logicals and delete the following line
eventMask = logical(eventMask);

% remove ripples too close to edges
ripple_inds(ripple_inds<win) = [];
ripple_inds(ripple_inds>(numel(eventMask)-win)) = [];

% this doesn't change anything in this case, but maybe elseful otherwise?
ripple_inds = NREM_mask_ind(ripple_inds);

% binary mask for randoms (assumes 1000 Hz)
% rippleShuff = false(numel(ripple_inds),2*win+1,nIter);
eventShuff = false(numel(ripple_inds),2*win+1,nIter);

% create output structure
nullDistributionMat = zeros(nIter,length(edges)-1);

% assumes 1000 Hz
win_t = -win:1:win;
win_size = numel(win_t);


% parallelized for loops (n workers defined by parpool(n))
parfor iter = 1:nIter  
    evTmp = false(numel(ripple_inds),win_size);

    % loop through ripples
    for ri = 1:numel(ripple_inds)
        
        % extract and randomize events
        evTmp(ri,:) = eventMask(ripple_inds(ri)-win:ripple_inds(ri)+win);
        evTmp(ri,:) = evTmp(ri,randperm(win_size));
    end
    
    % add randomized ripples to structure (outside of above for loop)
    eventShuff(:,:,iter) = evTmp;
end

% adds event counts within predefined bins (outside of parfor, but fast)
for iter = 1:nIter
    for e =  1:numel(edges)-1     
        nullDistributionMat(iter,e) = sum(sum(eventShuff(:,win_t>=edges(e) & win_t<edges(e+1), iter),1));
    end
end

% toc


return


