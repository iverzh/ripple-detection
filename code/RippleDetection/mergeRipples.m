function mergedLocs = mergeRipples(ripple,thresh,fs)

thresh = thresh * fs; %convert from seconds to samples


if any(diff(ripple.ind) < thresh)
    while any(diff(ripple.ind) < thresh)
        interRipAll = diff(ripple.ind);
        mergeMask = interRipAll < thresh;
        mergedLocs = [];
        rip =  1;
        while rip < numel(ripple.ind)

            if mergeMask(rip)
                ii1 = ripple.ind(rip);
                ii2 = ripple.ind(rip+1);

                filtDat = ripple.rippleband(ii1:ii2);

                newInd = find(filtDat == max(filtDat)) - 1;

                mergedLocs = [mergedLocs, ii1 + newInd];

                rip = rip + 2;


            else
                mergedLocs = [mergedLocs, ripple.ind(rip)];
                rip = rip + 1;

            end

        end

        ripple.ind = mergedLocs;
    end
else
    mergedLocs = ripple.ind;
end
