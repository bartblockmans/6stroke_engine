function [Xc, Yc] = AddManifolds(Xc, Yc, params, IV, EV, SP)
%ADD_MANIFOLDS  Add intake, exhaust, and scavenge manifolds to cylinder cavity.
%
%   [Xc, Yc] = AddManifolds(Xc, Yc, params, IV, EV, SP)
%
%   Inputs:
%       Xc, Yc : cylinder cavity coordinates [mm]
%       params : structure with geometry parameters (from parameters.m)
%       IV : intake valve status (1=open, 0=closed)
%       EV : exhaust valve status (1=open, 0=closed)
%       SP : scavenge port status (1=open, 0=closed)
% -------------------------------------------------------------------------
% Bart Blockmans, 2024 - bart@blockmans.net
% -------------------------------------------------------------------------


% Switch between valve situation
if IV == 1

    % Construct intact manifold
    X_intake = [params.xmi1, params.xir, fliplr(params.xmi2)];
    Y_intake = [params.ymi1, params.yir, fliplr(params.ymi2)];

    % Merge
    [Xc, Yc] = MergeManifold(Xc, Yc, X_intake, Y_intake);

end
if EV == 1 

    % Construct exhaust manifold
    X_exhaust = fliplr([params.xme2, params.xer, fliplr(params.xme1)]);
    Y_exhaust = fliplr([params.yme2, params.yer, fliplr(params.yme1)]);

    % Merge
    [Xc, Yc] = MergeManifold(Xc, Yc, X_exhaust, Y_exhaust);

end
if SP == 1

    % Construct left scavenge port
    X_scav_l = [params.xmsl1, params.xsrl, fliplr(params.xmsl2)];
    Y_scav_l = [params.ymsl1, params.ysrl, fliplr(params.ymsl2)];

    % Merge
    [Xc, Yc] = MergeManifold(Xc, Yc, X_scav_l, Y_scav_l);

    % Construct right scavenge port
    X_scav_r = fliplr([params.xmsr1, params.xsrr, fliplr(params.xmsr2)]);
    Y_scav_r = fliplr([params.ymsr1, params.ysrr, fliplr(params.ymsr2)]);

    % Merge
    [Xc, Yc] = MergeManifold(Xc, Yc, X_scav_r, Y_scav_r);

end

% Make sure Xc & Yc represents a closed polygon
if Xc(1) ~= Xc(end) || Yc(1) ~= Yc(end)
    Xc = [Xc, Xc(1)];
    Yc = [Yc, Yc(1)];
end

% Nested function for merging
function [Xc, Yc] = MergeManifold(Xc, Yc, X_merge, Y_merge)

% Get the two endpoints of the manifold (connection points)
P1_manifold = [X_merge(1), Y_merge(1)];
P2_manifold = [X_merge(end), Y_merge(end)];

% Find closest points in cylinder cavity to manifold endpoints
dist1 = sqrt((Xc - P1_manifold(1)).^2 + (Yc - P1_manifold(2)).^2);
dist2 = sqrt((Xc - P2_manifold(1)).^2 + (Yc - P2_manifold(2)).^2);

[~, idx1] = min(dist1);
[~, idx2] = min(dist2);

% Handle edge case where both endpoints map to the same point
if idx1 == idx2
    % If both endpoints are closest to the same point, use adjacent points
    idx2 = idx1 + 1;
    if idx2 > length(Xc)
        idx2 = 1; % Wrap around for closed curve
    end
end

% Ensure idx1 < idx2 (for typical case where points are adjacent)
if idx1 > idx2
    % Swap indices
    temp = idx1;
    idx1 = idx2;
    idx2 = temp;
end

% Remove the section between idx1 and idx2 (inclusive)
% Insert manifold geometry in place of removed section
if idx1 == 1
    % Handle case where idx1 is at the start
    Xc_new = [X_merge, Xc(idx2+1:end)];
    Yc_new = [Y_merge, Yc(idx2+1:end)];
elseif idx2 == length(Xc)
    % Handle case where idx2 is at the end
    Xc_new = [Xc(1:idx1-1), X_merge];
    Yc_new = [Yc(1:idx1-1), Y_merge];
else
    % Normal case: remove section and insert manifold
    Xc_new = [Xc(1:idx1-1), X_merge, Xc(idx2+1:end)];
    Yc_new = [Yc(1:idx1-1), Y_merge, Yc(idx2+1:end)];
end

Xc = Xc_new;
Yc = Yc_new;

end

end

