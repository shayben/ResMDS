function [S,grad] = stressCrit(Y,disparities)
%STRESS Stress criterion for nonmetric multidimensional scaling.

% The Euclidean norm of the differences between the distances and the
% disparities, normalized by the Euclidean norm of the distances.
%
% Zero weights are ok as long as the corresponding dissimilarities are
% finite.  However, any zero distances will throw an error:  this criterion
% is not differentiable when two points are coincident - it involves, in
% effect, abs(Y(i,k)-Y(j,k)).  This does sometimes happen, even though a
% configuration that is a local minimizer cannot have coincident points.

distances = pdist(Y);
diffs = distances - disparities;
sumDiffSq = sum(diffs.^2);
sumDistSq = sum(distances.^2);
S = sqrt(sumDiffSq ./ sumDistSq);
if nargout > 1
    [n,p] = size(Y);
    grad = zeros(n,p);
    if sumDiffSq > 0
        if all(distances > 0)
            dS = squareform((diffs./sumDiffSq - distances./sumDistSq) ./ distances);
            repcols = zeros(1,n);
            for i = 1:p
                repcols = repcols + 1;
                dY = Y(:,repcols) - Y(:,repcols)';
                grad(:,i) = sum(dS.*dY,2) .* S;
            end
        else
            error('stats:mdscale:ColocatedPoints', ...
                  ['Points in the configuration have co-located.  Try a' ...
                   ' different\nstarting point, or use a different criterion.']);
        end
    end
end
