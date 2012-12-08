function [Y,stress]=resmds(M, d, resample, structure, showresult, embeddingfun)
% [Y,stress]=resmds(M, d, resample, structure)
% M - distance / similarity matrix.
% d - target space dimension
% structure - block matrix structure (separate objects in the distance
% matrix)
% ex: M=(tM./max(tM(:))); d=3; resample=5; structure=[ChromosomeBins(1) ; diff(ChromosomeBins)];
%% initialize
if ~exist('showresult','var')
    showresult=0;
end
if ~exist('embeddingfun','var')
    embeddingfun=@mdscale;
end

%% Resample
blockify=mat2cell(M,structure,structure);
interpblk=cell(size(blockify));

for i=1:length(blockify)
    for j=1:length(blockify)
        xsample=linspace(1,size(blockify{i,j},2),ceil(size(blockify{i,j},2)/resample));
        ysample=linspace(1,size(blockify{i,j},1),ceil(size(blockify{i,j},1)/resample));
        [XI,YI]=meshgrid(xsample,ysample);
        interpblk{i,j}=interp2(blockify{i,j},XI,YI);
    end
end
tM=cell2mat(interpblk); 
tM=(tM+tM')./2; tM=tM./max(tM(:)); tM(logical(eye(size(tM))))=1; %fix interpolation roundoff errors

%% Embedd
%[Y,stress]=mdscale(tM,d);
[Y,stress]=embeddingfun(tM,d);

%% Interpolate
structuresampled=cellfun(@(x) size(x,2),interpblk(1,:));
cY=mat2cell(Y,structuresampled,3);
interpY=cell(size(cY));
for i=1:length(structuresampled)
    xsample=linspace(1,structuresampled(i),structure(i));
    interpY{i}=interp1(cY{i},xsample,'spline');
end
iY=cell2mat(interpY);

if showresult
    colvec=[];
    for i=1:length(interpblk)
        colvec=[colvec repmat(i,1,structure(i))];
    end
    figure; 
    scatter3(iY(:,1),iY(:,2),iY(:,3),30,colvec,'filled'); axis vis3d;
    D=sqrt(1-(M./max(M(:))));
    title(num2str(stressCrit(iY,D(tril(true(size(D)),-1))')));
end

end