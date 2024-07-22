
function cc=bwconncomp_vector(x)

% Function that can replace bwconncomp if used in a logical vector.
% Does not require Image Processing Toolbox.

x = [0 x 0];
xdf = diff(x);
[m,pini] = find(xdf==1);
[~,pfin] = find(xdf==-1);

cc.NumObjects = length(m);
cc.PixelIdxList = cell(length(m),1);

for i=1:length(m)
    cc.PixelIdxList{i} = pini(i):pfin(i)-1;
end



