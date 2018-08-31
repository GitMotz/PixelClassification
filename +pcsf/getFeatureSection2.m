function cellImSections = getFeatureSection2(cellIms,PadSize)
cellImSections = cell(size(cellIms));
for i = 1:length(cellIms)
cellImSections{i} = cellIms{i}(1+PadSize:end-PadSize,1+PadSize:end-PadSize);
end