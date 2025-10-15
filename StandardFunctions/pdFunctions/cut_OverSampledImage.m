function cutImage=cut_OverSampledImage(OSimage,ImageSize)
if numel(ImageSize)==1
    cutImage=zeros(ImageSize,1);
else
    cutImage=zeros(ImageSize);
end
cutImage=OSimage((floor(end/2)-floor(size(cutImage,1)/2))+(1:size(cutImage,1)),(floor(end/2)-floor(size(cutImage,2)/2))+(1:size(cutImage,2)),(floor(end/2)-floor(size(cutImage,3)/2))+(1:size(cutImage,3)),(floor(end/2)-floor(size(cutImage,4)/2))+(1:size(cutImage,4)));
