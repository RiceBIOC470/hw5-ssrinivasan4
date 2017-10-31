%HW5

% Note. You can use the code readIlastikFile.m provided in the repository to read the output from
% ilastik into MATLAB.

%% Problem 1. Starting with Ilastik

% Part 1. Use Ilastik to perform a segmentation of the image stemcells.tif
% in this folder. Be conservative about what you call background - i.e.
% don't mark something as background unless you are sure it is background.
% Output your mask into your repository. What is the main problem with your segmentation?  

%done in Ilastik
%the main problem with being conservative about background is that there is
%poor segmentation across cells. You can see that while some nuclei are
%visible, there is very poor seperation of cells and several cells are
%clumped together in a bigger mass.

% Part 2. Read you segmentation mask from Part 1 into MATLAB and use
% whatever methods you can to try to improve it. 

q1pt2=h5read('Segmentation (Label 2).h5', '/exported_data');
q1pt2im=squeeze(q1pt2);
imshow(q1pt2im, []);



q1pt2im_bw = bwconncomp(q1pt2im);
stats_1 = regionprops(q1pt2im_bw,'Area');
area_1 = [stats_1.Area];
fused_im1 = area_1 > mean(area_1) + std(area_1);
sublist_1 = q1pt2im_bw.PixelIdxList(fused_im1);
sublist_1=cat(1, sublist_1{:});
fusedmask1=false(size(q1pt2im));
fusedmask1(sublist_1)=1;
imshow(fusedmask1, 'InitialMagnification', 'fit');

img1_sqrt = round(1.2*sqrt(mean(area_1))/pi);
img1_erode = imerode(fusedmask1,strel('disk',img1_sqrt));
img1_outside = ~imdilate(fusedmask1,strel('disk',1));
img1_basin = imcomplement(bwdist(img1_outside));
img1_basin = imimposemin(img1_basin,img1_erode|img1_outside);
img1_watershed = watershed(img1_basin);
fusedmask1=uint8(fusedmask1);
newNuclearmask1= img1_watershed > 10 | (q1pt2im - fusedmask1);
imshow(newNuclearmask1, []);

% Part 3. Redo part 1 but now be more aggresive in defining the background.
% Try your best to use ilastik to separate cells that are touching. Output
% the resulting mask into the repository. What is the problem now?

%done in Ilastik. 
%the problem now is that there are several small holes in the smaller cells
%that are not as bright, so the nuclei of these cells are not visible. Some
%of cells are still touching.

% Part 4. Read your mask from Part 3 into MATLAB and try to improve
% it as best you can.
q1pt4=h5read('Segmentation (Label 2)__2.h5', '/exported_data');
q1pt4im=squeeze(q1pt4);
imshow(q1pt4im, []);

q1p4_er = imopen(q1pt4im,strel('disk',5));
imshow(q1p4_er, []);

%image quality from Ilastik after segmentation shows that most cells are
%seperated and thus using imopen to erode and dilate identifies objects.
%However, several cells are stil unable to be seperated from the rest.
%% Problem 2. Segmentation problems.

% The folder segmentationData has 4 very different images. Use
% whatever tools you like to try to segement the objects the best you can. Put your code and
% output masks in the repository. If you use Ilastik as an intermediate
% step put the output from ilastik in your repository as well as an .h5
% file. Put code here that will allow for viewing of each image together
% with your final segmentation. 
addpath('segmentationData/');

reader1=bfGetReader('segmentationData/bacteria.tif');
iplane=reader1.getIndex(1-1,1-1,1-1)+1;
bact_im=bfGetPlane(reader1, iplane); 
imshow(bact_im, []);

bact_thres=bact_im > 100;
imshow(bact_thres, []);
bact_er=imerode(bact_im, strel('disk',5));
imshow(bact_er, []);

imshowpair(bact_im, bact_er);

reader2=bfGetReader('segmentationData/cellPhaseContrast.png');
iplane=reader2.getIndex(1-1,1-1,1-1)+1;
phase_im=bfGetPlane(reader2, iplane); 
imshow(phase_im, []);

phase_thres=phase_im > 150;
imshow(phase_thres, []);
phase_fill=imfill(phase_thres, 'holes');
imshowpair(phase_im, phase_fill);

reader3=bfGetReader('segmentationData/worms.tif');
iplane=reader3.getIndex(1-1,1-1,1-1)+1;
worms_im=bfGetPlane(reader3, iplane); 
imshow(worms_im, []);

worms_ilastik=h5read('Seg_worms.h5', '/exported_data');
worms_ilastik=squeeze(worms_ilastik);
imshow(worms_ilastik, []);

worms_ws=imopen(worms_ilastik, strel('disk', 1));
imshow(worms_ws, []);

subplot(1,2,1);
imshow(worms_im, []);
subplot(1,2,2);
imshow(worms_ws, []);

reader4=bfGetReader('segmentationData/yeast.tif');
iplane=reader4.getIndex(1-1,1-1,1-1)+1;
yeast_im=bfGetPlane(reader4, iplane); 
imshow(yeast_im, []);
yeast_sm=imfilter(yeast_im, fspecial('gaussian', 5, 3));
imshow(yeast_sm);
yeast_bg=imopen(yeast_im, strel('disk', 200));
yeast_bgsub=imsubtract(yeast_sm, yeast_bg);

imshowpair(yeast_im, yeast_bgsub);
