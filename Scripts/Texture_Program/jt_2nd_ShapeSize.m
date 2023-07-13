%% NOTE: Image and ROI should have isotropic voxel size
% I: input image
% T: tumor ROI

function fea=jt_2nd_ShapeSize(I,T)
  %Impletemented Feature List
  IFL={'Area';'Volume';'Compactness1';'Compactness2';'Max3dDiam';'SphericalDisprop';'Spherity';'Surf2VolRatio';};
  
  % the following will be calculated anyway
  Area=surfaceArea3D(T);
  Volume=sum(T(:));
  
  for i=1:length(IFL);
    cF=IFL{i};  %current feature name
    if strcmp(cF,'Area')
      fea.Area=Area;
    elseif strcmp(cF,'Volume')
      fea.Volume=Volume;
    elseif strcmp(cF,'Compactness1')
      fea.Compactness1=Volume/(sqrt(pi)*power(Area,2/3));
    elseif strcmp(cF,'Compactness2')
      fea.Compactness2=36*pi*Volume*Volume/power(Area,3);
    elseif strcmp(cF,'Max3dDiam')
      fea.Max3dDiam=Max3dDiam(T);
    elseif strcmp(cF,'SphericalDisprop')
      r=power(3*Volume/(4*pi),1/3);
      fea.SphericalDisprop=Area/(4*pi*r*r);
    elseif strcmp(cF,'Spherity')
      fea.Spherity=power(pi,1/3)*power(6*Volume,2/3)/Area;
    elseif strcmp(cF,'Surf2VolRatio')
      fea.Surf2VolRatio=Area/Volume;
    else
      disp(['The requested feature {' cF '} has NOT been implemented.' 'The available features include: ']);
      disp(IFL');
    end
  end
end

%%
%% maximum distance: the largest pairwise eucledean distance
% NOT recommended as MaxDistance is sensitive to outliers
function D=Max3dDiam(T)
  ind=find(T);  %voxels in ROI
  [i j k]=ind2sub(size(T),ind);
  coord=[i j k];
  cc=repmat(mean(coord),length(i),1); %centroid
  
  d=sqrt(sum(abs(coord-cc).^2,2));  %distance of each voxel to the centroid
  %d=pdist([I J K]);
  
  D=max(d)*2;
end

%%
function Area=surfaceArea3D(T,option)

%% Implementation 1
% % Construct kernel where we can count the number of 6-connected neighbor voxels.
% conn6Kernel = zeros([3,3,3]);
% conn6Kernel(2,2,1) = 1;
% conn6Kernel(1,2,2) = 1;
% conn6Kernel(2,1,2) = 1;
% conn6Kernel(2,3,2) = 1;
% conn6Kernel(3,2,2) = 1;
% conn6Kernel(2,2,3) = 1;
% 
% % For each voxel, determine how many 6-connected neighbors it has.
% sumOfFaces=convn(T, conn6Kernel, 'same');
% % Find number of exposed faces for each voxel.
% surfaceArea=6*T-sumOfFaces;
% % Mask out zero voxels that have negative exposed faces.
% surfaceArea(surfaceArea<0)=0;
% 
% % Now we simply label the volume and sum up the values of each region.
% binaryVolume=surfaceArea>0;
% cc=bwconncomp(binaryVolume,6);
% 
% % Now it's labeled, so now we measure the PixelValues in each blob.
% measurements=regionprops(cc,surfaceArea,'PixelValues');
% 
% % Go through the regions, listing each regions exposed surface area.
% numberOfRegions=length(measurements);
% 
% if numberOfRegions>1
%   disp(['Wanrning: there are ' num2str(numberOfRegions) ' tumors in the brain?']);
% end
% 
% Area=0;
% for k = 1:numberOfRegions
%   % For this region, find out the sum of exposed faces.
%   thesePixelValues=[measurements(k).PixelValues];
%   % Sum up the number of exposed faces for each voxel in the blob.
%   thisRegionsArea = sum(thesePixelValues);
% %   fprintf('The exposed surface area for region %d is %d\n',k, thisRegionsArea);
%   Area=Area+thisRegionsArea;
% end

%% Implementation 2
[f,v]=isosurface(T,0.5);

a=v(f(:,2),:)-v(f(:,1),:);
b=v(f(:,3),:)-v(f(:,1),:);
c=cross(a,b,2);
Area=1/2*sum(sqrt(sum(c.^2, 2)));

% visualize the tumore VOI
% patch('Faces',f,'Vertices',v,'facecolor','r');
% view(30,-15);
% axis vis3d;
% colormap copper

end