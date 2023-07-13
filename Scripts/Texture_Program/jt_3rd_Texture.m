% NOTE 1: Image and ROI should have isotropic voxel size
% NOTE 2: InverseVariance is unstable, suggest removing it
% NOTE 3: GLCM could be time consuming, parfor was used for parallel computing

%%
% Input:
%     I: input image
%     T: tumor ROI
%     nLGLCM: number of offsets of GLCM
%     nLGLRL: number of levels of GLRL
% Output:
%     fea: structure of extracted features
%  use fieldnames(fea) in command window to list all feature names

% Example: fea=jt_3rd_Texture(I,T,32,16);

function fea=jt_3rd_Texture(I,T,nLGLCM,nLGLRL)
% nLGLCM=8;nLGLRL=8; % for debugging
warning off %#ok<WNOFF>
  %% GLCM: 
  IT=I;
  IT(~T)=NaN;
  feaGLCM=GLCM3d(IT,nLGLCM);
  
  %% GLRL
  feaGLRL=GLRL3d(I,T,nLGLRL);

  %% combine GLCM and GLRL
  fea=catstruct(feaGLCM,feaGLRL);

end

%% GLCM
% Input: 
%     IT - ROI with intensity, background was set to 0
%     nL - number of levels
% Output:
%     G - mean GLCM (numLevel x numLevel)
function G=GLCM3d(IT,nL)
  %calculaition matrix in 2d
  % implemented 10 directions
 warning off %#ok<WNOFF>
 sz=size(IT);
  
  d=1;  %delta: distance
  
  dao=[0 d; -d d; -d 0; -d -d]; %offset: distance and orientation
  numMatrix=size(dao,1);  %number of DAOs
  glcmXY=zeros(nL,nL,numMatrix,sz(3));
  parfor i=1:sz(3)  %xy plane
    cIT=squeeze(IT(:,:,i)); %xy plane
     if length(find(~isnan(cIT)))>4
        [glcmXY(:,:,:,i),~] = graycomatrix(cIT,'NumLevels',nL,'Offset',dao,'G',[],'Symmetric',true); 
     end
  end
  dao=[0 d; -d d; -d -d]; %offset: distance and orientation
  numMatrix=size(dao,1);  %number of DAOs
  glcmYZ=zeros(nL,nL,numMatrix,sz(1));
  parfor i=1:sz(1)  %yz plane
    cIT=squeeze(IT(i,:,:)); 
     if length(find(~isnan(cIT)))>4
        [glcmYZ(:,:,:,i),~] = graycomatrix(cIT,'NumLevels',nL,'Offset',dao,'G',[],'Symmetric',true);  
     end
  end
  dao=[-d d; -d 0; -d -d]; %offset: distance and orientation
  numMatrix=size(dao,1);  %number of DAOs
  glcmXZ=zeros(nL,nL,numMatrix,sz(2));
  parfor i=1:sz(2)  %xz plane
    cIT=squeeze(IT(:,i,:)); %xy plane
     if length(find(~isnan(cIT)))>4
        [glcmXZ(:,:,:,i),~] = graycomatrix(cIT,'NumLevels',nL,'Offset',dao,'G',[],'Symmetric',true);  
     end
  end
%   glcmXZ(:,:,:,218)=0;
  meanG=cat(3,mean(glcmXY,4),mean(glcmYZ,4),mean(glcmXZ,4));%mean over all slices
  %% calculate features
  GLCMall = GLCM_Features4(meanG);
  
%   G.Autocorrelation=GLCMall.autoc;
%   G.ClusterProminence=GLCMall.cprom;
%   G.ClusterShade=GLCMall.cshad;
%   G.ClusterTendency=GLCMall.ctend;
%   G.Contrast=GLCMall.contr;
%   G.Correlation=GLCMall.corrp;
%   G.DifferenceEntropy=GLCMall.denth;
%   G.Dissimilarity=GLCMall.dissi;
%   G.Energy=GLCMall.energ;
%   G.Entropy=GLCMall.entro;
%   G.Homogeneity1=GLCMall.homom;
%   G.Homogeneity2=GLCMall.homop;
%   G.IMC1=GLCMall.inf1h;
%   G.IMC2=GLCMall.inf2h;
%   G.IDMN=GLCMall.idmnc;
%   G.IDN=GLCMall.indnc;
%   G.InverseVariance=GLCMall.invva;
%   G.MaximumProbability=GLCMall.maxpr;
%   G.SumAverage=GLCMall.savgh;
%   G.SumEntropy=GLCMall.senth;
%   G.SumVariance=GLCMall.svarh;
%   G.Variance=GLCMall.sosvh;

  G.Autocorrelation=mean(GLCMall.autoc); 
  G.ClusterProminence=mean(GLCMall.cprom);
  G.ClusterShade=mean(GLCMall.cshad);
  G.ClusterTendency=mean(GLCMall.ctend);
  G.Contrast=mean(GLCMall.contr);
  G.Correlation=mean(GLCMall.corrp);
  G.DifferenceEntropy=mean(GLCMall.denth);
  G.Dissimilarity=mean(GLCMall.dissi);
  G.Energy=mean(GLCMall.energ);
  G.Entropy=mean(GLCMall.entro);
  G.Homogeneity1=mean(GLCMall.homom);
  G.Homogeneity2=mean(GLCMall.homop);
  G.IMC1=mean(GLCMall.inf1h);
  G.IMC2=mean(GLCMall.inf2h);
  G.IDMN=mean(GLCMall.idmnc);
  G.IDN=mean(GLCMall.indnc);
  G.InverseVariance=mean(GLCMall.invva);
  G.MaximumProbability=mean(GLCMall.maxpr);
  G.SumAverage=mean(GLCMall.savgh);
  G.SumEntropy=mean(GLCMall.senth);
  G.SumVariance=mean(GLCMall.svarh);
  G.Variance=mean(GLCMall.sosvh);
end

%% GLRL
function meanG=GLRL3d(I,T,nL)
% Input: 
%     I: input image
%     T: tumor ROI
%     nL - number of levels
% Output:
%     G - GLRL matrix (numLevel x maxLength x numDirections)
warning off %#ok<WNOFF>
mask=T;volume=I;
[iV,jV,kV] = find3d(mask);
boxBound(1,1) = min(iV);
boxBound(1,2) = max(iV);
boxBound(2,1) = min(jV);
boxBound(2,2) = max(jV);
boxBound(3,1) = min(kV);
boxBound(3,2) = max(kV);
maskBox = mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ROIbox = volume(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

ROIbox = double(ROIbox);
IT = ROIbox;
IT(~maskBox) = NaN;

IT= double(IT);
temp = IT(~isnan(IT));
u = mean(temp);
sigma = std(temp);
ITtemp = IT;
ITtemp(IT>(u + 3*sigma)) = NaN;
ITtemp(IT<(u - 3*sigma)) = NaN;

maskBox(isnan(ITtemp)) = 0;

IT= double(IT);
maxVal = max(IT(:));
minVal = min(IT(:));
IT = round((nL-1)*(IT-minVal)/(maxVal-minVal))+1;

nL=1:nL;

nLevel=length(nL);
if nLevel > 100
    adjust=10000;
else
    adjust=1000;
end
levelTemp=max(nL)+1;
IT(isnan(IT))=levelTemp;
nL=[nL,levelTemp];


% discretization
uniqueVol=round(nL*adjust)/adjust;
IT=round(IT*adjust)/adjust;
NL=length(nL) - 1;
% disp([NL size(uniqueVol) size(IT)]);

%INITIALIZATION
directions=13;  %13 directions in 3D
sizeV=size(IT);
numInit=ceil(max(sizeV)*sqrt(3)); % Max run length
GLRLMatrix=zeros(NL+1,numInit,directions);


% START COMPUTATION
% Directions [1,0,0], [0 1 0], [1 1 0] and [-1 1 0] : 2D directions
% (x:right-left, y:top-bottom, z:3rd dimension)  
if numel(size(IT)) == 3
    nComp=sizeV(3); % We can add-up the GLRLMatrixs taken separately in every image in the x-y plane
else
    nComp=1;
end
for i=1:nComp
    image=IT(:,:,i);
    uniqueIm=unique(image);
    NLtemp=length(uniqueIm);
    indexRow=zeros(NLtemp,1);
    temp=image;
    
    for j=1:NLtemp
        indexRow(j)=find(uniqueIm(j)==uniqueVol);
        image(temp==uniqueIm(j))=j;
    end
    
    % [1,0,0]
    GLRLMatrixtemp=rle_0(image,NLtemp);
    nRun=size(GLRLMatrixtemp,2);
    GLRLMatrix(indexRow(1:NLtemp),1:nRun,1)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,1)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
    
    % [0 1 0]
    GLRLMatrixtemp=rle_0(image',NLtemp);
    nRun=size(GLRLMatrixtemp,2);
    GLRLMatrix(indexRow(1:NLtemp),1:nRun,2)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,2)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
    
    % [1 1 0]
    seq=zigzag(image);
    GLRLMatrixtemp=rle_45(seq,NLtemp);
    nRun=size(GLRLMatrixtemp,2);
    GLRLMatrix(indexRow(1:NLtemp),1:nRun,3)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,3)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
    
    % [-1 1 0]
    seq=zigzag(fliplr(image));
    GLRLMatrixtemp=rle_45(seq,NLtemp);
    nRun=size(GLRLMatrixtemp,2);
    GLRLMatrix(indexRow(1:NLtemp),1:nRun,4)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,4)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
end

if numel(size(IT)) == 3 % 3D DIRECTIONS
    % Directions [0,0,1], [1 0 1] and [-1 0 1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    nComp=sizeV(1); % We can add-up the GLRLMatrixs taken separately in every image in the x-z plane
    image=zeros(sizeV(3),sizeV(2));
    for i=1:nComp
        for j=1:sizeV(3)
            image(j,1:end)=IT(i,1:end,j);
        end
        uniqueIm=unique(image);
        NLtemp=length(uniqueIm);
        indexRow=zeros(NLtemp,1);
        temp=image;
        for j=1:NLtemp
            indexRow(j)=find(uniqueIm(j)==uniqueVol);
            image(temp==uniqueIm(j))=j;
        end
        
        % [0,0,1]
        GLRLMatrixtemp=rle_0(image',NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,5)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,5)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
        
        % [1 0 1]
        seq=zigzag(image);
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,6)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,6)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
        
        % [-1 0 1]
        seq=zigzag(fliplr(image));
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,7)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,7)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
    end

    % Directions [0,1,1] and [0 -1 1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    nComp=sizeV(2); % We can add-up the GLRLMatrixs taken separately in every image in the y-z plane
    image=zeros(sizeV(1),sizeV(3));
    for i=1:nComp
        for j=1:sizeV(3)
            image(1:end,j)=IT(1:end,i,j);
        end
        uniqueIm=unique(image);
        NLtemp=length(uniqueIm);
        indexRow=zeros(NLtemp,1);
        temp=image;
        for j=1:NLtemp
            indexRow(j)=find(uniqueIm(j)==uniqueVol);
            image(temp==uniqueIm(j))=j;
        end
        
        % [0,1,1]
        seq=zigzag(image);
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,8)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,8)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
        
        % [0 -1 1]
        seq=zigzag(fliplr(image));
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,9)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,9)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
    end

    % Four corners: [1,1,1], [-1,1,1], [-1,1,-1], [1,1,-1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    image=zeros(sizeV(3),sizeV(2));
    temp=rand(sizeV(3),sizeV(2));
    diagTemp=spdiags(temp);
    szDiag=size(diagTemp);
    diagMat1=zeros(szDiag(1),szDiag(2),sizeV(1));
    diagMat2=zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
    for i=1:sizeV(1)
        for j=1:sizeV(3)
            image(j,1:end)=IT(i,1:end,j);
        end
        try
            diagMat1(:,:,i)=spdiags(image);
        catch
            % Add a column at the beginning to prevent errors
            temp=spdiags(image);
            numberDiff=abs(size(temp,2)-size(diagMat1,2));
            if mod(numberDiff,2) % Odd difference number
                temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
                diagMat1(:,:,i)=temp(:,1:end-1);
            else
                diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
            end
        end
        try
            diagMat2(:,:,i)=spdiags(fliplr(image));
        catch
            % Add a column at the beginning to prevent errors
            temp=spdiags(fliplr(image));
            numberDiff=abs(size(temp,2)-size(diagMat2,2));
            if mod(numberDiff,2) % Odd difference number
                temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
                diagMat2(:,:,i)=temp(:,1:end-1);
            else
                diagMat2(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
            end
        end
    end
    for j=1:szDiag(2)
        index=(diagMat1(:,j,1)~=0);
        nTemp=sum(index);
        image1=zeros(sizeV(1),nTemp);
        image2=zeros(sizeV(1),nTemp);
        for k=1:sizeV(1)
            image1(k,1:nTemp)=diagMat1(index(1:end),j,k)';
            image2(k,1:nTemp)=diagMat1(index(1:end),j,k)';
        end
        
        % 2 first corners
        uniqueIm=unique(image1);
        NLtemp=length(uniqueIm);
        indexRow=zeros(NLtemp,1);
        temp=image1;
        for i=1:NLtemp
            indexRow(i)=find(uniqueIm(i)==uniqueVol);
            image1(temp==uniqueIm(i))=i;
        end
        seq=zigzag(image1);
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,10)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,10)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
        seq=zigzag(fliplr(image1));
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,11)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,11)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
        
        % 2 last corners
        uniqueIm=unique(image2);
        NLtemp=length(uniqueIm);
        indexRow=zeros(NLtemp,1);
        temp=image2;
        for i=1:NLtemp
            indexRow(i)=find(uniqueIm(i)==uniqueVol);
            image2(temp==uniqueIm(i))=i;
        end
        seq=zigzag(image2);
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,12)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,12)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
        seq=zigzag(fliplr(image2));
        GLRLMatrixtemp=rle_45(seq,NLtemp);
        nRun=size(GLRLMatrixtemp,2);
        GLRLMatrix(indexRow(1:NLtemp),1:nRun,13)=GLRLMatrix(indexRow(1:NLtemp),1:nRun,13)+GLRLMatrixtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLMatrix
    end
end

% REMOVE UNECESSARY COLUMNS
GLRLMatrix(end,:,:)=[];
ind=find(sum(GLRLMatrix),1,'last');
GLRLMatrix(:,(ind+1):end,:)=[];


%% Textures
sz=size(GLRLMatrix);
G.SRE=0;    % 45. Short Run Emphasis (SRE)
G.LRE=0;    % 46. Long Run Emphasis (LRE)
G.GLN=0;    % 47. Gray-Level Nonuniformity (GLN)
G.RLN=0;    % 48. Run-Length Nonuniformity (RLN)
G.RP=0;     % 49. Run Percentage (RP)
G.LGLRE=0;  % 50. Low Gray-Level Run Emphasis (LGLRE)
G.HGLRE=0;  % 51. High Gray-Level Run Emphasis (HGLRE)
G.SRLGLE=0; % 52. Short Run Low Gray-Level Emphasis (SRLGLE)
G.SRHGLE=0; % 53. Short Run High Gray-Level Emphasis (SRHGLE)
G.LRLGLE=0; % 54. Long Run Low Gray-Level Emphasis (LRLGLE)
G.LRHGLE=0; % 55. Long Run High Gray-Level Emphasis (LRHGLE)

for i=1:sz(3)
  GLRLM=GLRLMatrix(:,:,i);
  
  nRuns=sum(GLRLM(:));
  cVect=1:sz(2); rVect=1:sz(1);% Row and column vectors
  [cMat,rMat]=meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLRLM
  pg=sum(GLRLM,2)'; % Gray-Level Run-Number Vector
  pr=sum(GLRLM); % Run-Length Run-Number Vector

  G.SRE=G.SRE+(pr*(cVect.^(-2))')/nRuns;
  G.LRE=G.LRE+(pr*(cVect.^2)')/nRuns; % cu wenli  
  G.GLN=G.GLN+sum(pg.^2)/nRuns;    
  G.RLN=G.RLN+sum(pr.^2)/nRuns;  
  G.RP=G.RP+nRuns/(pr*cVect');  
	G.LGLRE=G.LGLRE+(pg*(rVect.^(-2))')/nRuns;
  G.HGLRE=G.HGLRE+(pg*(rVect.^2)')/nRuns;
  G.SRLGLE=G.SRLGLE+sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^(-2))))/nRuns;
  G.SRHGLE=G.SRHGLE+sum(sum(GLRLM.*(rMat.^2).*(cMat.^(-2))))/nRuns; 
  G.LRLGLE=G.LRLGLE+sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^2)))/nRuns;
  G.LRHGLE=G.LRHGLE+sum(sum(GLRLM.*(rMat.^2).*(cMat.^2)))/nRuns;
end

meanG.ShortRunEmphasis=G.SRE;    % 45. Short Run Emphasis (SRE)
meanG.LongRunEmphasis=G.LRE;    % 46. Long Run Emphasis (LRE)
meanG.GrayLevelNonuniformity=G.GLN;    % 47. Gray-Level Nonuniformity (GLN)
meanG.RunLengthNonuniformity=G.RLN;    % 48. Run-Length Nonuniformity (RLN)
meanG.RunPercentage=G.RP;     % 49. Run Percentage (RP)
meanG.LowGrayLevelRunEmphasis=G.LGLRE;  % 50. Low Gray-Level Run Emphasis (LGLRE)
meanG.HighGrayLevelRunEmphasis=G.HGLRE;  % 51. High Gray-Level Run Emphasis (HGLRE)
meanG.ShortRunLowGrayLevelEmphasis=G.SRLGLE; % 52. Short Run Low Gray-Level Emphasis (SRLGLE)
meanG.ShortRunHighGrayLevelEmphasis=G.SRHGLE; % 53. Short Run High Gray-Level Emphasis (SRHGLE)
meanG.LongRunLowGrayLevelEmphasis=G.LRLGLE; % 54. Long Run Low Gray-Level Emphasis (LRLGLE)
meanG.LongRunHighGrayLevelEmphasis=G.LRHGLE; % 55. Long Run High Gray-Level Emphasis (LRHGLE)

end

function [iV,jV,kV] = find3d(mask3M)
  indV = find(mask3M(:));
  [iV,jV,kV] = fastind2sub(size(mask3M),indV);
  iV = iV';
  jV = jV';
  kV = kV';
end

function varargout = fastind2sub(siz,ndx)
warning off %#ok<WNOFF>  
nout = max(nargout,1);
  if length(siz)<=nout,
    siz = [siz ones(1,nout-length(siz))];
  else
    siz = [siz(1:nout-1) prod(siz(nout:end))];
  end
  n = length(siz);
  k = [1 cumprod(siz(1:end-1))];
  ndx = ndx - 1;
  for i = n:-1:1,
    varargout{i} = floor(ndx/k(i)) + 1;
    ndx = ndx - (varargout{i}-1) * k(i);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine structured by Jos van der Geest
function A = catstruct(varargin)
% CATSTRUCT   Concatenate or merge structures with different fieldnames
%   X = CATSTRUCT(S1,S2,S3,...) merges the structures S1, S2, S3 ...
%   into one new structure X. X contains all fields present in the various
%   structures. An example:
%
%     A.name = 'Me' ;
%     B.income = 99999 ;
%     X = catstruct(A,B) 
%     % -> X.name = 'Me' ;
%     %    X.income = 99999 ;
%
%   If a fieldname is not unique among structures (i.e., a fieldname is
%   present in more than one structure), only the value from the last
%   structure with this field is used. In this case, the fields are 
%   alphabetically sorted. A warning is issued as well. An axample:
%
%     S1.name = 'Me' ;
%     S2.age  = 20 ; S3.age  = 30 ; S4.age  = 40 ;
%     S5.honest = false ;
%     Y = catstruct(S1,S2,S3,S4,S5) % use value from S4
%
%   The inputs can be array of structures. All structures should have the
%   same size. An example:
%
%     C(1).bb = 1 ; C(2).bb = 2 ;
%     D(1).aa = 3 ; D(2).aa = 4 ;
%     CD = catstruct(C,D) % CD is a 1x2 structure array with fields bb and aa
%
%   The last input can be the string 'sorted'. In this case,
%   CATSTRUCT(S1,S2, ..., 'sorted') will sort the fieldnames alphabetically. 
%   To sort the fieldnames of a structure A, you could use
%   CATSTRUCT(A,'sorted') but I recommend ORDERFIELDS for doing that.
%
%   When there is nothing to concatenate, the result will be an empty
%   struct (0x0 struct array with no fields).
%
%   NOTE: To concatenate similar arrays of structs, you can use simple
%   concatenation: 
%     A = dir('*.mat') ; B = dir('*.m') ; C = [A ; B] ;

%   NOTE: This function relies on unique. Matlab changed the behavior of
%   its set functions since 2013a, so this might cause some backward
%   compatibility issues when dulpicated fieldnames are found.
%
%   See also CAT, STRUCT, FIELDNAMES, STRUCT2CELL, ORDERFIELDS

% version 4.1 (feb 2015), tested in R2014a
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% Created in 2005
% Revisions
%   2.0 (sep 2007) removed bug when dealing with fields containing cell
%                  arrays (Thanks to Rene Willemink)
%   2.1 (sep 2008) added warning and error identifiers
%   2.2 (oct 2008) fixed error when dealing with empty structs (thanks to
%                  Lars Barring)
%   3.0 (mar 2013) fixed problem when the inputs were array of structures
%                  (thanks to Tor Inge Birkenes).
%                  Rephrased the help section as well.
%   4.0 (dec 2013) fixed problem with unique due to version differences in
%                  ML. Unique(...,'last') is no longer the deafult.
%                  (thanks to Isabel P)
%   4.1 (feb 2015) fixed warning with narginchk
warning off %#ok<WNOFF>
narginchk(1,Inf) ;
N = nargin ;

if ~isstruct(varargin{end}),
    if isequal(varargin{end},'sorted'),
        narginchk(2,Inf) ;
        sorted = 1 ;
        N = N-1 ;
    else
        error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted".') ;
    end
else
    sorted = 0 ;
end

sz0 = [] ; % used to check that all inputs have the same size

% used to check for a few trivial cases
NonEmptyInputs = false(N,1) ; 
NonEmptyInputsN = 0 ;

% used to collect the fieldnames and the inputs
FN = cell(N,1) ;
VAL = cell(N,1) ;

% parse the inputs
for ii=1:N,
    X = varargin{ii} ;
    if ~isstruct(X),
        error('catstruct:InvalidArgument',['Argument #' num2str(ii) ' is not a structure.']) ;
    end
    
    if ~isempty(X),
        % empty structs are ignored
        if ii > 1 && ~isempty(sz0)
            if ~isequal(size(X), sz0)
                error('catstruct:UnequalSizes','All structures should have the same size.') ;
            end
        else
            sz0 = size(X) ;
        end
        NonEmptyInputsN = NonEmptyInputsN + 1 ;
        NonEmptyInputs(ii) = true ;
        FN{ii} = fieldnames(X) ;
        VAL{ii} = struct2cell(X) ;
    end
end

if NonEmptyInputsN == 0
    % all structures were empty
    A = struct([]) ;
elseif NonEmptyInputsN == 1,
    % there was only one non-empty structure
    A = varargin{NonEmptyInputs} ;
    if sorted,
        A = orderfields(A) ;
    end
else
    % there is actually something to concatenate
    FN = cat(1,FN{:}) ;    
    VAL = cat(1,VAL{:}) ;    
    FN = squeeze(FN) ;
    VAL = squeeze(VAL) ;
    
    
    [UFN,ind] = unique(FN, 'last') ;
    % If this line errors, due to your matlab version not having UNIQUE
    % accept the 'last' input, use the following line instead
    % [UFN,ind] = unique(FN) ; % earlier ML versions, like 6.5
    
    if numel(UFN) ~= numel(FN),
        warning('catstruct:DuplicatesFound','Fieldnames are not unique between structures.') ;
        sorted = 1 ;
    end
    
    if sorted,
        VAL = VAL(ind,:) ;
        FN = FN(ind,:) ;
    end
    
    A = cell2struct(VAL, FN);
    A = reshape(A, sz0) ; % reshape into original format
end

end

%% GLRL by Xunkai Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function oneglrlm = rle_0(si,NL)
% RLE   image gray level Run Length matrix for 0degree
%    
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
%  -------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007

warning off %#ok<WNOFF>
% Assure row number is exactly the gray level
[m,n]=size(si);

oneglrlm=zeros(NL,n);

for i=1:m
    x=si(i,:);
    % run length Encode of each vector
    index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
    len = diff([ 0 index ]); % run lengths
    val = x(index);          % run values
    temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
    oneglrlm = temp + oneglrlm; % accumulate each contribution
end
end
 
function oneglrlm = rle_45(seq,NL)
% RLE   image gray level Run Length matrix for 45 and 135
% This file is to handle the zigzag scanned sequence for 45 or 135 degree
% direction. Note for 135, just swap the left and the right colum
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
%  -------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007

warning off %#ok<WNOFF>
% Assure row number is exactly the gray leve;
% number of seqence
m =length(seq);
% number to store the possible max coloums
n = findmaxnum(seq);
%

oneglrlm=zeros(NL,n);

for i=1:m
    x=seq{i};
    % run length Encode of each vector
    index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
    len = diff([ 0 index ]); % run lengths
    val = x(index);          % run values
    temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
    oneglrlm = temp + oneglrlm; % accumulate each contribution
end
end

function seq = zigzag(SI)
%
%  Description:
%  ------------
%  This function is used to build the corresponding sequences of a given
%  scaled gray level image matrix from 45' degree direction. The whole process is using zigzag method
%  It can handle nonsquare image matrix
%
% Author:
% -------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
%
% History:
%  -------
% Creation: beta  Date: 01/11/2007
% Revision: 1.0   Date: 12/11/2007
% 
% Trick: all the sequence starts or ends lie on the boundary.

% initializing the variables
%----------------------------------
c = 1; % initialize colum indicator
r = 1; % initialize row   indicator
warning off %#ok<WNOFF>
rmin = 1; % row   boundary checker
cmin = 1; % colum boundary checker

rmax = size(SI, 1); % get row numbers
cmax = size(SI, 2); % get colum numbers

%
i = 1; % counter for current ith element
j = 1; % indicator for determining sequence interval

% intialize sequence mark
sq_up_begin=1;

sq_down_begin=1;

% % Output contain value and its flag status
% the first row contain value
% the second row contain its flag
output = zeros(1, rmax * cmax);
% sequence counter
%
% % Position Matrix
% position =zeros(1, rmax * cmax);
%----------------------------------

while ((r <= rmax) && (c <= cmax))

    % for current point, judge its zigzag direction up 45, or down 45, or
    % 0,or down 90

    if (mod(c + r, 2) == 0)      % up 45 direction
        %  if we currently walk to the left first colum
        if (r == rmin)
            % First, record current point
            output(i) = SI(r, c);
            % if we walk to right last colum
            if (c == cmax)
                % add row number move straight down 90
                r   = r + 1;
                sq_up_end = i;
                sq_down_begin = i+1;
                seq{j}=output(sq_up_begin:sq_up_end);
                j = j + 1;
                %

            else
                % Continue to move to next (1,c+1) point
                % This next point should be the begin point of next sequence
                c = c + 1;
                sq_up_end = i;
                sq_down_begin = i+1;

                seq{j}=output(sq_up_begin:sq_up_end);

                j = j + 1;

            end;

            % add couter
            i = i + 1;
            % if we currently walk to the last column
        elseif ((c == cmax) && (r < rmax))
            % first record the point
            output(i) = SI(r, c);
            % then move straight down to next row
            r = r + 1;
            
            sq_up_end = i;
            seq{j}=output(sq_up_begin:sq_up_end);
            sq_down_begin =i+1;
            j=j+1;
                        
            % add counter
            i = i + 1;
            % all other cases i.e. nonboundary points
        elseif ((r > rmin) && (c < cmax))
            output(i) = SI(r, c);
            % move to next up 45 point
            r = r - 1;
            c = c + 1;
            % add counter
            i = i + 1;
        end;
        % down 45 direction
    else
        % if we walk to the last row
        if ((r == rmax) && (c <= cmax))
            % firstly record current point
            output(i) = SI(r, c);
            % move right to next point
            c = c + 1;
            sq_down_end = i;
            seq{j}=output(sq_down_begin:sq_down_end);
            sq_up_begin =i+1;
            j = j + 1;
            % add counter
            i = i + 1;
            % if we walk to the first column
        elseif (c == cmin)
            %first record current point
            output(i) = SI(r, c);
            %
            if (r == rmax)
                c = c + 1;
                
                sq_down_end = i;
                seq{j}=output(sq_down_begin:sq_down_end);
                sq_up_begin =i+1;
                j = j + 1;

            else
                r = r + 1;
                % record sequence end
                sq_down_end = i;
                seq{j}=output(sq_down_begin:sq_down_end);
                sq_up_begin =i+1;
                j = j + 1;

            end;

            i = i + 1;
            % all other cases without boundary point
        elseif ((r < rmax) && (c > cmin))
            %
            output(i) = SI(r, c);
            %             position(i) = sub2ind(SI,r,c);
            r = r + 1;
            c = c - 1;
            % keep down_info
            i = i + 1;
        end;

    end;

    if ((r == rmax) && (c == cmax))          % bottom right element
        output(i) = SI(r, c);
        sq_end = i;
        seq{j}=output(sq_end);
        %         position(i) = sub2ind(SI,r,c);
        break
    end;
end;
end

function maxnum=findmaxnum(seq)
%
%  this function is obtain the maximum numbers of the given sequence
%  note the sequence is stored in cell mode
%
%
% See also zigzag
%
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
% ---------------------------------------------
% Creation: beta  Date: 01/10/2007
% Revision: 1.0   Date: 12/11/2007
%
warning off %#ok<WNOFF>
if iscell(seq)

    numseq=length(seq);
    maxnum=1;
    for i=1:numseq
        temp = seq{i};
        numseq = length(temp);
        if numseq > maxnum
            maxnum =numseq;
        end
    end
else
    error('I was only designed to handle cell sequence')
end
end


%% GLCM toolbox by Avinash Uppuluri
% http://www.mathworks.com/matlabcentral/fileexchange/22187-glcm-texture-features/content/GLCM_Features1.m
function [out] = GLCM_Features4(glcmin,pairs)
% This is an update of GLCM_Features2 (vectorized) without ismember()
%
% GLCM_Features2 helps to calculate the features from the different GLCMs
% that are input to the function. The GLCMs are stored in a i x j x n
% matrix, where n is the number of GLCMs calculated usually due to the
% different orientation and displacements used in the algorithm. Usually
% the values i and j are equal to 'NumLevels' parameter of the GLCM
% computing function graycomatrix(). Note that matlab quantization values
% belong to the set {1,..., NumLevels} and not from {0,...,(NumLevels-1)}
% as provided in some references
% http://www.mathworks.com/access/helpdesk/help/toolbox/images/graycomatrix
% .html
% 
% This vectorized version of GLCM_FEatures1.m reduces the 19 'for' loops
% used in the earlier code to 5 'for' loops
% http://blogs.mathworks.com/loren/2006/07/12/what-are-you-really-measuring
% /
% Using tic toc and cputime as in above discussion
%
% Although there is a function graycoprops() in Matlab Image Processing
% Toolbox that computes four parameters Contrast, Correlation, Energy,
% and Homogeneity. The paper by Haralick suggests a few more parameters
% that are also computed here. The code is not fully vectorized and hence
% is not an efficient implementation but it is easy to add new features
% based on the GLCM using this code. Takes care of 3 dimensional glcms
% (multiple glcms in a single 3D array)
% 
% If you find that the values obtained are different from what you expect 
% or if you think there is a different formula that needs to be used 
% from the ones used in this code please let me know. 
% A few questions which I have are listed in the link 
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/239608
%
%
%
% Features computed 
% Autocorrelation: [2]                      (out.autoc)
% Contrast: matlab/[1,2]                    (out.contr)
% Correlation: matlab                       (out.corrm)
% Correlation: [1,2]                        (out.corrp)
% Cluster Prominence: [2]                   (out.cprom)
% Cluster Shade: [2]                        (out.cshad)
% Dissimilarity: [2]                        (out.dissi)
% Energy: matlab / [1,2]                    (out.energ)
% Entropy: [2]                              (out.entro)
% Homogeneity: matlab                       (out.homom)
% Homogeneity: [2]                          (out.homop)
% Maximum probability: [2]                  (out.maxpr)
% Sum of squares: Variance [1]              (out.sosvh)
% Sum average [1]                           (out.savgh)
% Sum variance [1]                          (out.svarh)
% Sum entropy [1]                           (out.senth)
% Difference variance [1]                   (out.dvarh)
% Difference entropy [1]                    (out.denth)
% Information measure of correlation1 [1]   (out.inf1h)
% Informaiton measure of correlation2 [1]   (out.inf2h)
% Inverse difference (INV) is homom [3]     (out.homom)
% Inverse difference normalized (INN) [3]   (out.indnc) 
% Inverse difference moment normalized [3]  (out.idmnc)
%
% The maximal correlation coefficient was not calculated due to
% computational instability 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
%
% Formulae from MATLAB site (some look different from
% the paper by Haralick but are equivalent and give same results)
% Example formulae: 
% Contrast = sum_i(sum_j(  (i-j)^2 * p(i,j) ) ) (same in matlab/paper)
% Correlation = sum_i( sum_j( (i - u_i)(j - u_j)p(i,j)/(s_i.s_j) ) ) (m)
% Correlation = sum_i( sum_j( ((ij)p(i,j) - u_x.u_y) / (s_x.s_y) ) ) (p[2])
% Energy = sum_i( sum_j( p(i,j)^2 ) )           (same in matlab/paper)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + |i-j|) ) ) (as in matlab)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + (i-j)^2) ) ) (as in paper)
% 
% Where:
% u_i = u_x = sum_i( sum_j( i.p(i,j) ) ) (in paper [2])
% u_j = u_y = sum_i( sum_j( j.p(i,j) ) ) (in paper [2])
% s_i = s_x = sum_i( sum_j( (i - u_x)^2.p(i,j) ) ) (in paper [2])
% s_j = s_y = sum_i( sum_j( (j - u_y)^2.p(i,j) ) ) (in paper [2])
%
% 
% Normalize the glcm:
% Compute the sum of all the values in each glcm in the array and divide 
% each element by it sum
%
% Haralick uses 'Symmetric' = true in computing the glcm
% There is no Symmetric flag in the Matlab version I use hence
% I add the diagonally opposite pairs to obtain the Haralick glcm
% Here it is assumed that the diagonally opposite orientations are paired
% one after the other in the matrix
% If the above assumption is true with respect to the input glcm then
% setting the flag 'pairs' to 1 will compute the final glcms that would result 
% by setting 'Symmetric' to true. If your glcm is computed using the
% Matlab version with 'Symmetric' flag you can set the flag 'pairs' to 0
%
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973
% 2. L. Soh and C. Tsatsoulis, Texture Analysis of SAR Sea Ice Imagery
% Using Gray Level Co-Occurrence Matrices, IEEE Transactions on Geoscience
% and Remote Sensing, vol. 37, no. 2, March 1999.
% 3. D A. Clausi, An analysis of co-occurrence texture statistics as a
% function of grey level quantization, Can. J. Remote Sensing, vol. 28, no.
% 1, pp. 45-62, 2002
% 4. http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
%
%
% Example:
%
% Usage is similar to graycoprops() but needs extra parameter 'pairs' apart
% from the GLCM as input
% I = imread('circuit.tif');
% GLCM2 = graycomatrix(I,'Offset',[2 0;0 2]);
% stats = GLCM_features4(GLCM2,0)
% The output is a structure containing all the parameters for the different
% GLCMs
%
% [Avinash Uppuluri: avinash_uv@yahoo.com: Last modified: 04/05/2010]
warning off %#ok<WNOFF>
% If 'pairs' not entered: set pairs to 0 
if ((nargin > 2) || (nargin == 0))
   error('Too many or too few input arguments. Enter GLCM and pairs.');
elseif ( (nargin == 2) ) 
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
elseif (nargin == 1) % only GLCM is entered
    pairs = 0; % default is numbers and input 1 for percentage
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
       error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
end


format long e
if (pairs == 1)
    newn = 1;
    for nglcm = 1:2:size(glcmin,3)
        glcm(:,:,newn)  = glcmin(:,:,nglcm) + glcmin(:,:,nglcm+1);
        newn = newn + 1;
    end
elseif (pairs == 0)
    glcm = glcmin;
end

size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);

% checked 
out.autoc = zeros(1,size_glcm_3); % Autocorrelation: [2] 
out.contr = zeros(1,size_glcm_3); % Contrast: matlab/[1,2]
out.corrm = zeros(1,size_glcm_3); % Correlation: matlab
out.corrp = zeros(1,size_glcm_3); % Correlation: [1,2]
out.cprom = zeros(1,size_glcm_3); % Cluster Prominence: [2]
out.cshad = zeros(1,size_glcm_3); % Cluster Shade: [2]
out.ctend = zeros(1,size_glcm_3); % Cluster Tendency: %%Steven
out.dissi = zeros(1,size_glcm_3); % Dissimilarity: [2]
out.energ = zeros(1,size_glcm_3); % Energy: matlab / [1,2]
out.entro = zeros(1,size_glcm_3); % Entropy: [2]
out.homom = zeros(1,size_glcm_3); % Homogeneity: matlab
out.homop = zeros(1,size_glcm_3); % Homogeneity: [2]
out.invva = zeros(1,size_glcm_3); % Inverse variance %%Steven
out.maxpr = zeros(1,size_glcm_3); % Maximum probability: [2]

out.sosvh = zeros(1,size_glcm_3); % Sum of sqaures: Variance [1]
out.savgh = zeros(1,size_glcm_3); % Sum average [1]
out.svarh = zeros(1,size_glcm_3); % Sum variance [1]
out.senth = zeros(1,size_glcm_3); % Sum entropy [1]
out.dvarh = zeros(1,size_glcm_3); % Difference variance [4]
%out.dvarh2 = zeros(1,size_glcm_3); % Difference variance [1]
out.denth = zeros(1,size_glcm_3); % Difference entropy [1]
out.inf1h = zeros(1,size_glcm_3); % Information measure of correlation1 [1]
out.inf2h = zeros(1,size_glcm_3); % Informaiton measure of correlation2 [1]
%out.mxcch = zeros(1,size_glcm_3);% maximal correlation coefficient [1]
%out.invdc = zeros(1,size_glcm_3);% Inverse difference (INV) is homom [3]
out.indnc = zeros(1,size_glcm_3); % Inverse difference normalized (INN) [3]
out.idmnc = zeros(1,size_glcm_3); % Inverse difference moment normalized [3]

glcm_sum  = zeros(size_glcm_3,1);
glcm_mean = zeros(size_glcm_3,1);
glcm_var  = zeros(size_glcm_3,1);

% http://www.fp.ucalgary.ca/mhallbey/glcm_mean.htm confuses the range of 
% i and j used in calculating the means and standard deviations.
% As of now I am not sure if the range of i and j should be [1:Ng] or
% [0:Ng-1]. I am working on obtaining the values of mean and std that get
% the values of correlation that are provided by matlab.
u_x = zeros(size_glcm_3,1);
u_y = zeros(size_glcm_3,1);
s_x = zeros(size_glcm_3,1);
s_y = zeros(size_glcm_3,1);

% checked p_x p_y p_xplusy p_xminusy
p_x = zeros(size_glcm_1,size_glcm_3); % Ng x #glcms[1]  
p_y = zeros(size_glcm_2,size_glcm_3); % Ng x #glcms[1]
p_xplusy = zeros((size_glcm_1*2 - 1),size_glcm_3); %[1]
p_xminusy = zeros((size_glcm_1),size_glcm_3); %[1]
% checked hxy hxy1 hxy2 hx hy
hxy  = zeros(size_glcm_3,1);
hxy1 = zeros(size_glcm_3,1);
hx   = zeros(size_glcm_3,1);
hy   = zeros(size_glcm_3,1);
hxy2 = zeros(size_glcm_3,1);

corm = zeros(size_glcm_3,1);
corp = zeros(size_glcm_3,1);

for k = 1:size_glcm_3
    
    glcm_sum(k) = sum(sum(glcm(:,:,k)));
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k); % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k)); % compute mean after norm
    glcm_var(k)  = (std2(glcm(:,:,k)))^2;
    
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            p_x(i,k) = p_x(i,k) + glcm(i,j,k); 
            p_y(i,k) = p_y(i,k) + glcm(j,i,k); % taking i for j and j for i
            %if (ismember((i + j),[2:2*size_glcm_1])) 
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcm(i,j,k);
            %end
            %if (ismember(abs(i-j),[0:(size_glcm_1-1)])) 
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) +...
                    glcm(i,j,k);
            %end
        end
    end
    
end

% marginal probabilities are now available [1]
% p_xminusy has +1 in index for matlab (no 0 index)
% computing sum average, sum variance and sum entropy:


%Q    = zeros(size(glcm));

i_matrix  = repmat([1:size_glcm_1]',1,size_glcm_2);
j_matrix  = repmat([1:size_glcm_2],size_glcm_1,1);
% i_index = [ 1 1 1 1 1 .... 2 2 2 2 2 ... ]
i_index   = j_matrix(:);
% j_index = [ 1 2 3 4 5 .... 1 2 3 4 5 ... ]
j_index   = i_matrix(:);
xplusy_index = [1:(2*(size_glcm_1)-1)]';
xminusy_index = [0:(size_glcm_1-1)]';
mul_contr = abs(i_matrix - j_matrix).^2;

%% for calculating inverse variance
mul_contr_i_notequal_j = mul_contr; %Steven
for ij=1:size(mul_contr_i_notequal_j,1)
  mul_contr_i_notequal_j(ij,ij)=1;
end
glcm_i_notequal_j=glcm;
for kk=1:size_glcm_3
  for ij=1:size(glcm,1)
    glcm_i_notequal_j(ij,ij,k)=0;
  end
end
%%

mul_dissi = abs(i_matrix - j_matrix);
%div_homop = ( 1 + mul_contr); % used from the above two formulae
%div_homom = ( 1 + mul_dissi);

for k = 1:size_glcm_3 % number glcms
    
  
    out.contr(k) = sum(sum(mul_contr.*glcm(:,:,k)));
    out.dissi(k) = sum(sum(mul_dissi.*glcm(:,:,k)));
    out.energ(k) = sum(sum(glcm(:,:,k).^2));
    out.entro(k) = - sum(sum((glcm(:,:,k).*log(glcm(:,:,k) + eps))));
    out.homom(k) = sum(sum((glcm(:,:,k)./( 1 + mul_dissi))));
    out.homop(k) = sum(sum((glcm(:,:,k)./( 1 + mul_contr))));

    %     out.invva(k) = sum(sum((glcm(:,:,k)./( mul_contr)))); %%Steven:
    %     This is incorrect. The line follows is correct.
    out.invva(k) = sum(sum((glcm_i_notequal_j(:,:,k)./(mul_contr_i_notequal_j)))); %%Steven
    
    % [1] explains sum of squares variance with a mean value;
    % the exact definition for mean has not been provided in 
    % the reference: I use the mean of the entire normalized glcm     
    out.sosvh(k) = sum(sum(glcm(:,:,k).*((i_matrix - glcm_mean(k)).^2)));
    out.indnc(k) = sum(sum(glcm(:,:,k)./( 1 + (mul_dissi./size_glcm_1) )));
    out.idmnc(k) = sum(sum(glcm(:,:,k)./( 1 + (mul_contr./(size_glcm_1^2)))));
    out.maxpr(k) = max(max(glcm(:,:,k)));
    
    u_x(k)       = sum(sum(i_matrix.*glcm(:,:,k))); 
    u_y(k)       = sum(sum(j_matrix.*glcm(:,:,k))); 
    % using http://www.fp.ucalgary.ca/mhallbey/glcm_variance.htm for s_x
    % s_y : This solves the difference in value of correlation and might be
    % the right value of standard deviations required 
    % According to this website there is a typo in [2] which provides
    % values of variance instead of the standard deviation hence a square
    % root is required as done below:
    s_x(k)  = (sum(sum( ((i_matrix - u_x(k)).^2).*glcm(:,:,k) )))^0.5;
    s_y(k)  = (sum(sum( ((j_matrix - u_y(k)).^2).*glcm(:,:,k) )))^0.5;
    
   corp(k) = sum(sum((i_matrix.*j_matrix.*glcm(:,:,k))));
   corm(k) = sum(sum(((i_matrix - u_x(k)).*(j_matrix - u_y(k)).*glcm(:,:,k)))); 
   
   out.autoc(k) = corp(k);
   out.corrp(k) = (corp(k) - u_x(k)*u_y(k))/(s_x(k)*s_y(k));
   out.corrm(k) = corm(k) / (s_x(k)*s_y(k)); 
   
   out.cprom(k) = sum(sum(((i_matrix + j_matrix - u_x(k) - u_y(k)).^4).*...
                glcm(:,:,k))); 
   out.cshad(k) = sum(sum(((i_matrix + j_matrix - u_x(k) - u_y(k)).^3).*...
                glcm(:,:,k)));        
   out.ctend(k) = sum(sum(((i_matrix + j_matrix - u_x(k) - u_y(k)).^2).*...
                glcm(:,:,k)));   %%Steven
            

    

   out.savgh(k) = sum((xplusy_index + 1).*p_xplusy(:,k));
   % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
   out.senth(k) =  - sum(p_xplusy(:,k).*...
            log(p_xplusy(:,k) + eps));
   
    % compute sum variance with the help of sum entropy
    out.svarh(k) = sum((((xplusy_index + 1) - out.senth(k)).^2).*...
        p_xplusy(:,k));
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)    
    
    % compute difference variance, difference entropy,
    % out.dvarh2(k) = var(p_xminusy(:,k));
    % but using the formula in 
    % http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
    % we have for dvarh
    out.denth(k) = - sum((p_xminusy(:,k)).*...
        log(p_xminusy(:,k) + eps));
    out.dvarh(k) = sum((xminusy_index.^2).*p_xminusy(:,k));
    
    % compute information measure of correlation(1,2) [1]
    hxy(k) = out.entro(k);
    glcmk  = glcm(:,:,k)';
    glcmkv = glcmk(:);
    
    hxy1(k) =  - sum(glcmkv.*log(p_x(i_index,k).*p_y(j_index,k) + eps));
    hxy2(k) =  - sum(p_x(i_index,k).*p_y(j_index,k).*...
        log(p_x(i_index,k).*p_y(j_index,k) + eps));
     hx(k) = - sum(p_x(:,k).*log(p_x(:,k) + eps));
     hy(k) = - sum(p_y(:,k).*log(p_y(:,k) + eps));   
    
    out.inf1h(k) = ( hxy(k) - hxy1(k) ) / ( max([hx(k),hy(k)]) );
    out.inf2h(k) = ( 1 - exp( -2*( hxy2(k) - hxy(k) ) ) )^0.5;
    
    %     eig_Q(k,:)   = eig(Q(:,:,k));
    %     sort_eig(k,:)= sort(eig_Q(k,:),'descend');
    %     out.mxcch(k) = sort_eig(k,2)^0.5;
    % The maximal correlation coefficient was not calculated due to
    % computational instability 
    % http://murphylab.web.cmu.edu/publications/boland/boland_node26.html    
end
end