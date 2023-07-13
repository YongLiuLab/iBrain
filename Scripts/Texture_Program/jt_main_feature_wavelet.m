function jt_main_feature_wavelet(input_txt,input_mask,out_dir)

warning off %#ok<WNOFF>
1
n = 1;                   % Decomposition Level
w = 'sym4';              % Near symmetric wavelet

% % t2imglist = textread('I:\MRI_Glioma_Tiantan\T13D_T2_FLAIR_data\t2list.txt','%s');
% % maskimglist =  textread('I:\MRI_Glioma_Tiantan\T13D_T2_FLAIR_data\masklist.txt','%s');
resdir = strcat(out_dir,'\WT\');
if isempty(dir(resdir))
    mkdir(resdir)
end
mkdir(resdir);
t2imglist = textread(input_txt,'%s');
maskimglist =  textread(input_mask,'%s');
fea_nm_1st ={'size';'energy';'kurtosis';'maximum';'mean';'mad';'median';'minimum';'range';...
    'rms';'skewness';'std';'entropy';'uniformity';'var'};

% fea_nm_1st ={'energy';'kurtosis';'maximum';'mean';'mad';'median';'minimum';'range';...
%     'rms';'skewness';'std';'entropy';'uniformity';'var'};
fea_nm_2nd ={'Area';'Volume';'Compactness1';'Compactness2';'Max3dDiam';'SphericalDisprop';'Spherity';'Surf2VolRatio';};

fea_nm_3rd = { 'Autocorrelation';'ClusterProminence';'ClusterShade';'ClusterTendency';...
    'Contrast';'Correlation';'DifferenceEntropy';'Dissimilarity';'Energy';'Entropy';'Homogeneity1';'Homogeneity2';'IMC1';'IMC2';...
    'IDMN';'IDN';'InverseVariance';'MaximumProbability';'SumAverage';'SumEntropy';'SumVariance';'Variance';'ShortRunEmphasis';...
    'LongRunEmphasis';'GrayLevelNonuniformity';'RunLengthNonuniformity';'RunPercentage';'LowGrayLevelRunEmphasis';'HighGrayLevelRunEmphasis';...
    'ShortRunLowGrayLevelEmphasis';'ShortRunHighGrayLevelEmphasis';'LongRunLowGrayLevelEmphasis';'LongRunHighGrayLevelEmphasis';...
    };

for i = 1:length(t2imglist)
    i
    [pth t2nm]= (fileparts((t2imglist{i})));
    [pth masknm] = (fileparts((maskimglist{i})));
%     if ~strcmp(t2nm,masknm)
%         fprintf('please double check the subject %d with name %s, for the t2 img and mask img not match\n',i,t2nm);
%     end
    t2_nm = t2imglist{i};
    mask_nm = maskimglist{i};
    
    [tu_img t2_img roi_img] = jt_read_mask(mask_nm,t2_nm);
    WT = wavedec3(roi_img,n,w);    % Multilevel 3D wavelet decomposition.
    Rec = cell(1,8);
    Rec{1} = waverec3(WT,'aaa');
    Rec{2} = waverec3(WT,'aad');
    Rec{3} = waverec3(WT,'ada');
    Rec{4} = waverec3(WT,'add');
    Rec{5} = waverec3(WT,'daa');
    Rec{6} = waverec3(WT,'dad');
    Rec{7} = waverec3(WT,'dda');
    Rec{8} = waverec3(WT,'ddd');
    for kk = 1:8
        
        fea_1st(i,kk) = jt_1st_feature(Rec{kk},tu_img);
        fea_2nd(i,kk) =jt_2nd_ShapeSize(Rec{kk},tu_img);
        fea_3rd(i,kk)=jt_3rd_Texture(Rec{kk},tu_img,32,16);
        
    end
    %
    %     fea_1st(i) = jt_1st_feature(t2_img,tu_img);
end



%%%% write feature to txt file

fea1 = fea_1st;
fea2=fea_2nd;
fea3=fea_3rd;
%%
for kk = 1:8
   
    filename= strcat([resdir,'wavelet_feature_',num2str(kk),'.txt']);
    fid = fopen(filename,'w+');
    
    if fid==-1
        err=strcat('can not open a txt file\n');
        error(err);
    end
    s = strcat('Subject');
    for i = 1:length(fea_nm_1st);
        s = strcat(s,'\t',fea_nm_1st{i});
    end
    
    for i = 1:length(fea_nm_2nd);
        s = strcat(s,'\t',fea_nm_2nd{i});
    end

    
    for i = 1:length(fea_nm_3rd);
        s = strcat(s,'\t',fea_nm_3rd{i});
    end
    s = strcat(s,'\n');
    
    
    fprintf(fid,s);
    
     fea_1st= fea1(:,kk);
    fea_2nd= fea2(:,kk);
    fea_3rd= fea3(:,kk);
    
    
    for i = 1:length(t2imglist)
        [pth t2nm]= fileparts(t2imglist{i});
        s = strcat(t2nm);
        
        for k = 1:length(fea_nm_1st);
            data = fea_1st(i).(fea_nm_1st{k});
            
            if data>1000000
                s = strcat(s,'\t',num2str(data,'%10.4e'));
            else
                s = strcat(s,'\t',num2str(data,'%10.4f'));
            end
        end
        
        for k = 1:length(fea_nm_2nd);
            data = fea_2nd(i).(fea_nm_2nd{k});
            
            if data>1000000
                s = strcat(s,'\t',num2str(data,'%10.4e'));
            else
                s = strcat(s,'\t',num2str(data,'%10.4f'));
            end
        end
        
        
        for k = 1:length(fea_nm_3rd);
            data = fea_3rd(i).(fea_nm_3rd{k});
            
            if isnan(data)
                s = strcat(s,'\t',num2str(-1.234,'%10.4f'));
            else
                if data>1000000
                    s = strcat(s,'\t',num2str(data,'%10.4e'));
                else
                    s = strcat(s,'\t',num2str(data,'%10.4f'));
                end
            end
        end
        
        s = strcat(s,'\n');
        
        fprintf(fid,s);
        
    end
    fclose(fid);
end
