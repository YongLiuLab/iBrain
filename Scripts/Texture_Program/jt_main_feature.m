function jt_main_feature(t2imglist,maskimglist,out_dir)
% pause(7200)
% warning off %#ok<WNOFF>
resdir = strcat(out_dir);
if isempty(dir(resdir))
    mkdir(resdir)
end
% t2imglist = textread(input_txt,'%s');
% maskimglist =  textread(input_mask,'%s');
fea_nm_1st ={'energy';'energy1';'kurtosis';'maximum';'mean';'mad';'median';'minimum';'range';...
    'rms';'skewness';'std';'entropy';'uniformity';'var'};
fea_nm_2nd ={'Area';'Volume';'Compactness1';'Compactness2';'Max3dDiam';'SphericalDisprop';'Spherity';'Surf2VolRatio';};

fea_nm_3rd = { 'Autocorrelation';'ClusterProminence';'ClusterShade';'ClusterTendency';...
    'Contrast';'Correlation';'DifferenceEntropy';'Dissimilarity';'Energy';'Entropy';'Homogeneity1';'Homogeneity2';'IMC1';'IMC2';...
    'IDMN';'IDN';'InverseVariance';'MaximumProbability';'SumAverage';'SumEntropy';'SumVariance';'Variance';'ShortRunEmphasis';...
    'LongRunEmphasis';'GrayLevelNonuniformity';'RunLengthNonuniformity';'RunPercentage';'LowGrayLevelRunEmphasis';'HighGrayLevelRunEmphasis';...
    'ShortRunLowGrayLevelEmphasis';'ShortRunHighGrayLevelEmphasis';'LongRunLowGrayLevelEmphasis';'LongRunHighGrayLevelEmphasis';...
    };



for i = 1:length(t2imglist)
    i
    [pth t2nm]= (fileparts(t2imglist{i}));
    [pth masknm] = (fileparts(maskimglist{i}));
    if ~strcmp(t2nm,masknm)
        fprintf('please double check the subject %d with name %s, for the t2 img and mask img not match\n',i,t2nm);
    end
    t2_nm = t2imglist{i};
    mask_nm = maskimglist{i};
    
    [tu_img t2_img roi_img] = jt_read_mask(mask_nm,t2_nm);
    [x_minax,y_minax,z_minax] = cut_img(tu_img);
    tu_img = tu_img(x_minax(1):x_minax(2),y_minax(1):y_minax(2),z_minax(1):z_minax(2));
    t2_img=t2_img(x_minax(1):x_minax(2),y_minax(1):y_minax(2),z_minax(1):z_minax(2));
    fea_1st(i) = jt_1st_feature(t2_img,tu_img);
    fea_2nd(i) =jt_2nd_ShapeSize(t2_img,tu_img);
    fea_3rd(i)=jt_3rd_Texture(t2_img,tu_img,32,16);
end

%%%% write feature to txt file
%%
fid = fopen([resdir,'Feature_1_2.txt'],'w+');

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
s = strcat(s,'\n');

fprintf(fid,s);


for i = 1:length(t2imglist)
    [pth t2nm]= (fileparts(t2imglist{i}));
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
    
    s = strcat(s,'\n');
    
    fprintf(fid,s);
    
end
fclose(fid);
%%
fid = fopen([resdir,'Feature_3.txt'],'w+');

if fid==-1
    err=strcat('can not open a txt file\n');
    error(err);
end
s = strcat('Subject');


for i = 1:length(fea_nm_3rd);
    s = strcat(s,'\t',fea_nm_3rd{i});
end
s = strcat(s,'\n');

fprintf(fid,s);


for i = 1:length(t2imglist)
    [pth t2nm]= (fileparts(t2imglist{i}));
    s = strcat(t2nm);
    
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
