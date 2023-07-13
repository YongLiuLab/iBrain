from docxtpl import DocxTemplate,InlineImage,RichText
from docx.shared import Mm, Inches, Pt
import pandas as pd
import numpy as np
import scipy.io as scio
import csv
import nibabel as nib
import os
import argparse
import datetime
import scipy.io as sio
import matplotlib
import matplotlib.colors as mcolors
import re
from nilearn import plotting
from matplotlib import pyplot as plt
from PIL import Image
from scipy.stats import gaussian_kde
from collections import Counter
from openpyxl import load_workbook

mm_to_cm_scale_rate=1000
reference = {
                "wmVolume": [422142,547709],
                "gmVolume": [549124,657656],
                "totalVolume": [979473,1197158],
                "wmPer": [42.27, 46.65],
                "gmPer": [53.35, 57.73],
                "wmADper": [43.16, 47.93],
                "gmADper": [52.07, 56.84],

                "hippoLeft1Volume": [1888.83,2414.38],
                "hippoRight1Volume": [1955.12,2437.3],
                "hippoLeft2Volume": [2026.42,2601.32],
                "hippoRight2Volume": [1526.9,1590.35],
                "hippoAllVolume": [7565.58,9234.23],
}

stdDistribution = {
    "NCgmMean": 603389.93,
    "NCgmStd": 54265.81,
    "ADgmMean": 556724.19,
    "ADgmStd": 59083.45,

    "NChippoVolumeMean": 8399.91,
    "NChippoVolumeStd": 834.33,
    "ADhippoVolumeMean": 6873.39,
    "ADhippoVolumeStd": 1027.67,

    "NCgmMeanPer": 55.54,
    "NCgmStdPer": 2.19,
    "ADgmMeanPer": 54.45,
    "ADgmStdPer": 2.38,
}


#generate M by N slice images
def draw2DSlices(origimgPath, savePath,n_row, n_col):
    total_img_num = int(n_row*n_col)
    orig_img = nib.load(origimgPath).get_fdata()
    x_resolution = orig_img.shape[0]
    y_resolution = orig_img.shape[1]
    slice_num = orig_img.shape[2]
    start_slice = 0
    end_slice = slice_num-total_img_num
    slice_gap = int((end_slice-start_slice)/total_img_num)
    transversal_orig = np.rot90(orig_img, 1)
    plt.figure(figsize=(x_resolution/10*n_col, y_resolution/10*n_row))
    num = 0
    for ax in range(total_img_num):
        num = num+1
        current_slice=int(start_slice+num*slice_gap)
        plt.subplot(n_row, n_col, ax + 1)
        plt.axis("off")
        plt.imshow(transversal_orig[:, :, current_slice], cmap='gray', aspect='auto')
    plt.subplots_adjust(0, 0, 1, 1,hspace=0, wspace=0)
    plt.savefig(savePath,  bbox_inches='tight', pad_inches=0, dpi=72)
    plt.close()


#generate white matter(in green) and grey matter(in red) marked segment image
def drawSegment(origimgPath, gmimgPath, wmimgPath, savePath):
    orig_img = nib.load(origimgPath)
    wm_img_orig = nib.load(wmimgPath)
    wm_img = np.array(nib.load(wmimgPath).get_fdata())
    gm_img = np.array(nib.load(gmimgPath).get_fdata())
    zero_data = np.zeros_like(gm_img)
    combined_roi = np.array([zero_data, gm_img, wm_img])
    combined_roi = np.argmax(combined_roi,axis=0)
    combined_roi = nib.Nifti1Image(combined_roi, wm_img_orig.affine, wm_img_orig.header)
    light_green = (11 / 255, 254 / 255, 4 / 255)
    light_red = (247 / 255, 4 / 255, 3 / 255)
    discreteCmap = mcolors.ListedColormap([light_red, light_green])
    plotting.plot_roi(roi_img=combined_roi, bg_img=orig_img, draw_cross=False, dim=0, cmap=discreteCmap,output_file=savePath)



def drawRisk(riskimgPath,savePath):
    risk_img=nib.load(riskimgPath)
    plotting.plot_glass_brain(risk_img, threshold=2, colorbar=True, plot_abs=False,output_file=savePath)


#generate white matter(in green) and grey matter(in red) marked segment image
def drawRoiSlice(roiPath, bgPath, savePath, n_row, n_col):
    roiImg=nib.load(roiPath)
    transversal_roiImg = np.rot90(roiImg.get_fdata(), 1)
    bgImg=nib.load(bgPath)
    transversal_bgImg = np.rot90(bgImg.get_fdata(), 1)
    total_img_num = int(n_row * n_col)
    bgImg_data = nib.load(bgPath).get_fdata()
    #x_resolution = bgImg_data.shape[0]
    #y_resolution = bgImg_data.shape[1]
    x_resolution = 369
    y_resolution = 520
    slice_num = bgImg_data.shape[2]
    start_slice = 0
    end_slice = slice_num - total_img_num
    slice_gap = int((end_slice-start_slice)/total_img_num)
    fig, axes = plt.subplots(nrows=n_row, ncols=n_col, figsize=(x_resolution/10*n_col, y_resolution/10*n_row))
    for i in range(n_row):
        for j in range(n_col):
            slice_to_plot = start_slice + (i * n_col + j)*slice_gap
            axes[i][j].axis("off")
            axes[i][j].imshow(transversal_bgImg[:, :, slice_to_plot], cmap='gray', aspect='auto')
            temp_slice_mask = (transversal_roiImg[:, :, slice_to_plot]<-1.96)
            rgba = np.zeros((*transversal_roiImg[:, :, slice_to_plot].shape, 4), dtype=np.uint8)
            rgba[temp_slice_mask, :3] = [255, 0, 0]
            rgba[..., 3] = np.where(temp_slice_mask, 255, 0).astype(np.uint8)
            axes[i][j].imshow(rgba, aspect='auto',cmap='Reds')
    plt.subplots_adjust(0, 0, 1, 1, hspace=0, wspace=0)
    plt.savefig(savePath,  bbox_inches='tight', pad_inches=0, dpi=72)
    plt.close()

def drawRoi(imgPath,roiPath,savePath):
    plotting.plot_roi(roiPath, bg_img=imgPath, draw_cross=False, dim=0, cmap=plt.cm.autumn,output_file=savePath)


def drawDis(u,sig,subject,savePath,u2=None,sig2=None,show_max = 'Subject GM volume', xlabel="Value",ylabel="Probability"):
    import math
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    x = np.linspace(u - 4 * sig, u + 4 * sig, 50)
    stdG = lambda x, u, sig: np.exp(-(x - u) ** 2 / (2 * sig ** 2)) / (math.sqrt(2 * math.pi) * sig)
    y_sig = stdG(x, u, sig)
    figure = plt.figure(figsize=[8,5])
    plt.plot(x, y_sig, "r-", linewidth=2, label="Normal Population")
    plt.fill_between(x, 0, y_sig, alpha=0.2, color="orange")
    plt.axvline(x=subject, color='r', linestyle='-')
    plt.plot(subject, stdG(subject, u, sig), 'ks')
    loc = (subject + subject / 100, stdG(subject, u, sig))
    plt.annotate(show_max, xytext=loc, xy=loc)

    if u2 is not None:
        x2 = np.linspace(u2 - 4 * sig2, u2 + 4 * sig2, 50)
        y_sig2 = stdG(x2, u2, sig2)
        plt.plot(x2, y_sig2, "b-", linewidth=2, label="AD Population")
        plt.fill_between(x2, 0, y_sig2, alpha=0.2,  color="lightblue")
        plt.plot(subject, stdG(subject, u2, sig2), 'ks')
        loc2 = (subject + subject/100, stdG(subject, u2, sig2))
        plt.annotate(show_max, xytext=loc2, xy=loc2)

    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.xlabel(xlabel,fontdict={'size':14})
    plt.ylabel(ylabel,fontdict={'size':14})
    plt.legend(loc='upper right',prop={'size':14})
    plt.savefig(savePath)
    plt.close()

def drawHis(ref,subject,savePath,ref2=None,show_max = 'subject healthy score',xlabel="value",ylabel="probability"):
    arr = ref.ravel()
    hist, bin_edges = np.histogram(arr, density=True)
    kde = gaussian_kde(arr)
    x = np.linspace(arr.min(), arr.max(), 100)
    plt.plot(x, kde(x), "r-", linewidth=2, label="Normal Population")
    plt.fill_between(x, 0, kde(x), alpha=0.2, color='orange')
    plt.axvline(x=subject, color='r', linestyle='-')
    plt.plot(subject, float(kde.evaluate(subject)), 'ks')
    subject_loc = (subject, float(kde.evaluate(subject)))
    plt.annotate(show_max, xytext=subject_loc, xy=subject_loc)
    if ref2 is not None:
        arr2= ref2.ravel()
        hist2, bin_edges2 = np.histogram(arr2, density=True)
        kde2 = gaussian_kde(arr2)
        x2 = np.linspace(arr2.min(), arr.max(), 100)
        plt.plot(x2, kde2(x2), "b-", linewidth=2, label="AD Population")
        plt.fill_between(x2, 0, kde2(x2), alpha=0.2, color='lightblue')
        plt.plot(subject, float(kde2.evaluate(subject)), 'ks')
        subject_loc2 = (subject, float(kde2.evaluate(subject)))
        plt.annotate(show_max, xytext=subject_loc2, xy=subject_loc2)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlabel(xlabel,fontdict={'size':14})
    plt.ylabel(ylabel,fontdict={'size':14})
    plt.legend(loc='upper left', prop={'size': 14})
    plt.savefig(savePath)
    plt.close()

def combineBilateralRegionName(inputReport):
    if inputReport.strip():
        delimiters = ["、"]
        allregionNames = re.split("|".join(map(re.escape, delimiters)), inputReport)
        no_left_allregionNames=[s.replace("_L", "") for s in allregionNames]
        no_leftright_allregionNames = [s.replace("_R", "") for s in no_left_allregionNames]
        counter = Counter(no_leftright_allregionNames)
        duplicates = [k for k, v in counter.items() if v > 1]
        left_list = ["_L"+s for s in duplicates]
        right_list = ["_R"+s for s in duplicates]
        bilateral_list = [s for s in duplicates]
        left_removed_list = [s for s in allregionNames if s not in left_list]
        leftright_removed_list = [s for s in left_removed_list if s not in right_list]
        leftright_removed_list_with_comma = [s+"、" for s in leftright_removed_list]
        bilateral_list_comma = [s+"、" for s in bilateral_list]
        output_list = bilateral_list_comma + leftright_removed_list_with_comma
        output_report = "".join(output_list)
    else:
        output_report = ""
    return output_report

def extractSameRegionInReport(inputReport1,inputReport2):
    if inputReport1.strip() and inputReport2.strip():
        delimiters = ["、"]
        allregionNames1 = re.split("|".join(map(re.escape, delimiters)), inputReport1)
        allregionNames2 = re.split("|".join(map(re.escape, delimiters)), inputReport2)
        output_list = set(allregionNames1) & set(allregionNames2)
        output = "".join(list(output_list))
    else:
        output = ""
    return output

def reportBuild(iBrainPath, WorkingDir, SubjectID, ReportWholeBrain, ReportInIndividualSpace, MiniOutput,savePath,writeMode,tempFCsession):
    #generate report related file path based on report mode#

    if MiniOutput=='0':#convert matlab value to python boolean
        MiniOutput=False
    else:
        MiniOutput=True

    if ReportWholeBrain=='0':#convert matlab value to python boolean
        ReportWholeBrain=False
    else:
        ReportWholeBrain=True

    if ReportInIndividualSpace=='0':
        ReportInIndividualSpace=False
    else:
        ReportInIndividualSpace=True

    templatePath = os.path.join(iBrainPath, 'Report_template')
    modelPath = os.path.join(iBrainPath,'model_data')
    atlasPath = os.path.join(iBrainPath,'Atlas')
    hippoPath = os.path.join(atlasPath, 'hippo_mask.nii')
    if writeMode == 'T1':#only report T1
        if MiniOutput:
            gmPath = os.path.join(WorkingDir, 'T1Img', SubjectID, 'mri', 'mwp1'+SubjectID+'.nii')
            wmPath = os.path.join(WorkingDir, 'T1Img', SubjectID, 'mri', 'mwp2'+SubjectID+'.nii')
            originImagePath = os.path.join(WorkingDir, 'T1Img', SubjectID,'DNR_'+SubjectID+'.nii')

        else:
            gmPath = os.path.join(WorkingDir, 'T1ImgS', SubjectID, 'mri', 'mwp1'+SubjectID+'.nii')
            wmPath = os.path.join(WorkingDir, 'T1ImgS', SubjectID, 'mri', 'mwp2'+SubjectID+'.nii')
            originImagePath = os.path.join(WorkingDir, 'T1ImgC', SubjectID, 'DNR_'+SubjectID+'.nii')

        featurePath = os.path.join(WorkingDir, 'Results','T1ImgR2SN', SubjectID, 'features.csv')
        riskimgPath = os.path.join(WorkingDir, 'Results', 'T1ImgR2SN', SubjectID, 'RiskBrainMap.nii')
        T1riskCorrWithCogPath = os.path.join(WorkingDir, 'Results', 'T1ImgR2SN', SubjectID, 'T1_risk_corr_with_cognition.jpg')
        T1riskNetworkBelongings= os.path.join(WorkingDir, 'Results', 'T1ImgR2SN', SubjectID, 'T1_risk_map_to_Yeo7.jpg')
        atlas_GM_referencePath = os.path.join(templatePath, "BN246_volume_reference.csv")

        if ReportInIndividualSpace:
            originImagePath = os.path.join(WorkingDir, 'T1Img', SubjectID,  SubjectID + '.nii')
            gmPath = os.path.join(WorkingDir, 'T1Img', SubjectID, 'reg_mwp1' + SubjectID + '.nii')
            wmPath = os.path.join(WorkingDir, 'T1Img', SubjectID, 'reg_mwp2' + SubjectID + '.nii')
            hippoPath = os.path.join(WorkingDir, 'T1Img', SubjectID, 'reg_hippo_mask.nii')
            riskimgPath = os.path.join(WorkingDir, 'T1Img', SubjectID, 'reg_RiskBrainMap.nii')

        doc = DocxTemplate(os.path.join(templatePath, "Report_template_T1.docx"))



    #generate T1 report if required
    segmentReport = ""
    if (originImagePath is not None):
        SliceImagePath = os.path.join(savePath, 'slice.png')
        draw2DSlices(originImagePath, SliceImagePath, 5, 4)

    if(gmPath is not None) and (wmPath is not None):
        wmImage = nib.load(wmPath)
        gmImage = nib.load(gmPath)

        wmVolume = np.sum(wmImage.get_fdata().flatten())
        wmVolume = round(wmVolume/mm_to_cm_scale_rate, 2)
        gmVolume = np.sum(gmImage.get_fdata().flatten())
        gmVolume = round(gmVolume/mm_to_cm_scale_rate, 2)
        wmVolumeRef = str(round(reference["wmVolume"][0]/mm_to_cm_scale_rate,2)) + "-" + str(round(reference["wmVolume"][1]/mm_to_cm_scale_rate,2))
        gmVolumeRef = str(round(reference["gmVolume"][0]/mm_to_cm_scale_rate,2)) + "-" + str(round(reference["gmVolume"][1]/mm_to_cm_scale_rate,2))
        totalVolume = wmVolume + gmVolume
        totalVolume = round(totalVolume,2)

        wmVperData = float(int(wmVolume / totalVolume * 1000)) / 10
        gmVperData = float(int(gmVolume / totalVolume * 1000)) / 10
        wmVper = str(wmVperData) + "%"
        gmVper = str(gmVperData) + "%"
        wmVperRef = str(round(reference["wmPer"][0],2)) + "%-" + str(round(reference["wmPer"][1],2)) +"%"
        gmVperRef = str(round(reference["gmPer"][0],2)) + "%-" + str(round(reference["gmPer"][1],2)) + "%"

        if(wmVperData < reference["wmPer"][0]):
            segmentReport += "Subject whiter matter proportion reduced."
        elif(wmVperData > reference["wmPer"][1]):
            segmentReport += "Subject whiter matter proportion increased."
        if (gmVperData < reference["gmPer"][0]):
            segmentReport += "Subject grey matter proportion reduced."
        elif (gmVperData > reference["gmPer"][1]):
            segmentReport += "Subject grey matter proportion increased."
        if(segmentReport == ""):
            segmentReport += "Subject whiter matter and grey matter within normal range"

        SegmentImagePath = os.path.join(savePath, 'segment.png')
        drawSegment(originImagePath, gmPath, wmPath, SegmentImagePath)
        VolumeDistributionImagePath = os.path.join(savePath, "volume.png")
        drawDis(u=stdDistribution["NCgmMeanPer"],subject=gmVperData,xlabel="Grey matter proportion / % ",show_max="Subject grey matter proportion", sig=stdDistribution["NCgmStdPer"],u2=None,sig2=None,savePath=VolumeDistributionImagePath)
    else:
        segmentReport += "Missing segmentation result！"


    hippoReport = ""
    hippoDecreaseReport = ""
    hippoIncreaseReport = ""
    if(featurePath is not None):
        with open(featurePath,'r') as f:
            featureData = csv.DictReader(f)
            featureDict = [row for row in featureData][0]

        hippoLeft1Volume = round(float(featureDict['Hipp_L_2_1']) / mm_to_cm_scale_rate, 2)
        hippoRight1Volume = round(float(featureDict['Hipp_R_2_1']) / mm_to_cm_scale_rate, 2)
        hippoLeft2Volume = round(float(featureDict['Hipp_L_2_2']) / mm_to_cm_scale_rate, 2)
        hippoRight2Volume = round(float(featureDict['Hipp_R_2_2']) / mm_to_cm_scale_rate, 2)
        hippoVolume = round((hippoLeft1Volume + hippoRight1Volume + hippoLeft2Volume + hippoRight2Volume), 2)

        hippoLeft1VolumeRef = str(round(reference["hippoLeft1Volume"][0] / mm_to_cm_scale_rate, 2)) + '-' + str(
            round(reference["hippoLeft1Volume"][1] / mm_to_cm_scale_rate, 2))
        hippoRight1VolumeRef = str(round(reference["hippoRight1Volume"][0] / mm_to_cm_scale_rate, 2)) + '-' + str(
            round(reference["hippoRight1Volume"][1] / mm_to_cm_scale_rate, 2))
        hippoLeft2VolumeRef = str(round(reference["hippoLeft2Volume"][0] / mm_to_cm_scale_rate, 2)) + '-' + str(
            round(reference["hippoLeft2Volume"][1] / mm_to_cm_scale_rate, 2))
        hippoRight2VolumeRef = str(round(reference["hippoRight2Volume"][0] / mm_to_cm_scale_rate, 2)) + '-' + str(
            round(reference["hippoRight2Volume"][1] / mm_to_cm_scale_rate, 2))

        hippoLeft1VolumeMark = "-"
        hippoRight1VolumeMark = "-"
        hippoLeft2VolumeMark = "-"
        hippoRight2VolumeMark = "-"

        if (hippoLeft1Volume < reference["hippoLeft1Volume"][0] / mm_to_cm_scale_rate):
            hippoDecreaseReport += "Hipp_L_2_1、"
            hippoLeft1Volume = RichText(hippoLeft1Volume, color='FF0000', bold=True)
            hippoLeft1VolumeMark = RichText("↓", color='FF0000', bold=True)
        elif (hippoLeft1Volume > reference["hippoLeft1Volume"][1] / mm_to_cm_scale_rate):
            hippoIncreaseReport += "Hipp_L_2_1、"
            hippoLeft1Volume = RichText(hippoLeft1Volume, color='000000', bold=True)
            hippoLeft1VolumeMark = RichText("↑", color='000000', bold=True)

        if (hippoLeft2Volume < reference["hippoLeft2Volume"][0] / mm_to_cm_scale_rate):
            hippoDecreaseReport += "Hipp_L_2_2、"
            hippoLeft2Volume = RichText(hippoLeft2Volume, color='FF0000', bold=True)
            hippoLeft2VolumeMark = RichText("↓", color='FF0000', bold=True)
        elif (hippoLeft2Volume > reference["hippoLeft2Volume"][1] / mm_to_cm_scale_rate):
            hippoIncreaseReport += "Hipp_L_2_2、"
            hippoLeft2Volume = RichText(hippoLeft2Volume, color='000000', bold=True)
            hippoLeft2VolumeMark = RichText("↑", color='000000', bold=True)

        if (hippoRight1Volume < reference["hippoRight1Volume"][0] / mm_to_cm_scale_rate):
            hippoDecreaseReport += "Hipp_R_2_1、"
            hippoRight1Volume = RichText(hippoRight1Volume, color='FF0000', bold=True)
            hippoRight1VolumeMark = RichText("↓", color='FF0000', bold=True)
        elif (hippoRight1Volume > reference["hippoRight1Volume"][1] / mm_to_cm_scale_rate):
            hippoIncreaseReport += "Hipp_R_2_1、"
            hippoRight1Volume = RichText(hippoRight1Volume, color='000000', bold=True)
            hippoRight1VolumeMark = RichText("↑", color='000000', bold=True)

        if (hippoRight2Volume < reference["hippoRight2Volume"][0] / mm_to_cm_scale_rate):
            hippoDecreaseReport += "Hipp_R_2_2、"
            hippoRight2Volume = RichText(hippoRight2Volume, color='FF0000', bold=True)
            hippoRight2VolumeMark = RichText("↓", color='FF0000', bold=True)
        elif (hippoRight2Volume > reference["hippoRight2Volume"][1] / mm_to_cm_scale_rate):
            hippoIncreaseReport += "Hipp_R_2_2、"
            hippoRight2Volume = RichText(hippoRight2Volume, color='000000', bold=True)
            hippoRight2VolumeMark = RichText("↑", color='000000', bold=True)

        hippoIncreaseReport = combineBilateralRegionName(hippoIncreaseReport)
        hippoDecreaseReport = combineBilateralRegionName(hippoDecreaseReport)

        if "、" not in hippoDecreaseReport and  "、" not in hippoIncreaseReport:
            hippoReport = "Subject hippocampus volume within normal range."
        elif "、" not in hippoDecreaseReport and  "、" in hippoIncreaseReport:
            hippoReport = "Subject grey volume has been increased in: " + hippoIncreaseReport[0:-2]
        elif "、" in hippoDecreaseReport and  "、" not in hippoIncreaseReport:
            hippoReport = "Subject grey volume has been reduced in: " + hippoDecreaseReport[0:-2]
        else:
            hippoReport = "Subject grey volume has been reduced in: " + hippoDecreaseReport[0:-2] + "Whereas grey volume has been increased in: " + hippoIncreaseReport[0:-2]


        hippoImagePngPath = os.path.join(savePath,"hippoimage.png")
        drawRoi(originImagePath,hippoPath,hippoImagePngPath)
        hippoVolumeImagePath = os.path.join(savePath, "hippoVolume.png")
        drawDis(u=stdDistribution["NChippoVolumeMean"]/mm_to_cm_scale_rate, subject=hippoVolume, xlabel="hippocampus volume / mm3",show_max="subject hippo volume",
                sig=stdDistribution["NChippoVolumeStd"]/mm_to_cm_scale_rate, u2=None,
                sig2=None, savePath=hippoVolumeImagePath)

        ibrainScore = round(float(featureDict["ibrainScore"]) * 100, 2)
        brainScoreRefPath = os.path.join(modelPath, 'train_data_ibrain.mat')
        ref_ibrain = sio.loadmat(brainScoreRefPath)["nc_train_subject_ibrain_score"]*100
        iBrainDistributionPath = os.path.join(savePath, 'iBrainDistribution.png')
        drawHis(ref=ref_ibrain, subject=ibrainScore, savePath=iBrainDistributionPath,
                ref2=None, xlabel="iBrain value", show_max="subject iBrain value", )
        if ibrainScore > 50:
            ibrainScore = RichText(ibrainScore, color='00FF00')
        else:
            ibrainScore = RichText(ibrainScore, color='FF0000')

    else:
        hippoReport += "Missing hippocampus information."


    riskReport = ""
    sigdecreaseReport = ""
    riskimgsavePath = os.path.join(savePath, 'RiskBrainMap.png')
    risk_features = []
    if os.path.exists(riskimgPath) and atlas_GM_referencePath is not None:
        drawRoiSlice(riskimgPath, originImagePath, riskimgsavePath, 5, 4)

        with open(atlas_GM_referencePath, 'r' ,encoding='GB2312') as f:
            risk_refData = csv.DictReader(f)
            risk_refDicts = [row for row in risk_refData]

        for ref in risk_refDicts:
            featureName = ref['feature']
            if featureName in featureDict:  #add risk table generation, search from excel reference, Only add if exist
                featureValue = float(featureDict[featureName])/mm_to_cm_scale_rate
                refDataF = str(round(float(ref['NC-mean'])/mm_to_cm_scale_rate - float(ref['NC-std'])/mm_to_cm_scale_rate, 2)) + "~" + str(
                        round(float(ref['NC-mean'])/mm_to_cm_scale_rate + float(ref['NC-std'])/mm_to_cm_scale_rate, 2))

                if featureValue<round(float(ref['NC-mean'])/mm_to_cm_scale_rate - 2*float(ref['NC-std'])/mm_to_cm_scale_rate, 2):
                    remark = "↓"
                    temp_region_name = featureName.split("_")
                    all_decreased_name = sigdecreaseReport.split("、")
                    if temp_region_name[0] not in all_decreased_name:
                        sigdecreaseReport += temp_region_name[0]+'、'

                    temp_data = {}
                    temp_data["feature_name"] = RichText(featureName, color='FF0000', bold=False)
                    temp_data["feature_value"] = RichText(str(round(featureValue, 2)), color='FF0000', bold=False)
                    temp_data["reference_range"] = RichText(str(refDataF), color='000000', bold=False)
                    temp_data["mark"] = RichText(remark, color='FF0000', bold=True)
                    risk_features.append(temp_data)

        sigdecreaseReport = combineBilateralRegionName(sigdecreaseReport)
        riskReport = "Subject exhibited altered lower grey matter volume mainly in: "+sigdecreaseReport[0:-2]
    else:
        empty_folder_path = atlasPath
        emptyimgPath=os.path.join(empty_folder_path,'empty_image.nii')
        drawRisk(emptyimgPath, riskimgsavePath)
        T1riskCorrWithCogPath = os.path.join(savePath,'T1_risk_corr_with_cognition.png')
        drawRisk(emptyimgPath, T1riskCorrWithCogPath)
        T1riskNetworkBelongings = os.path.join(savePath,'T1_risk_map_to_Yeo7.png')
        drawRisk(emptyimgPath, T1riskNetworkBelongings)
        riskReport += "No obvious altered lower grey matter volume detected."


    GMReport = ""
    sigincreaseReport = ""
    GMdecreaseReport = ""
    GMincreaseReport = ""
    whole_brain_features = []

    if(featurePath is not None):
        Add_folder_path = os.path.dirname(os.path.abspath(featurePath))
        WholebrainfeaturePath = os.path.join(Add_folder_path, 'features_wholebrain.csv')
        with open(WholebrainfeaturePath,'r') as f:
            AddfeatureData = csv.DictReader(f)
            AddfeatureDict = [row for row in AddfeatureData][0]

        with open(atlas_GM_referencePath, 'r' ,encoding='GB2312') as f:
            GM_refData = csv.DictReader(f)
            GM_refDicts = [row for row in GM_refData]

        for ref in GM_refDicts:
            featureName = ref['feature']
            featureValue = float(AddfeatureDict[featureName]) / mm_to_cm_scale_rate
            refDataF = str(
                round(float(ref['NC-mean']) / mm_to_cm_scale_rate - float(ref['NC-std']) / mm_to_cm_scale_rate,
                      2)) + "~" + str(
                round(float(ref['NC-mean']) / mm_to_cm_scale_rate + float(ref['NC-std']) / mm_to_cm_scale_rate, 2))
            remark = "-"
            if featureValue < (
                    float(ref['NC-mean']) / mm_to_cm_scale_rate - float(ref['NC-std']) / mm_to_cm_scale_rate):
                remark = "↓"
                temp_region_name = featureName.split("_")
                all_decreased_name = sigdecreaseReport.split("、")
                all_GMdecreased_name = GMdecreaseReport.split("、")
                if temp_region_name[0] not in all_decreased_name and temp_region_name[0] not in all_GMdecreased_name:
                    GMdecreaseReport += temp_region_name[0] + '、'
            elif featureValue > (
                    float(ref['NC-mean']) / mm_to_cm_scale_rate + 2 * float(ref['NC-std']) / mm_to_cm_scale_rate):
                remark = "↑"
                temp_region_name = featureName.split("_")
                all_increased_name = sigincreaseReport.split("、")
                if temp_region_name[0] not in all_increased_name:
                    sigincreaseReport += temp_region_name[0] + '、'
            elif featureValue > (float(ref['NC-mean']) / mm_to_cm_scale_rate + float(
                    ref['NC-std']) / mm_to_cm_scale_rate) and featureValue <= (
                    float(ref['NC-mean']) / mm_to_cm_scale_rate + 2 * float(ref['NC-std']) / mm_to_cm_scale_rate):
                remark = "↑"
                temp_region_name = featureName.split("_")
                all_GMincreased_name = GMincreaseReport.split("、")
                if temp_region_name[0] not in all_GMincreased_name:
                    GMincreaseReport += temp_region_name[0] + '、'

            temp_data = {}
            if remark == "↓":
                temp_data["feature_name"] = RichText(featureName, color='FF0000', bold=False)
                temp_data["feature_value"] = RichText(str(round(featureValue, 2))+'*', color='FF0000', bold=False)
                temp_data["reference_range"] = RichText(str(refDataF), color='000000', bold=False)
                temp_data["mark"] = RichText(remark, color='FF0000', bold=True)
            elif remark == "↑":
                temp_data["feature_name"] = RichText(featureName, color='000000', bold=False)
                temp_data["feature_value"] = RichText(str(round(featureValue, 2))+'*', color='000000', bold=False)
                temp_data["reference_range"] = RichText(str(refDataF), color='000000', bold=False)
                temp_data["mark"] = RichText(remark, color='000000', bold=True)
            else:
                temp_data["feature_name"] = RichText(featureName, color='000000', bold=False)
                temp_data["feature_value"] = RichText(str(round(featureValue, 2)), color='000000', bold=False)
                temp_data["reference_range"] = RichText(str(refDataF), color='000000', bold=False)
                temp_data["mark"] = RichText(remark, color='000000', bold=True)
            whole_brain_features.append(temp_data)

        GMdecreaseReport = combineBilateralRegionName(GMdecreaseReport)
        sigincreaseReport = combineBilateralRegionName(sigincreaseReport)
        GMincreaseReport = combineBilateralRegionName(GMincreaseReport)


        if "、" in GMdecreaseReport:
            GMReport += "Subject contains potential lower grey matter volume mainly in: " + GMdecreaseReport[0:-2] #normalized T value larger than -2 but lower than mean-std
        if "、" in sigincreaseReport:
            GMReport += "Subject exhibited altered higher grey matter volume mainly in: " + sigincreaseReport[0:-2]
        if "、" in GMincreaseReport:
            GMReport += "Subject contains potential higher grey matter volume mainly in: " + GMincreaseReport[0:-2]


    finalReport = ''
    if len(segmentReport) != 0:
        finalReport += 'After segmentation, overall result of white matter and grey matter in subject brain are as follows：' + segmentReport + '\n'
    if len(hippoReport) != 0:
        finalReport += 'Summarized gray matter volume of hippocampus in subject brain are as follows：' + hippoReport + '\n'
    if len(riskReport) != 0:
        finalReport += 'Regions with altered lower gray matter volume in subject brain are as follows：' + riskReport

    if ReportWholeBrain:
        finalReport += '\n'
        if len(GMReport) != 0:
            finalReport += 'Summarized whole brain gray matter volume in subject brain are as follows：' + GMReport + '\n'


    context = {'SubjectName': os.path.split(savePath.rstrip(os.path.sep))[-1],
               'ReportDate': datetime.date.today(),
               "Age": '',
               "Sex": '',
               "ScanDate": '',
               "finalReport": finalReport}

    if writeMode == 'T1' or writeMode == 'Combined':
        T1_attributes = {"BrainOrigSliceImage": InlineImage(doc, SliceImagePath, width=Mm(130), height=Mm(200)),
                     "BrainSegmentImage": InlineImage(doc, SegmentImagePath, width=Mm(130)),
                     "wmVolume": wmVolume,
                     "gmVolume": gmVolume,
                     "wmVolumeRef": wmVolumeRef,
                     "gmVolumeRef": gmVolumeRef,
                     "totalVolume": totalVolume,
                     "wmVper": wmVper,
                     "gmVper": gmVper,
                     "wmVperRef": wmVperRef,
                     "gmVperRef": gmVperRef,
                     "segmentReport": segmentReport,
                     "VolumeDistributionImage": InlineImage(doc, VolumeDistributionImagePath, width=Mm(100)),
                     "SubjectHippoImage": InlineImage(doc, hippoImagePngPath, width=Mm(130)),
                     "left1Volume": hippoLeft1Volume,
                     "right1Volume": hippoRight1Volume,
                     "left2Volume": hippoLeft2Volume,
                     "right2Volume": hippoRight2Volume,
                     "left1VolumeRef": hippoLeft1VolumeRef,
                     "right1VolumeRef": hippoRight1VolumeRef,
                     "left2VolumeRef": hippoLeft2VolumeRef,
                     "right2VolumeRef": hippoRight2VolumeRef,
                     "left1VolumeMark": hippoLeft1VolumeMark,
                     "right1VolumeMark": hippoRight1VolumeMark,
                     "left2VolumeMark": hippoLeft2VolumeMark,
                     "right2VolumeMark": hippoRight2VolumeMark,
                     "hippoReport": hippoReport,
                     "HippoDistributionImage": InlineImage(doc, hippoVolumeImagePath, width=Mm(100)),
                     "iBrain": ibrainScore,
                     "iBrainDistributionImage": InlineImage(doc, iBrainDistributionPath, width=Mm(100)),
                     "BrainRiskImage": InlineImage(doc, riskimgsavePath, width=Mm(130)),
                     "T1BrainMultiDimCognition": InlineImage(doc, T1riskCorrWithCogPath, width=Mm(100)),
                     "BrainRiskNetworkBelongings": InlineImage(doc, T1riskNetworkBelongings, width=Mm(100)),
                     "risk_features": risk_features,
                     "riskReport": riskReport
                     }
        if ReportWholeBrain:
            T1_addattributes={"T1WholeBrainReportTitle":'Whole brain grey matter volume report',
                              "whole_brain_features": whole_brain_features,
                              "GMReport": GMReport}
        else:
            T1_addattributes = {}


    if writeMode == 'T1' or writeMode == 'Combined':
        for key, value in T1_attributes.items():
            context[key] = value
        for key, value in T1_addattributes.items():
            context[key] = value

    doc.render(context)
    saveFileName = os.path.realpath(os.path.join(savePath,SubjectID+"_report.docx"))
    if os.path.exists(saveFileName):
        os.remove(saveFileName)
    doc.save(saveFileName)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--iBrainPath', type=str, default=None)
    parser.add_argument('--WorkingDir', type=str, default=None)
    parser.add_argument('--SubjectID', type=str, default=None)
    parser.add_argument('--ReportWholeBrain', type=str, default=None)
    parser.add_argument('--ReportInIndividualSpace', type=str, default=None)
    parser.add_argument('--MiniOutput', type=str, default=None)
    parser.add_argument('--savePath', type=str, default=None)
    parser.add_argument('--writeMode', type=str, default='T1')  # T1 stands for only write T1 report, BOLD stands for only write BOLD report, whereas Combined stands for write T1 and BOLD report
    parser.add_argument('--tempFCsession', type=str, default='1')
    args = parser.parse_args()
    reportBuild(args.iBrainPath, args.WorkingDir, args.SubjectID, args.ReportWholeBrain, args.ReportInIndividualSpace, args.MiniOutput, args.savePath, args.writeMode, args.tempFCsession)

    # Example Inputs:
    # iBrainPath = r"E:\MATLABtoolboxes\iBrain"
    # WorkingDir = r"E:\MATLABtoolboxes\iBrain_test_data\T1AndBOLD"
    # SubjectID = r"An_Su_Fang"
    # ReportWholeBrain = "1"
    # ReportInIndividualSpace = "1"
    # MiniOutput = "0"
    # savePath = r"E:\MATLABtoolboxes\iBrain_test_data\T1AndBOLD\\Results\Report\An_Su_Fang"
    # writeMode = r"Combined"
    # tempFCsession = "1"
    # reportBuild(iBrainPath, WorkingDir, SubjectID, ReportWholeBrain, ReportInIndividualSpace, MiniOutput, savePath, writeMode, tempFCsession)