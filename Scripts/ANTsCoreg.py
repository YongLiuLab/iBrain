import ants
import argparse
import os

def coregister_Atlas_to_indispace(ref_img,input_img,input_atlas_mask,savePath):
    ref_img_file = ants.image_read(ref_img)
    input_img_file = ants.image_read(input_img)
    input_atlas_mask_file = ants.image_read(input_atlas_mask)
    transform = ants.registration(ref_img_file, input_img_file, type_of_transform='SyN', metric='MI')
    if input_atlas_mask is not None:
        registered = ants.apply_transforms(ref_img_file, input_atlas_mask_file, transform['fwdtransforms'])
        saveName = os.path.basename(input_atlas_mask)
    else:
        registered = ants.apply_transforms(ref_img_file, input_img_file, transform['fwdtransforms'])
        saveName = os.path.basename(input_img)
    ants.image_write(registered, os.path.join(savePath, "reg_"+saveName))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--ref_img', type=str, default=None)
    parser.add_argument('--input_img', type=str, default=None)
    parser.add_argument('--input_atlas_mask', type=str, default=None)
    parser.add_argument('--savePath', type=str, default=None)
    args = parser.parse_args()
    coregister_Atlas_to_indispace(args.ref_img, args.input_img, args.input_atlas_mask, args.savePath)

    # ref_img = r"E:\MATLABtoolboxes\iBrain_test_data\T1AndBOLD\T1Img\An_Su_Fang\An_Su_Fang.nii"
    # input_img = r"E:\MATLABtoolboxes\iBrain\Template\cat12_MNI152_T1_1mm.nii"
    # input_atlas_mask = r"E:\MATLABtoolboxes\iBrain_test_data\T1AndBOLD\T1ImgS\An_Su_Fang\mri\mwp1An_Su_Fang.nii"
    # #input_atlas_mask = r"E:\MATLABtoolboxes\iBrain\Atlas\Whole512_1mm.nii"
    # savePath=r"E:\MATLABtoolboxes\iBrain_test_data\T1AndBOLD\T1Img\An_Su_Fang"
    # coregister_Atlas_to_indispace(ref_img, input_img, input_atlas_mask, savePath)
