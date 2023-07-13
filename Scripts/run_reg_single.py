import ants
import os
import argparse

def reg_run(filename,template_path,output_dir):
    image = ants.image_read(filename)
    file_path, ori_name = os.path.split(filename)
    imagedenoise = ants.denoise_image(image, ants.get_mask(image))
    image_n4 = ants.n4_bias_field_correction(imagedenoise)
    template = ants.image_read(template_path)
    out = ants.registration(template, image_n4, type_of_transform='SyN')
    reg_img = out['warpedmovout']
    output_name='DNR_'+ori_name
    ants.image_write(reg_img, os.path.join(output_dir,output_name))

if __name__== '__main__':
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--input_file', type=str, default=None)
    parser.add_argument('--template_path', type=str, default=None)
    parser.add_argument('--output_dir', type=str, default=None)
    args = parser.parse_args()
    reg_run(args.input_file,args.template_path,args.output_dir)

