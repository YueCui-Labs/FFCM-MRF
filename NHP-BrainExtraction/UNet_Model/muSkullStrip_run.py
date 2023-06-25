#!/usr/bin/env python
import torch
import torch.nn as nn
from function import predict_volumes
from model import UNet2d
import os, sys
import argparse
import nibabel as nib

def NHP_Skullstripping_run(input_path, model, out_dir, is_save_mask):
    NoneType=type(None)
    # Argument
    # parser=argparse.ArgumentParser(description='Skullstripping', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # optional=parser._action_groups.pop()
    # required=parser.add_argument_group('required arguments')
    # # Required Option
    # required.add_argument('-in', '--input_t1w', type=str, required=True, help='Input T1w Image for Skull Stripping')
    # required.add_argument('-model', '--predict_model', required=True, type=str, help='Predict Model')
    # # Optional Option
    # optional.add_argument('-out', '--out_dir', type=str, help='Output Dir')
    # optional.add_argument('-suffix', '--mask_suffix', type=str, default="pre_mask", help='Suffix of Mask')
    # optional.add_argument('-slice', '--input_slice', type=int, default=3, help='Number of Slice for Model Input')
    # optional.add_argument('-conv', '--conv_block', type=int, default=5, help='Number of UNet Block')
    # optional.add_argument('-kernel', '--kernel_root', type=int, default=16, help='Number of the Root of Kernel')
    # optional.add_argument('-rescale', '--rescale_dim', type=int, default=256, help='Number of the Root of Kernel')
    # optional.add_argument('-ed_iter', '--erosion_dilation_iteration', type=int, default=0, help='Number of Iteration for Erosion and Dilation')
    # parser._action_groups.append(optional)
    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(1)
    # args = parser.parse_args()

    print("=============================Permforming Skullstripping=============================")
    
    train_model=UNet2d(dim_in=3, num_conv_block=5, kernel_root=16)
    checkpoint=torch.load(model, map_location={'cuda:0':'cpu'})
    train_model.load_state_dict(checkpoint['state_dict'])
    model=nn.Sequential(train_model, nn.Softmax2d())

    predict_volumes(model, cimg_in=input_path, bmsk_in=None, rescale_dim=256, save_dice=False,
            save_nii=True, nii_outdir=out_dir, suffix="brain_mask", ed_iter=0)
    brain_mask_file=os.path.join(out_dir,'MRA_brain_mask.nii.gz')
    brain_mask_=nib.load(brain_mask_file)
    brain_mask=brain_mask_.get_fdata()
    input_image=nib.load(input_path).get_fdata()
    brain=input_image*brain_mask
    save_brain_path = os.path.join(out_dir,'MRA_brain.nii.gz') # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(brain, brain_mask_.affine,
                             brain_mask_.header), save_brain_path)
    if not is_save_mask:
        os.remove(brain_mask_file)