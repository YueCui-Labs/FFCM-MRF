import numpy as np
import SimpleITK as sitk
import os
import sys
import glob
import nibabel as nib
from tqdm import tqdm
import multiprocessing

from joblib import Parallel, delayed

def Resample(volume_path, is_label=True):

    itk_image = sitk.ReadImage(volume_path)
    # brain_mask_image = sitk.ReadImage(brain_mask_path, sitk.sitkFloat32)
    # itk_image=itk_image*brain_mask_image

    original_spacing = itk_image.GetSpacing()
    original_size = itk_image.GetSize()
    out_spacing = [np.min(original_spacing)]*3
    # out_spacing = [np.mean(original_spacing)] * 3
    # print(np.min(original_spacing))
    # print(out_spacing)
    out_size = [
        int(np.round(original_size[0] * (original_spacing[0] / out_spacing[0]))),
        int(np.round(original_size[1] * (original_spacing[1] / out_spacing[1]))),
        int(np.round(original_size[2] * (original_spacing[2] / out_spacing[2])))
    ]


    # 上述也可以直接用下面这句简写
    # out_size = [int(round(osz*ospc/nspc)) for osz, ospc, nspc in zip(original_size, original_spacing, out_spacing)]
    # print(out_size)

    resample = sitk.ResampleImageFilter()
    resample.SetOutputSpacing(out_spacing)
    resample.SetSize(out_size)
    resample.SetOutputDirection(itk_image.GetDirection())
    resample.SetOutputOrigin(itk_image.GetOrigin())
    resample.SetTransform(sitk.Transform())
    resample.SetDefaultPixelValue(itk_image.GetPixelIDValue())

    if is_label:  # 如果是mask图像，就选择sitkNearestNeighbor这种插值
        resample.SetInterpolator(sitk.sitkNearestNeighbor)
    else:  # 如果是普通图像，就采用sitkBSpline插值法
        resample.SetInterpolator(sitk.sitkBSpline)

    return resample.Execute(itk_image)


if __name__ == "__main__":

    data_path='C:/Users/86133/Documents/pre-graduate-tasks/M097/M097'
    for ds in os.listdir(data_path):
        # subs = sorted([sub for sub in os.listdir(os.path.join(data_path, ds))
        #                if os.path.isdir(os.path.join(data_path, ds, sub))])
        # for sub in subs:
        # res_path = os.path.join(data_path, 'resample')
        # # print(res_path)
        # if not os.path.exists(res_path):
        #     os.makedirs(res_path)
        # if not os.path.isdir(os.path.join(data_path, ds)):

            # ds='MIDAS';
            # sub='M001'
            # brain_path = 'C:/Users/86133/Documents/pre-graduate-tasks/TOF-MRA/T1_to_TOF/FieldCorrection/Fin-N4Cor001-MRA.nii.gz'
        seg_path = os.path.join(data_path, 'MRA_brain.nii.gz')
        data = nib.load(seg_path)
        affine, header = data.affine, data.header

        #MSC专用
        # brain_mask_path=os.path.join(data_path, ds, sub, 'diameters', 'brain_mask_resample.nii.gz')
        # brain_mask_image = sitk.ReadImage(brain_mask_path, sitk.sitkFloat32)
        # brain = sitk.GetArrayFromImage(brain_mask_image)

        mask_seg_resample = Resample(seg_path)
        # original_origin = mask_seg_resample.GetOrigin()  # 获取图像的原点坐标
        # original_spacing = mask_seg_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        # original_direction = mask_seg_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵
        #
        # mask_seg = sitk.GetArrayFromImage(mask_seg_resample)
        # mask_seg_ = mask_seg * brain  # 消除突变的slice
        #
        # # 保存原本的属性值到新图片
        # mask_seg_resample = sitk.GetImageFromArray(mask_seg_)
        # mask_seg_resample.SetSpacing(original_spacing)
        # # mask_seg_resample.SetSize(original_size)
        # mask_seg_resample.SetDirection(original_direction)
        # mask_seg_resample.SetOrigin(original_origin)

        # sitk.WriteImage(mask_seg_resample, mask_resample_path)


        # seg_resampled[seg_resampled>1]=0
        out_name = os.path.join(data_path, 'brain_resample.nii.gz')
        # nib.Nifti1Image(seg_resampled,affine).to_filename(out_name)
        sitk.WriteImage(mask_seg_resample,out_name)
        print(' finished.')
    
