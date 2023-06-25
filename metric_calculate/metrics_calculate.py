# import argparse
import math
import SimpleITK as sitk
import nibabel as nib
import numpy as np
from scipy import ndimage
import numpy as np

import os
import sys
import pandas as pd
import nibabel as nib
from collections import OrderedDict
import multiprocessing

from resample_img import Resample
from ExtractCenterline import *
from ExtractDiameter import *


def whole_brain_feature_extraction(seg_path,brain_mask_path):
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-d", "--dataset", type=str)
    # parser.add_argument("-s","--sub",type=int)
    #
    # args = parser.parse_args()

    seg_=os.path.split(seg_path)
    # print(seg_)
    # print(prefix)
    res_path = os.path.join(seg_[0], 'stats')
    print(res_path)
    if not os.path.exists(res_path):
        os.makedirs(res_path)



    resample_path = os.path.join(res_path, 'vessel_resample.nii.gz')

    save_brain_path=os.path.join(res_path, 'brain_mask_resample.nii.gz')

    brain_resample_path=os.path.join(res_path, 'brain_resample.nii.gz')

    if not os.path.exists(save_brain_path):
    # if True:# 重制brain_mask
        brain_mask_image=Resample(brain_mask_path)
        # 读取原本的属性值
        # original_size=brain_mask_image.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin=brain_mask_image.GetOrigin()  # 获取图像的原点坐标
        original_spacing=brain_mask_image.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction=brain_mask_image.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        sitk.WriteImage(brain_mask_image, save_brain_path)# 把brain_mask做重采样，先保存
        brain_mask_image = sitk.ReadImage(save_brain_path, sitk.sitkFloat32)

        brain = sitk.GetArrayFromImage(brain_mask_image)
        brain[brain>1]=0 # 通过阈值消除多余的slice
        # 保存原本的属性值到新图片
        brain_mask_image = sitk.GetImageFromArray(brain)
        brain_mask_image.SetSpacing(original_spacing)
        # brain_mask_image.SetSize(original_size)
        brain_mask_image.SetDirection(original_direction)
        brain_mask_image.SetOrigin(original_origin)

        sitk.WriteImage(brain_mask_image, save_brain_path)# 保存真正的重采样brain_mask
    # 用重制brain_mask去mask重采样的数据
    brain_mask_image = sitk.ReadImage(save_brain_path, sitk.sitkFloat32)
    brain = sitk.GetArrayFromImage(brain_mask_image)
    if not os.path.exists(resample_path):
    # if True:
        seg_resample = Resample(seg_path)
        # 读取原本的属性值
        # original_size = seg_resample.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = seg_resample.GetOrigin()  # 获取图像的原点坐标
        original_spacing = seg_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = seg_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        seg = sitk.GetArrayFromImage(seg_resample)
        seg_ = seg * brain  # 消除突变的slice

        # 保存原本的属性值到新图片
        seg_resample = sitk.GetImageFromArray(seg_)
        seg_resample.SetSpacing(original_spacing)
        # seg_resample.SetSize(original_size)
        seg_resample.SetDirection(original_direction)
        seg_resample.SetOrigin(original_origin)

        sitk.WriteImage(seg_resample, resample_path)

    if not os.path.exists(brain_resample_path):# 大脑图像重采样
    # if True:
        brain_path=brain_mask_path[:-12]+'.nii.gz'
        print(brain_path)
        seg_resample = Resample(brain_path)
        # 读取原本的属性值
        # original_size = seg_resample.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = seg_resample.GetOrigin()  # 获取图像的原点坐标
        original_spacing = seg_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = seg_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        seg = sitk.GetArrayFromImage(seg_resample)
        seg_ = seg * brain  # 消除突变的slice

        # 保存原本的属性值到新图片
        seg_resample = sitk.GetImageFromArray(seg_)
        seg_resample.SetSpacing(original_spacing)
        # seg_resample.SetSize(original_size)
        seg_resample.SetDirection(original_direction)
        seg_resample.SetOrigin(original_origin)

        sitk.WriteImage(seg_resample, brain_resample_path)




    # Second compute centerline and diameter files of resampled image
    tof_file = os.path.join(res_path, 'vessel_resample.nii.gz')
    tof = nib.load(tof_file)
    affine: np.ndarray = tof.affine
    header: nib.Nifti1Header = tof.header

    cl_file = os.path.join(res_path, 'centerline.nii.gz')
    diameter_file = os.path.join(res_path, 'diameter.nii.gz')

    if not os.path.exists(cl_file):
    # if True:
        extract_cl = ExtractCenterline(tof_file, cl_file)
        extract_cl.run()
    if not os.path.exists(diameter_file):
    # if True:
        extract_dia = ExtractDiameter(tof_file, diameter_file, 0)
        extract_dia.run()

    # Finally compute diameters related stat
    cl = nib.load(cl_file).get_fdata()
    dia = nib.load(diameter_file).get_fdata()
    cl_diameter = cl*dia

    # 总图的summary
    stat = OrderedDict(
        (
            ("Volume", np.nan),
            ("Volume normalization", np.nan),
            ('Vessel length density', np.nan),
            ('Vessel length', np.nan),
            ("Diameter mean", np.nan),
            ("Diameter median", np.nan),
            ("Diameter max", np.nan),
            ("Diameter min", np.nan),
            ("Diameter STD", np.nan),
            ("Diameter N", np.nan),

        )
    )
    mean_zoom_cube = np.mean(header.get_zooms()) ** 3
    seg_diameter: np.ndarray = cl_diameter[cl_diameter != 0].flatten()
    stat["Volume"] = np.count_nonzero(tof.get_fdata()) * mean_zoom_cube
    stat["Volume normalization"] = np.count_nonzero(tof.get_fdata()) / np.count_nonzero(brain)
    stat['Vessel length density'] = np.count_nonzero(cl) * np.mean(header.get_zooms()) / (
            np.count_nonzero(brain) * mean_zoom_cube)
    stat['Vessel length'] = np.count_nonzero(cl) * np.mean(header.get_zooms())
    stat["Diameter mean"] = seg_diameter.mean()
    stat["Diameter median"] = np.median(seg_diameter)
    stat["Diameter max"] = seg_diameter.max()
    stat["Diameter min"] = seg_diameter.min()
    stat["Diameter STD"] = seg_diameter.std()
    stat["Diameter N"] = seg_diameter.size
    # pd.DataFrame.from_dict(stat, orient='index').to_csv(os.path.join(res_path, "SUMMARY.csv"))
    pd.DataFrame.from_dict(stat, orient='index').to_csv(os.path.join(res_path, "feature_summary.csv"))






def atlas_region_features_extraction(seg_path, brain_mask_path, atlas_path,atlas_csv_path,prefix):
    seg_ = os.path.split(seg_path)
    res_path = os.path.join(seg_[0], prefix+'_features')
    print(res_path)
    if not os.path.exists(res_path):
        os.makedirs(res_path)

    resample_path = os.path.join(res_path, 'vessel_resample.nii.gz')

    save_brain_path = os.path.join(res_path, 'brain_mask_resample.nii.gz')

    if not os.path.exists(save_brain_path):
        # if True:# 重制brain_mask
        brain_mask_image = Resample(brain_mask_path)
        # 读取原本的属性值
        # original_size=brain_mask_image.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = brain_mask_image.GetOrigin()  # 获取图像的原点坐标
        original_spacing = brain_mask_image.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = brain_mask_image.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        sitk.WriteImage(brain_mask_image, save_brain_path)  # 把brain_mask做重采样，先保存
        brain_mask_image = sitk.ReadImage(save_brain_path, sitk.sitkFloat32)

        brain = sitk.GetArrayFromImage(brain_mask_image)
        brain[brain > 1] = 0  # 通过阈值消除多余的slice
        # 保存原本的属性值到新图片
        brain_mask_image = sitk.GetImageFromArray(brain)
        brain_mask_image.SetSpacing(original_spacing)
        # brain_mask_image.SetSize(original_size)
        brain_mask_image.SetDirection(original_direction)
        brain_mask_image.SetOrigin(original_origin)

        sitk.WriteImage(brain_mask_image, save_brain_path)  # 保存真正的重采样brain_mask
    # 用重制brain_mask去mask重采样的数据
    brain_mask_image = sitk.ReadImage(save_brain_path, sitk.sitkFloat32)
    brain = sitk.GetArrayFromImage(brain_mask_image)
    if not os.path.exists(resample_path):
        # if True:
        seg_resample = Resample(seg_path)
        # 读取原本的属性值
        # original_size = seg_resample.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = seg_resample.GetOrigin()  # 获取图像的原点坐标
        original_spacing = seg_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = seg_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        seg = sitk.GetArrayFromImage(seg_resample)
        seg_ = seg * brain  # 消除突变的slice

        # 保存原本的属性值到新图片
        seg_resample = sitk.GetImageFromArray(seg_)
        seg_resample.SetSpacing(original_spacing)
        # seg_resample.SetSize(original_size)
        seg_resample.SetDirection(original_direction)
        seg_resample.SetOrigin(original_origin)

        sitk.WriteImage(seg_resample, resample_path)

    # Second compute centerline and diameter files of resampled image
    tof_file = os.path.join(res_path, 'vessel_resample.nii.gz')
    tof = nib.load(tof_file)
    tof_data=tof.get_fdata()
    affine: np.ndarray = tof.affine
    header: nib.Nifti1Header = tof.header

    cl_file = os.path.join(res_path, 'centerline.nii.gz')
    diameter_file = os.path.join(res_path, 'diameter.nii.gz')

    if not os.path.exists(cl_file):
        # if True:
        extract_cl = ExtractCenterline(tof_file, cl_file)
        extract_cl.run()
    if not os.path.exists(diameter_file):
        # if True:
        extract_dia = ExtractDiameter(tof_file, diameter_file, 0)
        extract_dia.run()

    # Finally compute diameters related stat
    cl = nib.load(cl_file).get_fdata()
    dia = nib.load(diameter_file).get_fdata()
    cl_diameter = cl * dia


    atlas_param=pd.read_csv(atlas_csv_path, names=['min','max'],index_col=0) #index为对应脑区名，min为最小阈值（不包括），max为最大阈值（包括）

    mask = nib.load(atlas_path)
    print(1)
    atlas_data = mask.get_fdata()
    # #Gyrus分区
    Gyrus_list = []
    for index in atlas_param.index:
        gyrus_ = np.zeros(atlas_data.shape)
        gyrus_[np.where((atlas_data > atlas_param.loc[index, 'min']) & (atlas_data <= atlas_param.loc[index, 'max']))] = 1
        Gyrus_list.append(gyrus_)

    length_densitys = np.zeros(tof_data.shape)
    volumes_normal = np.zeros(tof_data.shape)
    diameters = np.zeros(tof_data.shape)
    volumes=np.zeros(tof_data.shape)
    lengths=np.zeros(tof_data.shape)


    # #Lobe level stat
    print("start")
    # #没有文件夹要先建文件夹！
    gyrus_path=os.path.join(res_path, 'Atlas_features')
    if not os.path.exists(gyrus_path):
        os.makedirs(gyrus_path)

    for index, gyrus_ in zip(atlas_param.index, Gyrus_list):
        stat = OrderedDict(
            (
                ("Volume", np.nan),
                ("Volume density", np.nan),
                ('Vessel length', np.nan),
                ('Vessel length density', np.nan),
                ("Diameter mean", np.nan),
                ("Diameter median", np.nan),
                ("Diameter max", np.nan),
                ("Diameter min", np.nan),
                ("Diameter STD", np.nan),
                ("Diameter N", np.nan)
            )
        )
        mean_zoom_cube = np.mean(header.get_zooms()) ** 3  # 单个体素分辨率

        cl_gyrus_dia=cl_diameter*gyrus_ # 中心线*脑区mask
        seg_diameter: np.ndarray = cl_gyrus_dia[cl_gyrus_dia != 0].flatten()
        pict=tof.get_fdata()*gyrus_# 总血管*脑区mask
        # 体积
        stat["Volume"] = np.count_nonzero(pict) * mean_zoom_cube
        stat["Volume density"] = np.count_nonzero(pict) /np.sum(gyrus_ == 1)
    #血管长度密度
        stat['Vessel length density'] = np.count_nonzero(cl * gyrus_) * np.mean(header.get_zooms()) / (
                    np.sum(gyrus_ == 1) * mean_zoom_cube)
        stat['Vessel length'] =np.count_nonzero(cl * gyrus_) * np.mean(header.get_zooms())
    # 血管直径
        stat["Diameter N"] = seg_diameter.size
        if seg_diameter.size == 0:
            print(0)
            stat["Diameter max"] = 0
            stat["Diameter min"] = 0
            stat["Diameter mean"] = 0
            stat["Diameter median"] = 0
            stat["Diameter STD"] = 0
        else:
            stat["Diameter mean"] = seg_diameter.mean()
            stat["Diameter max"] = seg_diameter.max()
            stat["Diameter min"] = seg_diameter.min()
            stat["Diameter median"] = np.median(seg_diameter)
            stat["Diameter STD"] = seg_diameter.std()

        length_densitys[np.where(gyrus_ == 1)] = stat['Vessel length density']
        volumes_normal[np.where(gyrus_ == 1)] = stat["Volume density"]
        lengths[np.where(gyrus_ == 1)] = stat['Vessel length']
        volumes[np.where(gyrus_ == 1)] = stat["Volume"]
        diameters[np.where(gyrus_ == 1)] = stat["Diameter mean"]


        csv_name=index+"_summary.csv" #脑网络组图谱
        pd.DataFrame.from_dict(stat, orient='index').to_csv(os.path.join(gyrus_path, csv_name))


        print(f" {index} finished.")

    save_path = os.path.join(res_path, 'length_density.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(length_densitys, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'volume_density.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(volumes_normal, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'length.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(lengths, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'volume.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(volumes, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'diameter_mean.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(diameters, tof.affine, tof.header), save_path)


def atlas_region_features_extraction_no_csv(seg_path, brain_mask_path, atlas_path):
    seg_ = os.path.split(seg_path)
    res_path = os.path.join(seg_[0], 'stats')
    print(res_path)
    if not os.path.exists(res_path):
        os.makedirs(res_path)

    resample_path = os.path.join(res_path, 'vessel_resample.nii.gz')

    save_brain_path = os.path.join(res_path, 'brain_mask_resample.nii.gz')

    brain_resample_path = os.path.join(res_path, 'brain_resample.nii.gz')

    atlas_resample_path=os.path.join(res_path, 'Atlas_resample.nii.gz')

    if not os.path.exists(save_brain_path):
        # if True:# 重制brain_mask
        brain_mask_image = Resample(brain_mask_path)
        # 读取原本的属性值
        # original_size=brain_mask_image.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = brain_mask_image.GetOrigin()  # 获取图像的原点坐标
        original_spacing = brain_mask_image.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = brain_mask_image.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        sitk.WriteImage(brain_mask_image, save_brain_path)  # 把brain_mask做重采样，先保存
        brain_mask_image = sitk.ReadImage(save_brain_path, sitk.sitkFloat32)

        brain = sitk.GetArrayFromImage(brain_mask_image)
        brain[brain > 1] = 0  # 通过阈值消除多余的slice
        # 保存原本的属性值到新图片
        brain_mask_image = sitk.GetImageFromArray(brain)
        brain_mask_image.SetSpacing(original_spacing)
        # brain_mask_image.SetSize(original_size)
        brain_mask_image.SetDirection(original_direction)
        brain_mask_image.SetOrigin(original_origin)

        sitk.WriteImage(brain_mask_image, save_brain_path)  # 保存真正的重采样brain_mask
    # 用重制brain_mask去mask重采样的数据
    brain_mask_image = sitk.ReadImage(save_brain_path, sitk.sitkFloat32)
    brain = sitk.GetArrayFromImage(brain_mask_image)
    if not os.path.exists(resample_path):
        # if True:
        seg_resample = Resample(seg_path)
        # 读取原本的属性值
        # original_size = seg_resample.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = seg_resample.GetOrigin()  # 获取图像的原点坐标
        original_spacing = seg_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = seg_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        seg = sitk.GetArrayFromImage(seg_resample)
        seg_ = seg * brain  # 消除突变的slice

        # 保存原本的属性值到新图片
        seg_resample = sitk.GetImageFromArray(seg_)
        seg_resample.SetSpacing(original_spacing)
        # seg_resample.SetSize(original_size)
        seg_resample.SetDirection(original_direction)
        seg_resample.SetOrigin(original_origin)

        sitk.WriteImage(seg_resample, resample_path)
    #图谱重采样
    if not os.path.exists(atlas_resample_path):
        # if True:
        atlas_resample = Resample(atlas_path)
        # 读取原本的属性值
        # original_size = mask_seg_resample.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = atlas_resample.GetOrigin()  # 获取图像的原点坐标
        original_spacing = atlas_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = atlas_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        mask_seg = sitk.GetArrayFromImage(atlas_resample)
        mask_seg_ = mask_seg * brain  # 消除突变的slice

        # 保存原本的属性值到新图片
        atlas_resample = sitk.GetImageFromArray(mask_seg_)
        atlas_resample.SetSpacing(original_spacing)
        # mask_seg_resample.SetSize(original_size)
        atlas_resample.SetDirection(original_direction)
        atlas_resample.SetOrigin(original_origin)

        sitk.WriteImage(atlas_resample, atlas_resample_path)
    if not os.path.exists(brain_resample_path):# 大脑图像重采样
    # if True:
        brain_path=brain_mask_path[:-12]+'.nii.gz'
        print(brain_path)
        seg_resample = Resample(brain_path)
        # 读取原本的属性值
        # original_size = seg_resample.GetSize()  # 获取图像的大小，size为图像的每一个维度的长度，即每个维度像素点的个数
        original_origin = seg_resample.GetOrigin()  # 获取图像的原点坐标
        original_spacing = seg_resample.GetSpacing()  # 获取每个维度上像素或体素之间的间距，单位mm，其中对于二维图像是像素，三维图像是体素
        original_direction = seg_resample.GetDirection()  # 获取图像的方向，即图像坐标系相对世界坐标系的角度，角度采用的是方向余弦矩阵

        seg = sitk.GetArrayFromImage(seg_resample)
        seg_ = seg * brain  # 消除突变的slice

        # 保存原本的属性值到新图片
        seg_resample = sitk.GetImageFromArray(seg_)
        seg_resample.SetSpacing(original_spacing)
        # seg_resample.SetSize(original_size)
        seg_resample.SetDirection(original_direction)
        seg_resample.SetOrigin(original_origin)

        sitk.WriteImage(seg_resample, brain_resample_path)

    # Second compute centerline and diameter files of resampled image
    tof_file = os.path.join(res_path, 'vessel_resample.nii.gz')
    tof = nib.load(tof_file)
    tof_data=tof.get_fdata()
    affine: np.ndarray = tof.affine
    header: nib.Nifti1Header = tof.header

    cl_file = os.path.join(res_path, 'centerline.nii.gz')
    diameter_file = os.path.join(res_path, 'diameter.nii.gz')

    if not os.path.exists(cl_file):
        # if True:
        extract_cl = ExtractCenterline(tof_file, cl_file)
        extract_cl.run()
    if not os.path.exists(diameter_file):
        # if True:
        extract_dia = ExtractDiameter(tof_file, diameter_file, 0)
        extract_dia.run()

    # Finally compute diameters related stat
    cl = nib.load(cl_file).get_fdata()
    dia_file = nib.load(diameter_file)
    dia=dia_file.get_fdata()
    cl_diameter = cl * dia




    mask = nib.load(atlas_resample_path)
    print(1)
    atlas_data = mask.get_fdata()
    atlas_param=np.unique(atlas_data)
    # atlas_param.remove(0)
    atlas_param=np.delete(atlas_param, 0)

    # #Gyrus分区
    Gyrus_list = []
    for label in atlas_param:
        if label==0:
            continue
        gyrus_ = np.zeros(atlas_data.shape)
        gyrus_[np.where(atlas_data==label)] = 1
        Gyrus_list.append(gyrus_)

    length_densitys = np.zeros(tof_data.shape)
    volumes_normal = np.zeros(tof_data.shape)
    diameters = np.zeros(tof_data.shape)
    volumes=np.zeros(tof_data.shape)
    lengths=np.zeros(tof_data.shape)


    # #Lobe level stat
    # print("start")
    # # #没有文件夹要先建文件夹！
    # gyrus_path=os.path.join(res_path, 'Atlas_features')
    # if not os.path.exists(gyrus_path):
    #     os.makedirs(gyrus_path)

    SUMMARY=pd.DataFrame()

    for index, gyrus_ in zip(atlas_param, Gyrus_list):
        # stat = OrderedDict(
        #     (
        #         ("Volume", np.nan),
        #         ("Volume density", np.nan),
        #         ('Vessel length', np.nan),
        #         ('Vessel length density', np.nan),
        #         ("Diameter mean", np.nan),
        #         ("Diameter median", np.nan),
        #         ("Diameter max", np.nan),
        #         ("Diameter min", np.nan),
        #         ("Diameter STD", np.nan),
        #         ("Diameter N", np.nan)
        #     )
        # )
        Region=pd.Series()
        mean_zoom_cube = np.mean(header.get_zooms()) ** 3  # 单个体素分辨率

        cl_gyrus_dia=cl_diameter*gyrus_ # 中心线*脑区mask
        seg_diameter: np.ndarray = cl_gyrus_dia[cl_gyrus_dia != 0].flatten()
        pict=tof.get_fdata()*gyrus_# 总血管*脑区mask
        # 体积
        # stat["Volume"] = np.count_nonzero(pict) * mean_zoom_cube
        # stat["Volume density"] = np.count_nonzero(pict) /np.sum(gyrus_ == 1)
        Region["Volume"] = np.count_nonzero(pict) * mean_zoom_cube
        Region["Volume density"] = np.count_nonzero(pict) /np.sum(gyrus_ == 1)
    #血管长度密度
        # stat['Vessel length density'] = np.count_nonzero(cl * gyrus_) * np.mean(header.get_zooms()) / (
        #             np.sum(gyrus_ == 1) * mean_zoom_cube)
        # stat['Vessel length'] =np.count_nonzero(cl * gyrus_) * np.mean(header.get_zooms())
        Region['Vessel length density'] = np.count_nonzero(cl * gyrus_) * np.mean(header.get_zooms()) / (
                np.sum(gyrus_ == 1) * mean_zoom_cube)
        Region['Vessel length'] = np.count_nonzero(cl * gyrus_) * np.mean(header.get_zooms())
    # 血管直径
    #     stat["Diameter N"] = seg_diameter.size
        Region["Diameter N"] = seg_diameter.size
        if seg_diameter.size == 0:
            print(0)
            # stat["Diameter max"] = 0
            # stat["Diameter min"] = 0
            # stat["Diameter mean"] = 0
            # stat["Diameter median"] = 0
            # stat["Diameter STD"] = 0
            Region["Diameter max"] = 0
            Region["Diameter min"] = 0
            Region["Diameter mean"] = 0
            Region["Diameter median"] = 0
            Region["Diameter STD"] = 0
        else:
            # stat["Diameter mean"] = seg_diameter.mean()
            # stat["Diameter max"] = seg_diameter.max()
            # stat["Diameter min"] = seg_diameter.min()
            # stat["Diameter median"] = np.median(seg_diameter)
            # stat["Diameter STD"] = seg_diameter.std()
            Region["Diameter mean"] = seg_diameter.mean()
            Region["Diameter max"] = seg_diameter.max()
            Region["Diameter min"] = seg_diameter.min()
            Region["Diameter median"] = np.median(seg_diameter)
            Region["Diameter STD"] = seg_diameter.std()

        length_densitys[np.where(gyrus_ == 1)] = Region['Vessel length density']
        volumes_normal[np.where(gyrus_ == 1)] = Region["Volume density"]
        lengths[np.where(gyrus_ == 1)] = Region['Vessel length']
        volumes[np.where(gyrus_ == 1)] = Region["Volume"]
        diameters[np.where(gyrus_ == 1)] = Region["Diameter mean"]


        Region.name=str(index)
        SUMMARY=SUMMARY.append(Region)

        # csv_name=str(index)+"_summary.csv" #脑网络组图谱
        # pd.DataFrame.from_dict(stat, orient='index').to_csv(os.path.join(gyrus_path, csv_name))


        print(f" {index} finished.")

    csv_name = "summary.csv"  # 脑网络组图谱
    SUMMARY.to_csv(os.path.join(res_path, csv_name), mode='w', index=True)

    save_path = os.path.join(res_path, 'length_density.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(length_densitys, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'volume_density.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(volumes_normal, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'length.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(lengths, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'volume.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(volumes, tof.affine, tof.header), save_path)
    save_path = os.path.join(res_path, 'diameter_mean.nii.gz')  # 保存脑区对应血管图像
    nib.save(nib.Nifti1Image(diameters, tof.affine, tof.header), save_path)


if __name__ == '__main__':
    dir='C:/Users/86133/Documents/pre-graduate-tasks/app_test/MSC0303'
    seg_=os.path.join(dir,'MRA_brain_vessel.nii.gz')
    mask_=os.path.join(dir,'MRA_brain_mask.nii.gz')
    atlas_=os.path.join(dir,'Atlas.nii.gz')
    atlas_region_features_extraction_no_csv(seg_, mask_, atlas_)
    # atlas_=os.path.join(dir,'Atlas_resample.nii.gz')
    # csv_=os.path.join(dir,'atlas.csv')
    # atlas_region_features_extraction(seg_,mask_,atlas_,csv_)






