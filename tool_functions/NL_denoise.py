import os
import skimage.restoration as sk_res
import nibabel as nib
import SimpleITK as sitk


def NL_denoise(image_path):
    input_image=nib.load(image_path)
    image=input_image.get_fdata()
    # output_image=sk_res.denoise_nl_means(image,fast_mode=True)
    output_image=sk_res.denoise_nl_means(image, patch_size=1, patch_distance=1, h=0.1, sigma=0.1,fast_mode=True)
    input_ = os.path.split(image_path)
    save_path = os.path.join(input_[0], 'MRA_brain_denoise.nii.gz')
    nib.save(nib.Nifti1Image(output_image, input_image.affine, input_image.header), save_path)


def bias_field_correct(input_image_path, mask_image_path):
    input_image=sitk.ReadImage(input_image_path, sitk.sitkFloat32)
    mask_image=sitk.ReadImage(mask_image_path, sitk.sitkUInt8)
    corrector=sitk.N4BiasFieldCorrectionImageFilter()
    output_image = corrector.Execute(input_image, mask_image)
    output_image = sitk.Cast(output_image, sitk.sitkFloat32)
    input_ = os.path.split(input_image_path)
    save_path = os.path.join(input_[0], 'MRA_brain_biasfield_correct.nii.gz')
    sitk.WriteImage(output_image, save_path)

