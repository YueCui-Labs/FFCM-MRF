function BiasFieldCorrection(app,input_image_path)
    [file_dir,~,~]=fileparts(input_image_path);
    mask_img=fullfile(file_dir,'MRA_brain_mask.nii.gz');
    input_image_path=py.str(input_image_path);
    mask_img=py.str(mask_img);
    
    py.sys.path().append(fullfile(pwd,'tool_functions'));
    NL_denoise = py.importlib.import_module('NL_denoise');
    NL_denoise.bias_field_correct(input_image_path,mask_img);
end