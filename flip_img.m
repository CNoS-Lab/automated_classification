function [to_match,file] = flip_img(to_match,filepath,file)

nii_img = load_nii([filepath file '.nii']);
to_match = to_match*(-1);
file = [file '_flipped.nii'];
nii_img.img = to_match;
save_nii(nii_img, [filepath,file]);

end