function new_image = match_sizes(input_image,target_img)

% resize and reslice the input image to the size of the target image
% (relevant to template matching because the slice #s are relative to 
% this template's dimensions) 

[path,filename,ext] = fileparts(input_image);

target_img = load_nii(target_img);

% use the target image's values as target dimension and position

target_size = target_img.hdr.dime.dim(2:4);
target_vox_size = min(target_img.hdr.dime.pixdim(2:4));
target_origin = target_img.hdr.hist.originator(1:3);

% change the input .img to match the voxel dimensions of the target
reslice_nii(input_image, [path filesep filename '_resampled.nii'], target_vox_size,0);
new_image = load_nii([path filesep filename '_resampled.nii']);
new_image.img = double(new_image.img);
% need to match size and origins of input_image and target_image
image_size = new_image.hdr.dime.dim(2:4);
image_origin = new_image.hdr.hist.originator(1:3);

new_image.img = crop_image(new_image, image_size, image_origin, target_size, target_origin);

to_save = new_image;
to_save.img = double(permute(to_save.img,[1 2 3]));
to_save.hdr.dime.dim(2:4) = size(to_save.img);
to_save.hdr.hist.originator(1:3) = target_origin;
to_save.hdr.dime.datatype = 64;
to_save.hdr.dime.bitpix = 64;

save_nii(to_save,[path filesep filename '_resampled.nii']);

new_image = new_image.img;

end

function new_image = crop_image(input_image, input_size, input_origin, target_size, target_origin)

diff = input_size-target_size; %difference between old and new size
diff_origin = input_origin-target_origin;
image_size = input_size; % copy the original_size
image = input_image.img;

start = [1,1,1];

for i = 1:3
   if diff_origin(i)>=0 && diff(i)-diff_origin(i)>=0 
        start(i) = diff_origin(i)+1;
   elseif diff_origin(i)<0 && diff(i)-diff_origin(i)>=0
        pad_size_begin = image_size; % copy the image size
        pad_size_begin(i) = -diff_origin(i); % padding must be the negative of the difference between the origin
        pad_zeros = zeros(pad_size_begin); % create padding matrix
        image = cat(i, pad_zeros, image); % pad image
   elseif diff_origin(i)>=0 && diff(i) - diff_origin(i)<0
        start(i) = diff_origin(i)+1; % since the diff_origin is larger, we have to shave off the lower indices to make the diff_origin = component origin
        pad_size_end = image_size;
        pad_size_end(i) = -(diff(i) - diff_origin(i)); % trust me, this is a thing
        pad_zeros = zeros(pad_size_end); % create padding matrix
        image = cat(i, image, pad_zeros); % pad image at the end this time
   elseif diff_origin(i)<0 && diff(i)-diff_origin(i)<0
        pad_size_begin = image_size;
        pad_size_begin(i) = -diff_origin(i);
        pad_zeros = zeros(pad_size_begin);
        image = cat(i, pad_zeros, image);
        pad_size_end = image_size;
        pad_size_end(i) = -(diff(i) - diff_origin(i));
        pad_zeros = zeros(pad_size_end);
        image = cat(i, image, pad_zeros);
   else
       error('Something went wrong');
   end
   image_size = size(image);
end

new_image = image(start(1):start(1)+target_size(1)-1, start(2):start(2)+target_size(2)-1, start(3):start(3)+target_size(3)-1);
end