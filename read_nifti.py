import nibabel as nib

# Load the NIfTI file
nifti_img = nib.load('/mnt/nfs/incisive3/uoa/prostate/002-000190/002-000190_MR_BL/Series-501/Image-1.nii.gz')  # Replace with your file path

# Get the image data and header
img_data = nifti_img.get_fdata()
img_header = nifti_img.header

# Work with the image data as needed
# For example, you can access dimensions like this:
print("Image shape:", img_data.shape)
print("Voxel dimensions:", img_header.get_zooms())