import numpy as np
import nibabel as nib

# Input and output paths
v1_img_path = "path/to/V1.nii.gz"
theta_img_path = "path/to/theta.nii.gz"

# Load the V1 image
v1_img = nib.load(v1_img_path)

# Extract the V1 data as a NumPy array
v1_data = v1_img.get_fdata()

# Calculate sin(theta)^2
numerator = np.sum(v1_data[:, :, :, 0:2] ** 2, axis=3)
denominator = np.sum(v1_data ** 2, axis=3) ** 2
sin_theta_squared = numerator / denominator

# Calculate theta in degrees
theta_sqrt = np.sqrt(sin_theta_squared)
theta = np.arcsin(theta_sqrt) * (180 / np.pi)

# Replace NaN and Inf with 0
theta = np.nan_to_num(theta, nan=0.0, posinf=0.0, neginf=0.0)

# Create and save the NIfTI image
theta_img = nib.Nifti1Image(theta, v1_img.affine, v1_img.header)
nib.save(theta_img, theta_img_path)

print(f"Theta image saved as {theta_img_path}")
