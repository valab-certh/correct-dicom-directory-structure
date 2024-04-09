import os
import pydicom as dicom
import nibabel as nib

def get_series_attributes(serie, file):
    """
    Get attributes of a series including the number of DICOM files, NIFTI files, NIFTI slices, sequences, and DICOM sequences.

    :param serie: The directory path of the series.
    :param file: The file to write potential issues to.
    :return: A tuple containing (num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences).
    """
    # Get the filenames of all files inside the Series directory
    filenames = os.listdir(serie)

    # Initialize Series attributes
    num_dicom_files = 0
    num_nifti_files = 0
    num_nifti_slices = 0
    num_sequences = 0
    nifti_slices = []
    dicom_sequences = []
    
    # Loop over each file in the Series
    for filename in filenames:
        file_path = os.path.join(serie, filename)

        if file_path.endswith(".dcm"): # Check if DICOM
            # Increase the total DICOM files
            num_dicom_files += 1
            
            # Read DICOM file. It might be corrupted.
            try:
                img = dicom.dcmread(file_path)
            except:
                file.write(f"\nDICOM file {file_path} is missing header information and cannot be read.\n")
                continue

            # Get the modality attribute DICOM tag  
            image_modality = img.Modality
            
            if image_modality == "MR": # Check if MR, then get the SequenceName DICOM tag
                try:
                    tag = img.SequenceName
                    dicom_sequences.append(tag)
                except AttributeError:
                    continue
            elif image_modality == "MG": # Check if MR, then get the SeriesTime DICOM tag
                try:
                    tag = img.SeriesTime
                    dicom_sequences.append(tag)
                except AttributeError:
                    continue
        elif file_path.endswith(".nii") or file_path.endswith(".nii.gz"): # Check if NIFTI
            # Increase the total NIFTI files
            num_nifti_files += 1

            # Read NIFTI file. It might be corrupted.
            try:
                img = nib.load(file_path)
            except:
                file.write(f"\nNIFTI file {file_path} is missing tags and cannot be read.\n")
                continue

            header = img.header
            nifti_files_shape = header.get_data_shape()

            # Append NIFTI slices. NIFTI might have invalid dimensions
            try:
                nifti_slices.append(nifti_files_shape[2])
            except IndexError:
                file.write(f"\nNIFTI File: {file_path} has no valid dimensions.\n")
        else: # In case neither DICOM nor NIFTI file.
            continue

    # Sum total NIFTI slices
    num_nifti_slices = sum(nifti_slices)

    # Get the number of the total DICOM sequences
    if len(dicom_sequences) > 0:
        num_sequences = len(list(set(dicom_sequences)))
    else:
        num_sequences = num_nifti_files

    return num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences
