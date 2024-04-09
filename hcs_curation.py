import os
import pydicom as dicom
import nibabel as nib

import shutil
from tqdm import tqdm

def write_nifti_outside_series(path, file):
    file.write(f"\nFound NIFTI file: {path} outside of Series folder.\n")
    file.write("The medical expert should examine the annotation file and place it in the corresponding Series folder.\n")

def write_file_path(series_path, file):
    file.write(f"Issue in file: {series_path}\n")

def write_num_dicom(num_dicom_files, file):
    file.write(f"Number of DICOM files: {num_dicom_files}\n")

def write_num_nifti(num_nifti_files, file):
    file.write(f"Number of NIFTI files: {num_nifti_files}\n")

def write_num_nifti_slices(num_nifti_slices, file):
    file.write(f"Number of NIFTI slices: {num_nifti_slices}\n")

def write_num_sequences(num_sequences, file):
    file.write(f"Number of different DICOM sequences: {num_sequences}\n")

def write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)
    write_num_sequences(num_sequences, file)

def write_empty_series(series_path, file):
    file.write(f"\nFound empty Series: {series_path}.\n")

def write_found_only_nifti(series_path, num_nifti_files, num_nifti_slices, file):
    file.write(f"\nFound only NIFTI files\n")
    write_file_path(series_path, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)

def write_found_only_dicom(series_path, num_dicom_files, num_sequences, file):
    file.write(f"\nFound only DICOM files. Annotation is required.\n")
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_sequences(num_sequences, file)

def write_split_is_required(series_path, num_sequences, file):
    file.write(f"\nSeries folder: {series_path} needs to be split in {num_sequences} different Series folders.\n")

def write_split_auto(file):
    file.write(f"Series folder can be split automatically.\n")

def write_split_manual(file):
    file.write(f"Series folder cannot be split automatically. Manual split should be performed by a medical expert.\n")

def write_dicom_nifti_slices_incompatible(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    file.write(f"\nFound incompatible number of DICOM files and NIFTI slices.\n")
    file.write("A medical expert should examine the specific Series folder and correct it.\n")
    write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)

def write_sequences_nifti_files_incompatible(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    file.write(f"\nFound incompatible number of DICOM sequences and NIFTI files.\n")
    file.write("A medical expert should examine the specific Series folder and correct it.\n")
    write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)

def write_more_nifti_than_dicom(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    file.write(f"\nFound more NIFTI files than DICOM files.\n")
    write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)

def get_series_attributes(serie, file):

    # Get the filenames of all files inside Series directory
    filenames = os.listdir(serie)

    # Initialize Series attributes
    num_dicom_files = 0
    num_nifti_files = 0
    num_nifti_slices = 0
    num_sequences = 0
    nifti_slices = []
    dicom_sequences = []
    
    # Loop over each file in Series
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

def get_crc_filename(filename):
    parent_dir = os.path.abspath(os.path.join(filename, os.pardir))
    basename = os.path.basename(filename)
    crc_basename = f".{basename}.crc"
    crc_path = os.path.join(parent_dir, crc_basename)

    return crc_path, crc_basename

def check_if_crc_exists(filename):
    
    crc_path, _ = get_crc_filename(filename)
    if os.path.exists(crc_path):
        return True
    else:
        return False

def get_filenames_for_each_sequence(current_series_path):

    dicom_filenames = [os.path.join(current_series_path, filename) for filename in os.listdir(current_series_path) if filename.endswith('.dcm')]
    nifti_filenames = [os.path.join(current_series_path, filename) for filename in os.listdir(current_series_path) if filename.endswith('.nii.gz') or filename.endswith('.nii')]
    
    basenames  = [os.path.basename(filename).split('.')[0] for filename in dicom_filenames if filename.endswith('.dcm')]

    # Initialize the dictionary containing which filenames correspond to which sequence
    filenames_for_sequence = {}
    for basename in basenames:
        filenames_for_sequence[f"{basename}"] = []

    # print(filenames_for_sequence)
    # print_filenames(dicom_filenames)

    for filename in dicom_filenames:
        base = os.path.basename(filename).split('.')[0]
        filenames_for_sequence[f"{base}"].append(filename)
    
    # Add nifti filenames in the corresponding list
    if len(nifti_filenames) > 0:
        for filename in nifti_filenames:
            # print(os.path.basename(filename).split('.')[0])
            for filenames_list in filenames_for_sequence.values():
                basenames_list = [os.path.basename(path).split('.')[0] for path in filenames_list]
                # print(basenames_list)
                if os.path.basename(filename).split('.')[0] in basenames_list:
                    filenames_list.append(filename)
                    break

    # for filenames_list in filenames_for_sequence.values():
    #     print(filenames_list)

    return filenames_for_sequence

def find_next_available_series_name(series_names_in_parent_dir):

    counter = 1
    new_series_name = f"Series-{counter}"
    while (f"Series-{counter}" in series_names_in_parent_dir):
        counter += 1
        new_series_name = f"Series-{counter}"

    return new_series_name

def split_series(current_series_path):

    parent_dir_path = os.path.abspath(os.path.join(current_series_path, os.pardir))
    series_names_in_parent_dir = [os.path.basename(filename) for filename in os.listdir(parent_dir_path)]

    # print(f"Parent directory path: {parent_dir_path}")

    filenames_for_sequence = get_filenames_for_each_sequence(current_series_path)

    # print_dict_values(filenames_for_sequence)

    for key in filenames_for_sequence:
        try:
            next_name = find_next_available_series_name(series_names_in_parent_dir)
        except:
            print(f"Error at {current_series_path}")
        series_names_in_parent_dir.append(next_name)
        # print(f"Next avalailable Series name: {next_name}")

        # Normally I would mkdir with the new name
        new_dir_path = os.path.join(parent_dir_path, next_name)
        if not os.path.exists(new_dir_path):
            os.mkdir(new_dir_path)

        # Then copy all cooresponding images per key to the new direcory
        for image_path in filenames_for_sequence[key]:
            basename = os.path.basename(image_path)
            # print(basename)
            move_path = os.path.join(new_dir_path, basename)
            # print(f"Dst dicom path: {move_path}")
            shutil.move(image_path, move_path)

            if check_if_crc_exists(image_path):
                crc_path, crc_basename = get_crc_filename(image_path)
                move_crc_path = os.path.join(new_dir_path, crc_basename)
                # print(f"Dst dicom path: {move_crc_path}")
                shutil.move(crc_path, move_crc_path)

    # print(f"Series names in parent directory: {series_names_in_parent_dir}")

    # When finishing, maybe I need to delete the empty series directory
    if os.path.isdir(current_series_path):
        if not os.listdir(current_series_path):
            shutil.rmtree(current_series_path)
            print(f"Deleted the empty directory: {current_series_path}")

def write_issues_to_report(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences, file):

    # If there are no DICOM files
    if num_dicom_files == 0: 
        if num_nifti_files == 0: # Empty folder -> Issue
            # write_empty_series(serie, file)
            os.rmdir(serie)
        # Only NIFTI files found
        else: # Issue
            write_found_only_nifti(serie, num_nifti_files, num_nifti_slices, file)
    # DICOM files exist
    else: 
        # Series not annotated
        if num_nifti_files == 0: # Issue
                # write_found_only_dicom(serie, num_dicom_files, num_sequences, file)
                if num_sequences >= 1: # Multiple sequence - needs to be split -> Issue
                    # We can handle this with split_data
                    # write_split_is_required(serie, num_sequences, file)
                    # write_split_auto(file)
                    split_series(serie)
                    # write_series_info(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
        # More NIFTI than DICOM files 
        elif num_dicom_files < num_nifti_files: # Excessive NIFTI files -> Issue
            write_more_nifti_than_dicom(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
        # More DICOM than NIFTI
        elif num_dicom_files > num_nifti_files: 
            if num_sequences == num_nifti_files: # 1 sequence per NIFTI file
                if num_dicom_files == num_nifti_slices: # DICOM match to NIFTI slices
                    if num_sequences > 1: # Multiple sequence - needs to be split -> Issue
                        write_split_is_required(serie, num_sequences, file)
                        write_split_manual(file)
                        write_series_info(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
                else: # DICOM don't match to NIFTI slices -> Issue
                    # write_dicom_nifti_slices_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
                    split_series(serie)
            else: # Sequences don't match the number of NIFTI files
                # write_sequences_nifti_files_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices,  num_sequences, file)
                split_series(serie)
        # num_dicom_files = num_nifti_files
        else: 
            if num_sequences == num_nifti_files: # 1 sequence per NIFTI file
                if num_dicom_files == num_nifti_slices: # DICOM match to NIFTI slices
                    if num_sequences > 1: # Multiple sequence - needs to be split -> Issue
                        # We can handle this with split_annotated_data
                        # write_split_is_required(serie, num_sequences, file)
                        # write_split_auto(file)
                        split_series(serie)
                        # write_series_info(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
                else: # DICOM don't match to NIFTI slices -> Issue
                    write_dicom_nifti_slices_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
            else: # Sequences don't match the number of NIFTI files
                if num_sequences >= 1:
                    # write_sequences_nifti_files_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices,  num_sequences, file)
                    split_series(serie)

if __name__ == "__main__":

    cancer_type = "breast"
    data_provider = "hcs"

    # Define the report file path
    file_path = f"reports/{data_provider}/{cancer_type}/report_issues.txt"

    # Remove the existing report file if it exists
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
            print(f"{file_path} has been deleted successfully.")
        except OSError as e:
            print(f"Error deleting {file_path}: {e}")

    print(f"\nCreating report for {data_provider}-{cancer_type}.\n")

    database_path = "/mnt/nfs/incisive"
    # working_path = f"hcs/breastCopy"
    # working_path = f"{database_path}/hcs/breast/"
    working_path = "/mnt/nfs/incisive3/hcs/breast"

    with open(file_path, "a") as file:

        # Get all patients paths
        patients_paths = [os.path.join(working_path, filename) for filename in os.listdir(working_path) if os.path.isdir(os.path.join(working_path, filename))]

        # print(patients_paths)

        # Get all studies paths
        studies_paths = []
        for patient_path in patients_paths:
            studies = os.listdir(patient_path)
            for study in studies:
                if os.path.join(patient_path, study).endswith(".nii") or os.path.join(patient_path, study).endswith(".nii.gz"): # Check if NIFTI outside of Series -> Issue
                    nifti_path = os.path.join(patient_path, study)
                    write_nifti_outside_series(nifti_path, file)
                elif os.path.isdir(os.path.join(patient_path, study)):
                    studies_paths.append(os.path.join(patient_path, study))

        # print(studies_paths)

        # Get all series paths
        series_paths = []
        for study_path in studies_paths:
            series = os.listdir(study_path)
            for serie in series:
                if os.path.join(study_path, serie).endswith(".nii") or os.path.join(study_path, serie).endswith(".nii.gz"): # Check if NIFTI outside of Series -> Issue
                    nifti_path = os.path.join(study_path, serie)
                    write_nifti_outside_series(nifti_path, file)
                elif os.path.isdir(os.path.join(study_path, serie)):
                    series_paths.append(os.path.join(study_path, serie))
                
        # print(len(series_paths))

        # For each series get the attributes and write existing issues
        for series_path in tqdm(series_paths, total=len(series_paths)):
            # For each serie get attributes
            num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences = get_series_attributes(series_path, file)

            # Write the issues regarding the contents of the series the .txt file
            write_issues_to_report(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences, file)