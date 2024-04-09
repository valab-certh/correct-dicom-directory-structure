import os
import pydicom as dicom
import nibabel as nib
import shutil

from tqdm import tqdm

def delete_crc_files(root_directory):
    """
    Recursively deletes files with the '.crc' extension in the specified root directory and its subdirectories.

    :param root_directory: The root directory to start searching for '.crc' files.
    """
    for patient_folder in os.listdir(root_directory):
        patient_path = os.path.join(root_directory, patient_folder)
        if os.path.isdir(patient_path):
            delete_crc_files_in_patient(patient_path)

def delete_crc_files_in_patient(patient_path):
    """
    Recursively deletes '.crc' files in a patient's directory.

    :param patient_path: The path to the patient's directory.
    """
    for study_folder in os.listdir(patient_path):
        study_path = os.path.join(patient_path, study_folder)
        if os.path.isdir(study_path):
            delete_crc_files_in_study(study_path)

def delete_crc_files_in_study(study_path):
    """
    Recursively deletes '.crc' files in a study's directory.

    :param study_path: The path to the study's directory.
    """
    for series_folder in os.listdir(study_path):
        series_path = os.path.join(study_path, series_folder)
        if os.path.isdir(series_path):
            delete_crc_files_in_series(series_path)

def delete_crc_files_in_series(series_path):
    """
    Deletes '.crc' files in a series' directory.

    :param series_path: The path to the series' directory.
    """
    for file_name in os.listdir(series_path):
        if file_name.endswith('.crc'):
            file_path = os.path.join(series_path, file_name)
            try:
                os.remove(file_path)
                print(f"Deleted file: {file_path}")
            except PermissionError as e:
                print(f"Permission Error: {e}")

def replace_folders(destination_directory, source_directory):
    """
    Replace directories in the destination directory with those from the source directory.

    :param destination_directory: The destination directory where directories will be replaced.
    :param source_directory: The source directory containing directories to be copied and replaced.
    """
    print("Replacing common directories")
    
    # Get a list of all subdirectories in the source directory
    source_subdirectories = [os.path.join(source_directory, d) for d in os.listdir(source_directory) if os.path.isdir(os.path.join(source_directory, d))]

    for source_subdir in tqdm(source_subdirectories, total=len(source_subdirectories)):
        subdir_name = os.path.basename(source_subdir)
        destination_subdir = os.path.join(destination_directory, subdir_name)

        if os.path.exists(destination_subdir):
            # Directory with the same name already exists in destination, so replace it
            shutil.rmtree(destination_subdir)
        shutil.copytree(source_subdir, destination_subdir)

    print("Directories copied and replaced successfully.")

def get_directory_filenames(directory):
    """
    Get the list of filenames in the specified directory.

    :param directory: The directory path.
    :return: A list of filenames in the directory.
    """
    return os.listdir(directory)

def print_filenames(filenames):
    """
    Print the base names of filenames in a list.

    :param filenames: A list of filenames.
    """
    for filename in filenames:
        print(os.path.basename(filename))

def print_dict_values(dictionary):
    """
    Print the values of a dictionary.

    :param dictionary: A dictionary.
    """
    for key in dictionary:
        print(dictionary[key])

def get_crc_filename(filename):
    """
    Get the CRC filename and path corresponding to a given filename.

    :param filename: The input filename.
    :return: The CRC path and CRC basename.
    """
    parent_dir = os.path.abspath(os.path.join(filename, os.pardir))
    basename = os.path.basename(filename)
    crc_basename = f".{basename}.crc"
    crc_path = os.path.join(parent_dir, crc_basename)

    return crc_path, crc_basename

def check_if_crc_exists(filename):
    """
    Check if a CRC file exists for a given filename.

    :param filename: The input filename.
    :return: True if the CRC file exists, False otherwise.
    """
    crc_path, _ = get_crc_filename(filename)
    return os.path.exists(crc_path)

def get_filenames_for_each_sequence(current_series_path, dicom_sequences):
    """
    Get a dictionary that maps sequence names to corresponding filenames in a series path.

    :param current_series_path: The path to the series directory.
    :param dicom_sequences: List of sequence names to be searched for.
    :return: A dictionary where keys are sequence names and values are lists of filenames.
    """
    dicom_filenames = [os.path.join(current_series_path, filename) for filename in os.listdir(current_series_path) if filename.endswith('.dcm')]
    nifti_filenames = [os.path.join(current_series_path, filename) for filename in os.listdir(current_series_path) if filename.endswith('.nii.gz') or filename.endswith('.nii')]
    
    filenames_for_sequence = {seq_name: [] for seq_name in dicom_sequences}

    for filename in dicom_filenames:
        try:
            img = dicom.dcmread(filename)
            image_modality = img.Modality
            
            if image_modality == "MR":
                try:
                    tag = img.SequenceName
                    filenames_for_sequence[tag].append(filename)
                except AttributeError:
                    continue
            elif image_modality == "MG":
                try:
                    tag = img.SeriesTime
                    filenames_for_sequence[tag].append(filename)
                except AttributeError:
                    continue
        except AttributeError:
            break
    
    if len(nifti_filenames) > 0:
        for filename in nifti_filenames:
            for filenames_list in filenames_for_sequence.values():
                basenames_list = [os.path.basename(path).split('.')[0] for path in filenames_list]
                # print(basenames_list)
                if os.path.basename(filename).split('.')[0] in basenames_list:
                    filenames_list.append(filename)
                    break

    return filenames_for_sequence

def find_next_available_series_name(series_names_in_parent_dir):
    """
    Find the next available series name not present in the parent directory.

    :param series_names_in_parent_dir: List of series names in the parent directory.
    :return: The next available series name.
    """
    counter = 1
    new_series_name = f"Series-{counter}"
    while new_series_name in series_names_in_parent_dir:
        counter += 1
        new_series_name = f"Series-{counter}"
    return new_series_name

def split_series(current_series_path, names_of_sequences):
    """
    Split a series into multiple series based on sequence names.

    :param current_series_path: The path to the current series directory.
    :param names_of_sequences: List of sequence names to split the series.
    """
    parent_dir_path = os.path.abspath(os.path.join(current_series_path, os.pardir))
    series_names_in_parent_dir = [os.path.basename(filename) for filename in os.listdir(parent_dir_path)]

    filenames_for_sequence = get_filenames_for_each_sequence(current_series_path, names_of_sequences)

    for key in filenames_for_sequence:
        try:
            next_name = find_next_available_series_name(series_names_in_parent_dir)
        except Exception as e:
            print(f"Error at {current_series_path}: {e}")
        series_names_in_parent_dir.append(next_name)

        new_dir_path = os.path.join(parent_dir_path, next_name)
        if not os.path.exists(new_dir_path):
            os.mkdir(new_dir_path)

        for image_path in filenames_for_sequence[key]:
            basename = os.path.basename(image_path)
            move_path = os.path.join(new_dir_path, basename)
            shutil.move(image_path, move_path)

            if check_if_crc_exists(image_path):
                crc_path, crc_basename = get_crc_filename(image_path)
                move_crc_path = os.path.join(new_dir_path, crc_basename)
                shutil.move(crc_path, move_crc_path)

def write_nifti_outside_series(path, file):
    """
    Write a message about a NIFTI file found outside of a Series folder.

    :param file: The file to write to.
    :param path: The path of the NIFTI file.
    """
    file.write(f"\nFound NIFTI file: {path} outside of Series folder.\n")
    file.write("The medical expert should examine the annotation file and place it in the corresponding Series folder.\n")

def write_patient_not_pass_qc(path, file):
    """
    Write a message about a patient folder that did not pass the Quality Check tool.

    :param file: The file to write to.
    :param path: The path of the Study.
    """
    file.write(f"\nFound Study folder: {path} that did not pass the Quality Check.\n")
    file.write("The Data Provider should pass this patient throught the Quality Check.\n")

def write_file_path(series_path, file):
    """
    Write the path of a file to the output file.

    :param file: The file to write to.
    :param series_path: The path of the file to be written.
    """
    file.write(f"Issue in file: {series_path}\n")

def write_num_dicom(num_dicom_files, file):
    """
    Write the number of DICOM files to the output file.

    :param file: The file to write to.
    :param num_dicom_files: The number of DICOM files.
    """
    file.write(f"Number of DICOM files: {num_dicom_files}\n")

def write_num_nifti(num_nifti_files, file):
    """
    Write the number of NIFTI files to the output file.

    :param file: The file to write to.
    :param num_nifti_files: The number of NIFTI files.
    """
    file.write(f"Number of NIFTI files: {num_nifti_files}\n")

def write_num_nifti_slices(num_nifti_slices, file):
    """
    Write the number of NIFTI slices to the output file.

    :param file: The file to write to.
    :param num_nifti_slices: The number of NIFTI slices.
    """
    file.write(f"Number of NIFTI slices: {num_nifti_slices}\n")

def write_num_sequences(num_sequences, file):
    """
    Write the number of different DICOM sequences to the output file.

    :param file: The file to write to.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write(f"Number of different DICOM sequences: {num_sequences}\n")

def write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    """
    Write series information to the output file.

    :param file: The file to write to.
    :param series_path: The path of the series.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)
    write_num_sequences(num_sequences, file)

def write_empty_series(series_path, file):
    """
    Write a message about an empty Series folder to the output file.

    :param file: The file to write to.
    :param series_path: The path of the empty Series folder.
    """
    file.write(f"\nFound empty Series: {series_path}.\n")

def write_found_only_nifti(series_path, num_nifti_files, num_nifti_slices, file):
    """
    Write a message about a Series folder found with only NIFTI files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    """
    file.write(f"\nFound only NIFTI files\n")
    write_file_path(series_path, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)

def write_found_only_dicom(series_path, num_dicom_files, num_sequences, file):
    """
    Write a message about a Series folder found with only DICOM files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write(f"\nFound only DICOM files. Annotation is required.\n")
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_sequences(num_sequences, file)

def write_split_is_required(series_path, num_sequences, file):
    """
    Write a message about a Series folder requiring splitting to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_sequences: The number of sequences to split.
    """
    file.write(f"\nSeries folder: {series_path} needs to be split in {num_sequences} different Series folders.\n")

def write_split_auto(file):
    """
    Write a message that automatic splitting of a Series folder is possible to the output file.

    :param file: The file to write to.
    """
    file.write(f"Series folder can be split automatically.\n")

def write_split_manual(file):
    """
    Write a message that manual splitting of a Series folder is required to the output file.

    :param file: The file to write to.
    """
    file.write(f"Series folder cannot be split automatically. Manual split should be performed by a medical expert.\n")

def write_dicom_nifti_slices_incompatible(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    """
    Write a message about incompatible number of DICOM files and NIFTI slices to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write(f"\nFound incompatible number of DICOM files and NIFTI slices.\n")
    file.write("A medical expert should examine the specific Series folder and correct it.\n")
    write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)

def write_sequences_nifti_files_incompatible(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    """
    Write a message about incompatible number of DICOM sequences and NIFTI files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write(f"\nFound incompatible number of DICOM sequences and NIFTI files.\n")
    file.write("A medical expert should examine the specific Series folder and correct it.\n")
    write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)

def write_more_nifti_than_dicom(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file):
    """
    Write a message about finding more NIFTI files than DICOM files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write(f"\nFound more NIFTI files than DICOM files.\n")
    write_series_info(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)

def write_issues_to_report(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences, file):
    """
    Write issues related to a series to a report file.

    :param serie: The series being analyzed.
    :param num_dicom_files: The number of DICOM files in the series.
    :param num_nifti_files: The number of NIFTI files in the series.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    :param dicom_sequences: List of DICOM sequences.
    :param file: The file to write the report to.
    """
    if num_dicom_files == 0: 
        if num_nifti_files == 0:
            # Issue: Empty series
            if os.path.isdir(serie):
                if not os.listdir(serie):
                    shutil.rmtree(serie)
                    print(f"Deleted the empty directory: {serie}")
            # write_empty_series(serie, file)
        else:
            # Issue: Only NIFTI files found
            write_found_only_nifti(serie, num_nifti_files, num_nifti_slices, file)
    else: 
        if num_nifti_files == 0:
                if num_sequences > 1:
                    # Issue: Series with multiple sequences requires splitting
                    split_series(serie, dicom_sequences) 
        elif num_dicom_files < num_nifti_files:
            # Issue: More NIFTI files than DICOM files
            write_more_nifti_than_dicom(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
        elif num_dicom_files > num_nifti_files: 
            if num_sequences == num_nifti_files: 
                if num_dicom_files == num_nifti_slices: 
                    if num_sequences > 1:
                        # Issue: Series with multiple sequences requires splitting 
                        write_split_is_required(serie, num_sequences, file)
                        write_split_manual(file)
                        write_series_info(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
                else:
                    # Issue: Incompatible number of DICOM files and NIFTI slices
                    write_dicom_nifti_slices_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
            else:
                # Issue: Incompatible number of DICOM sequences and NIFTI files 
                write_sequences_nifti_files_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices,  num_sequences, file)
        else: 
            if num_sequences == num_nifti_files: 
                if num_dicom_files == num_nifti_slices: 
                    if num_sequences > 1:
                        # Issue: Series with multiple sequences requires splitting
                        split_series(serie, dicom_sequences)
                else:
                    # Issue: Incompatible number of DICOM files and NIFTI slices
                    write_dicom_nifti_slices_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, file)
            else:
                if num_sequences >= 1:
                    # Issue: Incompatible number of DICOM sequences and NIFTI files
                    write_sequences_nifti_files_incompatible(serie, num_dicom_files, num_nifti_files, num_nifti_slices,  num_sequences, file)



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

def main():
    # Provide a list of the cancer types to check
    cancer_types = ["breast", "colorectal", "lung", "prostate"]

    # Select the data provider
    data_providers = {"breast": ["dp1", "dp2"],
                     "colorectal": ["dp1", "dp2"],
                     "lung": ["dp1", "dp2"],
                     "prostate": ["dp1", "dp2"]}

    for cancer_type in cancer_types:
        for data_provider in data_providers[cancer_type]:

            # Define the target directory to examine based on cancer type and data provider name
            database_path = "incisive2"
            working_path = f"{database_path}/{cancer_type}/{data_provider}/data"

            # Define the report file path
            txt_file_path = f"reports/{data_provider}/{cancer_type}/report_issues.txt"

            # Remove the existing report file if it exists
            if os.path.exists(txt_file_path):
                try:
                    os.remove(txt_file_path)
                    print(f"{txt_file_path} has been deleted successfully.")
                except OSError as e:
                    print(f"Error deleting {txt_file_path}: {e}")

            print(f"\nCreating report for {data_provider}-{cancer_type}.\n")

            with open(txt_file_path, "a") as file:
                # Get all patients paths
                patients_paths = [os.path.join(working_path, filename) for filename in os.listdir(working_path) if os.path.isdir(os.path.join(working_path, filename))]

                # Get all studies and check for NIFTI files outside of Series
                studies_paths = []
                for patient_path in patients_paths:
                    studies = os.listdir(patient_path)
                    for study in studies:
                        study_path = os.path.join(patient_path, study)
                        if study_path.endswith((".nii", ".nii.gz")):
                            # Check if NIFTI is outside of Series -> Issue
                            write_nifti_outside_series(study_path, file)
                        # TODO Include check about QC
                        elif "Study" in study_path:
                            write_patient_not_pass_qc(study_path, file)
                        elif os.path.isdir(study_path):
                            studies_paths.append(study_path)

                # Get all series paths
                series_paths = []
                for study_path in studies_paths:
                    series = os.listdir(study_path)
                    for serie in series:
                        serie_path = os.path.join(study_path, serie)
                        if serie_path.endswith((".nii", ".nii.gz")):
                            # Check if NIFTI outside of Series -> Issue
                            write_nifti_outside_series(serie_path, file)
                        elif os.path.isdir(serie_path):
                            series_paths.append(serie_path)
                        
                # For each series, get the attributes and write existing issues
                for series_path in tqdm(series_paths, total=len(series_paths)):
                    # Get attributes for the series
                    num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences = get_series_attributes(series_path, file)

                    # Write the issues regarding the contents of the series to the .txt file
                    write_issues_to_report(series_path, num_dicom_files, num_nifti_files, num_nifti_slices, num_sequences, dicom_sequences, file)

if __name__ == '__main__':
    main()