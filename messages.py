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

