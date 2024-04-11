import contextlib
import os
import shutil
from pathlib import Path

import nibabel as nib
import pydicom as dicom
from tqdm import tqdm


def delete_crc_files_in_series(series_path: Path) -> None:
    """
    Deletes '.crc' files in a series' directory.

    :param series_path: The path to the series' directory.
    """
    for file_name in os.listdir(series_path):
        if file_name.suffix == ".crc":
            file_path = Path(series_path, file_name)
            with contextlib.suppress(PermissionError):
                Path.unlink(file_path)


def delete_crc_files_in_study(study_path: Path) -> None:
    """
    Recursively deletes '.crc' files in a study's directory.

    :param study_path: The path to the study's directory.
    """
    for series_folder in os.listdir(study_path):
        series_path = Path(study_path, series_folder)
        if Path.is_dir(series_path):
            delete_crc_files_in_series(series_path)


def delete_crc_files_in_patient(patient_path: Path) -> None:
    """
    Recursively deletes '.crc' files in a patient's directory.

    :param patient_path: The path to the patient's directory.
    """
    for study_folder in os.listdir(patient_path):
        study_path = Path(patient_path, study_folder)
        if Path.is_dir(study_path):
            delete_crc_files_in_study(study_path)


def delete_crc_files(root_directory: Path) -> None:
    """
    Recursively deletes files with the '.crc' extension in the specified root directory and its subdirectories.

    :param root_directory: The root directory to start searching for '.crc' files.
    """
    for patient_folder in os.listdir(root_directory):
        patient_path = Path(root_directory, patient_folder)
        if Path.is_dir(patient_path):
            delete_crc_files_in_patient(patient_path)


def replace_folders(destination_directory: Path, source_directory: Path) -> None:
    """
    Replace directories in the destination directory with those from the source directory.

    :param destination_directory: The destination directory where directories will be replaced.
    :param source_directory: The source directory containing directories to be copied and replaced.
    """

    # Get a list of all subdirectories in the source directory
    source_subdirectories = [
        Path(source_directory, d)
        for d in os.listdir(source_directory)
        if Path.is_dir(Path(source_directory, d))
    ]

    for source_subdir in tqdm(source_subdirectories, total=len(source_subdirectories)):
        subdir_name = Path.name(source_subdir)
        destination_subdir = Path(destination_directory, subdir_name)

        if Path.exists(destination_subdir):
            # Directory with the same name already exists in destination, so replace it
            shutil.rmtree(destination_subdir)
        shutil.copytree(source_subdir, destination_subdir)


def get_directory_filenames(directory: Path) -> list:
    """
    Get the list of filenames in the specified directory.

    :param directory: The directory path.
    :return: A list of filenames in the directory.
    """
    return os.listdir(directory)

def get_crc_filename(filename: Path) -> tuple:
    """
    Get the CRC filename and path corresponding to a given filename.

    :param filename: The input filename.
    :return: The CRC path and CRC basename.
    """
    parent_dir = Path.resolve(Path(filename, os.pardir))
    basename = Path.name(filename)
    crc_basename = f".{basename}.crc"
    crc_path = Path(parent_dir, crc_basename)

    return crc_path, crc_basename

def check_if_crc_exists(filename: Path) -> Path:
    """
    Check if a CRC file exists for a given filename.

    :param filename: The input filename.
    :return: True if the CRC file exists, False otherwise.
    """
    crc_path, _ = get_crc_filename(filename)
    return Path.exists(crc_path)

def get_filenames_for_each_sequence(
    current_series_path: Path,
    dicom_sequences: list,
) -> dict:
    """
    Get a dictionary that maps sequence names to corresponding filenames in a series path.

    :param current_series_path: The path to the series directory.
    :param dicom_sequences: List of sequence names to be searched for.
    :return: A dictionary where keys are sequence names and values are lists of filenames.
    """
    dicom_filenames = [
        Path(current_series_path, filename)
        for filename in os.listdir(current_series_path)
        if filename.suffix == ".dcm"
    ]
    nifti_filenames = [
        Path(current_series_path, filename)
        for filename in os.listdir(current_series_path)
        if filename.suffix in (".nii.gz", ".nii")
    ]

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
                basenames_list = [
                    Path.name(path).split(".")[0] for path in filenames_list
                ]
                if Path.name(filename).split(".")[0] in basenames_list:
                    filenames_list.append(filename)
                    break

    return filenames_for_sequence


def find_next_available_series_name(series_names_in_parent_dir: list) -> str:
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


def split_series(current_series_path: Path, names_of_sequences: list) -> None:
    """
    Split a series into multiple series based on sequence names.

    :param current_series_path: The path to the current series directory.
    :param names_of_sequences: List of sequence names to split the series.
    """
    parent_dir_path = Path.resolve(Path(current_series_path, os.pardir))
    series_names_in_parent_dir = [
        Path.name(filename) for filename in os.listdir(parent_dir_path)
    ]

    filenames_for_sequence = get_filenames_for_each_sequence(
        current_series_path,
        names_of_sequences,
    )

    for key in filenames_for_sequence:
        with contextlib.suppress(Exception):
            next_name = find_next_available_series_name(series_names_in_parent_dir)
        series_names_in_parent_dir.append(next_name)

        new_dir_path = Path(parent_dir_path, next_name)
        if not Path.exists(new_dir_path):
            Path.mkdir(new_dir_path)

        for image_path in filenames_for_sequence[key]:
            basename = Path.name(image_path)
            move_path = Path(new_dir_path, basename)
            shutil.move(image_path, move_path)

            if check_if_crc_exists(image_path):
                crc_path, crc_basename = get_crc_filename(image_path)
                move_crc_path = Path(new_dir_path, crc_basename)
                shutil.move(crc_path, move_crc_path)


def write_nifti_outside_series(path: Path, file: Path) -> None:
    """
    Write a message about a NIFTI file found outside of a Series folder.

    :param file: The file to write to.
    :param path: The path of the NIFTI file.
    """
    file.write(f"\nFound NIFTI file: {path} outside of Series folder.\n")
    file.write(
        "The medical expert should examine the annotation file and place it in the corresponding Series folder.\n",
    )


def write_patient_not_pass_qc(path: Path, file: Path) -> None:
    """
    Write a message about a patient folder that did not pass the Quality Check tool.

    :param file: The file to write to.
    :param path: The path of the Study.
    """
    file.write(f"\nFound Study folder: {path} that did not pass the Quality Check.\n")
    file.write(
        "The Data Provider should pass this patient throught the Quality Check.\n",
    )


def write_file_path(series_path: Path, file: Path) -> None:
    """
    Write the path of a file to the output file.

    :param file: The file to write to.
    :param series_path: The path of the file to be written.
    """
    file.write(f"Issue in file: {series_path}\n")


def write_num_dicom(num_dicom_files: int, file: Path) -> None:
    """
    Write the number of DICOM files to the output file.

    :param file: The file to write to.
    :param num_dicom_files: The number of DICOM files.
    """
    file.write(f"Number of DICOM files: {num_dicom_files}\n")


def write_num_nifti(num_nifti_files: int, file: Path) -> None:
    """
    Write the number of NIFTI files to the output file.

    :param file: The file to write to.
    :param num_nifti_files: The number of NIFTI files.
    """
    file.write(f"Number of NIFTI files: {num_nifti_files}\n")


def write_num_nifti_slices(num_nifti_slices: int, file: Path) -> None:
    """
    Write the number of NIFTI slices to the output file.

    :param file: The file to write to.
    :param num_nifti_slices: The number of NIFTI slices.
    """
    file.write(f"Number of NIFTI slices: {num_nifti_slices}\n")


def write_num_sequences(num_sequences: int, file: Path) -> None:
    """
    Write the number of different DICOM sequences to the output file.

    :param file: The file to write to.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write(f"Number of different DICOM sequences: {num_sequences}\n")


def write_series_info(
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: Path,
) -> None:
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


def write_empty_series(series_path: Path, file: Path) -> None:
    """
    Write a message about an empty Series folder to the output file.

    :param file: The file to write to.
    :param series_path: The path of the empty Series folder.
    """
    file.write(f"\nFound empty Series: {series_path}.\n")


def write_found_only_nifti(
    series_path: Path,
    num_nifti_files: int,
    num_nifti_slices: int,
    file: Path,
) -> None:
    """
    Write a message about a Series folder found with only NIFTI files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    """
    file.write("\nFound only NIFTI files\n")
    write_file_path(series_path, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)


def write_found_only_dicom(
    series_path: Path,
    num_dicom_files: int,
    num_sequences: int,
    file: Path,
) -> None:
    """
    Write a message about a Series folder found with only DICOM files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write("\nFound only DICOM files. Annotation is required.\n")
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_sequences(num_sequences, file)


def write_split_is_required(series_path: Path, num_sequences: int, file: Path) -> None:
    """
    Write a message about a Series folder requiring splitting to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_sequences: The number of sequences to split.
    """
    file.write(
        f"\nSeries folder: {series_path} needs to be split in {num_sequences} different Series folders.\n",
    )


def write_split_auto(file: Path) -> None:
    """
    Write a message that automatic splitting of a Series folder is possible to the output file.

    :param file: The file to write to.
    """
    file.write("Series folder can be split automatically.\n")


def write_split_manual(file: Path) -> None:
    """
    Write a message that manual splitting of a Series folder is required to the output file.

    :param file: The file to write to.
    """
    file.write(
        "Series folder cannot be split automatically. Manual split should be performed by a medical expert.\n",
    )


def write_dicom_nifti_slices_incompatible(
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: Path,
) -> None:
    """
    Write a message about incompatible number of DICOM files and NIFTI slices to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write("\nFound incompatible number of DICOM files and NIFTI slices.\n")
    file.write(
        "A medical expert should examine the specific Series folder and correct it.\n",
    )
    write_series_info(
        series_path,
        num_dicom_files,
        num_nifti_files,
        num_nifti_slices,
        num_sequences,
        file,
    )


def write_sequences_nifti_files_incompatible(
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: Path,
) -> None:
    """
    Write a message about incompatible number of DICOM sequences and NIFTI files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write("\nFound incompatible number of DICOM sequences and NIFTI files.\n")
    file.write(
        "A medical expert should examine the specific Series folder and correct it.\n",
    )
    write_series_info(
        series_path,
        num_dicom_files,
        num_nifti_files,
        num_nifti_slices,
        num_sequences,
        file,
    )


def write_more_nifti_than_dicom(
    series_path: int,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: Path,
) -> None:
    """
    Write a message about finding more NIFTI files than DICOM files to the output file.

    :param file: The file to write to.
    :param series_path: The path of the Series folder.
    :param num_dicom_files: The number of DICOM files.
    :param num_nifti_files: The number of NIFTI files.
    :param num_nifti_slices: The number of NIFTI slices.
    :param num_sequences: The number of DICOM sequences.
    """
    file.write("\nFound more NIFTI files than DICOM files.\n")
    write_series_info(
        series_path,
        num_dicom_files,
        num_nifti_files,
        num_nifti_slices,
        num_sequences,
        file,
    )


def write_issues_to_report(
    serie: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    dicom_sequences: list,
    file: Path,
) -> None:
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
            if Path.is_dir(serie) and not os.listdir(serie):
                shutil.rmtree(serie)
        else:
            # Issue: Only NIFTI files found
            write_found_only_nifti(serie, num_nifti_files, num_nifti_slices, file)
    elif num_nifti_files == 0:
        if num_sequences > 1:
            # Issue: Series with multiple sequences requires splitting
            split_series(serie, dicom_sequences)
    elif num_dicom_files < num_nifti_files:
        # Issue: More NIFTI files than DICOM files
        write_more_nifti_than_dicom(
            serie,
            num_dicom_files,
            num_nifti_files,
            num_nifti_slices,
            num_sequences,
            file,
        )
    elif num_dicom_files > num_nifti_files:
        if num_sequences == num_nifti_files:
            if num_dicom_files == num_nifti_slices:
                if num_sequences > 1:
                    # Issue: Series with multiple sequences requires splitting
                    write_split_is_required(serie, num_sequences, file)
                    write_split_manual(file)
                    write_series_info(
                        serie,
                        num_dicom_files,
                        num_nifti_files,
                        num_nifti_slices,
                        num_sequences,
                        file,
                    )
            else:
                # Issue: Incompatible number of DICOM files and NIFTI slices
                write_dicom_nifti_slices_incompatible(
                    serie,
                    num_dicom_files,
                    num_nifti_files,
                    num_nifti_slices,
                    num_sequences,
                    file,
                )
        else:
            # Issue: Incompatible number of DICOM sequences and NIFTI files
            write_sequences_nifti_files_incompatible(
                serie,
                num_dicom_files,
                num_nifti_files,
                num_nifti_slices,
                num_sequences,
                file,
            )
    elif num_sequences == num_nifti_files:
        if num_dicom_files == num_nifti_slices:
            if num_sequences > 1:
                # Issue: Series with multiple sequences requires splitting
                split_series(serie, dicom_sequences)
        else:
            # Issue: Incompatible number of DICOM files and NIFTI slices
            write_dicom_nifti_slices_incompatible(
                serie,
                num_dicom_files,
                num_nifti_files,
                num_nifti_slices,
                num_sequences,
                file,
            )
    elif num_sequences >= 1:
        # Issue: Incompatible number of DICOM sequences and NIFTI files
        write_sequences_nifti_files_incompatible(
            serie,
            num_dicom_files,
            num_nifti_files,
            num_nifti_slices,
            num_sequences,
            file,
        )


def get_series_attributes(serie: Path, file: Path) -> tuple:
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
        file_path = Path(serie, filename)

        if file_path.suffix == ".dcm":  # Check if DICOM
            # Increase the total DICOM files
            num_dicom_files += 1

            # Read DICOM file. It might be corrupted.
            try:
                img = dicom.dcmread(file_path)
            except TypeError:
                file.write(
                    f"\nDICOM file {file_path} is missing header information and cannot be read.\n",
                )
                continue

            # Get the modality attribute DICOM tag
            image_modality = img.Modality

            if (
                image_modality == "MR"
            ):  # Check if MR, then get the SequenceName DICOM tag
                try:
                    tag = img.SequenceName
                    dicom_sequences.append(tag)
                except AttributeError:
                    continue
            elif (
                image_modality == "MG"
            ):  # Check if MR, then get the SeriesTime DICOM tag
                try:
                    tag = img.SeriesTime
                    dicom_sequences.append(tag)
                except AttributeError:
                    continue
        elif (
            file_path.suffix in (".nii", ".nii.gz")
        ):  # Check if NIFTI
            # Increase the total NIFTI files
            num_nifti_files += 1

            # Read NIFTI file. It might be corrupted.
            try:
                img = nib.load(file_path)
            except TypeError:
                file.write(
                    f"\nNIFTI file {file_path} is missing tags and cannot be read.\n",
                )
                continue

            header = img.header
            nifti_files_shape = header.get_data_shape()

            # Append NIFTI slices. NIFTI might have invalid dimensions
            try:
                nifti_slices.append(nifti_files_shape[2])
            except IndexError:
                file.write(f"\nNIFTI File: {file_path} has no valid dimensions.\n")
        else:  # In case neither DICOM nor NIFTI file.
            continue

    # Sum total NIFTI slices
    num_nifti_slices = sum(nifti_slices)

    # Get the number of the total DICOM sequences
    if len(dicom_sequences) > 0:
        num_sequences = len(list(set(dicom_sequences)))
    else:
        num_sequences = num_nifti_files

    return (
        num_dicom_files,
        num_nifti_files,
        num_nifti_slices,
        num_sequences,
        dicom_sequences,
    )


def find_average_resolution() -> None:
    # Define a dictionary to store histograms for each image modality
    resolutions = {}

    # Provide a list of the cancer types to check
    cancer_types = ["breast", "colorectal", "lung", "prostate"]

    # Select the data provider
    dps = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    def find_modality_resolution(serie_path: Path) -> tuple:
        dicom_filenames = [
            Path(serie_path, filename)
            for filename in os.listdir(serie_path)
            if filename.suffix == ".dcm"
        ]
        if len(dicom_filenames):
            dicom_name = dicom_filenames[0]
            try:
                image = dicom.dcmread(
                    dicom_name,
                    specific_tags=[
                        (0x0008, 0x0060),
                        (0x0028, 0x0010),
                        (0x0028, 0x0011),
                    ],
                )
                modality = image.Modality
                resolution = (image.Rows, image.Columns)
            except TypeError:
                return None, None
        else:
            return None, None

        return modality, resolution

    # Function to update the resolutions
    def update_resolutions(resolution: tuple, image_modality: str) -> None:
        if resolution and image_modality:
            if image_modality not in resolutions:
                resolutions[image_modality] = []
            resolutions[image_modality].append(resolution)

    for cancer_type in cancer_types:
        data_providers = dps[cancer_type]
        for data_provider in data_providers:
            # Define the target directory to examine based on cancer type and data provider name
            database_path = r"prm/incisive2"
            working_path = f"{database_path}/{cancer_type}/{data_provider}/data"

            # Get all patients paths
            patients_paths = [
                Path(working_path, filename)
                for filename in os.listdir(working_path)
                if Path.is_dir(Path(working_path, filename))
            ]

            # Get all studies and check for NIFTI files outside of Series
            studies_paths = []
            for patient_path in patients_paths:
                studies = os.listdir(patient_path)
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)

            # Get all series paths
            series_paths = []
            for study_path in studies_paths:
                series = os.listdir(study_path)
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)

            for serie_path in tqdm(series_paths, total=len(series_paths)):
                if len(os.listdir(serie_path)) > 0:
                    modality, resolution = find_modality_resolution(serie_path)
                    update_resolutions(resolution, modality)

    # Calculate the average resolution for each image modality
    average_resolutions = {}
    for image_modality, resolutions_list in resolutions.items():
        total_rows = 0
        total_columns = 0
        for resolution in resolutions_list:
            total_rows += resolution[0]
            total_columns += resolution[1]
        average_rows = total_rows / len(resolutions_list)
        average_columns = total_columns / len(resolutions_list)
        average_resolutions[image_modality] = (average_rows, average_columns)


def find_devices() -> None:
    # Define a dictionary to store histograms for each image modality
    devices = {}

    # Provide a list of the cancer types to check
    cancer_types = ["breast", "colorectal", "lung", "prostate"]

    # Select the data provider
    dps = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    def find_manufacturer_model(serie_path: Path) -> tuple:
        dicom_filenames = [
            Path(serie_path, filename)
            for filename in os.listdir(serie_path)
            if filename.suffix == ".dcm"
        ]
        if len(dicom_filenames):
            dicom_name = dicom_filenames[0]
            try:
                image = dicom.dcmread(
                    dicom_name,
                    specific_tags=[(0x0008, 0x0070), (0x0008, 0x1090)],
                )
                manufacturer = image.Manufacturer
                model = image.ManufacturerModelName
            except TypeError:
                return None, None
        else:
            return None, None

        return manufacturer, model

    # Function to update the resolutions
    def update_devices(manufacturer: str, model: str) -> None:
        if manufacturer and model:
            device = f"{manufacturer} - {model}"
            if device not in devices:
                devices[device] = 1
            else:
                devices[device] += 1

    for cancer_type in cancer_types:
        data_providers = dps[cancer_type]
        for data_provider in data_providers:
            # Define the target directory to examine based on cancer type and data provider name
            database_path = r"prm/incisive2"
            working_path = f"{database_path}/{cancer_type}/{data_provider}/data"

            # Get all patients paths
            patients_paths = [
                Path(working_path, filename)
                for filename in os.listdir(working_path)
                if Path.is_dir(Path(working_path, filename))
            ]

            # Get all studies and check for NIFTI files outside of Series
            studies_paths = []
            for patient_path in patients_paths:
                studies = os.listdir(patient_path)
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)

            # Get all series paths
            series_paths = []
            for study_path in studies_paths:
                series = os.listdir(study_path)
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)

            for serie_path in tqdm(series_paths, total=len(series_paths)):
                for _ in range(len(os.listdir(serie_path))):
                    if len(os.listdir(serie_path)) > 0:
                        manufacturer, model = find_manufacturer_model(serie_path)
                        update_devices(manufacturer, model)

    sorted_devices = dict(
        sorted(devices.items(), key=lambda item: item[1], reverse=True),
    )

    # Remove the existing report file if it exists
    file_path = r"tmp/devices.txt"
    if Path.exists(file_path):
        with contextlib.suppress(OSError):
            Path.unlink(file_path)

    # Write the average resolutions
    with Path.open(file_path, "a") as f:
        for device, occ in sorted_devices.items():
            f.write(f"Manufacturer and Model: {device}, Occurences: {occ}\n")


def count_dicom_images() -> None:
    def count_dicom_files(folder_path: Path) -> int:
        dicom_count = 0

        for root, _, files in os.walk(folder_path):
            for filename in files:
                file_path = Path(root, filename)
                if file_path.suffix == ".dcm":
                    dicom_count += 1

        return dicom_count

    # Provide a list of the cancer types to check
    cancer_types = ["breast", "colorectal", "lung", "prostate"]

    # Select the data provider
    dps = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    for cancer_type in cancer_types:
        dicom_per_cancer_type = 0
        data_providers = dps[cancer_type]
        for data_provider in data_providers:
            # Define the target directory to examine based on cancer type and data provider name
            database_path = r"prm/incisive2"
            working_path = f"{database_path}/{cancer_type}/{data_provider}/data"

            # Get all patients paths
            patients_paths = [
                Path(working_path, filename)
                for filename in os.listdir(working_path)
                if Path.is_dir(Path(working_path, filename))
            ]

            # Get all studies and check for NIFTI files outside of Series
            studies_paths = []
            for patient_path in patients_paths:
                studies = os.listdir(patient_path)
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)

            # Get all series paths
            series_paths = []
            for study_path in studies_paths:
                series = os.listdir(study_path)
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)

            total_dicom = 0
            for serie_path in tqdm(series_paths, total=len(series_paths)):
                serie_dicom_count = count_dicom_files(serie_path)
                total_dicom += serie_dicom_count

            dicom_per_cancer_type += total_dicom


def find_number_images_per_modality() -> None:
    # Define a dictionary to store histograms for each image modality
    images_per_modality = {}

    # Provide a list of the cancer types to check
    cancer_types = ["breast", "colorectal", "lung", "prostate"]

    # Select the data provider
    dps = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    def find_modality(serie_path: Path) -> str:
        dicom_filenames = [
            Path(serie_path, filename)
            for filename in os.listdir(serie_path)
            if filename.suffix == ".dcm"
        ]
        if len(dicom_filenames):
            dicom_name = dicom_filenames[0]
            try:
                image = dicom.dcmread(dicom_name, specific_tags=[(0x0008, 0x0060)])
                modality = image.Modality
            except TypeError:
                return None
        else:
            return None

        return modality

    # Function to update the resolutions
    def update_images_per_modality(modality: str) -> None:
        if modality:
            if modality not in images_per_modality:
                images_per_modality[modality] = 1
            else:
                images_per_modality[modality] += 1

    for cancer_type in cancer_types:
        data_providers = dps[cancer_type]
        for data_provider in data_providers:
            # Define the target directory to examine based on cancer type and data provider name
            database_path = r"prm/incisive2"
            working_path = f"{database_path}/{cancer_type}/{data_provider}/data"

            # Get all patients paths
            patients_paths = [
                Path(working_path, filename)
                for filename in os.listdir(working_path)
                if Path.is_dir(Path(working_path, filename))
            ]

            # Get all studies and check for NIFTI files outside of Series
            studies_paths = []
            for patient_path in patients_paths:
                studies = os.listdir(patient_path)
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)

            # Get all series paths
            series_paths = []
            for study_path in studies_paths:
                series = os.listdir(study_path)
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)

            for serie_path in tqdm(series_paths, total=len(series_paths)):
                if len(os.listdir(serie_path)) > 0:
                    for _ in range(len(os.listdir(serie_path))):
                        modality = find_modality(serie_path)
                        update_images_per_modality(modality)

    sorted_images_per_modality = dict(
        sorted(images_per_modality.items(), key=lambda item: item[1], reverse=True),
    )

    # Remove the existing report file if it exists
    file_path = r"tmp/images_per_modality.txt"
    if Path.exists(file_path):
        with contextlib.suppress(OSError):
            Path.unlink(file_path)

    # Write the average resolutions
    with Path.open(file_path, "a") as f:
        for modality, occ in sorted_images_per_modality.items():
            f.write(f"For Modality: {modality} - Number of Images: {occ}\n")


def create_directory_structure(base_dir: Path) -> None:
    """
    Create directory structure as specified.

    Args:
        base_dir (str): Base directory where the structure will be created.

    Returns:
        None
    """
    # Check if 'reports' directory exists
    reports_dir = Path(base_dir, "reports")
    if Path.exists(reports_dir):
        shutil.rmtree(reports_dir)

    # Create 'reports' directory
    Path.mkdir(reports_dir, parents=True)

    # Create data provider directories inside 'reports'
    for dp in ["dp1", "dp2"]:
        dp_dir = Path(reports_dir, dp)
        Path.mkdir(dp_dir, parents=True)

        # Create subdirectories inside each data provider directory
        for sub_dir in ["breast", "colorectal", "lung", "prostate"]:
            sub_dir_path = Path(dp_dir, sub_dir)
            Path.mkdir(sub_dir_path, parents=True)


def correct_dicom_directory_structure(database_path: Path = r"prm/incisive2") -> None:
    # Provide a list of the cancer types to check
    cancer_types = ["breast", "colorectal", "lung", "prostate"]

    # Select the data provider
    data_providers = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    create_directory_structure("tmp")

    for cancer_type in cancer_types:
        for data_provider in data_providers[cancer_type]:
            # Define the target directory to examine based on cancer type and data provider name

            working_path = Path(rf"{database_path}/{cancer_type}/{data_provider}/data")

            # Define the report file path
            txt_file_path = Path(
                rf"tmp/reports/{data_provider}/{cancer_type}/report_issues.txt",
            )
            print(f"Generating issues .txt file in {txt_file_path}")

            # Remove the existing report file if it exists
            if Path.exists(txt_file_path):
                with contextlib.suppress(OSError):
                    Path.unlink(txt_file_path)

            with Path.open(txt_file_path, "a") as file:
                # Get all patients paths
                patients_paths = [
                    Path(working_path, filename)
                    for filename in os.listdir(working_path)
                    if Path.is_dir(Path(working_path, filename))
                ]

                # Get all studies and check for NIFTI files outside of Series
                studies_paths = []
                for patient_path in patients_paths:
                    studies = os.listdir(patient_path)
                    for study in studies:
                        study_path = Path(patient_path, study)
                        if (
                            study_path.suffix in (".nii", ".nii.gz")
                        ):
                            # Check if NIFTI is outside of Series -> Issue
                            write_nifti_outside_series(study_path, file)
                        elif "Study" in str(study_path):
                            write_patient_not_pass_qc(study_path, file)
                        elif Path.is_dir(study_path):
                            studies_paths.append(study_path)

                # Get all series paths
                series_paths = []
                for study_path in studies_paths:
                    series = os.listdir(study_path)
                    for serie in series:
                        serie_path = Path(study_path, serie)
                        if (
                            study_path.suffix in (".nii", ".nii.gz")
                        ):
                            # Check if NIFTI outside of Series -> Issue
                            write_nifti_outside_series(serie_path, file)
                        elif Path.is_dir(serie_path):
                            series_paths.append(serie_path)

                # For each series, get the attributes and write existing issues
                for series_path in tqdm(series_paths, total=len(series_paths)):
                    # Get attributes for the series
                    (
                        num_dicom_files,
                        num_nifti_files,
                        num_nifti_slices,
                        num_sequences,
                        dicom_sequences,
                    ) = get_series_attributes(series_path, file)

                    # Write the issues regarding the contents of the series to the .txt file
                    write_issues_to_report(
                        series_path,
                        num_dicom_files,
                        num_nifti_files,
                        num_nifti_slices,
                        num_sequences,
                        dicom_sequences,
                        file,
                    )


def main_cli() -> None:
    import fire

    fire.Fire(correct_dicom_directory_structure)


if __name__ == "__main__":
    database_path = Path("prm/incisive2")
    correct_dicom_directory_structure(database_path)
