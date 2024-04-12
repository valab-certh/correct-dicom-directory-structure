from __future__ import annotations  # noqa: EXE002

import contextlib
import os
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

import nibabel
import pydicom as dicom
from tqdm import tqdm

if TYPE_CHECKING:
    import io


def delete_crc_files_in_series(series_path: Path) -> None:
    for file_name in list(series_path.iterdir()):
        if file_name.suffix == ".crc":
            file_path = Path(series_path, file_name)
            with contextlib.suppress(PermissionError):
                Path.unlink(file_path)


def delete_crc_files_in_study(study_path: Path) -> None:
    for series_folder in list(study_path.iterdir()):
        series_path = Path(study_path, series_folder)
        if Path.is_dir(series_path):
            delete_crc_files_in_series(series_path)


def delete_crc_files_in_patient(patient_path: Path) -> None:
    for study_folder in list(patient_path.iterdir()):
        study_path = Path(patient_path, study_folder)
        if Path.is_dir(study_path):
            delete_crc_files_in_study(study_path)


def delete_crc_files(root_directory: Path) -> None:
    for patient_folder in list(root_directory.iterdir()):
        patient_path = Path(root_directory, patient_folder)
        if Path.is_dir(patient_path):
            delete_crc_files_in_patient(patient_path)


def replace_folders(destination_directory: Path, source_directory: Path) -> None:
    source_subdirectories = [
        Path(source_directory, d)
        for d in list(source_directory.iterdir())
        if Path.is_dir(Path(source_directory, d))
    ]
    for source_subdir in tqdm(source_subdirectories, total=len(source_subdirectories)):
        subdir_name = Path.name(source_subdir)
        destination_subdir = Path(destination_directory, subdir_name)
        if Path.exists(destination_subdir):
            shutil.rmtree(destination_subdir)
        shutil.copytree(source_subdir, destination_subdir)


def get_directory_filenames(directory: Path) -> list[Path]:
    return list(directory.iterdir())


def get_crc_filename(filename: Path) -> tuple[Path, str]:
    parent_dir = Path.resolve(Path(filename, os.pardir))
    basename = Path.name(filename)
    crc_basename = f".{basename}.crc"
    crc_path = Path(parent_dir, crc_basename)
    return (crc_path, crc_basename)


def check_if_crc_exists(filename: Path) -> bool:
    crc_path, _ = get_crc_filename(filename)
    return Path.exists(crc_path)


def get_filenames_for_each_sequence(  # noqa: C901, PLR0912
    current_series_path: Path,
    dicom_sequences: list[str],
) -> dict[str, list[Path]]:
    dicom_filenames = [
        Path(current_series_path, filename)
        for filename in list(current_series_path.iterdir())
        if filename.suffix == ".dcm"
    ]
    nifti_filenames = [
        Path(current_series_path, filename)
        for filename in list(current_series_path.iterdir())
        if filename.suffix in (".nii.gz", ".nii")
    ]
    filenames_for_sequence: dict[str, list[Path]] = {}
    for seq_name in dicom_sequences:
        filenames_for_sequence[seq_name] = []
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


def find_next_available_series_name(series_names_in_parent_dir: list[str]) -> str:
    counter = 1
    new_series_name = f"Series-{counter}"
    while new_series_name in series_names_in_parent_dir:
        counter += 1
        new_series_name = f"Series-{counter}"
    return new_series_name


def split_series(current_series_path: Path, names_of_sequences: list[str]) -> None:
    parent_dir_path = Path.resolve(Path(current_series_path, os.pardir))
    series_names_in_parent_dir = [
        Path.name(filename) for filename in list(parent_dir_path.iterdir())
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
            basename = Path.name(Path(image_path))
            move_path = Path(new_dir_path, basename)
            shutil.move(image_path, move_path)
            if check_if_crc_exists(Path(image_path)):
                crc_path, crc_basename = get_crc_filename(Path(image_path))
                move_crc_path = Path(new_dir_path, crc_basename)
                shutil.move(crc_path, move_crc_path)


def write_nifti_outside_series(path: Path, file: io.TextIOWrapper) -> None:
    file.write(f"\nFound NIFTI file: {path} outside of Series folder.\n")
    file.write(
        "The medical expert should examine the annotation file and place it in the corresponding Series folder.\n",  # noqa: E501
    )


def write_patient_not_pass_qc(path: Path, file: io.TextIOWrapper) -> None:
    file.write(f"\nFound Study folder: {path} that did not pass the Quality Check.\n")
    file.write(
        "The Data Provider should pass this patient throught the Quality Check.\n",
    )


def write_file_path(series_path: Path, file: io.TextIOWrapper) -> None:
    file.write(f"Issue in file: {series_path}\n")


def write_num_dicom(num_dicom_files: int, file: io.TextIOWrapper) -> None:
    file.write(f"Number of DICOM files: {num_dicom_files}\n")


def write_num_nifti(num_nifti_files: int, file: io.TextIOWrapper) -> None:
    file.write(f"Number of NIFTI files: {num_nifti_files}\n")


def write_num_nifti_slices(num_nifti_slices: int, file: io.TextIOWrapper) -> None:
    file.write(f"Number of NIFTI slices: {num_nifti_slices}\n")


def write_num_sequences(num_sequences: int, file: io.TextIOWrapper) -> None:
    file.write(f"Number of different DICOM sequences: {num_sequences}\n")


def write_series_info(  # noqa: PLR0913
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: io.TextIOWrapper,
) -> None:
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)
    write_num_sequences(num_sequences, file)


def write_empty_series(series_path: Path, file: io.TextIOWrapper) -> None:
    file.write(f"\nFound empty Series: {series_path}.\n")


def write_found_only_nifti(
    series_path: Path,
    num_nifti_files: int,
    num_nifti_slices: int,
    file: io.TextIOWrapper,
) -> None:
    file.write("\nFound only NIFTI files\n")
    write_file_path(series_path, file)
    write_num_nifti(num_nifti_files, file)
    write_num_nifti_slices(num_nifti_slices, file)


def write_found_only_dicom(
    series_path: Path,
    num_dicom_files: int,
    num_sequences: int,
    file: io.TextIOWrapper,
) -> None:
    file.write("\nFound only DICOM files. Annotation is required.\n")
    write_file_path(series_path, file)
    write_num_dicom(num_dicom_files, file)
    write_num_sequences(num_sequences, file)


def write_split_is_required(
    series_path: Path,
    num_sequences: int,
    file: io.TextIOWrapper,
) -> None:
    file.write(
        f"\nSeries folder: {series_path} needs to be split in {num_sequences} different Series folders.\n",  # noqa: E501
    )


def write_split_auto(file: io.TextIOWrapper) -> None:
    file.write("Series folder can be split automatically.\n")


def write_split_manual(file: io.TextIOWrapper) -> None:
    file.write(
        "Series folder cannot be split automatically. Manual split should be performed by a medical expert.\n",  # noqa: E501
    )


def write_dicom_nifti_slices_incompatible(  # noqa: PLR0913
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: io.TextIOWrapper,
) -> None:
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


def write_sequences_nifti_files_incompatible(  # noqa: PLR0913
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: io.TextIOWrapper,
) -> None:
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


def write_more_nifti_than_dicom(  # noqa: PLR0913
    series_path: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    file: io.TextIOWrapper,
) -> None:
    file.write("\nFound more NIFTI files than DICOM files.\n")
    write_series_info(
        series_path,
        num_dicom_files,
        num_nifti_files,
        num_nifti_slices,
        num_sequences,
        file,
    )


def write_issues_to_report(  # noqa: PLR0913, PLR0912, C901
    serie: Path,
    num_dicom_files: int,
    num_nifti_files: int,
    num_nifti_slices: int,
    num_sequences: int,
    dicom_sequences: list[str],
    file: io.TextIOWrapper,
) -> None:
    if num_dicom_files == 0:
        if num_nifti_files == 0:
            if Path.is_dir(serie) and (not list(serie.iterdir())):
                shutil.rmtree(serie)
        else:
            write_found_only_nifti(serie, num_nifti_files, num_nifti_slices, file)
    elif num_nifti_files == 0:
        if num_sequences > 1:
            split_series(serie, dicom_sequences)
    elif num_dicom_files < num_nifti_files:
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
                write_dicom_nifti_slices_incompatible(
                    serie,
                    num_dicom_files,
                    num_nifti_files,
                    num_nifti_slices,
                    num_sequences,
                    file,
                )
        else:
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
                split_series(serie, dicom_sequences)
        else:
            write_dicom_nifti_slices_incompatible(
                serie,
                num_dicom_files,
                num_nifti_files,
                num_nifti_slices,
                num_sequences,
                file,
            )
    elif num_sequences >= 1:
        write_sequences_nifti_files_incompatible(
            serie,
            num_dicom_files,
            num_nifti_files,
            num_nifti_slices,
            num_sequences,
            file,
        )


def get_series_attributes(  # noqa: PLR0912, C901
    serie: Path,
    file: io.TextIOWrapper,
) -> tuple[int, int, int, int, list[str]]:
    filenames = list(serie.iterdir())
    num_dicom_files = 0
    num_nifti_files = 0
    num_nifti_slices = 0
    num_sequences = 0
    nifti_slices = []
    dicom_sequences = []
    for filename in filenames:
        file_path = Path(serie, filename)
        if file_path.suffix == ".dcm":
            num_dicom_files += 1
            try:
                img = dicom.dcmread(file_path)
            except TypeError:
                file.write(
                    f"\nDICOM file {file_path} is missing header information and cannot be read.\n",  # noqa: E501
                )
                continue
            image_modality = img.Modality
            if image_modality == "MR":
                try:
                    tag = img.SequenceName
                    dicom_sequences.append(tag)
                except AttributeError:
                    continue
            elif image_modality == "MG":
                try:
                    tag = img.SeriesTime
                    dicom_sequences.append(tag)
                except AttributeError:
                    continue
        elif file_path.suffix in (".nii", ".nii.gz"):
            num_nifti_files += 1
            try:
                mask = nibabel.loadsave.load(str(file_path))
            except TypeError:
                file.write(
                    f"\nNIFTI file {file_path} is missing tags and cannot be read.\n",
                )
                continue
            try:
                nifti_slices.append(mask.shape[2])
            except IndexError:
                file.write(f"\nNIFTI File: {file_path} has no valid dimensions.\n")
        else:
            continue
    num_nifti_slices = sum(nifti_slices)
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


def find_devices() -> None:  # noqa: C901
    devices = {}
    cancer_types = ["breast", "colorectal", "lung", "prostate"]
    dps = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    def find_manufacturer_model(serie_path: Path) -> tuple[str, str]:
        dicom_filenames = [
            Path(serie_path, filename)
            for filename in list(serie_path.iterdir())
            if filename.suffix == ".dcm"
        ]
        if len(dicom_filenames):
            dicom_name = dicom_filenames[0]
            try:
                image = dicom.dcmread(dicom_name, specific_tags=[(8, 112), (8, 4240)])
                manufacturer = image.Manufacturer
                model = image.ManufacturerModelName
            except TypeError:
                return ("", "")
        else:
            return ("", "")
        return (manufacturer, model)

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
            database_path = "prm/incisive2"
            working_path = Path(f"{database_path}/{cancer_type}/{data_provider}/data")
            patients_paths = [
                Path(working_path, filename)
                for filename in list(working_path.iterdir())
                if Path.is_dir(Path(working_path, filename))
            ]
            studies_paths = []
            for patient_path in patients_paths:
                studies = list(patient_path.iterdir())
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)
            series_paths = []
            for study_path in studies_paths:
                series = list(study_path.iterdir())
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)
            for serie_path in tqdm(series_paths, total=len(series_paths)):
                for _ in range(len(list(serie_path.iterdir()))):
                    if len(list(serie_path.iterdir())) > 0:
                        manufacturer, model = find_manufacturer_model(serie_path)
                        update_devices(manufacturer, model)
    sorted_devices = dict(
        sorted(devices.items(), key=lambda item: item[1], reverse=True),
    )
    file_path = Path("tmp/devices.txt")
    if Path.exists(file_path):
        with contextlib.suppress(OSError):
            Path.unlink(file_path)
    with Path.open(file_path, "a") as f:
        for device, occ in sorted_devices.items():
            f.write(f"Manufacturer and Model: {device}, Occurences: {occ}\n")


def count_dicom_images() -> None:  # noqa: C901
    def count_dicom_files(folder_path: Path) -> int:
        dicom_count = 0
        for root, _, files in os.walk(folder_path):
            for filename in files:
                file_path = Path(root, filename)
                if file_path.suffix == ".dcm":
                    dicom_count += 1
        return dicom_count

    cancer_types = ["breast", "colorectal", "lung", "prostate"]
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
            database_path = "prm/incisive2"
            working_path = Path(f"{database_path}/{cancer_type}/{data_provider}/data")
            patients_paths = [
                Path(working_path, filename)
                for filename in list(working_path.iterdir())
                if Path.is_dir(Path(working_path, filename))
            ]
            studies_paths = []
            for patient_path in patients_paths:
                studies = list(patient_path.iterdir())
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)
            series_paths = []
            for study_path in studies_paths:
                series = list(study_path.iterdir())
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)
            total_dicom = 0
            for serie_path in tqdm(series_paths, total=len(series_paths)):
                serie_dicom_count = count_dicom_files(serie_path)
                total_dicom += serie_dicom_count
            dicom_per_cancer_type += total_dicom


def find_number_images_per_modality() -> None:  # noqa: C901
    images_per_modality = {}
    cancer_types = ["breast", "colorectal", "lung", "prostate"]
    dps = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }

    def find_modality(serie_path: Path) -> str:
        modality = ""
        dicom_filenames = [
            Path(serie_path, filename)
            for filename in list(serie_path.iterdir())
            if filename.suffix == ".dcm"
        ]
        if len(dicom_filenames):
            dicom_name = dicom_filenames[0]
            try:
                image = dicom.dcmread(dicom_name, specific_tags=[(8, 96)])
                modality = image.Modality
            except TypeError:
                return ""
        else:
            return ""
        return modality

    def update_images_per_modality(modality: str) -> None:
        if modality:
            if modality not in images_per_modality:
                images_per_modality[modality] = 1
            else:
                images_per_modality[modality] += 1

    for cancer_type in cancer_types:
        data_providers = dps[cancer_type]
        for data_provider in data_providers:
            database_path = "prm/incisive2"
            working_path = Path(f"{database_path}/{cancer_type}/{data_provider}/data")
            patients_paths = [
                Path(working_path, filename)
                for filename in list(working_path.iterdir())
                if Path.is_dir(Path(working_path, filename))
            ]
            studies_paths = []
            for patient_path in patients_paths:
                studies = list(patient_path.iterdir())
                for study in studies:
                    study_path = Path(patient_path, study)
                    if Path.is_dir(study_path):
                        studies_paths.append(study_path)
            series_paths = []
            for study_path in studies_paths:
                series = list(study_path.iterdir())
                for serie in series:
                    serie_path = Path(study_path, serie)
                    if Path.is_dir(serie_path):
                        series_paths.append(serie_path)
            for serie_path in tqdm(series_paths, total=len(series_paths)):
                if len(list(serie_path.iterdir())) > 0:
                    for _ in range(len(list(serie_path.iterdir()))):
                        modality = find_modality(serie_path)
                        update_images_per_modality(modality)
    sorted_images_per_modality = dict(
        sorted(images_per_modality.items(), key=lambda item: item[1], reverse=True),
    )
    file_path = Path("tmp/images_per_modality.txt")
    if Path.exists(file_path):
        with contextlib.suppress(OSError):
            Path.unlink(file_path)
    with Path.open(file_path, "a") as f:
        for modality, occ in sorted_images_per_modality.items():
            f.write(f"For Modality: {modality} - Number of Images: {occ}\n")


def create_directory_structure(base_dir: Path) -> None:
    reports_dir = Path(base_dir, "reports")
    if Path.exists(reports_dir):
        shutil.rmtree(reports_dir)
    Path.mkdir(reports_dir, parents=True)
    for dp in ["dp1", "dp2"]:
        dp_dir = Path(reports_dir, dp)
        Path.mkdir(dp_dir, parents=True)
        for sub_dir in ["breast", "colorectal", "lung", "prostate"]:
            sub_dir_path = Path(dp_dir, sub_dir)
            Path.mkdir(sub_dir_path, parents=True)


def correct_dicom_directory_structure(  # noqa: C901
    database_path: Path = Path("prm/incisive2"),
) -> None:
    cancer_types = ["breast", "colorectal", "lung", "prostate"]
    data_providers = {
        "breast": ["dp1", "dp2"],
        "colorectal": ["dp1", "dp2"],
        "lung": ["dp1", "dp2"],
        "prostate": ["dp1", "dp2"],
    }
    create_directory_structure(Path("tmp"))
    for cancer_type in cancer_types:
        for data_provider in data_providers[cancer_type]:
            working_path = Path(f"{database_path}/{cancer_type}/{data_provider}/data")
            txt_file_path = Path(
                f"tmp/reports/{data_provider}/{cancer_type}/report_issues.txt",
            )
            if Path.exists(txt_file_path):
                with contextlib.suppress(OSError):
                    Path.unlink(txt_file_path)
            with Path.open(txt_file_path, "a") as file:
                patients_paths = [
                    Path(working_path, filename)
                    for filename in list(working_path.iterdir())
                    if Path.is_dir(Path(working_path, filename))
                ]
                studies_paths = []
                for patient_path in patients_paths:
                    studies = list(patient_path.iterdir())
                    for study in studies:
                        study_path = Path(patient_path, study)
                        if study_path.suffix in (".nii", ".nii.gz"):
                            write_nifti_outside_series(study_path, file)
                        elif "Study" in str(study_path):
                            write_patient_not_pass_qc(study_path, file)
                        elif Path.is_dir(study_path):
                            studies_paths.append(study_path)
                series_paths = []
                for study_path in studies_paths:
                    series = list(study_path.iterdir())
                    for serie in series:
                        serie_path = Path(study_path, serie)
                        if study_path.suffix in (".nii", ".nii.gz"):
                            write_nifti_outside_series(serie_path, file)
                        elif Path.is_dir(serie_path):
                            series_paths.append(serie_path)
                for series_path in tqdm(series_paths, total=len(series_paths)):
                    (
                        num_dicom_files,
                        num_nifti_files,
                        num_nifti_slices,
                        num_sequences,
                        dicom_sequences,
                    ) = get_series_attributes(series_path, file)
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
