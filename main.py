import os
import pydicom as dicom
from tqdm import tqdm
from get_series_attributes import *
from messages import *
from issues_cases import *

# Provide a list of the cancer types to check
cancer_types = ["prostate"]

# Select the data provider
data_provider = "goc"

for cancer_type in cancer_types:

    # Define the target directory to examine based on cancer type and data provider name
    # database_path = "/mnt/nfs/incisive2"
    # working_path = f"{database_path}/{cancer_type}/{data_provider}/data"

    # Define a custom target directory to examine
    working_path = "/mnt/nfs/incisive3/goc/prostate"

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

    with open(file_path, "a") as file:
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

