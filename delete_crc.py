import os

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

if __name__ == "__main__":
    # Provide the root directory path where the 'data' folder is located
    root_directory_path = '/mnt/nfs/incisive2/breast/auth'

    # Call the function to delete .crc files in the directory structure
    delete_crc_files(os.path.join(root_directory_path, 'data'))