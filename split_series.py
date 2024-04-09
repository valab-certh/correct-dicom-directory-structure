import os
import shutil
import pydicom as dicom

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
        except:
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

def main():

    current_series_path = "/mnt/nfs/incisive2/CERTH-curation/breast/uns/data/patient-001/study-001/Series-1"
    names_of_sequences = ['122947', '123239', '122641', '123112']

    split_series(current_series_path, names_of_sequences)

if __name__ == '__main__':
    main()