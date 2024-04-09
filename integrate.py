import os
import shutil
from tqdm import tqdm

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
    
if __name__ == "__main__":
    destination_directory = "/mnt/nfs/incisive2/lung/visaris/data"
    source_directory = "/mnt/nfs/incisive3/visaris/lung"

    replace_folders(destination_directory, source_directory)


