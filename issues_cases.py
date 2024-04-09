from messages import *
from split_series import split_series
import os
import shutil

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

