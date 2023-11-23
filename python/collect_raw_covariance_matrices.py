## Script to collect the raw covariance matrices from the output directory of the C++ code

import numpy as np
import sys,os
from glob import glob
from warnings import warn


def convert_suffix(suffix: str) -> int | str:
    if suffix.isdigit(): return int(suffix)
    if suffix != "full": raise ValueError("Unknown suffix")
    return suffix

def organize_filename(filename: str, output_groups: dict, jack: bool = False) -> None:
    # interpret the filename
    filename_no_ext = os.path.basename(filename) # remove the directory name
    filename_no_ext = ".".join(filename_no_ext.split(".")[:-1]) # remove the extension
    filename_parts = filename_no_ext.split("_") # split the rest by underscore
    term_name = filename_parts[0] # cN, RRN or EEN
    if term_name not in ("c2", "c3", "c4", "RR1", "RR2", "EE1", "EE2"): return # do not process other arrays
    if jack and term_name.startswith("c"): term_name += "j" # add "j" to cN if doing jack
    output_group_name = "_".join(filename_parts[1:-2]) # nN and mM or lL joined back
    indices = filename_parts[-2] # tracer numbers
    try:
        suffix = convert_suffix(filename_parts[-1]) # "full" or subsample number
    except ValueError:
        raise ValueError(f"Unrecognized {suffix = } in {filename = }")

    if output_group_name not in output_groups: # create an empty sub-dictionary
        output_groups[output_group_name] = {}
    
    matrix_name = term_name + "_" + indices
    if matrix_name not in output_groups[output_group_name]: # create an empty sub-sub-dictionary
        output_groups[output_group_name][matrix_name] = {}
    
    output_groups[output_group_name][matrix_name][suffix] = filename

def rmdir_safe(dirname: str) -> None:
    # remove directory if it exists and not empty, otherwise leave it
    if os.path.isdir(dirname) and len(os.listdir(dirname)) > 0:
        os.rmdir(dirname)

def collect_raw_covariance_matrices(cov_dir: str, print_function = print) -> None:
    cov_dir_all = os.path.join(cov_dir, 'CovMatricesAll/')
    cov_dir_jack = os.path.join(cov_dir, 'CovMatricesJack/')

    output_groups = {}

    # load the full matrices
    for input_filename in glob(cov_dir_all + "*.txt"):
        organize_filename(input_filename, output_groups, jack = False)

    # load the jack matrices if present
    for input_filename in glob(cov_dir_jack + "*.txt"):
        organize_filename(input_filename, output_groups, jack = True)
    
    print_function(f"Detected {len(output_groups)} output groups in {cov_dir}")

    return_dictionary = {}
    
    for output_group_name, output_group in output_groups.items():
        print_function(f"Processing output group {output_group_name}")
        # check that the different matrices have the same number of subsamples
        subsample_numbers = []
        for matrix_filenames_dictionary in output_group.values():
            subsample_number = 0
            for suffix in matrix_filenames_dictionary.keys():
                if isinstance(suffix, int): subsample_number += 1
                # subsamples have integer suffixes
            subsample_numbers.append(subsample_number)
        subsample_number = min(subsample_numbers)

        if any(this_subsample_number != subsample_number for this_subsample_number in subsample_numbers):
            warn(f"Some matrices in output group {output_group_name} have different number of subsamples. Will cut to the smallest number of subsamples.")
            # now cut all to the minimal number of subsamples, using the lowest numbers present
            for matrix_filenames_dictionary in output_group.values():
                subsample_suffixes_increasing = sorted([suffix for suffix in matrix_filenames_dictionary.keys() if isinstance(suffix, int)])
                if len(subsample_suffixes_increasing) == 0: continue # some arrays will not have subsamples
                subsample_suffix_max = subsample_suffixes_increasing[subsample_number - 1]
                for suffix in matrix_filenames_dictionary.keys():
                    if isinstance(suffix, int) and suffix > subsample_suffix_max:
                        matrix_filenames_dictionary.pop(suffix)
        
        # now create and fill the dictionary to be saved in the numpy file
        output_dictionary = {}
        for matrix_name, matrix_filenames_dictionary in output_group.items():
            output_dictionary[matrix_name] = dict()
            for suffix, input_filename in matrix_filenames_dictionary.items():
                output_dictionary[matrix_name][suffix] = np.loadtxt(input_filename)

            # special treatment for string suffixes (at the moment, only "full")
            for suffix in output_dictionary[matrix_name].keys():
                if isinstance(suffix, str):
                    output_dictionary[matrix_name + "_" + suffix] = output_dictionary[matrix_name].pop(suffix)
                    # this creates a separate array to be saved

            # now all the remaining suffixes must be integers so can be sorted easily
            output_dictionary[matrix_name] = np.array([output_dictionary[matrix_name][i_subsample] for i_subsample in sorted(output_dictionary[matrix_name].keys())])
            # this transformed the dictionary to numpy array, ordered by increasing subsample index

            # calculate the full as the mean of the subsamples
            full_matrix_computed = np.mean(output_dictionary[matrix_name], axis = 0)
            full_matrix_name = matrix_name + "_full"
            if full_matrix_name in output_dictionary:
                if not np.allclose(output_dictionary[full_matrix_name], full_matrix_computed):
                    warn(f"For {matrix_name} matrix, the loaded full is different from the average of subsamples. The latter will be saved.")
                    matrix_filenames_dictionary.pop("full") # remove the filename since it will be technically unused
            output_dictionary[full_matrix_name] = full_matrix_computed
        
        return_dictionary[output_group_name] = output_dictionary
        
        output_filename = os.path.join(cov_dir, f"Raw_Covariance_Matrices_{output_group_name}.npz")
        if os.path.exists(output_filename):
            warn(f"The default output filename for group {output_group_name}, {output_filename}, already exists. Will try to find a replacement.")
            i = 1
            while True:
                output_filename = os.path.join(cov_dir, f"Raw_Covariance_Matrices_{output_group_name}.{i}.npz")
                if not os.path.exists(output_filename):
                    warn(f"Found unused name {output_filename}, will save there.")
                    break
                i += 1
                
        np.savez_compressed(output_filename, **output_dictionary)

        # now that the file is saved (not any earlier to be sure), can remove all the text files
        # the list contains only the files that had their contents loaded and saved
        for matrix_filenames_dictionary in output_group.values():
            for input_filename in matrix_filenames_dictionary.values():
                os.remove(input_filename)
        
        print_function(f"Finished with output group {output_group_name}")

    # remove subdirectories too if they are empty
    rmdir_safe(cov_dir_all)
    rmdir_safe(cov_dir_jack)

    return return_dictionary

def load_raw_covariances(file_root: str, label: str, print_function = print) -> dict[str]:
    input_filename = os.path.join(file_root, f"Raw_Covariance_Matrices_{label}.npz")
    if not os.path.isfile(input_filename):
        print_function(f"Collecting the raw covariance matrices from {file_root}.")
        result = collect_raw_covariance_matrices(file_root, print_function)
        if label not in result:
            raise ValueError(f"Raw covariance matrices for {label} not produced. Check the n and m/max_l values.")
        return result[label]
    return np.load(input_filename)

def load_raw_covariances_smu(file_root: str, n: int, m: int, print_function = print) -> dict[str]:
    label = f"n{n}_m{m}"
    return load_raw_covariances(file_root, label, print_function)

def load_raw_covariances_legendre(file_root: str, n: int, max_l: int, print_function = print) -> dict[str]:
    label = f"n{n}_l{max_l}"
    return load_raw_covariances(file_root, label, print_function)

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv) != 2: # if too few
        print("Usage: python collect_raw_covariance_matrices.py {COVARIANCE_DIR}")
        sys.exit(1)

    cov_dir = str(sys.argv[1])

    collect_raw_covariance_matrices(cov_dir)