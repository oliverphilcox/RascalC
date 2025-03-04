## Script to collect the raw covariance matrices from the output directory of the C++ code

import numpy as np
import os
from glob import glob
from shutil import copy2, copytree
from warnings import warn
from collections.abc import Iterable
from .utils import rmdir_if_exists_and_empty


def convert_suffix(suffix: str) -> int | str:
    if suffix.isdigit(): return int(suffix)
    if suffix != "full": raise ValueError("Unknown suffix")
    return suffix


def organize_filename(filename: str, output_groups: dict, jack: bool = False) -> None:
    "Interpret the filename as saved by the C++ code and put it into output_groups array. Jack is for jackknife"
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


def save_safe(output_dir: str, output_group_name: str, output_dictionary: dict[str]):
    "Save the dictionary of Numpy arrays into directory avoiding name clashes"
    output_filename = os.path.join(output_dir, f"Raw_Covariance_Matrices_{output_group_name}.npz")
    if os.path.exists(output_filename):
        warn(f"The default output filename for group {output_group_name}, {output_filename}, already exists. Will try to find a replacement.")
        i = 1
        while True:
            output_filename = os.path.join(output_dir, f"Raw_Covariance_Matrices_{output_group_name}.{i}.npz")
            if not os.path.exists(output_filename):
                warn(f"Found unused name {output_filename}, will save there.")
                break
            i += 1
    
    os.makedirs(os.path.dirname(output_filename), exist_ok = True) # make sure the directory exists
    np.savez_compressed(output_filename, **output_dictionary)


def collect_raw_covariance_matrices(cov_dir: str, cleanup: bool = True, print_function = print) -> dict[str, dict[str, np.ndarray[float]]]:
    """
    Collect the covariance matrices from text files written by the C++ code and organize them into a Numpy .npz file.
    With cleanup enabled (default), deletes the text files after collection.
    """

    cov_dir_all = os.path.join(cov_dir, 'CovMatricesAll/')
    cov_dir_jack = os.path.join(cov_dir, 'CovMatricesJack/')

    output_groups = {}

    # load the full matrices
    for input_filename in glob(cov_dir_all + "*.txt"):
        organize_filename(input_filename, output_groups, jack = False)

    # load the jack matrices if present
    for input_filename in glob(cov_dir_jack + "*.txt"):
        organize_filename(input_filename, output_groups, jack = True)
    
    print_function(f"Detected {len(output_groups)} output group(s) in {cov_dir}")

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
                matrix = np.loadtxt(input_filename)
                if matrix_name.startswith("c2") and matrix.ndim == 1: matrix = np.diag(matrix) # convert 1D c2 to a 2D diagonal matrix
                output_dictionary[matrix_name][suffix] = matrix

            # special treatment for string suffixes (at the moment, only "full")
            tmp_keys = list(output_dictionary[matrix_name].keys())
            for suffix in tmp_keys:
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
        
        save_safe(cov_dir, output_group_name, output_dictionary)

        # now that the file is saved (not any earlier to be sure), can remove all the text files
        # the list contains only the files that had their contents loaded and saved
        if cleanup:
            for matrix_filenames_dictionary in output_group.values():
                for input_filename in matrix_filenames_dictionary.values():
                    os.remove(input_filename)
        
        print_function(f"Finished with output group {output_group_name}")

    # remove subdirectories too if they are empty
    if cleanup:
        rmdir_if_exists_and_empty(cov_dir_all)
        rmdir_if_exists_and_empty(cov_dir_jack)

    return return_dictionary


def load_raw_covariances(file_root: str, label: str, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function = print) -> dict[str]:
    """
    Load the raw covariance matrices as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.

    label specifies the number of radial bins and either the number of angular (mu) bins (nN_mM) or the maximum (even) multipole index (nN_lL).

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.
    """
    input_filename = os.path.join(file_root, f"Raw_Covariance_Matrices_{label}.npz")
    if os.path.isfile(input_filename): raw_cov = dict(np.load(input_filename))
    else:
        print_function(f"Collecting the raw covariance matrices from {file_root}")
        result = collect_raw_covariance_matrices(file_root, print_function = print_function)
        if label not in result:
            raise ValueError(f"Raw covariance matrices for {label} not produced. Check the n and m/max_l values.")
        raw_cov = result[label]
    if n_samples is None: return raw_cov # return the full set
    elif isinstance(n_samples, int):
        if n_samples <= 0: raise ValueError("Number of samples must be positive if integer")
        n_samples = np.arange(n_samples)
    elif isinstance(n_samples, Iterable):
        if all(isinstance(_, int) for _ in n_samples): n_samples = np.array(n_samples, dtype = int)
        elif all(isinstance(_, bool) for _ in n_samples): n_samples = np.array(n_samples, dtype = bool)
        else: raise TypeError("n_samples elements must be all either int (indices) or bool (mask)")
    else: raise TypeError("n_samples must be None, positive int or iterable of int or bool")
    # select the given samples and update the averages
    keys = [key for key in raw_cov.keys() if not key.endswith("_full")]
    for key in keys:
        raw_cov[key] = raw_cov[key][n_samples]
        raw_cov[key + "_full"] = np.mean(raw_cov[key], axis = 0)
    return raw_cov


def load_raw_covariances_smu(file_root: str, n: int, m: int, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function = print) -> dict[str]:
    """
    Load the raw covariance matrices from the s_mu mode as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.
    n and m are numbers of radial and angular (mu) bins respectively.

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.
    """
    label = f"n{n}_m{m}"
    return load_raw_covariances(file_root, label, n_samples, print_function)


def load_raw_covariances_legendre(file_root: str, n: int, max_l: int, n_samples: None | int | Iterable[int] | Iterable[bool] = None, print_function = print) -> dict[str]:
    """
    Load the raw covariance matrices from the Legendre modes as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.
    n and m are numbers of radial and angular (mu) bins respectively.

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.
    """
    label = f"n{n}_l{max_l}"
    return load_raw_covariances(file_root, label, n_samples, print_function)


def cat_raw_covariance_matrices(n: int, mstr: str, input_roots: list[str], ns_samples: list[None | int | Iterable[int] | Iterable[bool]], output_root: str, collapse_factor: int = 1, print_function = print) -> dict[str]:
    """
    Catenate the raw covariance matrices from one or more input directories, assuming they are from similar runs: the same number of input random points, ``N2``, ``N3``, ``N4`` and ``loops_per_sample``.

    Parameters
    ----------
    n : integer
        The number of radial bins.

    mstr : string
        The next part of the output configuration label: ``m{number of angular bins}`` in s_mu mode and ``l{max_multipole}`` in Legendre modes.

    input_roots : list of strings
        Input directories to catenate from.

    ns_samples : list of the same length as ``input_roots``
        Each element allows to select subsamples flexibly from the corresponding input directory:

            - if None, returns all samples;
            - if a positive integer, returns as many samples from the beginning;
            - if a sequence of integers, returns subsamples with indices from this sequence;
            - sequence of boolean values is interpreted as a boolean mask for subsamples.

    output_root : string
        The output directory, in which the NumPy file with the resulting raw covariance matrix arrays is created.

    collapse_factor : positive integer
        (Optional) allows to reduce the number of subsamples by averaging over a given number of adjacent subsamples.
        Default value is 1, meaning no reduction.
    
    print_function : Callable
        (Optional) custom function to use for printing. Default is ``print``.

    Returns
    -------
    catenation_results : dict[str, np.ndarray[float]]
        Catenated raw covariance matrices as a dictionary with string keys and Numpy array values. All this information is also saved in a ``Raw_Covariance_Matrices*.npz`` file in the ``output_root``.
    """
    if collapse_factor <= 0: raise ValueError("Collapsing factor must be positive")
    if len(input_roots) < 1: raise ValueError("Need at least one input directory")
    if len(ns_samples) != len(input_roots): raise ValueError("Number of input dirs and subsamples to use from them must be the same")

    label = f"n{n}_{mstr}"
    result = {}
    for index, (input_root, n_samples) in enumerate(zip(input_roots, ns_samples)):
        input_file = load_raw_covariances(input_root, label, n_samples, print_function)
        # ignore full arrays
        input_file = {key: value for (key, value) in input_file.items() if not key.endswith("_full")}
        # check that the keys are the same, unless the result is brand new
        if len(result) > 0:
            result_keys = set(result.keys())
            input_keys = set(input_file.keys())
            if result_keys != input_keys:
                warn("Different sets of matrices present among the input files, will only use the overlapping ones.")
                common_keys = result_keys & input_keys
                result = {key: result[key] for key in common_keys}
                input_file = {key: input_file[key] for key in common_keys}
        # finally, loop over all the arrays
        for matrix_name, matrices in input_file.items():
            if matrix_name.endswith("_full"): continue # ignore full arrays
            if index != 0: result[matrix_name] = np.append(result[matrix_name], matrices, axis = 0)
            else: result[matrix_name] = matrices
    
    # loop over all the matrix names
    for matrix_name in list(result.keys()): # the dictionary will be changed
        if collapse_factor > 1:
            matrix_shape = result[matrix_name].shape
            result[matrix_name] = np.mean(result[matrix_name].reshape(matrix_shape[0] // collapse_factor, collapse_factor, *matrix_shape[1:]), axis = 1) # average over adjacent collapse_factor samples
        # make full arrays by averaging the subsamples
        matrix_name_full = matrix_name + "_full"
        result[matrix_name_full] = np.mean(result[matrix_name], axis = 0)
    
    save_safe(output_root, label, result)

    # copy other useful files from the first input root, unless identical with the output root
    # assuming they are the same among output roots; otherwise catenation should not be sensible
    if not os.path.samefile(input_roots[0], output_root):
        for pattern in ("weights", "xi*", "radial_binning*.csv"):
            for filename in glob(pattern, root_dir = input_roots[0]):
                src = os.path.join(input_roots[0], filename)
                (copytree if os.path.isdir(src) else copy2)(src, os.path.join(output_root, filename)) # need different functions for dirs and files

    return result