"Functions to collect the raw covariance matrices from the output directory of the C++ code into a numpy file, load them, and catenate raw covariance matrices from one directory or multiple directories into another directory."

import numpy as np
import numpy.typing as npt
import os
from glob import glob
from shutil import copy2, copytree
from typing import Callable, Iterable
from .utils import rmdir_if_exists_and_empty


def convert_suffix(suffix: str) -> int | str:
    if suffix.isdigit(): return int(suffix)
    if suffix != "full": raise ValueError("Unknown suffix")
    return suffix


def organize_filename(filename: str, output_groups: dict, jack: bool = False, threepcf: bool = False) -> None:
    "Interpret the filename as saved by the C++ code and put it into output_groups array. Jack is for jackknife"
    filename_no_ext = os.path.basename(filename) # remove the directory name
    filename_no_ext = ".".join(filename_no_ext.split(".")[:-1]) # remove the extension
    filename_parts = filename_no_ext.split("_") # split the rest by underscore
    term_name = filename_parts[0] # cN, RRN or EEN
    term_names = [f"c{i}" for i in range(3, 7)] if threepcf else ("c2", "c3", "c4", "RR1", "RR2", "EE1", "EE2")
    if term_name not in term_names: return # do not process other arrays
    if jack and term_name.startswith("c"): term_name += "j" # add "j" to cN if doing jack
    output_group_name = "_".join(filename_parts[1:-2]) # nN and mM or lL joined back
    indices = filename_parts[-2] # tracer numbers for 2PCF and integral index for 3PCF
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


def save_safe(output_dir: str, output_group_name: str, output_dictionary: dict[str], threepcf: bool = False, print_function: Callable[[str], None] = print) -> None:
    "Save the dictionary of Numpy arrays into directory avoiding name clashes"
    output_filename = os.path.join(output_dir, "Raw_Covariance_Matrices" + "_3PCF" * threepcf + f"_{output_group_name}.npz")
    if os.path.exists(output_filename):
        print_function(f"The default output filename for group {output_group_name}, {output_filename}, already exists. Will try to find a replacement.")
        i = 1
        while True:
            output_filename = os.path.join(output_dir, "Raw_Covariance_Matrices" + "_3PCF" * threepcf + f"_{output_group_name}.{i}.npz")
            if not os.path.exists(output_filename):
                print_function(f"Found unused name {output_filename}, will save there.")
                break
            i += 1
    
    os.makedirs(os.path.dirname(output_filename), exist_ok = True) # make sure the directory exists
    np.savez_compressed(output_filename, **output_dictionary)


def collect_raw_covariance_matrices(cov_dir: str, dry_run: bool = False, cleanup: bool = True, check_finished: bool = True, two_tracers: bool | None = None, print_function: Callable[[str], None] = print) -> list[dict[str, dict[str, npt.NDArray[np.float64]]]]:
    """
    Collect the covariance matrices from text files written by the C++ code and organize them into a Numpy .npz file.
    Returns a list of two dictionaries, for 2PCF and 3PCF separately, with the structure {output_group_name: {matrix_name: array}}, with {matrix_name: array} replicating the structure of the .npz file.
    With dry_run enabled, this function will only return the output group names via the dictionary, and not perform the actual collection.
    With cleanup enabled (default), deletes the text files after collection.
    With check_finished enabled (default), performs a heuristic check whether the run seems finished by looking for the presence of the full covariance matrices. If they are not found, the collection does not proceed and a warning is issued. This is to prevent disrupting ongoing runs by collecting and/or deleting the text files. If you want to check convergence of timed-out run(s) that did not produce the full matrices, you can disable this check, but be careful not to disrupt ongoing runs.
    With check_finished enabled, the code additionally checks presence of all types of matrices for two tracers. Pass two_tracers = True to enable this check explicitly, or two_tracers = False to disable it. By default (two_tracers = None), the code will try to guess whether the run is for two tracers by looking for the presence of the xi_22.dat file, which is only produced for two tracers by :func:`RascalC.run_cov`.
    """

    # 2PCF
    cov_dir_all = os.path.join(cov_dir, 'CovMatricesAll/')
    cov_dir_jack = os.path.join(cov_dir, 'CovMatricesJack/')

    output_groups = {}

    # load the full matrices
    for input_filename in glob(cov_dir_all + "*.txt"):
        organize_filename(input_filename, output_groups, jack=False, threepcf=False)

    # load the jack matrices if present
    for input_filename in glob(cov_dir_jack + "*.txt"):
        organize_filename(input_filename, output_groups, jack=True, threepcf=False)
    
    print_function(f"Detected {len(output_groups)} 2PCF output group(s) in {cov_dir}")

    # 3PCF - should be done separately, because the output group name might match 2PCF and some arrays would be overwritten
    cov_dir_3pcf_all = os.path.join(cov_dir, '3PCFCovMatricesAll/')

    output_groups_3pcf = {}

    # load the full 3PCF matrices
    for input_filename in glob(cov_dir_3pcf_all + "*.txt"):
        organize_filename(input_filename, output_groups_3pcf, jack=False, threepcf=True)
    
    print_function(f"Detected {len(output_groups_3pcf)} 3PCF output group(s) in {cov_dir}")

    if dry_run: return [{output_group_name: {} for output_group_name in output_groups}, {output_group_name: {} for output_group_name in output_groups_3pcf}] # return the group names via the list of dictionaries, skipping any actual collection

    return_dictionaries = [{}, {}] # for 2PCF and 3PCF separately
    
    for threepcf, these_output_groups in enumerate([output_groups, output_groups_3pcf]):
        for output_group_name, output_group in these_output_groups.items():
            print_function(f"Processing {2+threepcf}PCF output group {output_group_name}")

            # check if the run seems finished by looking for the presence of the full matrices, if check_finished is enabled
            if check_finished:
                unfinished_names = [matrix_name for matrix_name, matrix_filenames_dictionary in output_group.items() if "full" not in matrix_filenames_dictionary]
                if unfinished_names:
                    print_function(f"WARNING: {unfinished_names} full matrices not found for {2+threepcf}PCF output group {output_group_name}. This may indicate that the C++ code is still running. Will skip this output group to avoid disrupting an ongoing run by collecting and/or deleting the text files. If you are sure that the run(s) in this directory finished (e.g., it/they timed out without producing the full matrices), you can disable this check by setting check_finished to False, but be careful not to disrupt ongoing runs.")
                    continue
                # additional checks that all types of matrices are present
                if threepcf: # different names for 3PCF, and only single-tracer for now
                    expected_matrix_names = [f'c{n}_{i}' for n in range(3, 6) for i in range(2)] + ['c6_1'] # exclude c6_0 term, because it should be small but is also hard to compute (see Section 5.2.3 and Appendix A of https://arxiv.org/abs/1910.04764)
                else: # set 2PCF names, excluding jackknife
                    if two_tracers is None: two_tracers = os.path.isfile(os.path.join(cov_dir, f"xi/xi_22.dat")) # simple heuristic if not provided explicitly, copied from post_process_auto
                    expected_matrix_names = ['c2_' + i + j for i in "12" for j in "12"] + ['c3_' + i + ',' + index2 for i in "12" for index2 in ("11", "12", "22")] + ['c4_' + indices for indices in ("11,11", "11,22", "12,11", "12,12", "12,21", "21,22", "22,22")] if two_tracers else ['c2_11', 'c3_1,11', 'c4_11,11']
                not_found_names = [matrix_name for matrix_name in expected_matrix_names if matrix_name not in output_group]
                if not_found_names:
                    print_function(f"WARNING: {not_found_names} matrices not found in {2+threepcf}PCF output group {output_group_name}. This may indicate that the C++ code is still running. Will skip this output group to avoid disrupting an ongoing run by collecting and/or deleting the text files. If you are sure that the run(s) in this directory finished (e.g., it/they timed out without producing the full matrices), you can disable this check by setting check_finished to False, but be careful not to disrupt ongoing runs.")
                    continue

            # check that the different matrices have the same subsamples by index
            subsample_indices = [{suffix for suffix in matrix_filenames_dictionary.keys() if isinstance(suffix, int)} for matrix_filenames_dictionary in output_group.values()]
            common_subsample_indices = subsample_indices[0].intersection(*subsample_indices)

            if any(this_subsample_indices != common_subsample_indices for this_subsample_indices in subsample_indices):
                print_function(f"WARNING: Some matrices in output group {output_group_name} have different subsample indices. " + ("This may indicate that the C++ code is still running. Will skip this output group to avoid disrupting an ongoing run by collecting and/or deleting the text files. If you are sure that the run(s) in this directory finished (e.g., it/they timed out without producing the full matrices), you can disable this check by setting check_finished to False, but be careful not to disrupt ongoing runs." if check_finished else f"Will cut to the common subsample indices, leaving {len(common_subsample_indices)} subsamples."))
                if check_finished: continue # skip the output group to avoid disrupting an ongoing run
                # otherwise, cut all to the common subsample indices
                output_group = {matrix_name: {suffix: input_filename for suffix, input_filename in matrix_filenames_dictionary.items() if suffix in common_subsample_indices or isinstance(suffix, str)} for matrix_name, matrix_filenames_dictionary in output_group.items()}
            
            # now create and fill the dictionary to be saved in the numpy file
            output_dictionary = {}
            for matrix_name, matrix_filenames_dictionary in output_group.items():
                output_dictionary[matrix_name] = {}
                for suffix, input_filename in matrix_filenames_dictionary.items():
                    matrix = np.loadtxt(input_filename)
                    if matrix_name.startswith("c2") and matrix.ndim == 1: matrix = np.diag(matrix) # convert 1D c2 to a 2D diagonal matrix
                    if isinstance(suffix, str): # special treatment for string suffixes (at the moment, only "full")
                        output_dictionary[matrix_name + "_" + suffix] = matrix # this creates a separate array to be saved, which seems more convenient
                    elif isinstance(suffix, int): # integer suffixes are the subsample indices
                        output_dictionary[matrix_name][suffix] = matrix # save the subsample matrices in a sub-dictionary, which will be converted to a numpy array later
                    else: raise TypeError(f"Unexpected suffix type {type(suffix)} for {matrix_name} in {output_group_name}; expected int or str")

                # now all the remaining suffixes must be integers so can be sorted easily
                output_dictionary[matrix_name] = np.array([output_dictionary[matrix_name][i_subsample] for i_subsample in sorted(output_dictionary[matrix_name].keys())])
                # this transformed the dictionary to numpy array, ordered by increasing subsample index

                # calculate the full as the mean of the subsamples
                full_matrix_computed = np.mean(output_dictionary[matrix_name], axis = 0)
                full_matrix_name = matrix_name + "_full"
                if full_matrix_name in output_dictionary:
                    if not np.allclose(output_dictionary[full_matrix_name], full_matrix_computed):
                        print_function(f"For {matrix_name} matrix from {2+threepcf}PCF ouput group {output_group_name}, the loaded full is different from the average of subsamples. The latter will be saved.")
                        matrix_filenames_dictionary.pop("full") # remove the filename since it will be technically unused
                output_dictionary[full_matrix_name] = full_matrix_computed
            
            return_dictionaries[threepcf][output_group_name] = output_dictionary
            
            save_safe(cov_dir, output_group_name, output_dictionary, threepcf=threepcf, print_function=print_function)

            # now that the file is saved (not any earlier to be sure), can remove all the text files
            # the list contains only the files that had their contents loaded and saved
            if cleanup:
                for matrix_filenames_dictionary in output_group.values():
                    for input_filename in matrix_filenames_dictionary.values():
                        os.remove(input_filename)
            
            print_function(f"Finished with {2+threepcf}PCF output group {output_group_name}")

    # remove subdirectories too if they are empty
    if cleanup:
        if len(output_groups) > 0: # also don't want to delete 2PCF directories if no 2PCF output groups were found. This can happen when the run is ongoing at an early stage and the text files are not yet produced, removing a directory will cause an error when the C++ code eventually tries to write them
            rmdir_if_exists_and_empty(cov_dir_all)
            rmdir_if_exists_and_empty(cov_dir_jack)
        if len(output_groups_3pcf) > 0: # similarly, don't want to delete the 3PCF directory if no 3PCF output groups were found
            rmdir_if_exists_and_empty(cov_dir_3pcf_all)

    return return_dictionaries


def load_raw_covariances(file_root: str, label: str, threepcf: bool = False, n_samples: None | int | Iterable[int] | Iterable[bool] = None, check_finished: bool = True, two_tracers: bool | None = None, print_function: Callable[[str], None] = print) -> dict[str]:
    """
    Load the raw covariance matrices as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.

    label specifies the number of radial bins and either the number of angular (mu) bins (nN_mM) or the maximum (even) multipole index (nN_lL).

    threepcf chooses between 2PCF (False, default) and 3PCF (True).

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.

    With check_finished enabled (default), before collecting raw covariance matrices (if needed), performs a heuristic check whether the run seems finished by looking for the presence of the full covariance matrices. If they are not found, the collection does not proceed and a warning is issued. This is to prevent disrupting ongoing runs by collecting and/or deleting the text files. If you want to check convergence of timed-out run(s) that did not produce the full matrices, you can disable this check, but be careful not to disrupt ongoing runs.
    With check_finished enabled, the code additionally checks presence of all types of matrices for two tracers. Pass two_tracers = True to enable this check explicitly, or two_tracers = False to disable it. By default (two_tracers = None), the code will try to guess whether the run is for two tracers by looking for the presence of the xi_22.dat file, which is only produced for two tracers by :func:`RascalC.run_cov`.
    """
    input_filename = os.path.join(file_root, "Raw_Covariance_Matrices" + "_3PCF" * threepcf + f"_{label}.npz")
    if os.path.isfile(input_filename): raw_cov = dict(np.load(input_filename))
    else:
        print_function(f"Collecting the raw covariance matrices from {file_root}")
        result = collect_raw_covariance_matrices(file_root, check_finished=check_finished, two_tracers=two_tracers, print_function=print_function)[threepcf]
        if label not in result:
            raise ValueError(f"Raw {2+threepcf}PCF covariance matrices for {label} not produced. Check the n and m/max_l values.")
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


def load_raw_covariances_smu(file_root: str, n: int, m: int, n_samples: None | int | Iterable[int] | Iterable[bool] = None, check_finished: bool = True, two_tracers: bool | None = None, print_function: Callable[[str], None] = print) -> dict[str]:
    """
    Load the raw covariance matrices from the s_mu mode as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.
    n and m are numbers of radial and angular (mu) bins respectively.

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.

    With check_finished enabled (default), before collecting raw covariance matrices (if needed), performs a heuristic check whether the run seems finished by looking for the presence of the full covariance matrices. If they are not found, the collection does not proceed and a warning is issued. This is to prevent disrupting ongoing runs by collecting and/or deleting the text files. If you want to check convergence of timed-out run(s) that did not produce the full matrices, you can disable this check, but be careful not to disrupt ongoing runs.
    With check_finished enabled, the code additionally checks presence of all types of matrices for two tracers. Pass two_tracers = True to enable this check explicitly, or two_tracers = False to disable it. By default (two_tracers = None), the code will try to guess whether the run is for two tracers by looking for the presence of the xi_22.dat file, which is only produced for two tracers by :func:`RascalC.run_cov`.
    """
    label = f"n{n}_m{m}"
    return load_raw_covariances(file_root, label, threepcf=False, n_samples=n_samples, check_finished=check_finished, two_tracers=two_tracers, print_function=print_function)


def load_raw_covariances_legendre(file_root: str, n: int, max_l: int, n_samples: None | int | Iterable[int] | Iterable[bool] = None, check_finished: bool = True, two_tracers: bool | None = None, print_function = print) -> dict[str]:
    """
    Load the raw covariance matrices from the Legendre modes as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.
    n and m are numbers of radial and angular (mu) bins respectively.

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.

    With check_finished enabled (default), before collecting raw covariance matrices (if needed), performs a heuristic check whether the run seems finished by looking for the presence of the full covariance matrices. If they are not found, the collection does not proceed and a warning is issued. This is to prevent disrupting ongoing runs by collecting and/or deleting the text files. If you want to check convergence of timed-out run(s) that did not produce the full matrices, you can disable this check, but be careful not to disrupt ongoing runs.
    With check_finished enabled, the code additionally checks presence of all types of matrices for two tracers. Pass two_tracers = True to enable this check explicitly, or two_tracers = False to disable it. By default (two_tracers = None), the code will try to guess whether the run is for two tracers by looking for the presence of the xi_22.dat file, which is only produced for two tracers by :func:`RascalC.run_cov`.
    """
    label = f"n{n}_l{max_l}"
    return load_raw_covariances(file_root, label, threepcf=False, n_samples=n_samples, check_finished=check_finished, two_tracers=two_tracers, print_function=print_function)


def load_raw_covariances_3pcf_legendre(file_root: str, n: int, max_l: int, n_samples: None | int | Iterable[int] | Iterable[bool] = None, check_finished: bool = True, print_function = print) -> dict[str]:
    """
    Load the raw covariance matrices from the 3PCF Legendre mode as a dictionary. Uses the Numpy file if it exists, otherwise tried to run the collection function.

    file_root is the directory to look in.
    n and m are numbers of radial and angular (mu) bins respectively.

    n_samples allows to select subsamples flexibly:
    - if None (default) returns all samples;
    - if a positive integer, returns as many samples from the beginning;
    - if a sequence of integers, returns subsamples with indices from this sequence;
    - sequence of boolean values is interpreted as a boolean mask for subsamples.

    With check_finished enabled (default), before collecting raw covariance matrices (if needed), performs a heuristic check whether the run seems finished by looking for the presence of the full covariance matrices. If they are not found, the collection does not proceed and a warning is issued. This is to prevent disrupting ongoing runs by collecting and/or deleting the text files. If you want to check convergence of timed-out run(s) that did not produce the full matrices, you can disable this check, but be careful not to disrupt ongoing runs.
    """
    label = f"n{n}_l{max_l}"
    return load_raw_covariances(file_root, label, threepcf=True, n_samples=n_samples, check_finished=check_finished, print_function=print_function)


def cat_raw_covariance_matrices(n: int, mstr: str, input_roots: list[str], ns_samples: list[None | int | Iterable[int] | Iterable[bool]], output_root: str, collapse_factor: int = 1, threepcf: bool = False, check_finished: bool = True, two_tracers: bool | None = None, print_function: Callable[[str], None] = print) -> dict[str]:
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
    
    threepcf : boolean
        (Optional) whether the runs are for 3PCF. Default is False, meaning 2PCF.

    check_finished : boolean
        (Optional) whether to check if the runs are finished before collecting raw covariance matrices (if that is needed). Default is True.

    two_tracers : boolean | None
        (Optional) whether the runs are for two tracers (for an additional check for finished runs). If None (default), the code will try to guess by looking for the presence of the xi_22.dat file.

    print_function : Callable[[str], None]
        (Optional) custom function to use for printing. Default is ``print``.

    Returns
    -------
    catenation_results : dict[str, npt.NDArray[np.float64]]
        Catenated raw covariance matrices as a dictionary with string keys and Numpy array values. All this information is also saved in a ``Raw_Covariance_Matrices*.npz`` file in the ``output_root``.
    """
    if collapse_factor <= 0: raise ValueError("Collapsing factor must be positive")
    if len(input_roots) < 1: raise ValueError("Need at least one input directory")
    if len(ns_samples) != len(input_roots): raise ValueError("Number of input dirs and subsamples to use from them must be the same")

    label = f"n{n}_{mstr}"
    result = {}
    for index, (input_root, n_samples) in enumerate(zip(input_roots, ns_samples)):
        input_file = load_raw_covariances(input_root, label, threepcf=threepcf, n_samples=n_samples, check_finished=check_finished, two_tracers=two_tracers, print_function=print_function)
        # ignore full arrays
        input_file = {key: value for (key, value) in input_file.items() if not key.endswith("_full")}
        # check that the keys are the same, unless the result is brand new
        if len(result) > 0:
            result_keys = set(result.keys())
            input_keys = set(input_file.keys())
            if result_keys != input_keys:
                print_function("Different sets of matrices present among the input files, will only use the overlapping ones.")
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
    
    save_safe(output_root, label, result, print_function)

    # copy other useful files from the first input root, unless identical with the output root
    # assuming they are the same among output roots; otherwise catenation should not be sensible
    if not os.path.samefile(input_roots[0], output_root):
        for pattern in ("weights", "xi*", "radial_binning*.csv"):
            for filename in glob(pattern, root_dir = input_roots[0]):
                src = os.path.join(input_roots[0], filename)
                (copytree if os.path.isdir(src) else copy2)(src, os.path.join(output_root, filename)) # need different functions for dirs and files

    return result