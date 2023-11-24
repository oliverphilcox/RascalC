## Script to catenate (take a subset of subsamples or resample) raw covariance matrix files produced by collect_raw_covariance_matrices.py
## More specifically, copy part of partial results to other directory and recompute totals by averaging
## Determines single-field vs multi-field and jackknife automatically
## Do not use if subsamples have different numbers of pairs/triples/quadruplets

import numpy as np
import sys
from warnings import warn
from collect_raw_covariance_matrices import load_raw_covariances, save_safe


def cat_raw_covariance_matrices(n: int, mstr: str, input_roots: list[str], ns_samples: list[int], output_root: str, collapse_factor: int = 1, print_function = print) -> dict[str]:
    if collapse_factor <= 0: raise ValueError("Collapsing factor must be positive")
    if len(input_roots) < 1: raise ValueError("Need at least one input directory")
    if len(ns_samples) != len(input_roots): raise ValueError("Number of input dirs and subsamples to use from them must be the same")

    label = f"n{n}_{mstr}"
    result = {}
    for index, (input_root, n_samples) in enumerate(zip(input_roots, ns_samples)):
        input_file = load_raw_covariances(input_root, label, print_function)
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
        for matrix_name, matrices in input_file.values():
            if matrix_name.endswith("_full"): continue # ignore full arrays
            if n_samples is None: n_samples = len(matrices)
            if index != 0: result[matrix_name] = np.append(result[matrix_name], matrices[:n_samples], axis = 0)
            else: result[matrix_name] = matrices[:n_samples]
    
    # loop over all the matrix names
    matrix_names = result.keys() # the dictionary will be changed
    for matrix_name in matrix_names:
        if collapse_factor > 1:
            matrix_shape = result[matrix_name].shape
            result[matrix_name] = np.mean(result[matrix_name].reshape(matrix_shape[0] // collapse_factor, collapse_factor, *matrix_shape[1:]), axis = 1) # average over adjacent collapse_factor samples
        # make full arrays by averaging the subsamples
        matrix_name_full = matrix_name + "_full"
        result[matrix_name_full] = np.mean(result[matrix_name], axis = 0)
    
    save_safe(output_root, label, result)

    return result

if __name__ == "__main__": # if invoked as a script
    # PARAMETERS
    if len(sys.argv)<6: # if too few
        print("Usage: python cat_raw_covariance_matrices.py {N_R_BINS} {mN_MU_BINS/lMAX_L} {COVARIANCE_INPUT_DIR1} {N_SUBSAMPLES_TO_USE1} [{COVARIANCE_INPUT_DIR2} {N_SUBSAMPLES_TO_USE2} ...] [{COLLAPSE_FACTOR}] {COVARIANCE_OUTPUT_DIR}")
        sys.exit(1)

    n = int(sys.argv[1])
    mstr = str(sys.argv[2])
    input_roots = [str(s) for s in sys.argv[3:-1:2]]
    ns_samples = [int(s) for s in sys.argv[4:-1:2]]
    output_root = str(sys.argv[-1])
    collapse_factor = int(input_roots.pop()) if len(input_roots) > len(ns_samples) else 1 # recover the collapse factor if present

    cat_raw_covariance_matrices(n, mstr, input_roots, ns_samples, output_root, collapse_factor)