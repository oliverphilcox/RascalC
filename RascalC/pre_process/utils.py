import numpy as np
from astropy.io import fits
from warnings import warn


def my_str_to_bool(s: str) -> bool:
    # naive conversion to bool for all non-empty strings is True, and one can't give an empty string as a command line argument, so need to make it more explicit
    return s not in ("0", "false")


def parse_FKP_arg(FKP_weights: str) -> bool | tuple[float, str]:
    if not my_str_to_bool(FKP_weights): return False
    # determine if it actually has P0,NZ_name format. Such strings should convert to True
    arg_FKP_split = FKP_weights.split(",")
    if len(arg_FKP_split) == 2:
        return (float(arg_FKP_split[0]), arg_FKP_split[1])
    if len(arg_FKP_split) == 1: return True
    raise ValueError("FKP parameter matched neither USE_FKP_WEIGHTS (true/false in any register or 0/1) nor P0,NZ_name (float and string without space).")


def read_particles_fits_file(input_file: str, FKP_weights: bool | tuple[float, str] = False, mask: int = 0, use_weights: bool = True):
    # Read FITS file with particles. Can apply mask filtering and compute FKP weights in different ways. Works for DESI setups
    with fits.open(input_file) as f:
        data = f[1].data
        all_ra = data["RA"]
        all_dec = data["DEC"]
        all_z = data["Z"]
        colnames = data.columns.names
        all_w = data["WEIGHT"] if "WEIGHT" in colnames and use_weights else np.ones_like(all_z)
        if FKP_weights:
            all_w *= 1/(1+FKP_weights[0]*data[FKP_weights[1]]) if FKP_weights != True else data["WEIGHT_FKP"]
        if "WEIGHT" not in colnames and not FKP_weights: warn("No weights found, assigned unit weight to each particle.")
        if mask: filt = (data["STATUS"] & mask == mask) # all 1-bits from mask have to be set in STATUS; skip if mask=0
    result = np.array((all_ra, all_dec, all_z, all_w)).T
    if mask: result = result[filt]
    return result