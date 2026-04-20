from pathlib import Path
from typing import Literal


def validate_allcounts_format(allcounts_format: Literal[None, "pycorr", "lsstypes"], allcounts_files: list[str | Path]) -> Literal["pycorr", "lsstypes"]:
    """
    If the format of the allcounts files is not specified (None), try to guess it from the file extensions. If the format is specified, check that it is valid.

    Parameters
    ----------
    allcounts_format : None, "pycorr" or "lsstypes"
        The format of the allcounts files. If None, the format will be guessed from the filename extensions.

    allcounts_files : list of strings (paths) or Path objects
        Filenames for the allcounts files. Only used to guess the format if allcounts_format is None. Expected to all have .npy extensions for pycorr format or all .h5/.hdf5/.txt extensions for lsstypes format.

    Returns
    -------
    allcounts_format : string
        The specific format of the allcounts files, either "pycorr" or "lsstypes".
    """
    if allcounts_format is not None:
        if allcounts_format in ["pycorr", "lsstypes"]:
            return allcounts_format
        raise ValueError("Invalid allcounts_format, must be either None, 'pycorr' or 'lsstypes'")
    allcounts_files : list[Path] = [Path(allcounts_file) for allcounts_file in allcounts_files] # convert to Path objects
    if allcounts_files[0].suffix == ".npy":
        if all(allcounts_file.suffix == ".npy" for allcounts_file in allcounts_files):
            return "pycorr"
        raise ValueError("Inconsistent allcounts file extensions, expected all .npy for pycorr format (or all .h5/.hdf5/.txt for lsstypes format). Specify the format via allcounts_format='pycorr' (or 'lsstypes') to avoid this error.")
    if allcounts_files[0].suffix in [".h5", ".hdf5", ".txt"]:
        if all(allcounts_file.suffix in [".h5", ".hdf5", ".txt"] for allcounts_file in allcounts_files):
            return "lsstypes"
        raise ValueError("Inconsistent allcounts file extensions, expected all .h5/.hdf5/.txt for lsstypes format (or all .npy for pycorr format). Specify the format via allcounts_format='lsstypes' (or 'pycorr') to avoid this error.")
    raise ValueError("Could not guess the format of allcounts files, please specify via allcounts_format='pycorr' (or 'lsstypes')")