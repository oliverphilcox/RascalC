import sys, os

def adjust_path():
    # insert the parent directory, override the one installed in the system
    sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "../..")))
