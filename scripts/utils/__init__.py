import sys, os

def adjust_path():
    # insert the parent directory, override the one installed in the system
    sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "../..")))


def get_arg_safe(index: int, type = str, default: object = None) -> object:
    # get argument by index from sys.argv and convert it to the requested type if there are enough elements there
    # otherwise return the default value
    return type(sys.argv[index]) if len(sys.argv) > index else default
