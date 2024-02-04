#!/usr/bin/env python3


import os
import os.path as path
import fnmatch
import subprocess
import sys
from sys import version_info

import RascalC

base_url = "https://github.com/misharash/RascalC"
projectname = 'RascalC'
version = RascalC.__version__

# Make sure we are running on posix (Linux, Unix, MAC OSX)
if os.name != 'posix':
    sys.exit("Sorry, Windows is not supported")

min_py_major, min_py_minor = 3, 8
min_np_major, min_np_minor = 1, 23

# Enforce minimum python version
if version_info[0] < min_py_major or \
   (version_info[0] == min_py_major and version_info[1] < min_py_minor):
    raise RuntimeError('Sorry. Found python {0}.{1} but minimum required \
    python version is {2}.{3}'.format(version_info[0],
                                      version_info[1],
                                      min_py_major, min_py_minor))

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


def run_command(command):
    # print("about to execute command `{0}`. sources = {1}"
    # .format(command, sources))
    proc = subprocess.Popen(command, stderr=subprocess.STDOUT, shell=True)
    output, stderr = proc.communicate(input)
    status = proc.wait()
    if status:
        raise Exception("command = {0} failed with output = {1} status {2:d}\n"
                        .format(command, output, status))


# Only python >= 3.5 supports the recursive glob, hence
# defining the function that works on all reasonable pythons
# http://stackoverflow.com/questions/2186525/use-a-glob-to-
# find-files-recursively-in-python
def recursive_glob(rootdir='.', patterns=['*']):
    return [path.join(looproot, filename)
            for looproot, _, filenames in os.walk(rootdir)
            for filename in filenames for p in patterns
            if fnmatch.fnmatch(filename, p)]

def rd(filename):
    with open(filename) as f:
        return f.read()

# Taken from numpy setup.py
def setup_packages():

    # protect the user in case they run python setup.py not from root directory
    src_path = path.dirname(path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # create a list of the python extensions
    extensions = []

    # run make
    command = "cd RascalC && make"
    run_command(command)

    # find all the data-files required.
    # Now the lib + associated header files have been generated
    # and put in lib/ and include/
    # This step must run after ``make install``
    dirs_patterns = {'RascalC/bin': ["cov.*"]}
    data_files = []
    for d in dirs_patterns:
        patterns = dirs_patterns[d]
        f = recursive_glob(d, patterns)
        data_files.extend(f)

    # change them to be relative to package dir rather than root
    data_files = ["../{0}".format(d) for d in data_files]

    # Fix long description for PyPI
    try:
        import pypandoc
        long_description = pypandoc.convert('README.md', 'rst')
    except(IOError, ImportError):
        long_description = rd('README.md')

    # All book-keeping is done.
    classifiers = ['Development Status :: 3 - Alpha',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: POSIX',
                   'Programming Language :: C++',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10',
                   'Programming Language :: Python :: 3.11']
    metadata = dict(
        name=projectname,
        version=version,
        author='Michael Rashkovetskyi',
        author_email='misharash@gmail.com',
        maintainer='Michael Rashkovetskyi',
        maintainer_email='misharash@gmail.com',
        url=base_url,
        #download_url='{0}/archive/{1}-{2}.tar.gz'.format(
        #    base_url, projectname, version),
        description='A Fast Code for Galaxy Covariance Matrix Estimation',
        long_description=long_description,
        classifiers=classifiers,
        # license='MIT',
        # Solaris might work, Windows will almost certainly not work
        platforms=["Linux", "Mac OSX", "Unix"],
        keywords=['correlation functions', 'simulations',
                  'surveys', 'galaxies'],
        packages=[projectname] + [projectname + "." + submodule for submodule in ("comb", "obj", "post_process", "pre_process", "pycorr_utils", "xi")],
        # ext_package=projectname,
        # ext_modules=extensions,
        package_data={'': data_files},
        include_package_data=True,
        install_requires=['setuptools',
                          'numpy>={0}.{1}'.format(min_np_major, min_np_minor),
                          'scipy',
                          'astropy',
                          'pycorr'],
        zip_safe=False)

    # Now the actual setup
    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_packages()
