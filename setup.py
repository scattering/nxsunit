from setuptools import setup

def find_version(package, filename="__init__.py"):
    """Read package version string from __init__.py"""
    import os
    with open(os.path.join(package, filename)) as fid:
        for line in fid.readlines():
            if line.startswith('__version__'):
                line = line[:line.find('#')]  # strip comment, if any
                version = line.split('=', 2)[1].strip()
                if version[0] != version[-1] or version[0] not in "\"'":
                    break
                return version[1:-1]
    raise RuntimeError("Could not read version from %s/__init__.py"%package)

with open('README.md') as fid:
    long_description = fid.read()

setup(
    name='nxsunit',
    version=find_version('', 'nxsunit.py'),
    description="NeXus data field loader with unit conversion",
    long_description=long_description,
    author="Paul Kienzle",
    author_email="paul.kienzle@nist.gov",
    url="https://github.com/scattering/nxsunit",
    keywords="x-ray neutron units",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Public Domain',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    py_modules=['nxsunit'],
    scripts=['nxsunit.py'],
)
