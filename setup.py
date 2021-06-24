#!/usr/bin/env python

import sys
import subprocess
from setuptools import setup, find_packages

def check_minimap2():
    p = subprocess.Popen(['minimap2', '--version'], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()
        major, minor = line.strip().split('.')
        if int(major) >= 2:
            return True
    return False


def check_samtools():
    p = subprocess.Popen(['samtools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()
        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_mafft():
    p = subprocess.Popen(['mafft', '--help'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()
        line = line.strip()
        if line.startswith('MAFFT'):
            return True
    return False


def check_python():
    return sys.version_info.major == 3 and sys.version_info.minor >= 6    


def check_exonerate():
    p = subprocess.Popen(['exonerate'], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()
        if line.startswith('exonerate from exonerate'):
            major, minor = line.strip().split()[-1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 2 and int(minor) >= 2:
                return True
    return False


if __name__ == '__main__':
    if not check_python(): sys.exit('Dependency problem: python >= 3.6 is required')
    if not check_minimap2(): sys.exit('Dependency problem: minimap2 >= 2.0 not found')
    if not check_samtools(): sys.exit('Dependency problem: samtools >= 1.2 not found')
    if not check_mafft(): sys.exit('Dependency problem: mafft not found')
    if not check_exonerate(): sys.exit('Dependency problem: exonerate not found')

setup(
    name='tldr',
    version='1.2.1',
    author='Adam Ewing',
    author_email='adam.ewing@gmail.com',
    description=("Insertion finder for long noisy reads"),
    license='MIT',
    url='https://github.com/adamewing/tldr',
    scripts=['tldr/tldr'],
    packages=find_packages(),
    install_requires = [
        'cython>=0.29.0',
        'pysam>=0.8.1',
        'bx-python>=0.5.0',
        'scikit-learn>=0.20.0',
        'numpy>=1.9.0'
    ]

)


