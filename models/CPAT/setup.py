import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for CPAT  -- Coding Potential Assessment Tool
"""

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: CPAT requires Python 2.7"
	sys.exit()
	

def main():
	setup(  name = "CPAT",
			version = "1.2.4",
			py_modules = [ 'psyco_full' ],
			packages = find_packages( 'lib' ),
			package_dir = { '': 'lib' },
			package_data = { '': ['*.ps'] },
			scripts = glob.glob( "bin/*.py"),
			ext_modules = [],
			test_suite = 'nose.collector',
			setup_requires = ['nose>=0.10.4'],
			author = "Liguo Wang, Jung Hyun Park",
			author_email ="wangliguo78@gmail.com, hjpark@bcm.edu",
			platforms = ['Linux','MacOS'],
			requires = ['cython (>=0.17)'],
			install_requires = ['numpy','pysam'], 
			description = "CPAT (Coding Potential Assessment Tool)",
			long_description = "CPAT is an alignment-free method to predict RNA coding potential using four sequence features.",
			license='GNU General Public License',
			url = "http://rna-cpat.sourceforge.net/",
			zip_safe = False,
			dependency_links = [],
			classifiers=[
				'Development Status :: 5 - Production/Stable',
				'Environment :: Console',
				'Intended Audience :: Science/Research',
				'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: POSIX',
				'Programming Language :: Python',
				'Programming Language :: Python :: 2.7',
				'Topic :: Scientific/Engineering :: Bio-Informatics',
			],
			
			keywords='RNA coding potential prediction, lncRNA, lincRNA, logsitic regression',
			python_requires = '>=2.6, <3'
             )


if __name__ == "__main__":
	main()
