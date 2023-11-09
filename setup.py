from setuptools import setup
from setuptools import find_packages
import pathlib

# Get the long description from the README file
long_description = (pathlib.Path(__file__).parent.resolve() / 'README.md').read_text(encoding='utf-8')


setup(name='vlbiplanobs',
        # Versions should comply with PEP 440:
        # https://www.python.org/dev/peps/pep-0440/
        # [N!]N(.N)*[{a|b|rc}N][.postN][.devN]
        version='4.1.2',
        # one-line description or tagline of what your project does
        description='Planner for VLBI observations',
        # "Description" metadata field
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/bmarcote/vlbi_calculator',
        download_url='https://github.com/bmarcote/vlbi_calculator/tarball/2.0.3',
        author='B. Marcote',
        author_email='marcote@jive.eu',
        # Classifiers help users find your project by categorizing it.
        # For a list of valid classifiers, see https://pypi.org/classifiers/
        classifiers=[  # Optional
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 5 - Production/Stable',
            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Astronomy',
            # Pick your license as you wish
            'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
            # Specify the Python versions you support here. In particular, ensure
            # that you indicate you support Python 3. These classifiers are *not*
            # checked by 'pip install'. See instead 'python_requires' below.
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3 :: Only',
        ],
        keywords='astronomy, astrophysics, cosmology, radio, science, observations, planner, coordinate',
        python_requires='>=3.7, <4',
        # When your source code is in a subdirectory under the project root, e.g.
        # `src/`, it is necessary to specify the `package_dir` argument.
        # package_dir={'': 'src'},  # Optional
        packages=find_packages(),
        scripts=['bin/vlbiplanobs'],
        # This field lists other packages that your project depends on to run.
        # Any package you put here will be installed by pip when your project is
        # installed, so they must be valid existing projects.
        #
        # For an analysis of "install_requires" vs pip's requirements files see:
        # https://packaging.python.org/en/latest/requirements.html
        install_requires=['numpy', 'astropy>=4.0.2', 'astroplan>=0.7'],

        # If there are data files included in your packages that need to be
        # installed, specify them here.
        package_data={
            'stations': ['data/stations_catalog.inp'],
            'docs': ['docs/*'],
            'assets': ['assets/*']
        },

        # Although 'package_data' is the preferred approach, in some case you may
        # need to place data files outside of your packages. See:
        # http://docs.python.org/distutils/setupscript.html#installing-additional-files
        #
        # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
        # data_files=[('data', ['data/stations_catalog.inp'])],
        include_package_data=True,
        # List additional URLs that are relevant to your project as a dict.
        #
        # This field corresponds to the "Project-URL" metadata fields:
        # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
        #
        # Examples listed include a pattern for specifying where the package tracks
        # issues, where the source is hosted, where to say thanks to the package
        # maintainers, and where to support the project financially. The key is
        # what's used to render the link text on PyPI.
        project_urls={  # Optional
            'Bug Reports': 'https://github.com/bmarcote/vlbi_calculator/issues',
            # 'Funding': 'https://donate.pypi.org',
            # 'Say Thanks!': 'http://saythanks.io/to/example',
            'Source': 'https://github.com/bmarcote/vlbi_calculator/',
        },
    )
