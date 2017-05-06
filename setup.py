from setuptools import setup, find_packages
setup(
    name = "GeneLearn",
    version = "0.1",
    packages = find_packages(),
    scripts = ['say_hello.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = ['docutils>=0.3', 'ete2', 'scikit-learn',
                        'matplotlib'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        'hello': ['*.msg'],
    },

    # metadata for upload to PyPI
    author = "Me",
    author_email = "me@example.com",
    description = "This is an Example Package",
    license = "PSF",
    keywords = "hello world example examples",
    url = "https://bitbucket.org/berkeleylab/jgi-genelearn",   # project home page, if any

    # could also include long_description, download_url, classifiers, etc.
    test_suite = "nose.collector",
)