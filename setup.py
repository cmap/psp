import setuptools

setuptools.setup(
    name="psp",
    version="1.0",
    description="Package for processing GCP and P100 proteomics signatures",
    author="Lev Litichevskiy",
    author_email="lev@broadinstitute.org",
    url="https://github.com/cmap/proteomics-signature-pipeline.git",
    packages=["in_out", "utils", "dry", "steep"])
