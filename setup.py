import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PurityReviewer",
    version="1.0.0",
    author="Claudia Chu",
    author_email="cchu@broadinstitute.org",
    description="Suite of purity review dashboards",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/getzlab/PurityReviewer",
    project_urls={
        "Bug Tracker": "https://github.com/getzlab/PurityReviewer/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.8",
    install_requires = [
                        'ipykernel==6.22.0',
                        'pandas==1.5.2',
                        'AnnoMate>=1.0.0',
                        'setuptools>=70.0.0', # security vulnerability with versions up to 69.1.1 
                        'natsort',
                        'rpy2==3.4.2',
                        'firecloud-dalmatian==0.0.18',
                        'cnv_suite',
                        'numpy==1.26.0',
                        'firecloud-dalmatian',
                        'jupyter',
                        'scipy',
                        'sklearn'
                        'jupyter'
                       ]
)   
