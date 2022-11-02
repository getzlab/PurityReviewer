import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PurityReviewers",
    version="0.0.1",
    author="Claudia Chu",
    author_email="cchu@broadinstitute.org",
    description="Suite of purity review dashboards",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/getzlab/PurityReviewers",
    project_urls={
        "Bug Tracker": "https://github.com/getzlab/PurityReviewers/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.6",
    install_requires = ['dash',
                        'dash-bootstrap-components',
                        'fsspec',
                        'gcsfs',
                        'google-auth',
                        'google-api-core',
                        'hound',
                        'ipykernel',
                        'ipython',
                        'jupyter-dash',
                        'jupyterlab',
                        'matplotlib',
                        'pandas',
                        'pickleshare',
                        'pillow',
                        'pip',
                        'plotly',
                        'setuptools']
)   
