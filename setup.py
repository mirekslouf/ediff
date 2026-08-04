import setuptools

def get_version():
    with(open("src/ediff/__init__.py", "r")) as fh:
        for line in fh:
            if line.startswith("__version__"):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
        else:
            raise RuntimeError("Unable to find version string.")
                
def get_long_description():
    with open("README.md", "r") as fh: description = fh.read()
    return(description)
    
setuptools.setup(
    name="ediff",
    version=get_version(),
    author="Mirek Slouf",
    author_email="mirek.slouf@gmail.com",
    description="Processing of electron diffraction patterns",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/mirekslouf/ediff/",
    project_urls={
        "Documentation": "https://mirekslouf.github.io/ediff/"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
    license='MIT',
    package_dir={"":"src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    include_package_data=True,
    
    # GUI extras
    # Keeps the core ediff package lightweight — these are only installed
    # with `pip install ediff[gui]`, not with a plain `pip install ediff`.
    extras_require={
        "gui": ["qtpy", "PyQt6", "matplotlib", "opencv-python"],
    },
    
    # GUI can be run directly as a command after installation,
    # equivalent to `python -m ediff gui`.
    entry_points={
        "console_scripts": [
            "ediff-gui=ediff.gui.app:main",
        ],
    },
    
    # The bundled logo/icon images get included in the installed package 
    # — without this, pip install would succeed but the GUI would show blank 
    # images, since find_packages() only discovers .py files on its own.
    package_data={
        "ediff.gui.resources": ["*.png", "*.ico"],
    },
    )
