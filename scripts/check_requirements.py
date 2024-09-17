import subprocess
import sys


def check_requirements(requirements):
    """Checks that the requirements for the Chromosearch pipeline are installed and available.
    Exits the program if anything is missing with code 2.

    Args:
        requirements (dic): Dictionary of requirements, "packages" key as a list of tuples ("package_name", "version"), with the other keys only having a version string to check.
    """

    print("Checking requirements for chromosearch...")
    # Check each python module for being installed
    for package, version in requirements["packages"]:
        if not check_package(package, version):
            system_exit_missing_dependencies()

    # Check Prodigal and BLAST
    if not check_blast(requirements["blast_version"]) or not check_prodigal(
        requirements["prodigal_version"]
    ):
        system_exit_missing_dependencies()

    # Check Python version, doesn't return
    check_python_version(requirements["python_version"])

    return


def system_exit_missing_dependencies():

    print(
        "Chromosearch is missing some required dependencies, please install them, ensure they are available in PATH,and then rerun the pipeline."
    )
    sys.exit(2)


def check_package(package_name, version=None):
    """Checks whether a package in the requirements is installed and available with the correct version.

    Args:
        package_name (str): The name of the package to import
        version (str, optional): _description_. The version to check for the package. Defaults to None.

    Returns:
        bool: True for the package being installed, False in the opposite case. Version mismatch still returns True, since it is not critical.
    """
    try:

        pkg = __import__(package_name)
        if version:
            installed_version = pkg.__version__
            if installed_version != version:
                print(
                    f"Warning: {package_name} version mismatch. Expected: {version}, Installed: {installed_version}. This might lead to errors or inconsistent behaviour."
                )
        return True

    except ImportError:
        print(f"{package_name} is not installed.")
        return False


def check_blast(version):
    """Checks whether BLAST is installed and available. Also checks for the correct version.

    Args:
        version (str): Expected BLAST version.

    Returns:
        bool: True for BLAST being installed, False in the opposite case. Version mismatch still returns True, since it is not critical.
    """
    try:
        result = subprocess.run(
            ["blastp", "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if version not in result.stdout:
            print(
                "Warning: BLAST version mismatch, this may lead to errors or inconsistent  behaviour."
            )

        return True

    except FileNotFoundError:
        print("BLAST is not installed.")

        return False


def check_prodigal(version):
    """Checks whether Prodital is installed and available. Also checks for the correct version.

    Args:
        version (str): Expected Prodigal version.

    Returns:
        bool: True for Prodital being installed, False in the opposite case. Version mismatch still returns True, since it is not critical.
    """

    try:
        result = subprocess.run(
            ["prodigal", "-v"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # subprocess returns output in stderr for some reason
        if version not in result.stderr:
            print(
                "Warning: Prodigal version mismatch, this may lead to errors or inconsistent  behaviour."
            )

        return True

    except FileNotFoundError:
        print("Prodigal is not installed or cannot be found in PATH.")

        return False


def check_python_version(version):
    """Gives warning for different version of Python than the one tested."""

    # Check Python version
    python_version = sys.version.split()[0]
    if python_version != version:
        print(
            f"Warning: Python version mismatch. Expected: 3.12.5, Installed: {python_version}. This may lead to errors, or unexpected behaviour."
        )
