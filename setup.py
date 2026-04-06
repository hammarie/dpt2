#!/usr/bin/env python3
"""
Cross-platform setup and installation script for DPT2.
Handles environment setup, dependency installation, and validation.
"""
import sys
import os
import subprocess
import platform
from pathlib import Path


def print_header(message):
    """Print a formatted header."""
    print(f"\n{'=' * 60}")
    print(f"  {message}")
    print(f"{'=' * 60}\n")


def check_python_version():
    """Check if Python version is adequate."""
    print("Checking Python version...")
    if sys.version_info < (3, 9):
        print(f"❌ Python 3.9+ required, but you have {sys.version_info.major}.{sys.version_info.minor}")
        sys.exit(1)
    print(f"✓ Python {sys.version_info.major}.{sys.version_info.minor} OK")


def check_pip():
    """Check if pip is available."""
    print("Checking pip...")
    try:
        import pip
        print(f"✓ pip is available")
        return True
    except ImportError:
        print("❌ pip not found. Please install pip: https://pip.pypa.io/en/latest/installation/")
        return False


def install_package(editable=True):
    """Install the DPT2 package."""
    print_header("Installing DPT2")

    cmd = [sys.executable, "-m", "pip", "install"]
    if editable:
        cmd.append("-e")
    cmd.append(".")

    print(f"Running: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, cwd=Path(__file__).parent)

    if result.returncode == 0:
        print("\n✓ DPT2 installation successful")
        return True
    else:
        print("\n❌ DPT2 installation failed")
        return False


def check_blast():
    """Check if BLAST is installed."""
    print("\nChecking for BLAST+...")
    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            print(f"✓ BLAST+ is available: {version_line}")
            return True
    except FileNotFoundError:
        pass

    print("⚠ BLAST+ not found (optional - required only for --blast-screen)")
    print("Installation instructions:")
    if platform.system() == "Darwin":  # macOS
        print("  macOS: brew install blast")
    elif platform.system() == "Linux":
        print("  Ubuntu/Debian: sudo apt-get install ncbi-blast+")
    elif platform.system() == "Windows":
        print("  Windows: Download from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
    return False


def validate_installation():
    """Validate that the package is installed correctly."""
    print_header("Validating Installation")

    try:
        import mpx_dpcr
        print(f"✓ mpx_dpcr package found at: {mpx_dpcr.__file__}")
    except ImportError:
        print("❌ Failed to import mpx_dpcr")
        return False

    # Check key modules
    modules = [
        ("mpx_dpcr.fetch", "Sequence fetching"),
        ("mpx_dpcr.design", "Primer design"),
        ("mpx_dpcr.evaluate", "Thermodynamic evaluation"),
        ("mpx_dpcr.blast_check", "BLAST screening"),
    ]

    all_ok = True
    for module_name, description in modules:
        try:
            __import__(module_name)
            print(f"✓ {description} module OK")
        except ImportError as e:
            print(f"❌ {description} module failed: {e}")
            all_ok = False

    return all_ok


def print_quick_start():
    """Print quick start guide."""
    print_header("Quick Start")

    print("To design primers for genes RPP30 and MRGPRX1:\n")
    print("  python scripts/design_and_score_primers.py --genes RPP30 MRGPRX1\n")

    print("To use BLAST specificity screening:")
    print("  (First, set up the BLAST database - see INSTALL.md)\n")
    print("  python scripts/design_and_score_primers.py \\")
    print("    --genes RPP30 MRGPRX1 \\")
    print("    --blast-screen \\")
    print("    --blast-db /path/to/blast/database/human_rna\n")

    print("For more help:")
    print("  python scripts/design_and_score_primers.py --help")
    print("  See INSTALL.md for detailed instructions\n")


def main():
    """Main setup routine."""
    print_header("DPT2 Setup and Installation")
    print(f"Platform: {platform.system()} {platform.release()}")
    print(f"Python: {sys.version}\n")

    # Check prerequisites
    check_python_version()
    if not check_pip():
        sys.exit(1)

    # Install package
    if not install_package(editable=True):
        sys.exit(1)

    # Validate
    if not validate_installation():
        print("\n⚠ Installation validation failed - some modules are missing")
        print("Please check the error messages above and reinstall")
        sys.exit(1)

    # Optional checks
    check_blast()

    # Print summary
    print_header("Setup Complete!")
    print("✓ All required dependencies are installed")
    print("✓ DPT2 is ready to use\n")

    print_quick_start()


if __name__ == "__main__":
    main()
