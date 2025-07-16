"""
Compatibility utilities for handling NumPy 2.0 and other dependency issues.

This module provides functions to check for compatibility issues and provide
helpful error messages and solutions for common dependency problems.
"""

import sys
import warnings
from typing import Optional, Tuple, Dict, Any


def check_numpy_compatibility() -> Tuple[bool, Optional[str]]:
    """
    Check NumPy version compatibility.
    
    Returns
    -------
    Tuple[bool, Optional[str]]
        (is_compatible, error_message)
        - is_compatible: True if NumPy version is compatible
        - error_message: Error message if not compatible, None if compatible
    """
    try:
        import numpy as np
        numpy_version = np.__version__
        
        # Parse version
        major, minor = map(int, numpy_version.split('.')[:2])
        
        if major >= 2:
            error_msg = f"""
NumPy 2.0+ compatibility issue detected (current version: {numpy_version})

PROBLEM:
Your environment has NumPy 2.0 or higher, which has compatibility issues with 
some dependencies like pandas and pyarrow that haven't been fully updated yet.

SOLUTIONS:
1. Downgrade NumPy to 1.x:
   pip install "numpy<2.0"

2. Create a new environment with compatible versions:
   conda create -n pymut python=3.10 "numpy<2.0" "pandas<2.1"
   conda activate pymut
   pip install -r requirements.txt

3. Use pip to install compatible versions:
   pip install "numpy>=1.20.0,<2.0.0" "pandas>=1.3.0,<2.1.0"

4. If using conda, specify versions:
   conda install "numpy<2.0" "pandas<2.1"

TECHNICAL DETAILS:
NumPy 2.0 introduced breaking changes in the C API that affect compiled
extensions like pyarrow, which pandas depends on. Many packages are still
updating their compatibility.
"""
            return False, error_msg
        
        return True, None
        
    except ImportError:
        error_msg = """
NumPy not found in your environment.

SOLUTION:
Install NumPy with: pip install "numpy>=1.20.0,<2.0.0"
"""
        return False, error_msg
    except Exception as e:
        error_msg = f"""
Error checking NumPy compatibility: {e}

SOLUTION:
Try reinstalling NumPy: pip install --force-reinstall "numpy>=1.20.0,<2.0.0"
"""
        return False, error_msg


def check_pandas_compatibility() -> Tuple[bool, Optional[str]]:
    """
    Check pandas compatibility with current NumPy version.
    
    Returns
    -------
    Tuple[bool, Optional[str]]
        (is_compatible, error_message)
    """
    try:
        import pandas as pd
        return True, None
    except ImportError as e:
        if "numpy.core.multiarray" in str(e) or "_ARRAY_API" in str(e):
            error_msg = """
Pandas import failed due to NumPy compatibility issues.

PROBLEM:
Pandas was compiled against an older version of NumPy and cannot run with
NumPy 2.0+. This is a common issue with NumPy 2.0 compatibility.

SOLUTIONS:
1. Downgrade NumPy and reinstall pandas:
   pip install "numpy<2.0"
   pip install --force-reinstall pandas

2. Use compatible versions:
   pip install "numpy>=1.20.0,<2.0.0" "pandas>=1.3.0,<2.1.0"

3. Create a fresh environment:
   conda create -n pymut python=3.10
   conda activate pymut
   conda install "numpy<2.0" "pandas<2.1"
"""
        else:
            error_msg = f"""
Pandas import failed: {e}

SOLUTION:
Install pandas with: pip install "pandas>=1.3.0,<2.1.0"
"""
        return False, error_msg
    except Exception as e:
        error_msg = f"""
Error checking pandas compatibility: {e}

SOLUTION:
Try reinstalling pandas: pip install --force-reinstall "pandas>=1.3.0,<2.1.0"
"""
        return False, error_msg


def check_pyensembl_compatibility() -> Tuple[bool, Optional[str]]:
    """
    Check pyensembl availability and compatibility.
    
    Returns
    -------
    Tuple[bool, Optional[str]]
        (is_compatible, error_message)
    """
    try:
        from pyensembl import EnsemblRelease
        return True, None
    except ImportError:
        error_msg = """
pyensembl not found in your environment.

SOLUTION:
Install pyensembl with: pip install pyensembl

NOTE: pyensembl is required for transcript regions generation in the 
OncodriveCLUSTL pipeline. Without it, you can only run data preparation.
"""
        return False, error_msg
    except Exception as e:
        error_msg = f"""
Error checking pyensembl compatibility: {e}

SOLUTION:
Try reinstalling pyensembl: pip install --force-reinstall pyensembl
"""
        return False, error_msg


def check_all_dependencies() -> Dict[str, Tuple[bool, Optional[str]]]:
    """
    Check all major dependencies for compatibility issues.
    
    Returns
    -------
    Dict[str, Tuple[bool, Optional[str]]]
        Dictionary with dependency names as keys and (is_compatible, error_message) as values
    """
    checks = {
        'numpy': check_numpy_compatibility(),
        'pandas': check_pandas_compatibility(),
        'pyensembl': check_pyensembl_compatibility()
    }
    
    return checks


def print_compatibility_report(verbose: bool = True) -> bool:
    """
    Print a comprehensive compatibility report.
    
    Parameters
    ----------
    verbose : bool, default True
        Whether to print detailed error messages
    
    Returns
    -------
    bool
        True if all dependencies are compatible, False otherwise
    """
    print("=" * 60)
    print("DEPENDENCY COMPATIBILITY REPORT")
    print("=" * 60)
    
    checks = check_all_dependencies()
    all_compatible = True
    
    for dep_name, (is_compatible, error_msg) in checks.items():
        status = "✓ COMPATIBLE" if is_compatible else "✗ INCOMPATIBLE"
        print(f"{dep_name:>12}: {status}")
        
        if not is_compatible:
            all_compatible = False
            if verbose and error_msg:
                print(error_msg)
                print("-" * 60)
    
    if all_compatible:
        print("\n🎉 All dependencies are compatible!")
    else:
        print(f"\n⚠️  Some dependencies have compatibility issues.")
        print("Please follow the solutions above to resolve them.")
    
    print("=" * 60)
    return all_compatible


def safe_import_with_fallback(module_name: str, fallback_value: Any = None) -> Tuple[Any, bool]:
    """
    Safely import a module with fallback handling.
    
    Parameters
    ----------
    module_name : str
        Name of the module to import
    fallback_value : Any, default None
        Value to return if import fails
    
    Returns
    -------
    Tuple[Any, bool]
        (imported_module_or_fallback, import_success)
    """
    try:
        if module_name == 'pandas':
            import pandas as pd
            return pd, True
        elif module_name == 'numpy':
            import numpy as np
            return np, True
        elif module_name == 'pyensembl':
            from pyensembl import EnsemblRelease
            return EnsemblRelease, True
        else:
            module = __import__(module_name)
            return module, True
    except ImportError:
        return fallback_value, False
    except Exception:
        return fallback_value, False


def get_environment_info() -> Dict[str, str]:
    """
    Get information about the current Python environment.
    
    Returns
    -------
    Dict[str, str]
        Dictionary with environment information
    """
    info = {
        'python_version': sys.version,
        'platform': sys.platform,
    }
    
    # Try to get package versions
    packages = ['numpy', 'pandas', 'matplotlib', 'seaborn']
    
    for package in packages:
        try:
            module = __import__(package)
            info[f'{package}_version'] = getattr(module, '__version__', 'Unknown')
        except ImportError:
            info[f'{package}_version'] = 'Not installed'
        except Exception:
            info[f'{package}_version'] = 'Error getting version'
    
    return info


def suggest_environment_setup() -> str:
    """
    Provide suggestions for setting up a compatible environment.
    
    Returns
    -------
    str
        Environment setup suggestions
    """
    suggestions = """
RECOMMENDED ENVIRONMENT SETUP:

Option 1: Using pip (recommended for most users)
-----------------------------------------------
pip install "numpy>=1.20.0,<2.0.0"
pip install "pandas>=1.3.0,<2.1.0"
pip install "matplotlib>=3.4.0"
pip install "seaborn>=0.11.0"
pip install pyensembl

Option 2: Using conda
--------------------
conda create -n pymut python=3.10
conda activate pymut
conda install "numpy<2.0" "pandas<2.1" matplotlib seaborn
pip install pyensembl

Option 3: Using requirements.txt
-------------------------------
pip install -r requirements.txt

Option 4: If you have NumPy 2.0 issues
-------------------------------------
pip uninstall numpy pandas
pip install "numpy>=1.20.0,<2.0.0"
pip install "pandas>=1.3.0,<2.1.0"

VERIFICATION:
After installation, run:
python -c "import numpy, pandas; print('NumPy:', numpy.__version__); print('Pandas:', pandas.__version__)"
"""
    return suggestions


if __name__ == "__main__":
    # Run compatibility check when module is executed directly
    print_compatibility_report(verbose=True)
    
    if not check_all_dependencies()['numpy'][0]:
        print("\nENVIRONMENT SETUP SUGGESTIONS:")
        print(suggest_environment_setup())