'''
ediff.__main__ = top-level execution of the EDIFF package
=========================================================

Technical notes
---------------
* Modules use ``__main__``; packages use ``__main__.py``.
* For more information, see:
  https://docs.python.org/3/library/__main__.html

Purpose of __main__.py
----------------------
* ``python -m ediff``
    - print the (required) acknowledgement to stdout
    - the acknowledgement function is defined in the sibling ``__init__.py``

* ``python -m ediff gui``
    - launch the EDIFF GUI directly
    - the GUI interface is implemented in the ``ediff.gui`` subpackage
'''

import sys
from ediff import acknowledgement, acknowledgementGUI


def _print_dedicated_env_advice(err):
    print('The reliable fix is usually a dedicated environment, since')
    print('PyQt6 version conflicts with other installed tools (Spyder,')
    print('Qt Designer helpers, etc.) can\'t be solved by picking one')
    print('"correct" PyQt6 version — there may not be one that satisfies')
    print('everything installed in this environment at once.')
    print()
    print('    conda create -n ediff-gui python=3.11')
    print('    conda activate ediff-gui')
    print('    pip install ediff[gui]')
    print()
    print(f'(original error: {type(err).__name__}: {err})')
    
    
# Check whether the package is being executed as a script.
if __name__ == '__main__':

    # No command-line arguments: print the acknowledgement.
    if len(sys.argv) == 1:
        acknowledgement()

    # Process command-line arguments.
    else:
        # python -m ediff gui
        if sys.argv[1] == 'gui':
            # Imported here, not at module level, so that plain
            # `python -m ediff` (no arguments) never requires the GUI's
            # dependencies (PyQt6, matplotlib, opencv) to be installed —
            # only actually requesting the GUI does.
            try:
                from ediff.gui import main as gui_main

                acknowledgementGUI()
                
            except ImportError as err:
                print('Could not start the EDIFF GUI.')
                print('------')
                print('The GUI dependencies may not be installed. Run:')
                print('    pip install ediff[gui]')
                print(f'(original error: {err})')
                sys.exit(1)
            except Exception as err:
                # PyQt6 (or a related package) is present but not working
                # correctly — this is what a version conflict with another
                # PyQt6-based tool typically looks like: not a clean
                # ImportError, but some other failure part-way through
                # loading Qt bindings.
                print('Could not start the EDIFF GUI.')
                print('------')
                _print_dedicated_env_advice(err)
                sys.exit(1)

            try:
                gui_main()
            except Exception as err:
                # Same idea, but for failures that only surface once the
                # GUI actually tries to start (e.g. QApplication creation
                # itself failing due to a broken Qt platform plugin).
                print('The EDIFF GUI failed to start.')
                print('------')
                _print_dedicated_env_advice(err)
                sys.exit(1)