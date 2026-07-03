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
from ediff import acknowledgement

# Check whether the package is being executed as a script.
if __name__ == '__main__':

    # No command-line arguments: print the acknowledgement.
    if len(sys.argv) == 1:
        acknowledgement()

    # Process command-line arguments.
    else:

        # python -m ediff gui
        if sys.argv[1] == 'gui':
            print('TODO: Launch the EDIFF GUI.')
            print('------')
            print('* Call a function from the ediff.gui subpackage.')
            print('* Under development; to be implemented by Pavlina.')
