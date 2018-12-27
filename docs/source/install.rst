.. _install:

Installation
============

MAnorm is written in Python and supports Python 2.7 or Python 3.4+. It can be obtained easily from PyPI_ or Bioconda_,
the commands below show how to install the latest release of MAnorm.

.. warning::
   Starting from v1.2.0, MAnorm no longer support Window platform.

Install from PyPI
-----------------
The latest release of MAnorm is available at PyPI_, you can install via ``pip``:

.. code-block:: shell

    $ pip install -U manorm

Install with conda
------------------

You can also install MAnorm with conda_ through Bioconda_ channel:

.. code-block:: shell

   $ conda install -c bioconda manorm

Install from source code
------------------------

It's highly recommended to install MAnorm from PyPI_  or Bioconda_.
If you prefer to install it from source code, please read the following steps:

The source code of MAnorm is hosted on GitHub_, and ``pip`` is required for installation.

First, use ``git`` to clone the repository of MAnorm:

.. code-block:: shell

   $ git clone https://github.com/shao-lab/MAnorm.git

Then, install MAnorm in the source directory:

.. code-block:: shell

   $ cd MAnorm
   $ pip install .

Alternatively, you can download the release from the `GitHub release page`_. Once extracted the source code,
MAnorm can be installed with ``pip`` as before.

Galaxy Installation
-------------------
MAnorm is available on Galaxy_, you can incorporate MAnorm into your own Galaxy instance.

Please search and install MAnorm via the `Galaxy Tool Shed`_.

.. _PyPI: https://pypi.python.org/pypi/MAnorm
.. _Bioconda: https://bioconda.github.io
.. _conda: https://conda.io
.. _GitHub: https://github.com/shao-lab/MAnorm
.. _GitHub release page: https://github.com/shao-lab/MAnorm/releases
.. _Galaxy: https://galaxyproject.org
.. _`Galaxy Tool Shed`: https://toolshed.g2.bx.psu.edu/view/haydensun/manorm
