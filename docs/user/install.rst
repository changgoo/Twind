.. _install:

Installation
============

Since **Twind** is a pure Python module, it should be pretty easy to install.
In addition to widely used packages like `scipy <https://www.scipy.org>`_, `numpy <https://numpy.org/>`_,
`matplotlib <https://matplotlib.org>`_,
you'll need `xarray <http://xarray.pydata.org/en/stable/>`_
and `astropy <https://www.astropy.org>`_. Also, :ref:`paper-figures` needs
`seaborn <https://seaborn.pydata.org>`_,
`CMasher <https://cmasher.readthedocs.io>`_, and
`cmocean <https://matplotlib.org/cmocean/>`_.

Github sources
--------------
`pip <http://www.pip-installer.org/>`_ or `conda <https://conda.io>`_ can be used to install ``xarray``
(see `their installation instruction <http://xarray.pydata.org/en/stable/installing.html>`_) and ``astropy``.
Simply, the following command would work.

.. code-block:: bash

    conda install -c conda-forge xarray astropy

Then, you are ready to clone the repository.

.. code-block:: bash

    git clone https://github.com/changgoo/Twind.git your-twind-path

Add path to the source directory using ``sys.path``; e.g.,

.. code-block:: python

    import sys
    sys.path.insert('your-twind-path')

Package managers
----------------

**(Not packaged yet)**

The recommended way to install the stable version of **Twind** is using
`pip <http://www.pip-installer.org/>`_

.. code-block:: bash

    pip install -U twind
..
    or `conda <https://conda.io>`_

    .. code-block:: bash

        conda install -c conda-forge twind
