*******************
Quickstart
*******************
After installation, ``nettingi`` can be used in a python environment. The most basic use has default settings for RFI mitigation. Since `pub <>`_ found that XXXXXX mitigation method is the most universal for Green Bank Telescope data, it will default use that method with the optimized settings. 

IQRM
====================
If you wanted to use the IQRM mitigation method this is how you'd use it. Note that ``filename`` is the file path to the data with RFI you want to mitigate and ``replacement`` is either ``nans`` or ``noise`` depending on what the RFI should be replaced with.

.. code::

    import nettingi
    nettingi.rfi_iqrm(filename, replacement)

