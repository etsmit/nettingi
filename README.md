# nettingi
## Installing ``nettingi``

``nettingi`` works with Python 3.9.




### From GitHub

To install from github:

```bash
    $ git clone git@github.com:kaitlynchen1/nettingi.git
    $ cd nettingi
    $ pip install -e .
```
### Dependencies
``nettingi`` requires multiple packages that will be installed automatically if using ``pip`` or manually after cloning the repo and being in the repo directory with:

```bash
    $ pip install -r requirements.txt
```

## Quickstart
After installation, ``nettingi`` can be used in terminal or a python environment. The most basic use has default settings for RFI mitigation. Since [pub]() found that XXXXXX mitigation method is the most universal for Green Bank Telescope data, it will default use that method with the optimized settings. 

### In terminal
If using terminal, this is how you would use it. Note that ``filename`` is the file path to the data with RFI you want to mitigate and ``replacement`` is either ``nans`` or ``noise`` depending on what the RFI should be replaced with.

```bash
    $ python mitigateRFI -i filename -r replacement
```
