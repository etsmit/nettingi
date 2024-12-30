# NettingI
NettingI is a toolbox of radio frequency interference (RFI) excision algorithms, and was designed to work on raw pulsar data from the VEGAS: Versatile Green Bank Telescope Astronomical Spectrometer. This repository is still under development, with the intention of having two more mitigation algorithms implemented: [median absolute deviation](https://openscholarship.wustl.edu/cgi/viewcontent.cgi?article=1005&context=undergrad_etd) and [Cyclostationary Processing](https://cyclostationary.blog/). Read [Chen & Smith 2024](https://www.nrao.edu/students/2024/Reports/ChenKaitlyn.pdf) for more details. 

[Chen & Smith 2024](https://www.nrao.edu/students/2024/Reports/ChenKaitlyn.pdf) found that certain algorithms are more effective for flagging certain types of RFI than others. [Interquartile range mitigation (IQRM)](https://academic.oup.com/mnras/article/510/1/1393/6449380) was good at flagging the Iridium satellite, but not the FAA radar or GPS. [Spectral kurtosis estimator (SK)](https://academic.oup.com/mnrasl/article/406/1/L60/1041152) mitigation was good at flagging the FAA radar but not the Iridium satellite or GPS. 

While these algorithms are intended for real-time RFI mitigation use, some (like IQRM) can also be used on averaged data. 

## Installing ``nettingi``

``nettingi`` works with Python 3.9.

### From GitHub

To install from github:

```bash
    $ git clone git@github.com:kaitlynchen1/nettingi.git
    $ cd nettingi
    $ pip install -e .
```

### `pip` installation: coming soon

``nettingi`` will be installable with ``pip``.  The packaged code will be hosted on [PyPi](https://pypi.org/user/kaitlynchen).


```bash
    $ pip install nettingi
```


### Dependencies
``nettingi`` requires multiple packages that will be installed automatically if using ``pip`` or manually after cloning the repo and being in the repo directory with:

```bash
    $ pip install -r requirements.txt
```

## Quickstart
After installation, ``nettingi`` can be used in a python environment. The most basic use has default settings for RFI mitigation. There are more complex example uses in the [notebooks folder](https://github.com/kaitlynchen1/nettingi/blob/kait-dev/notebooks/running.ipynb).

### In Jupyter Notebooks
Note that ``filename`` is the file path to the data with RFI you want to mitigate and ``replacement`` is either ``nans`` or ``noise`` depending on what the RFI should be replaced with.
#### IQRM
```
    rfimit = importlib.import_module('rfi-mitigation')
    running = rfimit.rfi_iqrm(filename)
    running.run_all()
```

#### SK
```
    rfimit = importlib.import_module('rfi-mitigation')
    running = rfimit.rfi_sk(filename)
    running.run_all()
```

#### AOF
```
    rfimit = importlib.import_module('rfi-mitigation')
    running = rfimit.rfi_aof(filename)
    running.run_all()
```
