# NettingI
NettingI is a toolbox of radio frequency interference (RFI) excision algorithms, and was designed to work on raw pulsar data from the VEGAS: Versatile Green Bank Telescope Astronomical Spectrometer. This repository is still under development, with the intention of having two more mitigation algorithms implemented: [median absolute deviation](https://openscholarship.wustl.edu/cgi/viewcontent.cgi?article=1005&context=undergrad_etd) and [Cyclostationary Processing](https://cyclostationary.blog/). Read [Chen & Smith 2024](https://www.nrao.edu/students/2024/Reports/ChenKaitlyn.pdf) for more details. 

[Chen & Smith 2024](https://www.nrao.edu/students/2024/Reports/ChenKaitlyn.pdf) found that certain algorithms are more effective for flagging certain types of RFI than others. [Interquartile range mitigation (IQRM)](https://academic.oup.com/mnras/article/510/1/1393/6449380) was good at flagging the Iridium satellite, but not the FAA radar or GPS. [Spectral kurtosis estimator (SK)](https://academic.oup.com/mnrasl/article/406/1/L60/1041152) mitigation was good at flagging the FAA radar but not the Iridium satellite or GPS. 

While these algorithms are intended for real-time RFI mitigation use, some (like IQRM) can also be used on averaged data. 

## Installing ``nettingi``

``nettingi`` works with Python 3.9. To start a conda environment and enable it for Jupyter notebooks, type in bash:

```bash
$ conda create --name [env name] -python=3.9
$ conda activate [env name]
$ pip install ipython
$ pip install jupyter
$ python -m ipykernel install --user --name [env name] --display-name [env name]
```

### From GitHub

To install from github:

```bash
    $ git clone git@github.com:etsmit/nettingi.git
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
After installation, ``nettingi`` can be used in a python environment. The most basic use has default settings for RFI mitigation. There are more complex example uses in the [notebooks folder](https://github.com/etsmit/nettingi/tree/main/notebooks).

### In Jupyter Notebooks
Note that ``filename`` is the file path to the data with RFI you want to mitigate and ``replacement`` is either ``nans`` or ``noise`` depending on what the RFI should be replaced with.
#### IQRM
```
    running = rfimit.rfi_iqrm(filename, ...)
    running.run_all()
```

#### SK
```
    running = rfimit.rfi_sk(filename, ...)
    running.run_all()
```

#### AOF
```
    running = rfimit.rfi_aof(filename, ...)
    running.run_all()
```

## Github Etiquette

This code is maintained and improved by creating your own branch for each feature/bug fix/etc. and opening a Pull Request back to `main`. To keep your code and your PR's up to date, here are the steps:

1. Go to the base directory of the repository and update your main branch
```bash
$ cd nettingi
$ git pull origin main
```

2. Make a new branch and switch to it. For example, I recently added 2-dimensional flagging to the IQRM algorithm. Branch names should be short but informational, with spaces between words as "-".
```bash
$ git branch iqrm-2d
$ git checkout iqrm-2d
```

3. At this point, make any changes to the code you wish. Once you are satisfied with your progress (with any level of granularity) you can make "commits" by adding your changes to the staging area with `git add` and sending a commit with `git commit`. Make sure you are still in the base directory. You can preview your changes before adding them to the staging area with `git diff`.

```bash
$ git diff
$ git add .
$ git commit -m "2D IQRM implemented"
```
4. Once you are satisfied with your local changes, push them to the Github repository.

```bash
$ git push origin iqrm-2d
```

5. Finally, on Github, create a pull request to the `main` branch. Github will automatically check for conflicts, which shouldn't really happen. After this check, you will technically be allowed to merge your branch but it is reccommended that someone reviews the PR before that happens.










