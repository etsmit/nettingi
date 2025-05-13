# Code Breakdown

This is a breakdown of how nettingi is structured, with a guide on how to add or maintain RFI detection algorithms.

## Basic usage and Summary

RFI mitigation on VEGAS pulsar mode baseband files is done by instantiating a class specific to your example algorithm with all of the required inputs. This includes the input file, any parameters used by the algorithm, and other, more global arguments. Your algorithm class will contain all the code it needs to perform RFI detection, as well as the base class which has the functionality to run the full RFI mitigation and manage the output files.

The input data will be taken from `/jetstor/scratch/rfimit/unmitigated/rawdata/` by default. It will output the mitigated data to
`/jetstor/scratch/rfimit/mitigated/rawdata/`, and the "SRDPs" to `/jetstor/scratch/rfimit/mitigated/reduced/`. These SRDPs will contain the averaged data pre- and post- mitigation, the flag file(s), and any other intermediate arrays that are pertinent to the specific algorithm. For example, Spectral Kurtosis also saves the kurtosis values in addition to what's mentioned above.

### VPM Baseband data description

## Core


`core.py` contains the Parent class `mitigateRFI`, which is used under the hood of each of the separate tools. When you instantiate a class to perform a certain RFI algorithm, it will also have the functionality of `mitigateRFI` attached. This includes `run_all()`, which loads the data, performs RFI detection and mitigation, and writes out the mitigated data along with the SRDPs. When you add your own algorithm, you have to add functionality to `mitigateRFI.run_all()` so that it understands how to process data how you want. This is done at (currently) lines 107-137 in `core.py`. An example is shown below, for Spectral Kurtosis:


```python
            if self.det_method == 'SK':
                print('SK mitigation')
                flags_block, ss_sk_block, ms_sk_block = self.SK_detection(data)
                if bi == 0:
                    self.ss_sk_all = ss_sk_block
                    self.ms_sk_all = ms_sk_block
                else:
                    self.ss_sk_all = np.concatenate((self.ss_sk_all, ss_sk_block),axis=1)
                    self.ms_sk_all = np.concatenate((self.ms_sk_all, ms_sk_block),axis=1)
```

At this point in the code, you have the raw channelized voltages corresponding to one block of data, called `data`. If you want SK detection, this code will run `self.SK_detection()` and collect the flags corresponding to that block, as well as two intermediate arrays for the single- and multi-scale SK. Since we make no assumptions about how many blocks there are, we concatenate the results time-wise for all blocks after the first one. The minimum end result of this code is the flagging array, which will be used to perform data replacement/mitigation.

After flagging and mitigation, you do need to add any extra `numpy.save` statements starting on line 216. The code will already save the pre- and post- mitigation averaged spectra and the flagging data, but any intermediate arrays should also be added. For example, SK has to save the single- and multi-scale SK results. There is also a printout to a logfile, which details the output SRDP filenames and serves as an extra log of what parameters have been used.



```python
if self.det_method == 'SK':
            print(f'SS-SK: {self._ss_sk_filename}')
            np.save(self._ss_sk_filename, self.ss_sk_all)
            print(f'MS-SK: {self._ms_sk_filename}')

            np.save(self._ms_sk_filename, self.ms_sk_all)
            #need to add the logging thing
            log = '/data/scratch/SKresults/SK_log.txt'
            os.system(f"""echo "'{self._spect_filename}'...'\n===============================" >> {log}""")
```

`core.py` also contains functions to fine channelize the data (if you have a spectral line target) and load the SRDPs. This is useful if you want to just get the SRDPs and examine them without rerunning the whole mitigation. This analysis code is forthcoming. There is also an entry to run pulsar data reduction, but this function will `pass` for the time being as it requires a different Conda environment. Relevant pulsar scripts are in `pulsar_reduction/`.



## Example algorithm


You will place all the code relating to your algorithm's RFI detection in its own file. Note that this repo has `sk.py`,`iqrmrfi.py`, etc. associated with the other algorithms already in the toolbox. Ideally, you would write an `__init__` function that initializes every argument and parameter the algorithm needs, as well as more general arguments for the RFI mitigation process. You can use one of these `__init__` functions as a springboard to modify for your own uses. 

Notice how `sk.py` is structured. As part of the class initialization, I assign all the arguments to the `self` namespace, so that they are globally available no matter what function we are inside at any given moment. At this step I also derive the required detection thresholds using `SK_thresholds`, which uses `upperRoot` and `lowerRoot`. For the actual RFI detection, I run `SK_detection`, which makes use of `single_scale_SK_EST` and `multi_scale_SK_EST` to find the RFI, and the end result is a flagging array.







