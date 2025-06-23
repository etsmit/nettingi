# Pulsar reduction scripts


These scripts are meant for analysis of pulsar data (and some stuff for spectral line is in here too) after RFI mitigation. These scripts will NOT run while using nettingi, you must instead activate the GBO pulsar environment. The pulsar software is not maintained by any nettingi authors and is not installed as nettingi dependencies.


## Preamble

In order to activate the pulsar environment, add the following lines to your .bash_profile:

```bash
psrenv() {
    if [[ $(uname -r) =~ .*el6.* ]]
    then
	PSRHOME=/home/pulsar64
	source $PSRHOME/pulsar.bash
    elif [[ $(uname -r) =~ .*el7.* ]]
    then
	PSRHOME=/home/pulsar_rhel7
	source $PSRHOME/pulsar.presto3.bash
    elif [[ $(uname -r) =~ .*el8.* ]]
    then
	PSRHOME=/home/pulsar_rhel8/opt
	source $PSRHOME/pulsar.bash
    else
	PSRHOME=/home/pulsar_rhel8/opt
	source $PSRHOME/pulsar.bash
    fi
    export PATH=$HOME/opt/bin:$PATH
}
```
And then start a new terminal window and type

```bash
psrenv
```

It is not recommended to run your local nettingi conda environment simultaneously with the GBO pulsar environment; there will be clashing dependencies.

## Data Reduction

On the `leibniz` machine, the SRDPs corresponding to the mitigated data can be found at `/jetstor/scratch/rfimit/mitigated/reduced/`. From there, navigate to the folder that matches the output filename pattern you specified in your RFI mitigation. These filenames and paths are all printed out as part of the RFI mitigation, so if you still have that terminal up, you can use those for reference. Keep in mind that the directory structure takes the scan first, then the file number, so you will have to go two directories down to see e.g. "scan 0004, file 0000". 

Once you are in the SRDP directory, you should see a set of files:

 * `.npy`: These correspond to the intermediate files made during RFI mitigation, and are the pre- and post- mitigation spectrograms, the flagging results, and any other intermediate files.
 
 * `.py`: Pulsar reduction files. These will be used to do pulsar analysis. The method prescribed below only requires running a few of these at the top level, but they are mostly all used under-the-hood and can be ran directly with the correct inputs.

 * `.par`: A parfile which is soft-linked from `pulsar_reduction/parfiles/`. This contains all the details we need to accurately dedisperse, fold, and calculate TOA residuals.

 * `.raw`: This is the RFI-mitigated baseband file which undergoes analysis. It is linked from `../../../rawdata/`.

 Once you have activated the `psrenv` and NOT your own `nettingi`-enabled environment, you can run the following script:

 `$ python reduce_raw_data.py -o [output file basename] -d [DM] -p [parfile] [raw baseband file]`

The `output file basename` is the original file's full name without the extension - so this includes everything up to and including the file number. For example, `vegas_60299_76099_B0329+54_0004.0000`. The directory name already has this, so you can copy/paste it from the present working directory in your terminal. The `DM` is the pulsar's dispersion measure. For each of the pulsars we have data on, the dispersion measures are:

* B0329+54: 26.833

* B0355+54: 57.142

* J1713+0747: 15.992196

The parfile is simply the pulsar's name with a `.par` at the end, e.g. `B0329+54.par`. This should automatically be soft-linked into the directory. The final name is the raw baseband data file. This should also already by softlinked, and ends in `.raw`.


The code above will create an `rfifind` mask, search for candidates, and then dedisperse and fold the data.

# Data Analysis

The next step is to record data quality numbers to see how well RFI mitigation improved the data. Run the following:

`python get_results.py -f [base filename] -p [parfile] -d [DM] -t [period]`

The `base filename` is the same as above, all the way out the the file number. The `parfile` and `DM` are the same as above as well. The `period` is the pulsar period, which take on the following values:

* B0329+54: 0.71452

* B0355+54: 0.156384

* J1713+0747: 0.0045701365

This code will run through the following steps. After each one, it will wait for the user to log any results and press enter to continue.

## compare_rfi_masks.py

This will compare the `rfifind` masks on the unmitigated and mitigated data, show a few plots, and print some numbers. Our metric for success is if the percentage flagged in the mitigated data is smaller - this means that the RFI algorithm has removed RFI that `rfifind` would have otherwise found. 

## compare_profiles.py

This will compare the `dspsr`-folded pulse profiles and plot them over each other. The two numbers to record are the relative max height between the mit/unmit datasets, and the relative signal-to-noise ratio. Ideally, we want to see that the pulse profile strength has not been diminished by overflagging (~100%), and that the signal-to-noise ratio has been increased. (>100%).

## get_TOAS.py

This will finish calculating TOA residuals and uncertainties. The residuals and uncertainties should hopefully be lower, although overflagging has been known to artificially lower them. They are useful to record, but don't stand as a metric of success on their own.

## get_paz_TOAS.py

This will run the RFI-zapping command `paz` on the data before re-running TOA residual calculation. On the unmitigated side, this provides a nice goal to aim for our RFI mitigation. On the already-mitigated data, this serves as an extra test to whether or not it's useful to run both our RFI mitigation and `paz` in series.

## do_cands.py

This sifts and cleans the `accelcands` files given by the candidate search, in both the unmitigated/mitigated datasets and the `rfifind` mask/nomask datasets. Record the number of candidates in all 4 cases, and the SNR of the pulsar in all 4 cases. Ideally the number of pulsar candidates is lower in the mitigated dataset, and the SNR of our actual pulsar (at the right period and DM) is higher.












