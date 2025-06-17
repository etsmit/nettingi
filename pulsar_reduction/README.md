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

It is not recommended to run your local nettingi conda environment simultaneously with the GBO pulsar environment; there will be clashing python packages.

## Data Reduction and Analysis

bogos binted

On the `leibniz` machine, the SRDPs corresponding to the mitigated data can be found at `/jetstor/scratch/rfimit/mitigated/reduced/`. From there, navigate to the folder that matches the output filename pattern you specified in your RFI mitigation. These filenames and paths are all printed out as part of the RFI mitigation, so if you still have that terminal up, you can use those for reference. Keep in mind that the directory structure takes the scan first, then the file number, so you will have to go two directories down to see e.g. "scan 0004, file 0000". 

Once you are in the SRDP directory, you should see a set of files:

 * `.npy`: These correspond to the intermediate files made during RFI mitigation, and are the pre- and post- mitigation spectrograms, the flagging results, and any other intermediate files.
 
 *`.py`: Pulsar reduction files. These will be used to do pulsar analysis. The method presribed below only requires running a few of these at the top level, but they are mostly all used under-the-hood and can be ran directly with the correct inputs.

 *`.raw`: This is the RFI-mitigated baseband file which undergoes analysis. It is linked from `../../../rawdata/`.

 Once you have activate the `psrenv` and NOT your own `nettingi`-enabled environment, you can run the following script:

 `$ python pu









