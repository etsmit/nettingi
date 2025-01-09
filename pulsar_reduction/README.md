# Pulsar reduction scripts


These scripts are meant for analysis of pulsar data (and some stuff for spectral line is in here too) after RFI mitigation. These scripts will NOT run while using nettingi, you must instead activate the GBO pulsar environment. The pulsar software is not maintained by any nettingi authors and is not installed as nettingi dependencies.

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

