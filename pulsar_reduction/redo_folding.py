import os
import sys
import copy
import subprocess
import glob
import tempfile
import argparse
import numpy as np
import matplotlib.pyplot as plt
import time
from operator import attrgetter
from pathlib import Path
from presto import sifting
from presto import residuals
from presto.presto import get_baryv
from presto.infodata import infodata


class FullPath(argparse.Action):
    "Extension of argparse.Action class to yield full paths to single files"
    def __call__(self,parser,namespace,values,option_string=None):
        setattr(namespace,self.dest,os.path.abspath(os.path.expanduser(values)))


class FullPathAppend(argparse.Action):
    "Extension to argparse.Action class to yield full paths to groups of files"
    def __call__(self,parser,namespace,values,option_string=None):
        items = [os.path.abspath(os.path.expanduser(value)) for value in values]
        setattr(namespace,self.dest,items)


def get_args():
    "Parse arguments using the excellent package arparse"
    # Create parser and argument groups
    parser = argparse.ArgumentParser()
    general_args = parser.add_argument_group("general argument")
    search_args = parser.add_argument_group("search-mode arguments")
    fold_args = parser.add_argument_group("fold-mode arguments")
    
    # Add general arguments
    general_args.add_argument("-o", "--basenm", 
                              help=("Base name to use for output "
                                    "(default=base name of input files)"))
    general_args.add_argument("-O", "--outdir", default=".",
                              help="Output directory (default=%(default)s")
    general_args.add_argument("--outfile", default=None, 
                              help=("File to which output is directed "
                                    "(default=STDOUT)"))
    general_args.add_argument("--errfile", default=None, 
                              help=("File to which errors are directed "
                                    "(default=STDERR)"))
    general_args.add_argument("-n", "--nchan", type=int, 
                              help="Number of channels to form "
                              "(default=number of coarse voltage channels)")
    general_args.add_argument("-d", "--dm", type=float, default=0.0,
                              help="DM (pc/cc; default=%(default)s)")
    general_args.add_argument("--nbits", type=int, 
                              help="Number of bits in output PSRFITS files "
                              "(default=native number of bits used in input)")

    # Add arguments specific to search-mode data products
    search_args.add_argument("--no-search", action="store_true", 
                             help="Do not produce search-mode data")
    search_args.add_argument("-t", "--tint", type=float, default=40.96e-6,
                             help="Integration time (s; default=%(default)s)")
    search_args.add_argument("--low-dm", 
                             help=("The lowest DM to use when searching "
                                   "(pc/cc; default=dm-20.0)"))
    search_args.add_argument("--dm-step", type=float, default=0.05, 
                             help=("DM step size to use when searching (pc/cc; "
                                   "default=%(default)s)"))
    search_args.add_argument("--ndms", type=int, 
                             help=("Number of DMs to use when searching "
                                   "(default=(dm-low_dm)/dm_step)"))
    search_args.add_argument("--max-harms", type=int, default=16,
                             help=("Maximum number of harmonics to sum when "
                                   "searching (default=%(default)s)"))
    search_args.add_argument("--zmax", type=int, default=50,
                             help=("Maximum acceleration to use when "
                                   "searching (Fourier bins; "
                                   "default=%(default)s)"))
    search_args.add_argument("--no-zero-dm", action="store_true",
                             help="Do not produce zero-DM'd data products")
    search_args.add_argument("--no-topocentric", action="store_true",
                             help="Do not produce topocentric time series")

    # Add arguments specific to fold-mode data products
    fold_args.add_argument("--no-fold", action="store_true", 
                           help="Do not produce fold-mode data")
    fold_args.add_argument("-p", "--parfile", action=FullPath, 
                           help=("TEMPO-compatible parfile to use for "
                                 "folding"))
    fold_args.add_argument("-b", "--nbin", type=int, default=2048,
                           help=("Number of pulse profile bins "
                                 "(default=%(default)s)"))
    fold_args.add_argument("-c", "--cal-psr", action=FullPath, 
                           help="Pulsar calibration data")
    fold_args.add_argument("--fluxcal", action=FullPath, 
                           help="PSRCHIVE fluxcal data")
    fold_args.add_argument("--fluxcal-on", action=FullPath,  
                           help="On-source fluxcal data")
    fold_args.add_argument("--fluxcal-off", action=FullPath, 
                           help="Off-source fluxcal data")
    parser.add_argument("rawfiles", action=FullPathAppend, nargs="+", 
                        help="RAW input files to process")

    return parser.parse_args()


def execute(command,out=sys.stdout,err=sys.stderr,quiet=False,dry_run=False):
    """
    Execute a given command in the shell

    Parameters
    ----------
    command : str
        The command to execute.
    out : {str || file}
        Destination for standard output.  If a string, a file will be opened 
        with that name.
    err : {str || file}
        Destination for standard error.  If a string, a file will be opened
        with that name.
    quiet : {bool}
        Suppress all output.
    dry_run : {bool}
        Do not actually execute the command.

    Returns
    -------
    out : int
        subprocess return code after executing the command
    """
    if not dry_run: 
        process = subprocess.run(command, shell=True, stdout=out, stderr=err,
                                 universal_newlines=True)
        if not quiet:
            print("#####  "+command, file=out)
            out.flush()
            print(process.stdout, file=out)
            out.flush()
            print("\n\n", file=out)
            out.flush()
            # If stderr was redirected the same place as stdout, then 
            # no need to print it as well
            if out is not err: 
                print("#####  "+command, file=err)
                err.flush()
                print(process.stderr, file=err)
                err.flush()
                print("\n\n", file=err)
                err.flush()
        return process.returncode
    else:
        return None


def read_raw_header(infile, line_length=80):
    """
    Create a header dictionary from a GUPPI RAW file.
    
    Parameters
    ----------
    infile : str
        A GUPPI RAW file.
    line_length : {int}
        The length of a GUPPI RAW header line.
    
    Returns
    out : dict
        The header as a python dictionary.
    
    Notes
    -----
    The header is created from the first block of data only.
    """
    header = {}
    with open(infile, "rb") as f:
        while True:
            header_line = f.read(line_length).strip().decode("utf-8")
            if header_line.startswith("END"): 
                break
            else:
                key,value = header_line.split("=")
                header[key.strip()] = value.strip().replace("'","")
    return header


def do_search(infiles,basenm,dm,low_dm,dm_step,ndms,nbits=8,nchan=None,
              tint=40.96e-6,max_harms=16,zmax=50,parfile=None,zero_dm=True,
              make_topocentric=True,outfile=sys.stdout,errfile=sys.stderr):
    """
    Search for pulsars in GUPPI RAW data.

    Parameters
    ----------
    infiles : str || list
        GUPPI RAW file names.
    basenm : str
        The basename to use for output.
    dm : float
        The dispersion measure around which to search (pc/cc).
    low_dm : float
        The low DM to use when searching (pc/cc).
    dm_step : float
        The DM step-size to use when searching (pc/cc).
    ndms : int
        The number of trial DMs to generate
    nbits : {int}
        The number of bits of the output data.
    nchan : {int}
        The number of channels to form from the GUPPI RAW data.  If None, use
        the native number of coarse channels.
    tint : {float}
        The sampling time of the output generated from the GUPPI RAW data.
    max_harms : {int}
        Maximum number of harmonics to sum during searching.
    zmax : {int}
        Maximum Fourier acceleration over which to search (Fourier bins)
    parfile : {str}
        tempo1-compatible parfile to use when folding data.  If None, 
        no data will be folded.
    zero_dm : {bool}
        Analyze the data with zero-DM'ing turned on.
    make_topocentric : {bool}
        Analyze topocentric time series.

    Returns
    -------
    None

    Notes
    -----
    digifits is used to create PSRFITS-format data, which is then searched
    using PRESTO.  Barycentered data is always produced and searched.
    """
    from astropy.io import fits
    # infiles will be assumed to be a list later on
    if type(infiles) is str: infiles = [infiles]
    # Make a dictionary of the function arguments.  This will make formatting
    # the commands a lot simpler
    kwargs = {"infiles":" ".join(infiles),"basenm":basenm,"dm":dm,
              "low_dm":low_dm,"dm_step":dm_step,"ndms":int(ndms),"nbits":nbits,
              "nchan":nchan,"tint":tint,"max_harms":max_harms,"zmax":zmax,
              "parfile":parfile,"zero_dm":zero_dm,
              "make_topocentric":make_topocentric}

    ### DIGIFTS LEAVES A BUNCH OF KEYWORDS OUT AND PRESTO CANNOT PROCESS
    ### THE OUTPUT.  CONVERTING TO FILTERBANK FOR NOW.
    ### # Make search-mode PSRFITS files
    cmd = ("digifits -c -L 599.95324416 -t 40.96e-6 "
           " -p 4 -b {nbits} -D {dm} {infiles}".format(**kwargs))
    ret = execute(cmd, out=outfile, err=errfile)
    tmpfilenms = glob.glob("*.sf")
    tmpfilenms.sort()
    for ii,filenm in enumerate(tmpfilenms):
        newfilenm = "{basenm}_search_{num:04d}.fits".format(
            basenm=kwargs["basenm"],num=ii+1)
        print(newfilenm)
        cmd = "mv {0} {1}".format(filenm,newfilenm)
        ret = execute(cmd, out=outfile, err=errfile)
        # NCHNOFFS is set to "*" by digifits, so we need to set it to zero 
        # before making filterbank file
        with fits.open(newfilenm,"update") as f:
            f[-1].header["NCHNOFFS"] = 0
        # Make a filterbank file
        cmd = ("/data/rfimit/unmitigated/reduced/psrfits2fil.py --sumpols "
               "{0}".format(newfilenm))
        ret = execute(cmd, out=outfile, err=errfile)

    ### cmd = ("digifil -b {nbits} -d 4 -F {nchan}:D -D {dm} -K "
    ###      "-o {basenm}.fil -t 32 {infiles}".format(**kwargs))
    ### ret = execute(cmd, out=outfile, err=errfile)
        
    # Run rfifind
    ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
    cmd = "rfifind -o {basenm} -time 2 *search*.fil".format(**kwargs)
    ret = execute(cmd, out=outfile, err=errfile)

    # Generate bary-centric time series
    ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
    cmd = ("prepsubband -o {basenm}_bary_search_mask "
           "-mask {basenm}_rfifind.mask -lodm {low_dm} "
           "-dmstep {dm_step} -numdms {ndms:d} *search*.fil".format(**kwargs))
    ret = execute(cmd, out=outfile, err=errfile)
    cmd = ("prepsubband -o {basenm}_bary_search_nomask "
           "-lodm {low_dm} "
           "-dmstep {dm_step} -numdms {ndms:d} *search*.fil".format(**kwargs))
    ret = execute(cmd, out=outfile, err=errfile)

    if zero_dm:
        # Also generate barycentric time series using zero-DM'ing
        ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
        cmd = ("prepsubband -zerodm -o {basenm}_bary_zerodm_search_mask "
               "-mask {basenm}_rfifind.mask -lodm {low_dm} "
               "-dmstep {dm_step} -numdms {ndms:d} *search*.fil".format(
                   **kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
        # Also generate barycentric time series using zero-DM'ing
        ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
        cmd = ("prepsubband -zerodm -o {basenm}_bary_zerodm_search_nomask "
               "-lodm {low_dm} "
               "-dmstep {dm_step} -numdms {ndms:d} *search*.fil".format(
                   **kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
    
    if make_topocentric:
        # Also generate topocentric time series
        ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
        cmd = ("prepsubband -nobary -o {basenm}_topo_search_mask "
               "-mask {basenm}_rfifind.mask -lodm {low_dm} "
               "-dmstep {dm_step} -numdms {ndms:d} *search*.fil".format(
                   **kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
        cmd = ("prepsubband -nobary -o {basenm}_topo_search_nomask "
               "-lodm {low_dm} "
               "-dmstep {dm_step} -numdms {ndms:d} *search*.fil".format(
                   **kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
        
        if zero_dm:
            # Also generate topocentric time series using zero-DM'ing
            ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
            cmd = ("prepsubband -nobary -zerodm -o {basenm}_topo_zerodm_search_mask "
                   "-mask {basenm}_rfifind.mask -lodm {low_dm} "
                   "-dmstep {dm_step} -numdms {ndms:d} "
                   "*search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("prepsubband -nobary -zerodm -o {basenm}_topo_zerodm_search_nomask "
                   "-lodm {low_dm} "
                   "-dmstep {dm_step} -numdms {ndms:d} "
                   "*search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)

    # FFT all the time series ending in .dat
    cmd = "ls *search*.dat | xargs -L 1 realfft"
    ret = execute(cmd, out=outfile, err=errfile)
    
    # Get the barycentric velocity during the observation and then run 
    # accelsearch on all Fourier spectra ending in .fft
    try:
        i = infodata("{basenm}_rfifind.inf".format(**kwargs))
        baryv = get_baryv(i.RA, i.DEC, i.epoch, i.N*i.dt, obs=i.telescope)
    except FileNotFoundError:
        baryv = 0.0
    cmd = ("ls *search*.fft | xargs -L 1 accelsearch -numharm 16 -zmax 50 "
           "-baryv %f "%baryv)
    ret = execute(cmd, out=outfile, err=errfile)
    
    # Run single-pulse search on all time series ending in .dat
    cmd = "single_pulse_search.py -t 2.0 *search_mask*.dat"
    ret = execute(cmd, out=outfile, err=errfile)
    cmd = "single_pulse_search.py -t 2.0 *search_nomask*.dat"
    ret = execute(cmd, out=outfile, err=errfile)
    
    # Since we want to sift candidates from the bary/topocentric time series
    # with/without zero-DM'ing separately, make some glob strings that will be
    # used when sifting
    glob_strings = ["bary_search_mask","bary_search_nomask"]
    if zero_dm: 
        glob_strings.append("bary_zerodm_search_mask","bary_zerodm_search_nomask")
    if make_topocentric: 
        glob_strings.append("topo_search_mask","topo_search_nomask")
        if zero_dm: 
            glob_strings.append("topo_zerodm_search_mask","topo_zerodm_search_nomask")
    
    # Now sift through the candidates
    for gs in glob_strings:
        cands = sifting.read_candidates(
            glob.glob("*"+gs+"_search*ACCEL_50"))
        if len(cands): 
            cands = sifting.remove_duplicate_candidates(cands)
            cands.sort(key=attrgetter("sigma"), reverse=True)
            sifting.write_candlist(
                cands,"{}_{}.accelcands".format(basenm,gs))
        else:
            Path("{}_{}.accelcands".format(basenm,gs)).touch()
        
        # Create single-pulse plots for each group of time series
        cmd = "single_pulse_search.py -t 5.0 *%s_search*.singlepulse"%gs
        ret = execute(cmd, out=outfile, err=errfile)

    if parfile is not None:
        # Fold the full-resolution data
        ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
        cmd = ("prepfold -noxwin -o {basenm} -timing {parfile} "
               "*search*.fil".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)

        # Create and fold barycentric time series at the exact DM of the 
        # pulsar
        cmd = ("prepdata -o {basenm}_bary_mask_DM{dm:.2f} -dm {dm} "
               "-mask {basenm}_rfifind.mask *search*.fil".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
        cmd = ("prepdata -o {basenm}_bary_nomask_DM{dm:.2f} -dm {dm} "
               "*search*.fil".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
               
        cmd = ("prepfold -noxwin -timing {parfile} "
               "{basenm}_bary_DM{dm:.2f}.dat".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)

        if zero_dm:
            # Also fold the full-resolution data using zero-DM'ing
            ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
            cmd = ("prepfold -noxwin -zerodm -o {basenm}_zerodm -timing "
                   "{parfile} *search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            
            # Also create and fold barycentric time series at the exact DM
            # of the pulsar using zero-DM'ing
            cmd = ("prepdata -zerodm -o {basenm}_bary_zerodm_mask_DM{dm:.2f} "
                   "-dm {dm} -mask {basenm}_rfifind.mask "
                   "*search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("prepdata -zerodm -o {basenm}_bary_zerodm_nomask_DM{dm:.2f} "
                   "-dm {dm} "
                   "*search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)

            cmd = ("prepfold -noxwin -timing {parfile} "
                   "{basenm}_bary_zerodm_DM{dm:.2f}.dat".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)

        if make_topocentric:
            # Also create and fold topocentric time series at the exact DM
            # of the pulsar
            ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
            cmd = ("prepdata -nobary -o {basenm}_topo_mask_DM{dm:.2f} -dm {dm} "
                   "-mask {basenm}_rfifind.mask *search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("prepdata -nobary -o {basenm}_topo_nomask_DM{dm:.2f} -dm {dm} "
                   "*search*.fil".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            
            cmd = ("prepfold -noxwin -timing {parfile} "
                   "{basenm}_topo_DM{dm:.2f}.dat".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)

            if zero_dm:
                # Also create and fold topocentric time series at the exact
                # DM of the pulsar using zero-DM'ing
                ### TODO: CHANGE TO .fits WHEN DIGIFITS IS WORKING
                cmd = ("prepdata -zerodm -nobary "
                       "-o {basenm}_topo_zerodm_mask_DM{dm:.2f} -dm {dm} "
                       "-mask {basenm}_rfifind.mask "
                       "*search*.fil".format(**kwargs))
                ret = execute(cmd, out=outfile, err=errfile)
                cmd = ("prepdata -zerodm -nobary "
                       "-o {basenm}_topo_zerodm_nomask_DM{dm:.2f} -dm {dm} "
                       " "
                       "*search*.fil".format(**kwargs))
                ret = execute(cmd, out=outfile, err=errfile)

                cmd = ("prepfold -noxwin -timing {parfile} "
                       "{basenm}_topo_zerodm_DM{dm:.2f}.dat".format(**kwargs))
                ret = execute(cmd, out=outfile, err=errfile)


def do_fold(infiles,basenm,parfile,dm,nbin=2048,nbits=8,nchan=None,cal_psr=None,
            fluxcal=None,fluxcal_on=None,fluxcal_off=None,outfile=sys.stdout,
            errfile=sys.stderr):
    """
    Fold and, optionally, calibrate GUPPI RAW data.

    Parameters
    ----------
    infiles : str || list
        File name(s) of GUPPI RAW data.
    basenm : str
        The basename to use for output.
    parfile : str
        tempo1-compatible parfile to use when folding data.
    dm : float
        The dispersion measure around which to search (pc/cc).
    nbin : {int}
        Number of pulse profile bins to use when folding
    nbits : {int}
        The number of bits of the output data.
    nchan : {int}
        The number of channels to form from the GUPPI RAW data.  If None, use
        the native number of coarse channels.
    cal_psr : {str}
        File name of GUPPI RAW pulsar calibration data.
    fluxcal : {str}
        File name of PSRCHIVE fluxcal data.
    fluxcal_on : {str}
        File name of GUPPI RAW on-source fluxcal data.
    fluxcal_off : {str}
        File name of GUPPI RAW off-source fluxcal data.

    Returns
    -------
    None
    """
    # infiles will be assumed to be a list later on
    if type(infiles) is str: infiles = [infiles]
    # Make a dictionary of the function arguments.  This will make formatting
    # the commands a lot simpler
    kwargs = {"infiles":" ".join(infiles),"basenm":basenm,"dm":dm,
              "parfile":parfile,"nbin":nbin,"nbits":nbits,"nchan":nchan,
              "cal_psr":cal_psr,"fluxcal":fluxcal,"fluxcal_on":fluxcal_on,
              "fluxcal_off":fluxcal_off}

    # dspsr requires a tempo2-formatted parfile, so make one from the 
    # tempo1-formatted parfile that was given
    t2parfile = parfile.split(".")[0] + ".t2.par"
    kwargs["t2parfile"] = t2parfile
    cmd = "tempo2 -gr transform {parfile} {t2parfile}".format(**kwargs)
    ret = execute(cmd, out=outfile, err=errfile)

    # Fold the raw data using dspsr
    if nchan is not None:
        cmd = ("dspsr -t 8 -O {basenm}_fold -a psrfits -e fits -F {nchan}:D "
               "-D {dm} -d 4 -b {nbin} -E {t2parfile} -A -nsub 42 -L 1.0 "
               "-U 8192 {infiles}".format(**kwargs))
    else:
        cmd = ("dspsr -t 8 -O {basenm}_fold -a psrfits -e fits "
               "-D {dm} -d 4 -b {nbin} -E {t2parfile} -A -nsub 42 -L 1.0 "
               "-U 8192 {infiles}".format(**kwargs))
    print(cmd)
    ret = execute(cmd, out=outfile, err=errfile)
    cmd = ("psredit -m -c freq=1500.0,rcvr:hand=-1,be:phase=-1,rcvr:sa=-45deg "
           "*fold*")
    ret = execute(cmd, out=outfile, err=errfile)
    cmd = ("psradd -o {basenm}_fold_sum.fits -E {t2parfile} "
           "*fold_0*.fits*".format(**kwargs))
    ret = execute(cmd, out=outfile, err=errfile)

    # If we have raw calibration data taken at the position of the pulsar, 
    # fold it as well
    if cal_psr is not None:
        cal_basenm = os.path.basename(cal_psr.split(".")[0])
        kwargs["cal_basenm"] = cal_basenm
        cmd = ("dspsr -t 8 -O {cal_basenm}_cal -a psrfits -e fits "
               "-F {nchan}:D -D 0.0 -K -d 4 -b {nbin} -c 0.04 -p 0 -A "
               "-nsub 32 -L 1.0 -U 8192 {cal_psr}".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
        cmd = "pam -m --type PolnCal {cal_basenm}_cal_0001.fits".format(
            **kwargs)
        ret = execute(cmd, out=outfile, err=errfile)
        cmd = ("psredit -m -c freq=1500.0,rcvr:hand=-1,be:phase=-1,"
               "rcvr:sa=-45deg {cal_basenm}_cal_0001.*".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)
            
        # If don't already have a fluxcal, but we have on/off source 
        # fluxcal data, fold it as well
        if (fluxcal is None and fluxcal_on is not None and 
            fluxcal_off is not None):
            fluxcal_on_basenm = os.path.basename(fluxcal_on.split(".")[0])
            kwargs["fluxcal_on_basenm"] = fluxcal_on_basenm
            cmd = ("dspsr -t 8 -O {fluxcal_on_basenm}_cal -a psrfits "
                   "-e fits -F {nchan}:D -D 0.0 -K -d 4 -b 2048 -c 0.04 "
                   "-p 0 -A -nsub 32 -L 1.0 -U 8192  "
                   "{fluxcal_on}".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("pam -m --type FluxCalOn "
                   "{fluxcal_on_basenm}_cal_0001.fits".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("psredit -m -c freq=1500.0,rcvr:hand=-1,be:phase=-1,"
                   "rcvr:sa=-45deg "
                   "{fluxcal_on_basenm}_cal_0001.fits".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)

            fluxcal_off_basenm = os.path.basename(fluxcal_off.split(".")[0])
            kwargs["fluxcal_off_basenm"] = fluxcal_off_basenm
            cmd = ("dspsr -t 8 -O {fluxcal_off_basenm}_cal -a psrfits "
                   "-e fits -F {nchan}:D -D 0.0 -K -d 4 -b {nbin} -c 0.04 "
                   "-p 0 -A -nsub 32 -L 1.0 -U 8192 "
                   "{fluxcal_off}".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("pam -m --type FluxCalOff "
                   "{fluxcal_off_basenm}_cal_0001.fits".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
            cmd = ("psredit -m -c freq=1500.0,rcvr:hand=-1,be:phase=-1,"
                   "rcvr:sa=-45deg "
                   "{fluxcal_off_basenm}_cal_0001.fits".format(**kwargs))
            ret = execute(cmd, out=outfile, err=errfile)
                
            # Create a pac database
            with open("calfiles.txt", "w") as f:
                f.write(cal_basenm+"_cal_0001.fits\n")
                f.write(fluxcal_on_basenm+"_cal_0001.fits\n")
                f.write(fluxcal_off_basenm+"_cal_0001.fits\n")
                for filenm in glob.glob("*fold*.fits"):
                    f.write(filenm+"\n")
            cmd = "pac -w -W calfiles.txt"
            ret = execute(cmd, out=outfile, err=errfile)
            
            # Create a fluxcal
            cmd = "fluxcal -d database.txt"
            ret = execute(cmd, out=outfile, err=errfile)
        
        # If we do have fluxcal data, just go straight into create the 
        # pac database
        elif fluxcal is not None:
            with open("calfiles.txt", "w") as f:
                f.write(cal_basenm+"_cal_0001.fits\n")
                f.write(fluxcal+"\n")
                for filenm in glob.glob("*fold*.fits"):
                    f.write(filenm+"\n")
            cmd = "pac -w -W calfiles.txt"
            ret = execute(cmd, out=outfile, err=errfile)

        # If we have no fluxcal data at all, then create a pac database
        # with only the pulsar and pulsar cal data
        else:
            with open("calfiles.txt", "w") as f:
                f.write(cal_basenm+"cal_0001.fits\n")
                for filenm in glob.glob("*fold*.fits"):
                    f.write(filenm+"\n")
            cmd = "pac -w -W calfiles.txt"
            ret = execute(cmd, out=outfile, err=errfile)
    
        # Calibrate using the database we made in the previous section
        cmd = "pac -FcTd database.txt *fold*.fits"
        ret = execute(cmd, out=outfile, err=errfile)
        
        # Sum the calibrated data
        cmd = ("psradd -o {basenm}_fold_sum.calib -E {t2parfile} "
               "*fold*.calib*".format(**kwargs))
        ret = execute(cmd, out=outfile, err=errfile)

        # Make a fully scrunched version for generating a standard template
        cmd = "pam -e scr -FTp %s_fold_sum.calib*"%basenm
        ret = execute(cmd, out=outfile, err=errfile)
        
        # Make a frequency-scunched version for generating TOAs
        cmd = "pam -e fscr --setnchn 4 *fold*.calib*"
        ret = execute(cmd, out=outfile, err=errfile)            
            
    # If we have no calibration data at all, then just process the uncalibrated
    # pulsar data
    else:
        # Sum the uncalibrated data
        #cmd = ("psradd -o {basenm}_fold_sum.fits -E {t2parfile} "
        #       "*fold*.fits".format(**kwargs))
        #ret = execute(cmd, out=outfile, err=errfile)

        # Make a fully scrunched version for generating a standard template
        cmd = "pam -e scr -FTp {basenm}_fold_sum.fits".format(**kwargs)
        ret = execute(cmd, out=outfile, err=errfile)
        
        # Make a frequency-scunched version for generating TOAs
        cmd = "pam -e fscr --setnchn 8 *fold*.fits"
        ret = execute(cmd, out=outfile, err=errfile)            
    
    # Make a smoothed template
    cmd = "psrsmooth -W -t UD8 {basenm}_fold_sum.scr".format(**kwargs)
    ret = execute(cmd, out=outfile, err=errfile)
    
    # Rename the smooth template
    cmd = "mv {basenm}_fold_sum.scr.sm {basenm}_fold_sum.std".format(**kwargs)
    ret =execute(cmd, out=outfile, err=errfile)
        
    # Generate TOAs from the frequency-scrunched data
    cmd = ("pat -F -f princeton -s {basenm}_fold_sum.std *fold*.fscr "
           "> toas.tim".format(**kwargs))
    ret = execute(cmd, out=outfile, err=errfile)
        
    # Run tempo
    cmd = "tempo -f {parfile} toas.tim".format(**kwargs)
    ret = execute(cmd, out=outfile, err=errfile)

    # Make a plot of the pre-fit residuals
    r = residuals.read_residuals()
    plt.errorbar(r.bary_TOA,r.prefit_sec*1e3,yerr=r.uncertainty*1e3,fmt=".")
    plt.xlabel("MJD")
    plt.ylabel("Residual (ms)")
    plt.title("Pre-fit Residuals for {basenm}".format(**kwargs))
    plt.savefig("{basenm}_residuals.png".format(**kwargs), dpi=300,
                bbox_inches="tight", pad_inches=0.05)
    

def main():
    "Main program"
    # Get the arguments

    start_time = time.time()
    
    args = get_args()
    
    # Switch to specified outdir
    if not os.path.isdir(args.outdir): os.mkdir(args.outdir)
    os.chdir(args.outdir)

    # Parse the arguments that dont't have simple defaults
    basenm = (args.basenm if args.basenm is not None 
              else os.path.basename(args.rawfiles[0].split(".")[0]))

    if args.outfile is not None:
        if args.outfile == "STDOUT" or args.outfile == "stdout":
            outfile = sys.stdout
        elif args.outfile == "STDERR" or args.outfile == "stderr":
            outfile = sys.stderr  # Why would anyone do this?
        else:
            outfile = open(args.outfile, "w")
    else:
        outfile = sys.stdout

    if args.errfile is not None:
        if args.errfile == "STDERR" or args.outfile == "stderr":
            outfile = sys.stderr
        elif args.outfile == "STDOUT" or args.outfile == "stdout":
            outfile = sys.stdout
        else:
            errfile = open(args.errfile, "w")
    else:
        errfile = sys.stdout

    if args.nchan is None:
        header = read_raw_header(args.rawfiles[0])
        nchan = int(header["OBSNCHAN"])
    else:
        nchan = args.nchan

    if args.nbits is None:
        header = read_raw_header(args.rawfiles[0])
        if "NBITSREQ" in header:
            nbits = int(header["NBITSREQ"])
        else:
            nbits = int(header["NBITS"])
    else:
        nbits = args.nbits

    low_dm = args.low_dm if args.low_dm is not None else max(0,args.dm-20.0)

    if args.ndms is None:
        ndms = 2*(args.dm-low_dm)//args.dm_step
    else:
        ndms = args.ndms

    if not args.no_fold and args.parfile is None:
        print("ERROR: A parfile is required to produce fold-mode data")
        sys.exit(-1)
    
    # Search the data
    #if not args.no_search:
    #    do_search(args.rawfiles,basenm,args.dm,low_dm,args.dm_step,ndms,nbits,
    #              nchan,args.tint,args.max_harms,args.zmax,args.parfile,
    #              not args.no_zero_dm, not args.no_topocentric,outfile,errfile)

    # Fold the data
    if not args.no_fold:
        print('Folding...')
        do_fold(args.rawfiles,basenm,args.parfile,args.dm,args.nbin,nbits,
                nchan,args.cal_psr,args.fluxcal,args.fluxcal_on,
                args.fluxcal_off,outfile,errfile)

    outfile.close()
    errfile.close()

    stop_time = time.time()
    main_dur = np.around((stop_time - start_time)/60)
    print(f'Duration: {main_dur} min')

if __name__ == "__main__":
    main()  # Do it!
