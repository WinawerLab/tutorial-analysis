# Experiment

Run the experiment on the scanner's Mac computer. You should not need
to install anything (assuming you're using the winawerlab account),
but you will need to make sure this directory is somewhere on the
computer that you can find and you'll need to download all of
the [stimuli files](https://osf.io/jm42p/) (the files are too large to
include in the Github repository). Place those files into the
`data/stimuli` folder.

Type

```
python psychopy_example.py -h
```

to see the help string and
 
```
python psychopy_example.py
```
 
to run the experiment. There are several arguments you can set, but
for the purposes of gathering the tutorial data, you can trust the
defaults. This will do two runs, waiting for the scanner trigger (the
"5" key) to start each run. It will save the outputs (which contain
the timing information, behavioral results, etc.) as hdf5 files in the
`data/raw_behavioral` directory, which you should not commit to this
repository.

Note that the `winawerlab` environment you set up does not include the
`pylink` package, which is necessary to get the data from the
eye-tracker in this script. This is because `pylink` is very annoying
to set up (see the documentation
from [pygaze](https://www.pygaze.org/installation/)
and [psychopy](http://www.psychopy.org/api/hardware/pylink.html) for
more info) and because you won't need it: you (presumably) will not be
gathering eye-tracker data using your personal computer and it's
already installed and working on Apple machine at the scanner.
