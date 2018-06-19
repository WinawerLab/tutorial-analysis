#!/usr/bin/python
"""example psychopy script, based on my spatial_frequency_preferences experiment

You must first create your stimuli (as 2d arrays) and save them as a 3d numpy array (e.g.,
stimuli.npy). We also assume you've created and saved separate numpy arrays that give the indices
into your stimulus array, specifying the order you want to present them in.

For example, let's assume stimuli.npy contains 10 stimuli, each 1080 by 1080 pixels. Then
stimuli.shape is (10, 1080, 1080). You have 5 runs in your experiment and your subject name is
sub-01, so you've created sub-01_run00.npy, sub-01_run01.npy, sub-01_run02.npy, sub-01_run03.npy,
sub-01_run04.npy, each of which contains an array with a permutation of the numbers 0 through
9. This script will load in the stimuli and each of those run indices, then present them to the
subject in 5 runs, presenting all stimuli in each run, in the order specified by the run indices.

We also add some blank stimuli to the beginning and end of each run.

As an example task, we'll present a digit stream at fixation (alternating black and white), where
the subject must perform a one-back task.

This will then save out the results as an hdf5 file.

A lot of this script wrapper stuff for the purpose of making your life easier (e.g., so you don't
have to start this many times during one scanning session, so you don't accidentally save over some
output). If you want to see the core of this, how to display things using pscyhopy, update them
correctly, and record any button presses, look into the `run` function.

"""

import argparse
import h5py
import glob
import os
import datetime
import numpy as np
from psychopy import visual, core, event
from psychopy.tools import imagetools
from scipy import misc as smisc


def _create_blanks(blank_sec_length, on_msec_length, off_msec_length, stimulus_shape, blank_loc):
    """create the blank stimuli

    stimulus_shape: 2-tuple of ints. specifies the shape of a single stimulus.
    """
    nblanks = blank_sec_length / ((on_msec_length + off_msec_length) / 1000.)
    if nblanks != int(nblanks):
        raise Exception("Because of your timing ({loc}_blank_sec_length, on_msec_length, and "
                        "off_msec_length), I can't show blanks for the {loc} {length:.02f} seconds"
                        ". {loc}_blank_sec_length must be a multiple of on_msec_length+"
                        "off_msec_length!".format(loc=blank_loc, length=blank_sec_length))
    nblanks = int(nblanks)
    blanks = smisc.bytescale(np.zeros((nblanks, stimulus_shape[0], stimulus_shape[1])),
                             cmin=-1, cmax=1)
    return nblanks, blanks


def _set_params(stim_path, idx_path, on_msec_length=300, off_msec_length=200,
                fix_button_prob=1/6., final_blank_sec_length=16, init_blank_sec_length=16,
                size=[1920, 1080], monitor='CBI-prisma-projector', units='pix', fullscr=True,
                screen=1, color=128, colorSpace='rgb255',  **monitor_kwargs):
    """set the various experiment parameters
    """
    stimuli = np.load(stim_path)
    idx = np.load(idx_path)
    stimuli = stimuli[idx]
    expt_params = {}
    # the first dimension of stimuli (retrieved by len) is how many stimuli we have. the next two
    # are the size of the stimuli
    expt_params['stim_size'] = stimuli.shape[1:]
    expt_params['non_blank_stimuli_num'] = stimuli.shape[0]
    # In order to get the right amount of blank time at the end of the run, we insert an
    # appropriate amount of blank stimuli.
    final_nblanks, final_blanks = _create_blanks(final_blank_sec_length, on_msec_length,
                                                 off_msec_length, stimuli.shape[1:], 'final')
    init_nblanks, init_blanks = _create_blanks(init_blank_sec_length, on_msec_length,
                                               off_msec_length, stimuli.shape[1:], 'init')
    stimuli = np.concatenate([init_blanks, stimuli, final_blanks])
    expt_params['final_nblanks'] = final_nblanks
    expt_params['init_nblanks'] = init_nblanks

    # we show a digit with every other stimuli.
    digit_num = stimuli.shape[0] / 2
    probs = np.ones(10)/9
    digits = [int(np.random.uniform(0, 10)), ""]
    for i in range(digit_num-1):
        # because -1 and -3 (and all negative odd indices) will be blanks, we check two and
        # four ago to make sure there are no repeats here
        if np.random.uniform() < fix_button_prob and (len(digits) == 2 or digits[-4] != digits[-2]):
            digits.append(digits[-2])
        else:
            probs_tmp = probs.copy()
            probs_tmp[digits[-2]] = 0
            digits.append(np.random.choice(range(10), p=probs_tmp))
        # after every digit, we show a blank
        digits.append("")
    expt_params['fixation_text'] = iter(digits)
    expt_params['fixation_color'] = iter(['white', 'white', 'black', 'black'] * int(np.ceil(digit_num/2.)))
    # these are all a variety of kwargs used by monitor
    monitor_kwargs.update({'size': size, 'monitor': monitor, 'units': units, 'fullscr': fullscr,
                           'screen': screen, 'color': color, 'colorSpace': colorSpace})
    return stimuli, idx, expt_params, monitor_kwargs


def run(stim_path, idx_path, on_msec_length=300, off_msec_length=200, final_blank_sec_length=16,
        init_blank_sec_length=16, fix_pix_size=10, fix_button_prob=1/6., **monitor_kwargs):
    """run one run of the experiment

    stim_path specifies the path of the unshuffled experiment stimuli, while idx_path specifies the
    path of the shuffled indices to use for this run. This function will load in the stimuli at
    stim_path and rearrange them using the indices found at idx_path, then simply go through those
    stimuli in order, showing each stimuli for `on_msec_length` msecs and then a blank screen for
    `off_msec_length` msecs (or as close as possible, given the monitor's refresh rate).

    For fixation, we show a stream of digits whose colors alternate between black and white, with a
    `fix_button_prob` chance of repeating. Digits are presented with alternating stimuli ON and OFF
    blocks, so that a digit will be shown for on_msec_length+off_msec_length msecs and then there
    will be nothing at fixation for the next on_msec_length+off_msec_length msecs. For now, you
    can't change this.

    All stimuli loaded in from stim_path will be shown.


    Arguments
    ============

    stim_path: string, path to .npy file where stimuli are stored (as 3d array)

    idx_path: string, path to .npy file where shuffled indices are stored (as 1d array)

    on_msec_length: int, length of the ON blocks in milliseconds; that is, the length of time to
    display each stimulus before moving on

    off_msec_length: int, length of the OFF blocks in milliseconds; that is, the length of time to
    between stimuli

    fix_pix_size: int, the size of the fixation digits, in pixels.

    fix_button_prob: float. the probability that the fixation digit will repeat or the fixation dot
    will change color (will never repeat more than once in a row). For fixation digit, this
    probability is relative to each stimulus presentation / ON block starting; for fixation dot,
    it's each stimulus change (stimulus ON or OFF block starting).
    """
    stimuli, idx, expt_params, monitor_kwargs = _set_params(
        stim_path, idx_path, on_msec_length, off_msec_length, fix_button_prob,
        final_blank_sec_length, init_blank_sec_length,
        **monitor_kwargs)

    win = visual.Window(**monitor_kwargs)
    win.mouseVisible = False
    # linear gamma ramp
    win.gammaRamp = np.tile(np.linspace(0, 1, 256), (3, 1))

    fixation = visual.TextStim(win, expt_params['fixation_text'].next(), height=fix_pix_size,
                               color=None)
    fixation.color = expt_params['fixation_color'].next()
    # first one is special: we preload it, but we still want to include it in the iterator so the
    # numbers all match up (we don't draw or wait during the on part of the first iteration)
    img = visual.ImageStim(win, image=imagetools.array2image(stimuli[0]),
                           size=expt_params['stim_size'])

    wait_text = visual.TextStim(win, ("Press 5 to start\nq will quit this run\nescape will quit "
                                      "this session"))
    wait_text.draw()
    win.flip()
    # preload these to save time
    img.draw()
    fixation.draw()

    clock = core.Clock()
    # wait until receive 5, which is the scanner trigger
    all_keys = event.waitKeys(keyList=['5', 'q', 'escape'], timeStamped=clock)
    if 'q' in [k[0] for k in all_keys] or 'escape' in [k[0] for k in all_keys]:
        win.close()
        return all_keys, [], [], expt_params, idx

    keys_pressed = [(key[0], key[1]) for key in all_keys]
    timings = [("start", "off", clock.getTime())]
    fixation_info = []
    for i, stim in enumerate(stimuli):
        if i > 0:
            # we don't wait the first time, and all these have been preloaded while we were waiting
            # for the scan trigger
            if "fixation_text" in expt_params:
                fixation.text = expt_params['fixation_text'].next()
            img.image = imagetools.array2image(stim)
            fixation.color = expt_params['fixation_color'].next()
            img.draw()
            fixation.draw()
            next_stim_time = (i*on_msec_length + i*off_msec_length - 2)/1000.
            core.wait(abs(clock.getTime() - timings[0][2] - next_stim_time))
        win.flip()
        timings.append(("stimulus_%d" % i, "on", clock.getTime()))
        fixation_info.append((fixation.text, clock.getTime()))
        fixation.draw()
        next_stim_time = ((i+1)*on_msec_length + i*off_msec_length - 1)/1000.
        core.wait(abs(clock.getTime() - timings[0][2] - next_stim_time))
        win.flip()
        timings.append(("stimulus_%d" % i, "off", clock.getTime()))
        all_keys = event.getKeys(timeStamped=clock)
        if all_keys:
            keys_pressed.extend([(key[0], key[1]) for key in all_keys])
        if 'q' in [k[0] for k in all_keys] or 'escape' in [k[0] for k in all_keys]:
            break
    win.close()
    return keys_pressed, fixation_info, timings, expt_params, idx


def expt(stimuli_path, number_of_runs, first_run, subj_name, output_dir="../data/raw_behavioral",
         input_dir="../data/stimuli", **kwargs):
    """run a full experiment

    this just loops through the specified stims_path, passing each one to the run function in
    turn. any other kwargs are sent directly to run as well. it then saves the returned
    keys_pressed and frame intervals
    """
    if output_dir[-1] != '/':
        output_dir += '/'
    if input_dir[-1] != '/':
        input_dir += '/'
    file_path = "%s%s_%s_sess{sess}.hdf5" % (output_dir, datetime.datetime.now().strftime("%Y-%b-%d"), subj_name)
    sess_num = 0
    while glob.glob(file_path.format(sess=sess_num)):
        sess_num += 1
    idx_paths = [input_dir + "%s_run%02d_idx.npy" % (subj_name, i) for i in range(first_run, first_run+number_of_runs)]
    for p in idx_paths:
        if not os.path.isfile(p):
            raise IOError("Unable to find array of stimulus indices %s!" % p)
    print("Running %d runs, with the following stimulus:" % number_of_runs)
    print("\t%s" % stimuli_path)
    print("Will use the following indices:")
    print("\t%s" % "\n\t".join(idx_paths))
    for i, path in enumerate(idx_paths):
        keys, fixation, timings, expt_params, idx = run(stimuli_path, path, **kwargs)
        with h5py.File(file_path.format(sess=sess_num), 'a') as f:
            f.create_dataset("run_%02d_button_presses" % i, data=np.array(keys))
            f.create_dataset("run_%02d_fixation_data" % i, data=np.array(fixation).astype(str))
            f.create_dataset("run_%02d_timing_data" % i, data=np.array(timings))
            f.create_dataset("run_%02d_stim_path" % i, data=stimuli_path)
            f.create_dataset("run_%02d_idx_path" % i, data=path)
            f.create_dataset("run_%02d_shuffled_indices" % i, data=idx)
            for k, v in expt_params.iteritems():
                if k in ['fixation_color', 'fixation_text']:
                    continue
                f.create_dataset("run_%02d_%s" % (i, k), data=v)
            # also note differences from default options
            for k, v in kwargs.iteritems():
                if v is None:
                    f.create_dataset("run_%02d_%s" % (i, k), data=str(v))
                else:
                    f.create_dataset("run_%02d_%s" % (i, k), data=v)
        if 'escape' in [k[0] for k in keys]:
            break


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=("Run an experiment! This takes in the path to your unshuffled stimuli, the "
                     "name of your subject, and the number of runs, and passes that to expt. This "
                     "will then assume that your run indices (which shuffle the stimuli) are saved"
                     "in the INPUT_DIR at SUBJ_NAME_runNUM_idx.npy, where NUM runs from FIRST_RUN "
                     "to FIRST_RUN+NUMBER_OF_RUNS-1 (because this is python, 0-based indexing), "
                     "with all single-digit numbers represented as 0#."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("stimuli_path",
                        help="path to your unshuffled stimuli.")
    parser.add_argument("number_of_runs", help="number of runs you want to run", type=int)
    parser.add_argument("subj_name", help="name of the subject")
    parser.add_argument("--input_dir", '-i', help=("path to directory that contains your shuffled"
                                                   " run indices"),
                        default="data/stimuli")
    parser.add_argument("--output_dir", '-o', help="directory to place output in",
                        default="data/raw_behavioral")
    parser.add_argument("--first_run", '-f', type=int, default=0,
                        help=("Which run to run first. Useful if, for instance, you ran the first "
                              "two runs without problem and then had to quit out in the third. You"
                              " should then set this to 2 (because they're 0-indexed)."))
    args = vars(parser.parse_args())
    expt(**args)