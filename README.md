# Tutorial analysis

This repository contains the code necessary to analyze the tutorial
data set (found on the Winawer
Lab's [OSF page](https://osf.io/avp83/)). It also contains some of the
setup materials (the python environment, code used to run the
experiment, etc.)

# Setup

## Python

To setup your python system,
install [miniconda](https://conda.io/miniconda.html) (this just
contains python and `conda`, a very nifty package manager; choose
python 3.6), and then run `conda env create -f environment.yml`. This
will create a virtual environment
(see
[here](https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c) for
some more details on what that means) called `winawerlab` and installs
a variety of python packages that will be useful. Note that this
initial install will take a while, as there's a lot to grab.

In order to activate this virtual environment, simply type 

```
conda activate winawerlab
```

Now, when you type `python`, you'll use the python from within your
virtual environment and have access to all the packages. For most
day-to-day work, testing and the like, you probably want to use
Jupyter Lab or a Jupyter Notebook.

If you're completely new to python, there are a large amount of
resources on the web available for free. I highly recommend the
tutorials written by Jake
VanderPlas:
[A Whirlwind Tour of Python](https://github.com/jakevdp/WhirlwindTourOfPython) is
a brief introduction and
the
[Python Data Science Handbook](https://github.com/jakevdp/PythonDataScienceHandbook) goes
into more details.

## Git

In order to interact with the Winawer lab Github, you'll need
git. This may already be installed on your system (type `which git`
and see if it gives you a path), but if not, you can install
it
[here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

You can interact with git from the command line or via GUIs. Here are
two good GUIs: [Github's official one](https://desktop.github.com/)
(I'm not sure if this will interact with other git hosting services,
so it may be slightly limited)
and
[GitKraken](https://www.gitkraken.com/). Many
[text editors](#text-editors) also include a handy interface with git,
so you should check that!

git is a little confusing,
so [here's](http://neuroplausible.com/github) a brief cheat sheet that
covers some of the basics (it gives the commands to use on the command
line, but you'll be interacting with the same concepts if you interact
with git via a GUI).

It's also very helpful and good practice to make good use of
your [.gitignore](https://git-scm.com/docs/gitignore) file so you
don't add a bunch of unnecessary files. The one included in this
repository includes stuff you want to ignore from python, matlab, and
that annoying OSX `.DS_Store` file.

## Text editors

You'll need a text editor to edit text files, such as python scripts
and the markdown files that make up the content on
the
[Winawer Lab wiki](https://github.com/WinawerLab/winawerlab.github.io)
(you can also use it for writing MATLAB scripts, but the built-in
editor is pretty good for this purpose). You can use one of the
default editors on your system (`nano`, Notepad, etc.), but you'll
probably want one that has a little more features. Much digital ink
has been spilled on what the best text editor is (largely debating the
merits of emacs vs. vi/vim), and there's plenty on the internet to
read on the topic. Noah and Billy both
use [emacs](https://www.gnu.org/software/emacs/), but there's a steep
learning curve (Billy is happy to evangelize as to why you should use
emacs if you're interested!); Serra
uses [Sublime](https://www.sublimetext.com/). If you have no specific
preferences and just want something straightforward, I've heard great
things baout [Atom](https://atom.io/). Since it's created by Github,
it automatically works well with git (and so you won't need a separate
GUI to interact with git). It also has
a [markdown preview mode](https://github.com/atom/markdown-preview),
which is very useful for editing the files on the lab wiki. Note that
emacs and Sublime also have these features, so it will ultimately come
down to personal preference.

# Recommended readings
