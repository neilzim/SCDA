================================================================================
How to prepare, launch, and harvest a coronagraph design survey on NCCS Discover
================================================================================

0. Install the SCDA source tree in your NCCS Discover home directory
---------------------------------------------------------------------

Network restrictions disallow syncing the SCDA source code directly to NCCS from github. Instead, from your local workstation create a pristine copy of the master branch:

``$ git clone https://github.com/neilzim/SCDA.git SCDA-master``

or if such a clone already exists, update it via

.. code-block::

  $ cd SCDA-master
  $ git pull origin master

Sync the SCDA source tree to your Discover home directory:

``$ rsync -turp --progress SCDA-master discover.nccs.nasa.gov:~``

In your home directory on Discover, edit the **.bashrc** file to include the following lines:

.. code-block:: bash

  export AMPLBIN=$HOME/AMPL
  export AMPLLIB=$HOME/AMPL
  export SCDA=$HOME/SCDA-master
  export PATH=${PATH}:${AMPLBIN}:$SCDA
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${AMPLLIB}
  export PYTHONPATH=${PYTHONPATH}:$SCDA
  
For these shell variables to take effect, you need to open a new login session or run ``source ~/.bashrc`` in every existing console.

This bash profile will need to be modified if the locations of either the SCDA source tree or the AMPL/Gurobi installation change.

1. Define the survey parameters and create the optimization files
-----------------------------------------------------------------

On your local worksation, the easiest way to define a survey is to start from the example notebook in the SCDA source tree, **./examples/small_hexap_survey.ipynb** (viewable in a web browser at https://github.com/neilzim/SCDA/blob/master/examples/small_hexap_survey.ipynb).

Copy the notebook to a location outside of your SCDA source tree. Open the notebook in jupyter via

``$ jupyter notebook small_hexap_survey.ipynb``

From the jupyter interface, rename the notebook in a way descriptive to your survey. Execute the cells at least up to the step "Store the design survey as a serialized python object". Make sure the cell "Check the status of input files needed to run the AMPL program" evaluates as ``True`` before moving on.

2. Sync the survey files and input masks to NCCS Discover
---------------------------------------------------------

Sync the input masks:

``$ rsync -turp --progress /astro/opticslab1/SCDA/Apertures/InputMasks /discover/nobackup/youruserid/``

Sync the new survey:

``$ rsync -turp --progress /path/to/new/surveydir discover.nccs.nasa.gov:/discover/nobackup/youruserid/``

3. Launch the survey
--------------------

From a Discover console, use ssh to switch to the special **discover-cron** node:

``$ ssh discover-cron``

Change to the survey directory via

``$ cd /discover/nobackup/youruserid/surveydir``

Launch the queue manager script, specifying the existing survey archive file (extension .pkl) as the first argument, and directing the output to **queuefill.log**.

``$ scda_queuefill.py awesomesurvey.pkl >> queuefill.log``

If you print the queue log, it should indicate that it wrote a new crontab file, in this case named ``crontab_awesomesurvey``. Set the crontab so that scda_queufill is run automatically at hourly intervals:

``$ crontab crontab_awesomesurvey``

You can edit the crontab directly, for example to modify the minute of execution, via

``$ crontab -e``
