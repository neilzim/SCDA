================================================================================
How to prepare, launch, and harvest a coronagraph design survey on NCCS Discover
================================================================================

1. Define the survey parameters and create the optimization files
-----------------------------------------------------------------

Start from the example notebook in the github repo, https://github.com/neilzim/SCDA/blob/master/examples/small_hexap_survey.ipynb. Copy this to a new file outside of your local copy of the SCDA repo. Execute the cells at least up to the step "Store the design survey as a serialized python object". Be sure that the input mask files exist.  

2. Copy the survey files to NCCS Discover
-----------------------------------------

``$ rsync -turp /path/to/new/surveydir discover.nccs.nasa.gov:/discover/nobackup/youruserid/``  

3. Install the latest version of the SCDA python package in your home directory on NCCS Discover
------------------------------------------------------------------------------------------------

Network restrictions disallow syncing the SCDA code directly to your home directory via ``git pull origin master``, etc. Instead, from your local workstation open the github page https://github.com/neilzim/SCDA in a web browser and click the ``Download ZIP`` button to download the current master branch.

Copy this to your Disocver home directory:

``$ scp -p ~/Downloads/SCDA-master.zip discover.nccs.nasa.gov:~``

Log in to Discover, and unzip it directly in the home directory:

``$ unzip SCDA-master.zip``

and overwrite any existing files if prompted. Finally, make sure your bash profile (``~/.bashrc``) sets the PYTHONPATH to include ~/SCDA-master, so that scripts can import the scda.py module:

.. code-block:: bash

  export PYTHONPATH=${PYTHONPATH}:$HOME/SCDA-master

And while we're looking at the bash profile, also be sure it sets the shell variables necessary to run AMPL:

.. code-block:: bash

  export AMPLBIN=$HOME/AMPL
  export AMPLLIB=$HOME/AMPL
  export PATH=${PATH}:${AMPLBIN}
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${AMPLLIB}



