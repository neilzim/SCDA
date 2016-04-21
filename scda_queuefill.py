#!/usr/bin/env python

'''
Queue manager script for SCDA mask optimization programs
4-20-2016

AUTHOR

Neil Zimmerman
Space Telescope Science Institute

USAGE

Execute directly from the command line, specifying the name of the design
survey file (extension .pkl) as an argument. For example:

$ ./scda_queuefill.py april_survey01_15bw_ntz_2016-04-20.pkl

Informative output will print to stdout. To instead append the output to a text file, run

$ ./scda_queuefill.py april_survey01_15bw_ntz_2016-04-20.pkl >> /discover/nobackup/nzimmerm/april_survey01_15bw/queuefill.log

OUTPUTS

* Indicates how many programs are currently in the queue, what programs were
newly submitted, and how many programs in the survey have been submitted so
far.

* Creates a text file containing a crontab command, if it doesn't already exist.
This crontab file will have the same name as the survey, minus the .pkl
extension, and prepended with "crontab_". For the above example, the cron file
is named crontab_april_survey01_15bw_ntz_2016-04-20, and stored in the
directory /discover/nobackup/nzimmerm/april_survey01_15bw/

CRONTAB CONFIGURATION

The crontab file is configured to execute this script once per hour, begninning
5 minutes after the time of creation. The first number in the crontab file
specify the "minutes hand" of the clock when the crontab is triggered.

In order to configure the crontab on NCCS Discover, you need to run the crontab
command from the discover-cron node. This can be accomplished in a single ssh
command, following the NCCS password prompt:

$ ssh discover-cron crontab /discover/nobackup/nzimmerm/april_survey01_15bw/crontab_april_survey01_15bw_ntz_2016-04-20

Note that it is necessary to specify the full path to the crontab file.

'''

import sys
import os
import subprocess 
import getpass
import datetime
SCDA_location = os.environ["SCDA"]
sys.path.append(os.path.expanduser(SCDA_location))
import scda

QUEUE_MAX = 25

assert len(sys.argv) >= 2, "Missing design survey file (.pkl) argument"
try:
    survey_fname = os.path.abspath(sys.argv[1])
    survey = scda.load_design_param_survey(survey_fname)
except:
    print("Could not load design survey file: {0}".format(survey_fname))
    sys.exit(1)

cwd = os.getcwd()
survey_dir = os.path.dirname(os.path.abspath(survey_fname))
os.chdir(survey_dir)

ocount_str = subprocess.check_output("squeue -u {0:s} | wc -l".format(getpass.getuser()), shell=True)
qcount = int(ocount_str.split('\n')[0]) - 1

print("{0:s}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
print("{0:d} jobs in {1:s}'s queue".format(qcount, getpass.getuser()))
sys.stdout.flush()

max_submission_count = QUEUE_MAX - qcount

# Fill the queue
new_submission_count = 0
if new_submission_count < max_submission_count:
    for coron in survey.coron_list:
        if coron.ampl_submission_status is not True:
            try:
                subprocess.check_call("sbatch {0:s}".format(coron.fileorg['slurm fname']), shell=True)
                coron.ampl_submission_status = True
                print("          {0:s}".format(coron.fileorg['slurm fname']))
                sys.stdout.flush()
                new_submission_count += 1
            except subprocess.CalledProcessError:
                coron.ampl_submission_status = False
                break
        if new_submission_count == max_submission_count:
            break

# Tally all submissions to date, check for the existence of solution files, and update the corresponding statuses
overall_submission_count = 0
solution_count = 0
for coron in survey.coron_list:
    if coron.ampl_submission_status is True:
        overall_submission_count += 1
        if os.path.exists(coron.fileorg['sol fname']):
            coron.solution_status = True
            solution_count += 1

print("{0:d} out of {1:d} optimization jobs in the survey have been submitted, and {2:d} have solutions.".format(overall_submission_count, survey.N_combos, solution_count))

if solution_count == survey.N_combos:
    print ("Done!")

# Store updated survey object
survey.write(survey_fname)

# Write crontab file if it does not exist
crontab_fname = "crontab_{0:s}".format(os.path.basename(survey_fname)[:-4])
log_fname = os.path.join(survey_dir, "queuefill.log")
if not os.path.exists(crontab_fname):
    cron_fobj = open(crontab_fname, "w")
    dt = datetime.datetime.now()
    run_minute = (dt.minute + 5) % 60
    cron_fobj.write("{0:d} * * * *  (. /etc/profile ; . $HOME/setcronenv.sh ; /usr/local/other/SLES11/SIVO-PyD/1.1.2/bin/python $SCDA/{1:s} {2:s} 1>> {3:s} 2>&1)".format(run_minute, os.path.basename(__file__), survey_fname, log_fname))
    cron_fobj.write("\n")
    cron_fobj.close()
    print("Wrote crontab file to {0:s}".format(crontab_fname))

print("")
