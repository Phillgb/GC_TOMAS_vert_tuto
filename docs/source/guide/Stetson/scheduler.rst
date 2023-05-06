.. _How to use qsub:

Scheduler use 
=============

To execute tasks through a scheduler, the user need to provide that scheduler with
some informations on the ressources needed to execute those tasks and a list of 
actions that need to be followed to lauch the desired programs. These informations
and steps are passed to the scheduler in the form of text in what we call a *job script*.
The scheduler interpret this script and will execute the list of command given when 
the ressources required by the tasks have been freed. 

Creating a job script
---------------------

Here is an example of a typical job script used to execute a simulation of GEOS-Chem
TOMAS on Stetson.

.. code-block::bash
    #!/bin/bash
    #$ -S /bin/bash
    ./etc/profile
    #$ -o job_output
    #$ -M phillipe.gauvin-bourdon@dal.ca
    #$ -m be
    #$ -pe threads 36

    #######################################################################################################################
    ### GEOS-Chem run script for Stetson
    #######################################################################################################################
    # Set proper # of threads for OpenMP
    export OMP_NUM_THREADS=36

    # Activate proper packages
    source /software/share/set-path-gchp-intel-conda-2022.0.2-intelmpi

    # Navigate to run directory
    cd /data12/pgbourdon/GCC_v14.1.1/rundirs/gc-4x5_47L_merra2_fullchem_TOMAS15

    # Run GEOS-Chem
    ./gcclassic > runlog.txt

This script have 2 main sections, the first part being the directives passed to 
the scheduler regarding the ressources that will be needed to execute the tasks 
we want to run (called shebang) and the second part being a list of actions to be 
executed to launch the simulation.

The most commonly used shebang options are the following:

.. option:: #$ -S /bin/bash 
    
    Specify which interpreter to use to read and process the following command.
    (Using the bash shell interpret in the current case.)

.. option:: #$ -o job_output
    
    Specify that where the path where the standard output of the job will be written.

.. option:: #$ -M phillipe.gauvin-bourdon@dal.ca
    
    Specify at which email to send notifications corresponding with the current job.

.. option:: #$ -m be
    
    Specify under which circumstances notifications need to be sent to the email 
    address provided. 
    
    * b: Begginning of the job
    * e: End  of the job
    * a: Job is aborted or rescheduled
    * s: Job is suspended
    * n: No mail sent

.. option:: #$ -pe threads 36
    
    Specify the parallel environments options. Here qsub is asked to use 36 threads when 
    possible.

Following the *shebang*, the rest of the script simply list actions to be executed 
in order and that will lead to the execution of our GEOS-Chem TOMAS simulation.
Those lines are corresponding to the actions that a user would enter in the terminal 
command lines if they were to execute a simulation manually.

The main commands used when launching GEOS-Chem simulation are the following:

.. option:: export OMP_NUM_THREADS=36

    Specify the number of computationnal cores that GEOS-Chem should use. Normally, 
    needs to match the number specified with the *-pe* option above.

.. option::  cd /data12/pgbourdon/GCC_v14.1.1/rundirs/gc-4x5_47L_merra2_fullchem_TOMAS15

    Navigate to the run directory of the simulation to execute.

.. option:: ./gcclassic > runlog.txt

    Launch the GEOS-Chem classic simulation and pipes the command line output to 
    a file named :file:`runlog.txt`.

Qsub opperation
---------------

After having created your *job script*, the opperations containned within the script 
are passed to the scheduler **qsub** via the command line command :code:`qsub`. 
It does not matter what your current directory is when you are passing the command 
to read the *job script* to **qsub** as long as you built your *job script* as if 
the user was starting back from the home directory of Stetson. 

In its most simple form, the :code:`qsub` command only demand that a *job script*
to be passed as parameter, like the following :code:`qsub <job_script.sh>`. Since
when using Stetson, we would like for the scheduler to use the newnodes of the server
to execute the simulations of GEOS-Chem we need to pass an extra option with the 
:code:`qsub` command. 

.. option:: qsub -q 2017-batch geoschem_run.sh

    Typical **qsub** command line call to execute the *job script* :file:`geoschem_run.sh`
    on the new nodes of the Stetson cluster. 

Once your job has been submitted the following commands can be usefull to monitor 
the job status and manage it. 

.. option:: qstat -u <username> 

    Command returning all pending or running job assosciated with the specified 
    username. :code:`qstat` can also be used without any option to return a list 
    of all the job that **qsub** is treating at the moment.

.. option:: qdel <job_id>

    Command terminating the job with the specified ID.

    