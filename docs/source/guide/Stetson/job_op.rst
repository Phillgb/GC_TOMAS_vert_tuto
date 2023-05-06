Job submission
==============

This section will describe the processes of lauching programs on the Stetson 
server and using the scheduler to launch lenghty tasks.

Logging in into computing nodes
-------------------------------

The execution of any programs on stetson are done by logging into compute nodes 
directly via :code:`ssh`. The compute nodes on Stetson are numbered from 1 to 32, 
with the most recent and powerfull nodes being the 9 last ones (23 to 32). To log 
into the computing node 24 for example, the user would have to input the following 
in the command line :code:`ssh newnode24`. This will redirect you into the home 
directory on node 24. To execute a programs, you would then need to navigate to the 
appropriate working directory and launch the program as you would on a local machine.

Before logging into a new node, it is a good practice to verify the current load 
on each node. The usage of each node can be monitored via the `Cluster status 
website <http://stetson.phys.dal.ca>`_.

Loading packages
----------------

Some programs might require you to have some packages installed on the system, 
prior to their execution. For example, the compilation of the GEOS-Chem model 
require the user to have specific modules like netcdf and fortran compiler to work. 

On stetson, most packages needed to opperate GEOS-Chem with or without TOMAS are 
already available as packaged modules in shared :ref:`packages libraries`.

To load the desired packages in a session, the user need to execute the desired
libraries as a source code. For most use cases, the necessary packages to execute 
GEOS-Chem are present in the :file:`set-path-intel-gchp-conda-2022.0.2-intelmpi`
and can be loaded into a session via the following command 
:code:`source /software/share/set-path-intel-gchp-conda-2022.0.2-intelmpi`.

Scheduler
---------

Additionnally, to manually launching programs on Stetson, it is also possible to 
execute series of command by using a scheduler. Stetson uses the grid engine scheduler
**qsub**.

For lenghty simulations to be run on Stetson, it is recommended to use the scheduler 
in order to execute those simulations. The use of the scheduler allow for an even
repartition of the server ressources when multiples tasks are launched at the same
time and is necessary to execute the model in parallel mode. An other advantage 
of the scheduler, is that it allows for a better management of time. For example, 
if all nodes are occupied when a user manually access the server, launching a job 
via the scheduler will allow it to start as soon as the necessary ressources get 
freed, instead of rellying on the user manually accessing the server and guessing 
when the nodes will get freed up.

Information on how to work with qsub are available in the next section :ref:`How to use qsub` 
and in the original `qsub documentation`_.


.. _qsub documentation: https://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

