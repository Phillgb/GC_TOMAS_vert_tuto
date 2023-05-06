File Transfert
==============

When working on Stetson, you will most likely want to transfert data and files from 
your personnal computer to the cluster or the inverse. To do so, there are 2 main 
approachs, using the `Secure File Transfert Protocol`_ (SFTP) or `Globus`_ transfer service.

SFTP
----

The SFTP protocol allow the user to transfer files between the Stetson cluster and 
a local machine while staying in the command line environment. It can be a good option
for transfer if you want to quickly transfert a single file or a limited amount 
of small files to or from the Stetson cluster. With the command line interface, 
you can connect to Stetson using your :code:`username`, with the command that follows:

.. code-block:: bash

    sftp username@stetson.phys.dal.ca

This command will return the :code:`sftp` command prompt where transfers command 
can be issued. You can use the command :code:`help` to get a list of all available 
command in SFTP.

.. note::

    There are also multiples graphical interface options to use SFTP, such as `WinSCP`_
    (Windows), `MobaXterm`_ (Windows), `filezilla`_ (Windows, Mac and Linux) or 
    `cyberduck`_ (Mac and Windows). 

GLOBUS
------

For most use case, and particullarly when transfering large input files for the 
simulations or the output netcdf files resulting from the simulations, we recommend
using the Globus transfer system to move files in and out of the Stetson cluster.

Globus is a file transfer service both fast and reliable, which as an easy to use 
graphical interface. Before using this service, you will need to make sure that you 
have an active account at the `Digital Research Alliance of Canada`_ (formerly known as
Compute Canada). Once you have your account, using globus is as simple as visiting 
the main `Globus portal`_ and select "Compute Canada" when asked which organisation
you are affiliated with. After supplying your username and password, you will be 
redirected to the transfer page. 

Transfers on Globus are made between collections. To access the Stetson cluster 
collection, simply click on the collection box and start typing "stetson". After 
a couple seconds, a list of suggested collection corresponding to Stetson should 
appear, select the one named :code:`Globus Transfert`. You will then see a list 
of the directories present in the :ref:`Globus transfer space` described earlier. 

Navigate to your personnal folder, where the files to transfer are located and select 
all files to be tranfered as you would in the explorer of your local machine. 
You can then click on Transfer or sync in right hand side option bar of the interface 
to make your destination collection pane appear. In the search bar of this second 
pane, you can search for the collection name of the destination you want to put the 
selected files in. If you want to transfer files between the Stetson cluster and 
your personnal machine make sure that the Globus Connect personnal software is 
running or that you follow the steps in :ref:`Personal Globus Endpoint` if it is 
your first time transfering data to a personnal endpoint.

Once you have navigated to the folder where you want to store the file to transfer 
in the receiving collection pane, you can click on the *Start* blue button at the 
top of the origin collection pane to start the transfer of the files. Once the transfer 
as started, you can monitor its progress in the *Activity button* of the left hand side
of the screen. 

For more information, on the usage of Globus, you can visit the Digital Research
Alliance of Canada `Globus wiki`_ page.

.. _Personal Globus Endpoint:

Personal Globus Endpoint:
--------------------------

To tranfer file using Globus to and from your personnal machine, you will have first
to create a personnal collection refering to your personnal machine. This process 
need only to be done the first time you use Globus with your personnal machine.
Every subsequent time, you only need to make sure that the Globus Connect personnal
software is running on your machine before attempting a transfer to or from your 
personnal machine.

To install **Globus Connect personnal** follow these steps:
    
#. From the `Globus Portal`_, go to the *Collection* page by clicking the icon on the left side of the screen. 
#. Click on the "+ Get Globus Connect Personnal" button in the top of the screen
#. Click on the download link corresponding to your operating system
#. Follow the instruction on screen to install Globus Connect personnal
#. You will now be able to use globus to and from your personnal computer.

The Globus Connect personnal software need to be running and not paused every 
time you will want to access the collection point on your personnal machine. You 
can make sure that the software is not paused by "right clicking" on the Globus
icon in the notification tray of your system.

.. _Secure File Transfert Protocol: https://en.wikipedia.org/wiki/SSH_File_Transfer_Protocol
.. _additionnal SCP examples: https://www.hypexr.org/linux_scp_help.php
.. _Globus: https://www.globus.org/
.. _WinSCP: https://winscp.net/eng/index.php
.. _MobaXterm: https://mobaxterm.mobatek.net/
.. _filezilla: https://filezilla-project.org/
.. _cyberduck: https://cyberduck.io/?l=en
.. _Digital Research Alliance of Canada: https://alliancecan.ca/en
.. _Globus Portal: https://app.globus.org/file-manager
.. _Globus wiki: https://docs.alliancecan.ca/wiki/Globus
