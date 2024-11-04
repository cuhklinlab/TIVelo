Welcome to TIVelo's documentation!
==================================

.. figure:: ../fig/workflow_1.pdf
   :align: center

   Workflow of TIVelo

RNA velocity inference is a valuable tool for understanding cell development, differentiation, and disease progression. However, existing RNA velocity inference methods typically rely on explicit assumptions of ordinary differential equations (ODE), which prohibits them to capture complex transcriptome expression patterns. In this study, we introduce TIVelo, a novel RNA velocity estimation approach that first determines the velocity direction at the cell cluster level based on trajectory inference, before estimating velocity for individual cells. TIVelo calculates an orientation score to infer the direction at the cluster level without an explicit ODE assumption, which effectively captures complex transcriptional patterns, avoiding potential inconsistencies in velocity estimation for genes that do not follow the simple ODE assumption. We validated the effectiveness of TIVelo by its application to 16 real datasets and the comparison with five benchmarking methods.



.. note::

   This project is under active development.

User Guide
----------

To run baseline, we need to install::

   pip install --upgrade pip               
   pip install tensorflow==2.10.0                 
   pip install unitvelo==0.2.5.2  multivelo==0.1.3  velovi==0.3.1

.. toctree::
   :maxdepth: 2

   install


tivelo Tutorials
----------------



Commonly-used RNA velocity datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   notebooks/notebooks/Pancreas
   notebooks/notebooks/DentateGyrus
   notebooks/notebooks/DentateGyrus2
   notebooks/notebooks/Erythroid
   notebooks/notebooks/Hindbrain
   notebooks/notebooks/Reprogramming
   notebooks/notebooks/Organoid
   notebooks/notebooks/Retina
   notebooks/notebooks/scNT
   notebooks/notebooks/Hindbrain2

Datasets with scATAC-seq
~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   notebooks/notebooks/HSPCs
   notebooks/notebooks/HumanBrain
   notebooks/notebooks/MouseBrain
   notebooks/notebooks/MouseSkin

FUCCI datasets
~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   notebooks/notebooks/U2OS
   notebooks/notebooks/RPE1
