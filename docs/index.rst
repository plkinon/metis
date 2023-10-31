.. MetisDocu documentation master file, created by
   sphinx-quickstart on Thu Aug 17 10:27:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to metis's documentation!
=====================================

.. image:: images/logo.png

**Metis** is an object-oriented MATLAB code package that allows to carry out dynamical simulations of mechanical systems. The main feature is the possibility to use and compare different integration schemes. In principle, two different types of analysis can be carried out:

a. Single analysis: allows the integration of equations of motion at a constant time step size. Output are quantities like energy, energy differences, the system's state variables and constraint evaluations (among others) as a function of time.


b. Error analysis: the equations of motion are solved at different time step sizes. The output is a convergence analysis with a specified error quantity as a function of the time step size.

Metis is licensed under the MIT License.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 4

   GettingStarted
   codestructure
   
   
.. _Citation:

Citation
---------

If you used this code (partially) in your project or found it somehow useful, please cite one or both of the following works:

     .. code-block:: console

        @article{kinon_ggl_2023,
        author  = {Kinon, Philipp L. and Betsch, P. and Schneider, S.},
        title   = {The {GGL} variational principle for constrained mechanical systems},
        journal = {Multibody Syst. Dyn.},
        volume  = {57},
        pages   = {211--236},
        year    = {2023},
        doi     = {10.1007/s11044-023-09889-6}
        }

        @article{kinon_structure_2023,
        author  = {Kinon, Philipp L. and Betsch, P. and Schneider, S.},
        title   = {Structure-preserving integrators based on a new variational principle for constrained mechanical systems},
        journal = {Nonlinear Dyn.},
        year    = {2023},
        doi     = {10.1007/s11071-023-08522-7},
        note    = {doi: 10.1007/s11071-023-08522-7}
        }
