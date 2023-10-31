Getting started
=====

.. _Requirements:

Requirements
------------

1. `MATLAB <https://www.mathworks.com/products/matlab.html>`_
2. `matlab2tikz <https://github.com/matlab2tikz/matlab2tikz>`_ (optional)


.. _First steps:

First steps
------------

1. Clone this directory or download the .zip folder
2. Get `matlab2tikz <https://github.com/matlab2tikz/matlab2tikz>`_ (optional)
3. Open the `MATLAB <https://www.mathworks.com/products/matlab.html>`_ editor or run it with the shell script metis.sh
4. Open start_metis_single_analysis.m
5. Adjust <input_file_name> corresponding to a file from /input, for more info look at README_input.md

    .. code-block:: console

           [simulation, system, integrator, solver] = Metis('input/<input_file_name>',1,1);
6. Adjust the path to the matlab2tikz directory in the chosen input-file
7. Execute start_metis_single_analysis.m for a first simulation
8. Edit or change input file or create a new one in /input
9. For error analyses run start_metis_error_analysis.m with corresponding input-file
10. Have fun!

