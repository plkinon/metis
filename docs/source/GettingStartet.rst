How to get startet
=====

.. _First steps:

First steps
------------

1. Clone this directory or download the .zip folder
2. Get matlab2tikz (optional)
3. Open the MATLAB editor or run it with the shell script metis.sh
4. Open start_metis_single_analysis.m
5. Adjust <input_file_name> corresponding to a file from /input, for more info look at README_inputz
  .. code-block:: console
6. Adjust the path to the matlab2tikz directory in the chosen input-file
7. Execute start_metis_single_analysis.m for a first simulation
8. Edit or change input file or create a new one in /input
9. For error analyses run start_metis_error_analysis.m with corresponding input-file
10. Have fun!

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']
