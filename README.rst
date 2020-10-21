pyrad
=====

.. image:: https://travis-ci.com/menzel-gfdl/pylbl.svg?branch=reorganize-take-2
   :target: https://travis-ci.com/menzel-gfdl/pylbl
   :alt: Build Status

.. image:: https://readthedocs.org/projects/pylbl/badge/?version=latest
   :target: https://pylbl.readthedocs.io/en/latest/
   :alt: Documentation Status

Pyrad is a simple (one-dimensional), pure python, all-sky atmospheric radiation package.


Gas absorption coefficients can be calculated by:

.. literalinclude:: ../tests/example-gas-optics.py
   :language: python

.. image:: gas-optics.png

By default, the above code will download the necessary molecular line and total partition
function data from the web.  This can take a significant amount of time, especially if the
Gas objects are created often.  To retify this, I recommend creating local SQLite databases,
then re-using when creating Gas objects:

.. literalinclude:: ../tests/example-gas-optics-database.py
   :language: python

Clouds are generated in a stochastic fashion (typically found in GCMs):

.. literalinclude:: ../tests/example-stochastic-clouds.py
   :language: python

.. image:: liquid-stochastic-clouds.png

.. image:: ice-stochastic-clouds.png

and their optics are calculated using standard look-up table parameterizations:

.. literalinclude:: ../tests/example-cloud-optics.py
   :language: python

Aerosol optics are also calculated using a GCM parameterization:

.. literalinclude:: ../tests/example-aerosol-optics.py
   :language: python

