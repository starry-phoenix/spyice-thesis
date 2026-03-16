Running Tests
=============

Run the test suite:

.. code-block:: console

   hatch test tests\test_userinput.py

Check code coverage using Coverage.py:

.. code-block:: console

   hatch test --cover


Static Analysis
===============

Hatch uses Ruff for linting:

.. code-block:: console

   hatch fmt


Type Hinting
============

Hatch uses mypy for type hinting:

.. code-block:: console

   hatch run types:check