.. _user_guide:

User Guide
==========

This document provides a brief description of the compilation and
operation of Enzo, a structured `Adaptive Mesh Refinement
<http://en.wikipedia.org/wiki/Adaptive_mesh_refinement>`_ (SAMR, or
more loosely AMR) code which is primarily intended for use in
astrophysics and cosmology. The User's Guide is intended to explain
how to compile and run Enzo, the initial conditions generation code
and the various analysis tools bundled with Enzo. The instructions on
actually running the code are not comprehensive in that they are not
machine or platform-specific.  Arguably the most useful and important
piece of this guide is :doc:`../EnzoParameters`, which contains
descriptions of all of the roughly 300 possible input parameters (as
of September 2008). For more detailed information on the Enzo
algorithms and on running Enzo on different platforms, you should
refer to the :doc:`../tutorials/index`. Detailed information on the
algorithms used in Enzo will be available in the method paper
(unreleased as of September 2008). In the meantime, see the
:ref:`EnzoPrimaryReferences` for more concrete Enzo information.

This guide (and Enzo itself) was originally written by Greg
Bryan. Since the original writing of both the simulation code and the
User's Guide, the maintenance of Enzo and its associated tools and
documentation was for some time largely driven by the `Laboratory for
Computational Astrophysics <http://lca.ucsd.edu>`_ at `The University
of California, San Diego <http://www.ucsd.edu>`_, but it is now a
fully open source community with developers from Stanford, Columbia,
Princeton, UCSD, University of Colorado, Michigan State, UC Berkeley,
and many other universities.  Your input in improving both the code
and the User's Guide is appreciated -- developement of the code is
driven by working researchers, and we encourage everyone who has made
useful changes to contribute those changes back to the community and
to participate in the collaborative development of the code.  Email
inquiries and comments should be directed to the `Enzo Users' List
<http://mailman.ucsd.edu/mailman/listinfo/enzo-users-l>`_. Thank you!

.. toctree::
   :maxdepth: 2

   ExecutablesArgumentsOutputs.rst
   RunningEnzo.rst
   EnzoTestSuite.rst
   CosmologicalInitialConditions.rst
   EnzoOutputFormat.rst
   AnalyzingWithYT.rst
   SimulationNamesAndIdentifiers.rst
   EmbeddedPython.rst
   HierarchyFile.rst
   FlowChart.rst
