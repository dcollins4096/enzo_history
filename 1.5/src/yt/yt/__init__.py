"""
YT is a package written primarily in Python designed to make the task of
running Enzo easier.  It contains facilities for creating Enzo data (currently
in prototype form) as well as runnning Enzo simulations, simulating the actions
of Enzo on various existing data, and analyzing output from Enzo in a
wide-variety of methods.

An ever-growing collection of documentation is also available at
http://yt.enzotools.org/doc/ . Additionally, there is a
project site at http://yt.enzotools.org/ with recipes, a wiki, a subversion
changelog and a bug-reporting system.

YT is divided into several packages, all named after characters from Snow
Crash.

Lagos
=====

Lagos deals with data structures. It defines things like EnzoGrid, EnzoData,
Enzo2DData, EnzoSphere, etc. If you want to handle actual data, use Lagos.

Raven
=====

Raven is the plotting interface.  All data plotting goes through
Raven.

Enki
====

Enki is the package used to create data, and instantiate runs. It supports
creating Enzo Problems, and then using SWIG-interfaced Enzo calls to
actually create the data for those problems. Additionally, facilities are
being developed to use Enki to directly execute runs.

Right now, Enki is still largely experimental.  It provides some primitive
methods for interacting with Enzo, but more work needs to be done before
it reaches its vision.

Fido
====

Fido is the messenger/protector of data.  It takes data outputs, puts them
wherever you want, and then calls a function handler to deal with that data.
Ultimately Fido will deal with all file-handling; submission of runs to a
central (user-specific) database is in the works, and Fido will be the
entity that moves the files in and out of storage.

Deliverator
===========

The Deliverator is a TurboGears-based system for querying and displaying
images. Images are dumped from Raven into local, web-accessible storage
space, and then metadata about those images is submitted to The Deliverator.
The user (you) then goes to the Deliverator website and views those plots.

The base package YT provides facilities for configuration files and logging (via
the Python logger.)

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

#from yt.logger import *
#from yt.config import *
#from yt.funcs import *
