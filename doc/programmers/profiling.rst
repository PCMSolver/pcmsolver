Profiling
---------

You should obtain profiling information before attempting any optimization of
the code. There are many ways of obtaining this information, but we have only
experimented with the following:

#. Using Linux ``perf`` and related `tools <http://www.brendangregg.com/perf.html>`_.
#. Using ``gperftools``.
#. Using Intel VTune.

Profiling should be done using the standalone executable ``run_pcm`` and any of
the input files gathered under the ``tests/benchmark`` directory. These files
are copied to the build directory. If you are lazy, you can run the profiling
from the build directory:

.. code-block:: bash

   >>> cd tests/benchmark

   >>> env PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH
          python <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

Using ``perf``
==============

``perf`` is a tool available on Linux. Though part of the kernel tools, it is
not usually preinstalled on most Linux distributions. For visualization
purposes we also need `additional tools <https://github.com/brendangregg/perf-tools>`_,
in particular the `flame graph generation scripts <https://github.com/brendangregg/FlameGraph>`_
Probably your distribution has them prepackaged already.
``perf`` will trace all CPU events on your system, hence you might need to
fiddle with some kernel set up files to get permissions to trace events.

.. note::
   ``perf`` **is NOT** available on ``stallo``. Even if it were, you would
   probably not have permissions to record kernel traces.

These are the instructions I used:

1. Trace execution. This will save CPU stack traces to a ``perf.data`` file.
   Successive runs do not overwrite this file.

   .. code-block:: bash

      >>> cd tests/benchmark

      >>> perf record -F 99 -g -- env PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH python
                    <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

2. Get reports. There are different ways of getting a report from the
   ``perf.data`` file. The following will generate a call tree.

   .. code-block:: bash

      >>> perf report --stdio

3. Generate an interactive flame graph.

   .. code-block:: bash

      >>> perf script | stackcollapse-perf.pl > out.perf-folded

      >>> cat out.perf-folded | flamegraph.pl > perf-run_pcm.svg

Using ``gperftools``
====================

This set of tools was previously known as Google Performance Tools. The
executable needs to be linked against the ``profiler``, ``tcmalloc``
and ``unwind`` libraries.
CMake will attempt to find them. If this fails, you will have to install them,
you should either check if they are available for your distribution or compile
from source.
In principle, one could use the ``LD_PRELOAD`` mechanism to skip the *ad hoc*
compilation of the executable.

.. note::
   ``gperftools`` **is** available on ``stallo``, but it's an ancient version.

1. Configure the code with the ``--gperf`` option enabled. CPU and heap
   profiling, together with heap-checking will be available.

2. CPU profiling can be done with the following command:

   .. code-block:: bash

      >>> env CPUPROFILE=run_pcm.cpu.prof PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH
              python <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

  This will save the data to the ``run_pcm.cpu.prof`` file. To analyze the gathered
  data we can use the ``pprof`` script:

  .. code-block:: bash

     >>> pprof --text <build_dir>/bin/run_pcm run_pcm.cpu.prof

  This will print a table. Any row will look like the following:

  .. code-block:: bash

     2228   7.2%  24.8%    28872  93.4% pcm::utils::splineInterpolation

  where the columns respectively report:

  #. Number of profiling samples in this function.
  #. Percentage of profiling samples in this function.
  #. Percentage of profiling samples in the functions printed so far.
  #. Number of profiling samples in this function and its callees.
  #. Percentage of profiling samples in this function and its callees.
  #. Function name.

  For more details look `here <https://gperftools.github.io/gperftools/cpuprofile.html>`_

3. Heap profiling can be done with the following command:

   .. code-block:: bash

      >>> env HEAPPROFILE=run_pcm.hprof PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH
              python <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

  This will output a series of datafiles ``run_pcm.hprof.0000.heap``,
  ``run_pcm.hprof.0001.heap`` and so forth. You will have to kill execution
  when enough samples have been collected.
  Analysis of the heap profiling data can be done using ``pprof``. `Read more
  here <https://gperftools.github.io/gperftools/heapprofile.html>`_


Using Intel VTune
=================

This is probably the easiest way to profile the code.
`VTune <https://software.intel.com/en-us/intel-vtune-amplifier-xe>`_ is Intel software, it might be possible to get a personal, free license.
The instructions will hold on any machine where VTune is installed and you can
look for more details on the `online documentation <https://software.intel.com/en-us/vtune-amplifier-help>`_
You can, in principle, use the GUI. I haven't managed to do that though.

On ``stallo``, start an interactive job and load the following modules:

.. code-block:: bash

   >>> module load intel/2018a

   >>> module load CMake

   >>> module load VTune

   >>> export BOOST_INCLUDEDIR=/home/roberto/Software/boost/include

   >>> export BOOST_LIBRARYDIR=/home/roberto/Software/boost/lib

You will need to compile with optimizations activated, *i.e.* release mode.
It is better to first parse the input file and then call ``run_pcm``:

.. code-block:: bash

   >>> cd <build_dir>/tests/benchmark

   >>> env PYTHONPATH=../../lib64/python:$PYTHONPATH
       python ../../bin/go_pcm.py --inp=standalone_bubble.pcm

To start collecting hotspots:

.. code-block:: bash

   >>> amplxe-cl -collect hotspots ../../bin/run_pcm @standalone_bubble.pcm

VTune will generate a folder ``r000hs`` with the collected results. A report
for the hotspots can be generated with:

.. code-block:: bash

   >>> amplxe-cl -report hotspots -r r000hs > report
