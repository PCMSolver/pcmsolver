Timer class
-----------

The ``Timer`` class enables timing of execution throughout the module.
The class uses the `Boost.Timer <http://www.boost.org/doc/libs/1_59_0/libs/timer/doc/index.html>`_
and `Boost.Chrono <http://www.boost.org/doc/libs/1_59_0/doc/html/chrono.html>`_ libraries.
When compiled with timer support, i.e. when ``-DENABLE_TIMER=ON`` is passed to the ``setup.py``
script, we have to link against those Boost libraries.
Timing macros are available by inclusion of the ``Config.hpp`` header file.

The class is basically a wrapper around an ordered map of strings and cpu timers.
To time a code snippet:

.. code-block:: cpp

   TIMER_ON("code-snippet");
   // code-snippet
   TIMER_OFF("code-snippet");

The timings are printed out to the ``pcmsolver.timer.dat`` by a call
to the ``TIMER_DONE`` macro. This should obviously happen at the very end
of the execution!

.. doxygenfile:: TimerInterface.hpp
   :project: PCMSolver
