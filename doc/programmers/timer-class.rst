Timer class
-----------

The ``Timer`` class enables timing of execution throughout the module.
Timer support is enabled by passing ``-DENABLE_TIMER=ON`` to the ``setup.py``
script.
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
