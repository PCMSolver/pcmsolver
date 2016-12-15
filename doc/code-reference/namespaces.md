# Namespaces


We use namespaces to delimit the visibility of functions and classes defined in
the various subdirectories of the project.
Namespaces provide a convenient layered structure to the project and we use
them as a convention to signal which functions and classes are supposed to be
used in any given layer.
The top-level namespace is called `pcm` and includes all functions and classes
that can be called from the outside world, i.e. a C++ API.
Each subdirectory introduces a new namespace of the same name, nested into `pcm`.
Code that can be used _outside_ of a given subdirectory is put directly in the
`pcm` namespace, i.e. the outermost layer.
Finally, the namespace `detail`, at the third level of nesting, is used for
functions and classes that are used exclusively within the code in a given
subdirectory.
