General Structure
=================

lblabmaaer

Coding standards
----------------

Including header files
......................

Do not include header files unnecessarily. 
Follow these guidelines to decide whether to include or forward declare:
    1. class A makes no reference to class B. Do not include B;
    2. class A refers to class B as a friend. Forward declare B;
    3. class A contains functions with a class B object (value/pointer/reference) as parameter/return value. include B;
    4. class A is derived from class B. include B;
    5. class A contains a class B object. include B.

.. code-block:: cpp

    #include <iostream>
    
    


