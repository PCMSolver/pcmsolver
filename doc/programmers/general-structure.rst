General Structure
=================

lblabmaaer

Coding standards
----------------

General Object-Oriented design principles you should try to follow:
        1. Identify the aspects of your application that vary and separate them from what stays the same;
        2. Program to an interface, not an implementation;
        3. Favor composition over inheritance;
        4. Strive for loosely coupled designs between objects that interact;
        5. Classes should be open for extension, but closed for modification;
        6. Depend upon abstractions. Do not depend upon concrete classes;
        7. Principle of Least Knowledge. Talk only to your immediate friends;

[Sutter2004]_, [Cline1998]_, [CppFAQs]_

Including header files
......................

Do not include header files unnecessarily. Even if PCMSolver is not a big project, unnecessary include directives and/or forward declarations
introduce nasty interdependencies among different parts of the code. 
This reflects mainly in longer compilation times, but also in uglier looking code.

Follow these guidelines to decide whether to include or forward declare:
        1. class A makes no reference to class B. Neither include nor forward declare B;
        2. class A refers to class B as a friend. Neither include nor forward declare B;
        3. class A contains a pointer/reference to a class B object. Forward declare B;
        4. class A contains functions with a class B object (value/pointer/reference) as parameter/return value. Forward declare B;
        5. class A is derived from class B. include B;
        6. class A contains a class B object. include B.

.. code-block:: cpp
    
    #ifndef MYCLASS_H
    #define MYCLASS_H

    //==============================
    // Forward declared dependencies
    class Foo;
    class Bar;
    
    //==============================
    // Included dependencies
    #include <vector>
    #include "Parent.h"

    //==============================
    // The actual class
    class MyClass : public Parent         // Parent object, so #include "Parent.h"
    {
        public:
                std::vector<int> avector; // vector object, so #include <vector>
                Foo * foo;                // Foo pointer, so forward declare
                void Func(Bar & bar);     // Bar reference as parameter, so forward declare

                friend class MyFriend;    // friend declaration is not a dependency
                                          //    don't do anything about MyFriend
    };                                  

    #endif // MYCLASS_H
    
    


