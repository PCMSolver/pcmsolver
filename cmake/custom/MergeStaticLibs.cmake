#    Copyright (C) 2012 Modelon AB

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the BSD style license.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    FMILIB_License.txt file for more details.

#    You should have received a copy of the FMILIB_License.txt file
#    along with this program. If not, contact Modelon AB <http://www.modelon.com>.


# heavily modified by bast at kth dot se
# adapted to PCMSolver by roberto dot d dot remigio at uit dot no

# merge_static_libs(outlib lib1 lib2 ... libn) merges a number of static
# libs into a single static library
function(merge_static_libs outlib )
    set(libs ${ARGV})
    list(REMOVE_AT libs 0)
    # Create a dummy file that the target will depend on
    set(dummyfile ${CMAKE_CURRENT_BINARY_DIR}/${outlib}_dummy.c)
    file(WRITE ${dummyfile} "const char * dummy = \"${dummyfile}\";")

    add_library(${outlib} STATIC ${dummyfile})

    # First get the file names of the libraries to be merged
    foreach(lib ${libs})
        get_target_property(libtype ${lib} TYPE)
        if(NOT libtype STREQUAL "STATIC_LIBRARY")
            message(FATAL_ERROR "merge_static_libs can only process static libraries")
        endif()
        set(libfile ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${lib}.a)
        list(APPEND libfiles "${libfile}")
    endforeach()
    list(REMOVE_DUPLICATES libfiles)

    set(outfile ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${outlib}.a)
    foreach(lib ${libfiles})
        # objlistfile will contain the list of object files for the library
        set(objlistfile ${lib}.objlist)
        set(objdir ${lib}.objdir)
        set(objlistcmake  ${objlistfile}.cmake)
        # we only need to extract files once
        if(${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.check_cache IS_NEWER_THAN ${objlistcmake})
            #---------------------------------
            FILE(WRITE ${objlistcmake}
                "# Extract object files from the library
                message(STATUS \"Extracting object files from ${lib}\")
                EXECUTE_PROCESS(COMMAND ${CMAKE_AR} -x ${lib}
                    WORKING_DIRECTORY ${objdir})
                # save the alphabetically sorted list of object files
                EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c  \"import os; print('\\\\n'.join(sorted(os.listdir(os.curdir), key=str.lower)))\"
                OUTPUT_FILE ${objlistfile}
                WORKING_DIRECTORY ${objdir})"
                )
            #---------------------------------
            file(MAKE_DIRECTORY ${objdir})
            add_custom_command(
                OUTPUT ${objlistfile}
                COMMAND ${CMAKE_COMMAND} -P ${objlistcmake}
                DEPENDS ${lib}
                )
        endif()
        list(APPEND extrafiles "${objlistfile}")
        # relative path is needed by ar under MSYS
        file(RELATIVE_PATH objlistfilerpath ${objdir} ${objlistfile})
        add_custom_command(TARGET ${outlib} POST_BUILD
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/MergeStaticLibs.py
            ${CMAKE_AR} ${outfile} ${objlistfilerpath}
            WORKING_DIRECTORY ${objdir}
            )
    endforeach()
    add_custom_command(TARGET ${outlib} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "Running: ${CMAKE_RANLIB} ${outfile}"
        COMMAND ${CMAKE_RANLIB} ${outfile}
        )
    file(WRITE ${dummyfile}.base "const char* ${outlib}_sublibs=\"${libs}\";")
    add_custom_command(
        OUTPUT  ${dummyfile}
        COMMAND ${CMAKE_COMMAND}  -E copy ${dummyfile}.base ${dummyfile}
        DEPENDS ${libs} ${extrafiles}
        )
endfunction()
