function(configure_script infile outfile)
    configure_file(${infile} ${outfile} @ONLY)
    execute_process(COMMAND
        chmod 755 ${outfile} OUTPUT_QUIET)
endfunction()
