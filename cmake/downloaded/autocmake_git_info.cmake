# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

#.rst:
#
# Creates git_info.h in the build directory.
# This file can be included in sources to print
# Git repository version and status information
# to the program output.
#
# autocmake.yml configuration::
#
#   url_root: https://github.com/coderefinery/autocmake/raw/master/
#   fetch:
#     - "%(url_root)modules/git_info/git_info.h.in"

get_filename_component(_current_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

function(generate_git_info_header _header_location _header_name)
  # _header_location: where the Git info header file should be generated
  # _header_name: the Git info header name, complete with extension (.h, .hpp, .hxx or whatever)
  find_package(Git)

  set(_git_last_commit_hash "unknown")
  set(_git_last_commit_author "unknown")
  set(_git_last_commit_date "unknown")
  set(_git_branch "unknown")

  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%h -n 1
      OUTPUT_VARIABLE _git_last_commit_hash
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%aN -n 1
      OUTPUT_VARIABLE _git_last_commit_author
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%ad -n 1
      OUTPUT_VARIABLE _git_last_commit_date
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      OUTPUT_VARIABLE _git_branch
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )
  endif()

  configure_file(
    ${_current_dir}/git_info.h.in
    ${_header_location}/${_header_name}
    @ONLY
    )

  unset(_git_last_commit_hash)
  unset(_git_last_commit_author)
  unset(_git_last_commit_date)
  unset(_git_branch)

  add_custom_target(
    git_info
    ALL DEPENDS ${_header_location}/${_header_name}
    )
endfunction()
