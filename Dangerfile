# Adapted from https://github.com/capistrano/danger/blob/master/Dangerfile
# Q: What is a Dangerfile, anyway? A: See http://danger.systems/

# ------------------------------------------------------------------------------
# Additional pull request data
# ------------------------------------------------------------------------------
pr_number = github.pr_json["number"]
pr_url = github.pr_json["_links"]["html"]["href"]
# Sometimes its a README fix, or something like that - which isn't relevant for
# including in a CHANGELOG for example
declared_trivial = (github.pr_title + github.pr_body).include?("#trivial")

# Just to let people know
warn("PR is classed as Work in Progress") if github.pr_title.include? "[WIP]"

# Ensure a clean commits history
if git.commits.any? { |c| c.message =~ /^Merge branch '#{github.branch_for_base}'/ }
  fail('Please rebase to get rid of the merge commits in this PR')
end
commit_lint.check disable: [:subject_cap, :subject_period]

# Check code style with clang-format
code_style_validation.check file_extensions: ['.hpp', '.cpp', '.h'], ignore_file_patterns: [/^external\//]

# ------------------------------------------------------------------------------
# What changed?
# ------------------------------------------------------------------------------
# Were source files added? Were source files deleted?
# We look for added and removed files in the api, include, src, tests and tools directories.
# Files whose license header might need to be updated are kept here.
has_added_files = !git.added_files.grep(/^(api|include|src|tests|tools)/).empty?
has_deleted_files = !git.deleted_files.grep(/^(api|include|src|tests|tools)/).empty?
# Was any code modified in the api, cmake, external, include, src and tools directories?
has_code_changes = !git.modified_files.grep(/^(api|cmake|external|include|src|tools)/).empty?
# Was any code modified in the tests directory?
has_test_changes = !git.modified_files.grep(/^tests/).empty?
# Was the CHANGELOG.md file modified?
has_changelog_changes = git.modified_files.include?("CHANGELOG.md")
# Was the .gitattributes file modified?
# .gitattributes is used to keep track of license headers and has to be updated
# if files are added or removed
has_gitattributes_changes = git.modified_files.include?(".gitattributes")
# Was documentation added?
has_doc_changes = !git.modified_files.grep(/^doc/).empty?

# ------------------------------------------------------------------------------
# You've made changes to api|cmake|external|include|src|tools,
# but didn't write any tests?
# ------------------------------------------------------------------------------
if has_code_changes && !has_test_changes
  if %w(tests).any? { |dir| Dir.exist?(dir) }
    warn("There are code changes, but no corresponding tests. "\
         "Please include tests if this PR introduces any modifications in "\
         "behavior.",
         :sticky => false)
  else
    markdown <<-MARKDOWN
Thanks for the PR! This project lacks automated tests, which makes reviewing and approving PRs somewhat difficult. Please make sure that your contribution has not broken backwards compatibility or introduced any risky changes.

MARKDOWN
  end
end

# ------------------------------------------------------------------------------
# You've made nontrivial changes to api|cmake|external|include|src|tools,
# but didn't write any docs?
# ------------------------------------------------------------------------------
doc_changes_recommended = git.insertions > 15
if has_code_changes && !has_doc_changes && doc_changes_recommended && !declared_trivial
  warn("Consider adding supporting documentation to this change. Documentation sources can be found in the `doc` directory.")
end

# ------------------------------------------------------------------------------
# Have you updated CHANGELOG.md?
# ------------------------------------------------------------------------------
if !has_changelog_changes && has_code_changes
  markdown <<-MARKDOWN
Here's an example of a CHANGELOG.md entry:

```markdown
* [##{pr_number}](#{pr_url}): #{github.pr_title} - [@#{github.pr_author}](https://github.com/#{github.pr_author})
```
MARKDOWN
  warn("Please update CHANGELOG.md with a description of your changes. "\
       "If this PR is not a user-facing change (e.g. just refactoring), "\
       "you can disregard this.", :sticky => false)
end

# ------------------------------------------------------------------------------
# You've made changes to api|include|src|tests|tools,
# but didn't update .gitattributes?
# ------------------------------------------------------------------------------
if has_added_files && has_deleted_files && !has_gitattributes_changes
  warn("You have added source files without updating `.gitattributes`.", :sticky => false)
end
if has_deleted_files && !has_gitattributes_changes
  warn("You have removed source files without updating `.gitattributes`.", :sticky => false)
end

lgtm.check_lgtm
