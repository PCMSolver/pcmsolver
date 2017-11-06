# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into PCMSolver.
Our contribution guide is based on [Psi4 contribution guide](https://github.com/psi4/psi4/blob/master/.github/CONTRIBUTING.md)

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) the
  [PCMSolver/pcmsolver](https://github.com/PCMSolver/pcmsolver) repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the PCMSolver repository.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of PCMSolver on GitHub and open a
  [pull request](https://help.github.com/articles/using-pull-requests/) (PR)
  __towards the `release/1.Y` branch__.
  Note that after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.
  Each commit added to the PR will be validated for mergeability, compilation
  and test suite compliance; the results of these tests will be visible on the
  PR page.
* The title of your pull request should be marked with `[WIP]` if it’s a work
  in progress and with `#trivial` if it is a set of trivial changes.
* If you're providing a new feature, you must add test cases, documentation and
  update the `CHANGELOG.md` file.
* When the code is ready to go, make sure you run the full or relevant portion
  of the test suite on your local machine to check that nothing is broken.
* When you're ready to be considered for merging, check the "Ready to go" box
  on the PR page to let the PCMSolver team know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration (Travis for Linux and Mac) returns checkmarks, and multiple core
  developers give "Approved" reviews.

## Pull Request Requirements

The project is integrated with [Danger.Systems](http://danger.systems/ruby/).
On each PR, one CI job will run the integration and a [bot](https://github.com/minazobot) will
report which requirements are **not met** in your PR.
These reports can be _warnings_ and _errors_. You will discuss and solve both
of them with the reviewers.
The automatic rules are laid out in the `Dangerfile` and are used to enforce an
adequate level of testing, documentation and code quality.

### Danger.Systems Warnings

- PRs classed as Work in Progress.
- Codebase was modified, but no tests were added.
- Nontrivial changes to the codebase, but no documentation added.
- Codebase was modified, but `CHANGELOG.md` was not updated.
- Source files were added or removed, but `.gitattributes` was not updated.

### Danger.Systems Errors

- Commit message linting, based on some of [these recommendations](https://chris.beams.io/posts/git-commit/):
  - Commit subject is more than one word.
  - Commit subject is no longer than 50 characters.
  - Commit subject and body are separated by an empty line.

- Clean commit history, without merge commits.

- Code style for `.hpp`, `.cpp`, `.h` files follows the conventions in
  `.clang-format`.

## Licensing

We do not require any formal copyright assignment or contributor license
agreement.
**Any contributions intentionally sent upstream are presumed to be offered under
terms of the OSI-approved LGPLv3 License.**

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
