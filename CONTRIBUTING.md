# Visualizations: How to Contribute

### Making an Issue
First - before you start working on a change to the repository, please make
sure that there are no exisiting
[issues](https://github.com/SciLifeLab/visualizations/issues) relating to
whatever change you intend to make. If there aren't, please create one
so that others know that you're working on something.

### Workflow
The workflow for adding to this repository should be as follows:

1. [Create an issue](https://github.com/SciLifeLab/visualizations/issues)
   describing what you intend to work on
2. Fork this repository to your GitHub account
3. If adding a new visualization:
    1. Create a new subdirectory named after your script
    2. Add all files relating to that code to that directory
    3. Create example plots in the region of _800 px_ by _600 px_ and add these
       to the `examples` directory
    4. Write a `README.md` file and place this in your code's subdirectory
       _(see below for instructions about what's required)_
    5. Add a link and a thumbnail to the root `README.md` file, linking to the
       subdirectory (and readme) of your code
5. Submit a Pull Request and wait for the person who is currently responsible for
    PRs to review and merge the code.

### Readme Requirements
The readme that you write for your script is almost as important as the code itself.
It is key to making your code useful to others. As a minimum, it should include the
following sections (in the following order):

1. Introduction
2. Example Output
   1. Remember to use _relative_ links when embedding images and content into your
      Readme file, so that it works on any fork. See other Readme files for examples
      or the [GitHub blog post](https://help.github.com/articles/relative-links-in-readmes/)
3. Usage
4. Parameters _(command line or argument)_
5. Dependencies
6. Credits _(generic, copy from another Readme)_

Please look at some of the other Readme files in the repository before you start
to ensure that your documentation is consistent. Readme files should use Markdown
formatting - for help, see
[examples of GitHub markdown syntax](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
and try using the excellent [http://dillinger.io/](http://dillinger.io/)
markdown editor.

### Adding to the root Readme file
Finally, please add a link to your code's subdirectory to the root `README.md`
file with a one-line description of what the plot does. Next, add an example
of the plot output to the tables below. Avoid multiple variants of the same
plot, only use more than one image if very different outputs are generated.
See the [current Readme](https://raw.githubusercontent.com/SciLifeLab/visualizations/master/README.md)
code for examples of how to add these - they should be in a table cell which 
is 50% wide (`<td width="50%">`) with an accompanying header cell (`<th>`). Both
image and header should link to your code's subdirectory.

### Getting Help
If you have any queries, please get in touch with [Phil Ewels](https://github.com/ewels).
Thanks for contributing to our code!

